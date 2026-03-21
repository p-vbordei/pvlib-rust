/// Solves the single-diode equation for current `I` given voltage `V` using the Newton-Raphson method.
///
/// `I = I_L - I_0 * (exp((V + I*R_s) / (n*N_s*V_th)) - 1) - (V + I*R_s) / R_sh`
///
/// # Arguments
/// * `v` - Given voltage constraint in Volts.
/// * `photocurr` - Light-generated current `I_L` in Amperes.
/// * `saturation_curr` - Diode saturation current `I_0` in Amperes.
/// * `resistance_series` - Series resistance `R_s` in Ohms.
/// * `resistance_shunt` - Shunt resistance `R_sh` in Ohms.
/// * `n_ns_vth` - Product of diode ideality factor `n`, series cells `N_s`, and thermal voltage `V_th` in Volts.
///
/// # Returns
/// Current `I` in Amperes.
pub fn i_from_v(
    v: f64,
    photocurr: f64,
    saturation_curr: f64,
    resistance_series: f64,
    resistance_shunt: f64,
    n_ns_vth: f64,
) -> f64 {
    let mut i = photocurr - v / resistance_shunt;

    for _ in 0..100 {
        let exp_term = ((v + i * resistance_series) / n_ns_vth).exp();
        let f = photocurr
            - saturation_curr * (exp_term - 1.0)
            - (v + i * resistance_series) / resistance_shunt
            - i;
        let df_di = -saturation_curr * (resistance_series / n_ns_vth) * exp_term
            - (resistance_series / resistance_shunt)
            - 1.0;
        let step = f / df_di;
        i -= step;
        if step.abs() < 1e-6 {
            break;
        }
    }

    i
}

/// Solves the single-diode equation for voltage `V` given current `I`.
///
/// # Returns
/// Voltage `V` in Volts.
pub fn v_from_i(
    i: f64,
    photocurr: f64,
    saturation_curr: f64,
    resistance_series: f64,
    resistance_shunt: f64,
    n_ns_vth: f64,
) -> f64 {
    let mut v = n_ns_vth * ((photocurr - i) / saturation_curr).ln() - i * resistance_series;
    if v.is_nan() {
        v = 0.0;
    }

    for _ in 0..100 {
        let exp_term = ((v + i * resistance_series) / n_ns_vth).exp();
        let f = photocurr
            - saturation_curr * (exp_term - 1.0)
            - (v + i * resistance_series) / resistance_shunt
            - i;
        let df_dv = -saturation_curr / n_ns_vth * exp_term - 1.0 / resistance_shunt;
        let step = f / df_dv;
        v -= step;
        if step.abs() < 1e-6 {
            break;
        }
    }

    v
}

// ---------------------------------------------------------------------------
// Bishop88 single diode model
// ---------------------------------------------------------------------------

/// Default built-in voltage per cell junction for a:Si, CdTe (Mertens et al.)
pub const VOLTAGE_BUILTIN: f64 = 0.9;

const NEWTON_TOL: f64 = 1e-6;
const NEWTON_MAXITER: usize = 100;

/// Rough estimate of open circuit voltage.
///
/// Assumes infinite shunt resistance and zero series resistance:
/// `Voc_est = nNsVth * ln(Iph / I0 + 1)`
pub fn estimate_voc(photocurrent: f64, saturation_current: f64, n_ns_vth: f64) -> f64 {
    n_ns_vth * (photocurrent / saturation_current + 1.0).ln()
}

/// Output of the Bishop88 model at one operating point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bishop88Output {
    pub current: f64,
    pub voltage: f64,
    pub power: f64,
}

/// Gradients from the Bishop88 model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bishop88Gradients {
    /// dI/dVd
    pub di_dvd: f64,
    /// dV/dVd
    pub dv_dvd: f64,
    /// dI/dV
    pub di_dv: f64,
    /// dP/dV
    pub dp_dv: f64,
    /// d2P/(dV dVd)
    pub d2p_dvd: f64,
}

/// Parameters for the Bishop88 model including optional breakdown terms.
#[derive(Debug, Clone, Copy)]
pub struct Bishop88Params {
    pub photocurrent: f64,
    pub saturation_current: f64,
    pub resistance_series: f64,
    pub resistance_shunt: f64,
    pub n_ns_vth: f64,
    /// PVSyst recombination parameter (d^2 / mu*tau). Default 0.
    pub d2mutau: f64,
    /// Ns * Vbi. Default f64::INFINITY.
    pub ns_vbi: f64,
    /// Fraction of ohmic current in avalanche breakdown. Default 0.
    pub breakdown_factor: f64,
    /// Reverse breakdown voltage. Default -5.5 V.
    pub breakdown_voltage: f64,
    /// Avalanche breakdown exponent. Default 3.28.
    pub breakdown_exp: f64,
}

impl Bishop88Params {
    /// Create params with no recombination or breakdown terms.
    pub fn new(
        photocurrent: f64,
        saturation_current: f64,
        resistance_series: f64,
        resistance_shunt: f64,
        n_ns_vth: f64,
    ) -> Self {
        Self {
            photocurrent,
            saturation_current,
            resistance_series,
            resistance_shunt,
            n_ns_vth,
            d2mutau: 0.0,
            ns_vbi: f64::INFINITY,
            breakdown_factor: 0.0,
            breakdown_voltage: -5.5,
            breakdown_exp: 3.28,
        }
    }
}

/// Explicit calculation of (I, V, P) on the IV curve using Bishop (1988).
///
/// `diode_voltage` is the voltage across the diode (Vd = V + I*Rs).
pub fn bishop88(diode_voltage: f64, p: &Bishop88Params) -> Bishop88Output {
    let (out, _) = bishop88_inner(diode_voltage, p, false);
    out
}

/// Bishop88 with gradients.
pub fn bishop88_with_gradients(
    diode_voltage: f64,
    p: &Bishop88Params,
) -> (Bishop88Output, Bishop88Gradients) {
    let (out, grad) = bishop88_inner(diode_voltage, p, true);
    (out, grad.unwrap())
}

fn bishop88_inner(
    vd: f64,
    p: &Bishop88Params,
    gradients: bool,
) -> (Bishop88Output, Option<Bishop88Gradients>) {
    // recombination loss current
    let (i_recomb, v_recomb) = if p.d2mutau > 0.0 {
        let vr = p.ns_vbi - vd;
        (p.photocurrent * p.d2mutau / vr, vr)
    } else {
        (0.0, f64::INFINITY)
    };

    let v_star = vd / p.n_ns_vth;
    let g_sh = 1.0 / p.resistance_shunt;

    // breakdown term (only compute when breakdown_factor is nonzero)
    let (brk_pwr, i_breakdown) = if p.breakdown_factor != 0.0 {
        let brk_term = 1.0 - vd / p.breakdown_voltage;
        if brk_term <= 0.0 {
            (f64::INFINITY, f64::INFINITY)
        } else {
            let bp = brk_term.powf(-p.breakdown_exp);
            (bp, p.breakdown_factor * vd * g_sh * bp)
        }
    } else {
        (1.0, 0.0)
    };

    let i = p.photocurrent - p.saturation_current * v_star.exp_m1() - vd * g_sh - i_recomb
        - i_breakdown;
    let v = vd - i * p.resistance_series;
    let out = Bishop88Output {
        current: i,
        voltage: v,
        power: i * v,
    };

    if !gradients {
        return (out, None);
    }

    let grad_i_recomb = if p.d2mutau > 0.0 {
        i_recomb / v_recomb
    } else {
        0.0
    };
    let grad_2i_recomb = if p.d2mutau > 0.0 {
        2.0 * grad_i_recomb / v_recomb
    } else {
        0.0
    };

    let g_diode = p.saturation_current * v_star.exp() / p.n_ns_vth;

    let (grad_i_brk, grad2i_brk) = if p.breakdown_factor != 0.0 {
        let brk_term = 1.0 - vd / p.breakdown_voltage;
        if brk_term <= 0.0 {
            (f64::INFINITY, f64::INFINITY)
        } else {
            let brk_pwr_1 = brk_term.powf(-p.breakdown_exp - 1.0);
            let brk_pwr_2 = brk_term.powf(-p.breakdown_exp - 2.0);
            let brk_fctr = p.breakdown_factor * g_sh;
            let gi = brk_fctr * (brk_pwr + vd * (-p.breakdown_exp) * brk_pwr_1);
            let g2i = brk_fctr
                * (-p.breakdown_exp)
                * (2.0 * brk_pwr_1 + vd * (-p.breakdown_exp - 1.0) * brk_pwr_2);
            (gi, g2i)
        }
    } else {
        (0.0, 0.0)
    };

    let di_dvd = -g_diode - g_sh - grad_i_recomb - grad_i_brk;
    let dv_dvd = 1.0 - di_dvd * p.resistance_series;
    let di_dv = di_dvd / dv_dvd;
    let dp_dv = v * di_dv + i;

    let d2i_dvd = -g_diode / p.n_ns_vth - grad_2i_recomb - grad2i_brk;
    let d2v_dvd = -d2i_dvd * p.resistance_series;
    let d2p_dvd = dv_dvd * di_dv + v * (d2i_dvd / dv_dvd - di_dvd * d2v_dvd / (dv_dvd * dv_dvd))
        + di_dvd;

    let grad = Bishop88Gradients {
        di_dvd,
        dv_dvd,
        di_dv,
        dp_dv,
        d2p_dvd,
    };

    (out, Some(grad))
}

/// Find current at a given voltage using Bishop88 + Newton-Raphson.
///
/// Searches for the diode voltage `vd` such that `bishop88(vd).voltage == voltage`.
pub fn bishop88_i_from_v(voltage: f64, p: &Bishop88Params) -> f64 {
    // Start from voltage as initial guess (matching pvlib-python)
    let mut vd = voltage;
    // Only clamp to ns_vbi if it is finite, positive, and vd exceeds it
    if p.ns_vbi.is_finite() && p.ns_vbi > 0.0 && vd >= p.ns_vbi {
        vd = 0.9999 * p.ns_vbi;
    }

    for _ in 0..NEWTON_MAXITER {
        let (out, grad) = bishop88_inner(vd, p, true);
        let grad = grad.unwrap();
        let residual = out.voltage - voltage;
        let step = residual / grad.dv_dvd;
        vd -= step;
        if step.abs() < NEWTON_TOL {
            break;
        }
    }

    bishop88(vd, p).current
}

/// Find voltage at a given current using Bishop88 + Newton-Raphson.
///
/// Searches for the diode voltage `vd` such that `bishop88(vd).current == current`.
pub fn bishop88_v_from_i(current: f64, p: &Bishop88Params) -> f64 {
    let voc_est = estimate_voc(p.photocurrent, p.saturation_current, p.n_ns_vth);
    let mut vd = if p.ns_vbi.is_finite() && p.ns_vbi > 0.0 && voc_est >= p.ns_vbi {
        0.9999 * p.ns_vbi
    } else {
        voc_est
    };

    for _ in 0..NEWTON_MAXITER {
        let (out, grad) = bishop88_inner(vd, p, true);
        let grad = grad.unwrap();
        let residual = out.current - current;
        // di/dvd is the derivative of current w.r.t. diode voltage
        let step = residual / grad.di_dvd;
        vd -= step;
        if step.abs() < NEWTON_TOL {
            break;
        }
    }

    bishop88(vd, p).voltage
}

/// Find the maximum power point using Bishop88 + Newton-Raphson.
///
/// Searches for diode voltage where dP/dV = 0.
pub fn bishop88_mpp(p: &Bishop88Params) -> Bishop88Output {
    let voc_est = estimate_voc(p.photocurrent, p.saturation_current, p.n_ns_vth);
    let mut vd = if p.ns_vbi.is_finite() && p.ns_vbi > 0.0 && voc_est >= p.ns_vbi {
        0.9999 * p.ns_vbi
    } else {
        voc_est
    };

    for _ in 0..NEWTON_MAXITER {
        let (_, grad) = bishop88_inner(vd, p, true);
        let grad = grad.unwrap();
        if grad.d2p_dvd.abs() < 1e-30 {
            break;
        }
        // root-find dP/dV = 0, using d2P/(dV*dVd) as the derivative
        let step = grad.dp_dv / grad.d2p_dvd;
        vd -= step;
        // clamp vd to reasonable range
        if vd < 0.0 {
            vd = 0.0;
        }
        if step.abs() < NEWTON_TOL {
            break;
        }
    }

    bishop88(vd, p)
}
