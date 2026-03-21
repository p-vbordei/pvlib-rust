/// IV curve tools: utilities and parameter fitting for the single diode model.
///
/// Ported from pvlib-python `pvlib.ivtools`.

/// Result of fitting the single diode model to an IV curve or module specs.
#[derive(Debug, Clone)]
pub struct SdmFitResult {
    /// Light-generated (photo) current [A]
    pub photocurrent: f64,
    /// Diode saturation (dark) current [A]
    pub saturation_current: f64,
    /// Series resistance [ohm]
    pub resistance_series: f64,
    /// Shunt (parallel) resistance [ohm]
    pub resistance_shunt: f64,
    /// Product of diode ideality factor, cells in series, and thermal voltage [V]
    pub n_ns_vth: f64,
}

// ---------------------------------------------------------------------------
// rectify_iv_curve
// ---------------------------------------------------------------------------

/// Sort an IV curve by voltage, remove NaN/negative values, and merge duplicate voltages.
///
/// The returned vectors are sorted by voltage, contain only non-negative values,
/// and have no duplicate voltage entries (duplicates are averaged).
pub fn rectify_iv_curve(voltage: &[f64], current: &[f64]) -> (Vec<f64>, Vec<f64>) {
    // Pair up, filter NaN and negatives
    let mut pairs: Vec<(f64, f64)> = voltage
        .iter()
        .zip(current.iter())
        .filter(|(v, i)| v.is_finite() && i.is_finite() && **v >= 0.0 && **i >= 0.0)
        .map(|(v, i)| (*v, *i))
        .collect();

    // Sort by voltage, then descending current
    pairs.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal))
    });

    if pairs.is_empty() {
        return (vec![], vec![]);
    }

    // Merge duplicate voltages by averaging current
    let mut out_v: Vec<f64> = Vec::with_capacity(pairs.len());
    let mut out_i: Vec<f64> = Vec::with_capacity(pairs.len());

    let mut cur_v = pairs[0].0;
    let mut cur_sum = pairs[0].1;
    let mut cur_count = 1usize;

    for &(v, i) in &pairs[1..] {
        if (v - cur_v).abs() < f64::EPSILON * cur_v.abs().max(1.0) {
            cur_sum += i;
            cur_count += 1;
        } else {
            out_v.push(cur_v);
            out_i.push(cur_sum / cur_count as f64);
            cur_v = v;
            cur_sum = i;
            cur_count = 1;
        }
    }
    out_v.push(cur_v);
    out_i.push(cur_sum / cur_count as f64);

    (out_v, out_i)
}

// ---------------------------------------------------------------------------
// fit_sandia_simple
// ---------------------------------------------------------------------------

/// Fit the single diode equation to an IV curve using the Sandia simplified method.
///
/// The IV curve must be sorted by increasing voltage from 0 to `v_oc`, with
/// current decreasing from `i_sc` to 0.
///
/// `vlim` defines the fraction of `v_oc` below which the exponential term is
/// neglected (linear region). `ilim` defines the fraction of `i_sc` used to
/// identify the exponential region.
///
/// # Errors
/// Returns `Err` if parameter extraction fails (e.g. insufficient data,
/// negative slopes not found, or saturation current cannot be determined).
pub fn fit_sandia_simple(
    voltage: &[f64],
    current: &[f64],
    v_oc: Option<f64>,
    i_sc: Option<f64>,
    v_mp_i_mp: Option<(f64, f64)>,
    vlim: f64,
    ilim: f64,
) -> Result<SdmFitResult, String> {
    let n = voltage.len();
    if n < 6 || current.len() != n {
        return Err("Need at least 6 matching voltage/current points".into());
    }

    let v_oc = v_oc.unwrap_or(voltage[n - 1]);
    let i_sc = i_sc.unwrap_or(current[0]);

    let (v_mp, i_mp) = v_mp_i_mp.unwrap_or_else(|| {
        let mut best_idx = 0;
        let mut best_p = voltage[0] * current[0];
        for k in 1..n {
            let p = voltage[k] * current[k];
            if p > best_p {
                best_p = p;
                best_idx = k;
            }
        }
        (voltage[best_idx], current[best_idx])
    });

    // Step 1-3: linear fit of low-voltage region → beta0, beta1
    let (beta0, beta1) = sandia_beta0_beta1(voltage, current, vlim, v_oc)?;

    // Step 4-5: exponential fit → beta3, beta4
    let (beta3, beta4) = sandia_beta3_beta4(voltage, current, beta0, beta1, ilim, i_sc)?;

    // Step 6: calculate parameters
    sandia_simple_params(beta0, beta1, beta3, beta4, v_mp, i_mp, v_oc)
}

/// Linear fit on the low-voltage region of the IV curve.
fn sandia_beta0_beta1(
    v: &[f64],
    i: &[f64],
    vlim: f64,
    v_oc: f64,
) -> Result<(f64, f64), String> {
    let threshold = vlim * v_oc;
    // Find first index where v >= threshold
    let first_idx = v.iter().position(|&x| x >= threshold).unwrap_or(3).max(3);

    for idx in first_idx..=v.len() {
        let (slope, intercept) = polyfit1(&v[..idx], &i[..idx]);
        if slope < 0.0 {
            return Ok((intercept, -slope));
        }
    }
    Err("Parameter extraction failed: could not determine beta0, beta1 from linear region".into())
}

/// Exponential fit on the high-current-deficit region.
fn sandia_beta3_beta4(
    voltage: &[f64],
    current: &[f64],
    beta0: f64,
    beta1: f64,
    ilim: f64,
    i_sc: f64,
) -> Result<(f64, f64), String> {
    // y = beta0 - beta1*V - I; select points where y > ilim * i_sc
    let n = voltage.len();
    let mut xv: Vec<[f64; 3]> = Vec::new();
    let mut yv: Vec<f64> = Vec::new();

    for k in 0..n {
        let y = beta0 - beta1 * voltage[k] - current[k];
        if y > ilim * i_sc {
            xv.push([1.0, voltage[k], current[k]]);
            yv.push(y.ln());
        }
    }

    if xv.len() < 3 {
        return Err("Parameter extraction failed: insufficient points in exponential region".into());
    }

    // Least-squares: yv = xv * [beta2, beta3, beta4]^T
    let coef = lstsq_3(&xv, &yv)?;
    let beta3 = coef[1];
    let beta4 = coef[2];

    if beta3.is_nan() || beta4.is_nan() {
        return Err(format!(
            "Parameter extraction failed: beta3={}, beta4={}",
            beta3, beta4
        ));
    }
    Ok((beta3, beta4))
}

/// Calculate SDM parameters from regression coefficients.
fn sandia_simple_params(
    beta0: f64,
    beta1: f64,
    beta3: f64,
    beta4: f64,
    v_mp: f64,
    i_mp: f64,
    v_oc: f64,
) -> Result<SdmFitResult, String> {
    let n_ns_vth = 1.0 / beta3;
    let rs = beta4 / beta3;
    let gsh = beta1 / (1.0 - rs * beta1);
    let rsh = 1.0 / gsh;
    let iph = (1.0 + gsh * rs) * beta0;

    let io_vmp = calc_i0(v_mp, i_mp, iph, gsh, rs, n_ns_vth);
    let io_voc = calc_i0(v_oc, 0.0, iph, gsh, rs, n_ns_vth);

    let io = if io_vmp > 0.0 && io_voc > 0.0 {
        0.5 * (io_vmp + io_voc)
    } else if io_vmp > 0.0 {
        io_vmp
    } else if io_voc > 0.0 {
        io_voc
    } else {
        return Err("Parameter extraction failed: I0 is undetermined".into());
    };

    Ok(SdmFitResult {
        photocurrent: iph,
        saturation_current: io,
        resistance_series: rs,
        resistance_shunt: rsh,
        n_ns_vth,
    })
}

fn calc_i0(voltage: f64, current: f64, iph: f64, gsh: f64, rs: f64, n_ns_vth: f64) -> f64 {
    let x = (voltage + rs * current) / n_ns_vth;
    let denom = x.exp() - 1.0;
    if denom.abs() < 1e-30 {
        return f64::NAN;
    }
    (iph - current - gsh * (voltage + rs * current)) / denom
}

// ---------------------------------------------------------------------------
// fit_desoto
// ---------------------------------------------------------------------------

/// Fit the De Soto single diode model from module datasheet specifications.
///
/// Solves a system of 5 nonlinear equations using a hybrid Newton method
/// to determine the five SDM parameters at reference conditions.
///
/// # Arguments
/// * `v_mp` - Voltage at maximum power point [V]
/// * `i_mp` - Current at maximum power point [A]
/// * `v_oc` - Open-circuit voltage [V]
/// * `i_sc` - Short-circuit current [A]
/// * `alpha_sc` - Temperature coefficient of Isc [A/K]
/// * `beta_voc` - Temperature coefficient of Voc [V/K]
/// * `cells_in_series` - Number of cells in series
/// * `eg_ref` - Bandgap energy at reference [eV], default 1.121 for silicon
/// * `d_eg_dt` - Temperature dependence of bandgap [1/K], default -0.0002677
///
/// # Returns
/// `SdmFitResult` with parameters at reference conditions, where `n_ns_vth`
/// corresponds to `a_ref` (modified ideality factor).
///
/// # Errors
/// Returns `Err` if the Newton solver does not converge.
pub fn fit_desoto(
    v_mp: f64,
    i_mp: f64,
    v_oc: f64,
    i_sc: f64,
    alpha_sc: f64,
    beta_voc: f64,
    cells_in_series: i32,
    eg_ref: f64,
    d_eg_dt: f64,
) -> Result<SdmFitResult, String> {
    // Boltzmann constant in eV/K
    const K_EV: f64 = 8.617333262e-5;
    let t_ref = 25.0 + 273.15; // K

    // Initial guesses (Duffie & Beckman, p753)
    let a_0 = 1.5 * K_EV * t_ref * cells_in_series as f64;
    let il_0 = i_sc;
    let io_0 = i_sc * (-v_oc / a_0).exp();
    let rs_0 = {
        let ratio = (il_0 - i_mp) / io_0;
        if ratio > 0.0 {
            (a_0 * (1.0 + ratio).ln() - v_mp) / i_mp
        } else {
            0.1
        }
    };
    let rsh_0 = 100.0;

    let mut params = [il_0, io_0, rs_0, rsh_0, a_0];
    let specs = DesotoSpecs {
        i_sc,
        v_oc,
        i_mp,
        v_mp,
        beta_voc,
        alpha_sc,
        eg_ref,
        d_eg_dt,
        t_ref,
        k: K_EV,
    };

    // Newton-Raphson iteration with damping
    for _ in 0..500 {
        let f = desoto_equations(&params, &specs);
        let j = desoto_jacobian(&params, &specs);

        let delta = match solve_5x5(&j, &f) {
            Some(d) => d,
            None => {
                // Perturb parameters slightly and retry
                for p in params.iter_mut() {
                    *p *= 1.0 + 1e-6;
                }
                continue;
            }
        };

        // Damped step: limit each parameter change to avoid overshooting
        let mut alpha = 1.0_f64;
        for i in 0..5 {
            if params[i].abs() > 1e-30 {
                let rel = (delta[i] / params[i]).abs();
                if rel > 0.5 {
                    alpha = alpha.min(0.5 / rel);
                }
            }
        }

        let mut max_step = 0.0_f64;
        for i in 0..5 {
            let step = alpha * delta[i];
            params[i] -= step;
            let scale = params[i].abs().max(1e-30);
            max_step = max_step.max((step / scale).abs());
        }

        // Ensure I0 stays positive
        if params[1] <= 0.0 {
            params[1] = 1e-15;
        }

        if max_step < 1e-10 {
            return Ok(SdmFitResult {
                photocurrent: params[0],
                saturation_current: params[1],
                resistance_series: params[2],
                resistance_shunt: params[3],
                n_ns_vth: params[4],
            });
        }
    }

    Err("De Soto parameter estimation did not converge".into())
}

struct DesotoSpecs {
    i_sc: f64,
    v_oc: f64,
    i_mp: f64,
    v_mp: f64,
    beta_voc: f64,
    alpha_sc: f64,
    eg_ref: f64,
    d_eg_dt: f64,
    t_ref: f64,
    k: f64,
}

/// Evaluate the 5 De Soto equations. params = [IL, Io, Rs, Rsh, a].
fn desoto_equations(params: &[f64; 5], s: &DesotoSpecs) -> [f64; 5] {
    let (il, io, rs, rsh, a) = (params[0], params[1], params[2], params[3], params[4]);

    // eq1: short circuit
    let y0 = s.i_sc - il + io * ((s.i_sc * rs / a).exp() - 1.0) + s.i_sc * rs / rsh;

    // eq2: open circuit at Tref
    let y1 = -il + io * ((s.v_oc / a).exp() - 1.0) + s.v_oc / rsh;

    // eq3: max power point
    let vrs_mp = s.v_mp + s.i_mp * rs;
    let y2 = s.i_mp - il + io * ((vrs_mp / a).exp() - 1.0) + vrs_mp / rsh;

    // eq4: dP/dV = 0 at MPP
    let exp_mp = (vrs_mp / a).exp();
    let num = s.i_mp - s.v_mp * (io / a * exp_mp + 1.0 / rsh);
    let den = 1.0 + io * rs / a * exp_mp + rs / rsh;
    let y3 = num / den;

    // eq5: open circuit at T2
    let t2 = s.t_ref + 2.0;
    let voc2 = (t2 - s.t_ref) * s.beta_voc + s.v_oc;
    let a2 = a * t2 / s.t_ref;
    let il2 = il + s.alpha_sc * (t2 - s.t_ref);
    let eg2 = s.eg_ref * (1.0 + s.d_eg_dt * (t2 - s.t_ref));
    let io2 = io * (t2 / s.t_ref).powi(3) * ((1.0 / s.k) * (s.eg_ref / s.t_ref - eg2 / t2)).exp();
    let y4 = -il2 + io2 * ((voc2 / a2).exp() - 1.0) + voc2 / rsh;

    [y0, y1, y2, y3, y4]
}

/// Numerical Jacobian via finite differences for the De Soto system.
fn desoto_jacobian(params: &[f64; 5], s: &DesotoSpecs) -> [[f64; 5]; 5] {
    let f0 = desoto_equations(params, s);
    let mut jac = [[0.0; 5]; 5];

    for j in 0..5 {
        let mut p = *params;
        let h = (params[j].abs() * 1e-8).max(1e-14);
        p[j] += h;
        let f1 = desoto_equations(&p, s);
        for i in 0..5 {
            jac[i][j] = (f1[i] - f0[i]) / h;
        }
    }
    jac
}

// ---------------------------------------------------------------------------
// Linear algebra helpers (no external dep needed)
// ---------------------------------------------------------------------------

/// Simple 1st-degree polynomial fit (y = slope*x + intercept) via least squares.
fn polyfit1(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    let sx: f64 = x.iter().sum();
    let sy: f64 = y.iter().sum();
    let sxy: f64 = x.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
    let sx2: f64 = x.iter().map(|a| a * a).sum();

    let denom = n * sx2 - sx * sx;
    if denom.abs() < 1e-30 {
        return (0.0, sy / n);
    }
    let slope = (n * sxy - sx * sy) / denom;
    let intercept = (sy - slope * sx) / n;
    (slope, intercept)
}

/// Least-squares solve for 3 unknowns: A*x = b  where A is Mx3.
fn lstsq_3(a: &[[f64; 3]], b: &[f64]) -> Result<[f64; 3], String> {
    // Normal equations: A^T A x = A^T b
    let m = a.len();
    let mut ata = [[0.0; 3]; 3];
    let mut atb = [0.0; 3];

    for k in 0..m {
        for i in 0..3 {
            atb[i] += a[k][i] * b[k];
            for j in 0..3 {
                ata[i][j] += a[k][i] * a[k][j];
            }
        }
    }

    // Solve 3x3 via Cramer or Gauss elimination
    solve_3x3(&ata, &atb).ok_or_else(|| "Singular matrix in least-squares".to_string())
}

/// Solve 3x3 linear system via Gaussian elimination with partial pivoting.
fn solve_3x3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
    let mut aug = [[0.0; 4]; 3];
    for i in 0..3 {
        for j in 0..3 {
            aug[i][j] = a[i][j];
        }
        aug[i][3] = b[i];
    }

    for col in 0..3 {
        // partial pivot
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in (col + 1)..3 {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        if max_val < 1e-30 {
            return None;
        }
        aug.swap(col, max_row);

        let pivot = aug[col][col];
        for row in (col + 1)..3 {
            let factor = aug[row][col] / pivot;
            for j in col..4 {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    // Back substitution
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        x[i] = aug[i][3];
        for j in (i + 1)..3 {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }
    Some(x)
}

/// Solve 5x5 linear system via Gaussian elimination with partial pivoting.
fn solve_5x5(a: &[[f64; 5]; 5], b: &[f64; 5]) -> Option<[f64; 5]> {
    let mut aug = [[0.0; 6]; 5];
    for i in 0..5 {
        for j in 0..5 {
            aug[i][j] = a[i][j];
        }
        aug[i][5] = b[i];
    }

    // Compute a scale factor for relative pivot check
    let max_abs = aug
        .iter()
        .flat_map(|row| row.iter())
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max)
        .max(1e-300);

    for col in 0..5 {
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in (col + 1)..5 {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        if max_val < max_abs * 1e-15 {
            return None;
        }
        aug.swap(col, max_row);

        let pivot = aug[col][col];
        for row in (col + 1)..5 {
            let factor = aug[row][col] / pivot;
            for j in col..6 {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    let mut x = [0.0; 5];
    for i in (0..5).rev() {
        x[i] = aug[i][5];
        for j in (i + 1)..5 {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }
    Some(x)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Generate a synthetic IV curve from known SDM parameters using the
    /// single diode equation: I = IL - I0*(exp((V+I*Rs)/a) - 1) - (V+I*Rs)/Rsh
    fn generate_iv_curve(
        il: f64,
        i0: f64,
        rs: f64,
        rsh: f64,
        a: f64,
        n_points: usize,
    ) -> (Vec<f64>, Vec<f64>) {
        // Estimate Voc
        let v_oc_est = a * (il / i0).ln();
        let mut voltages = Vec::with_capacity(n_points);
        let mut currents = Vec::with_capacity(n_points);

        for k in 0..n_points {
            let v = v_oc_est * (k as f64) / (n_points as f64 - 1.0);
            // Newton-Raphson to solve for I
            let mut i = il - v / rsh;
            for _ in 0..200 {
                let exp_term = ((v + i * rs) / a).exp();
                let f = il - i0 * (exp_term - 1.0) - (v + i * rs) / rsh - i;
                let df = -i0 * rs / a * exp_term - rs / rsh - 1.0;
                let step = f / df;
                i -= step;
                if step.abs() < 1e-12 {
                    break;
                }
            }
            voltages.push(v);
            currents.push(i.max(0.0));
        }
        (voltages, currents)
    }

    #[test]
    fn test_rectify_iv_curve_basic() {
        let v = vec![3.0, 1.0, 2.0, 1.0, -0.5, f64::NAN];
        let i = vec![0.5, 2.0, 1.0, 1.5, 3.0, 1.0];
        let (rv, ri) = rectify_iv_curve(&v, &i);

        // Should remove NaN and negative voltage entries
        assert_eq!(rv.len(), ri.len());
        // Sorted by voltage
        for w in rv.windows(2) {
            assert!(w[0] <= w[1]);
        }
        // No negative values
        for &v in &rv {
            assert!(v >= 0.0);
        }
        for &i in &ri {
            assert!(i >= 0.0);
        }
        // Duplicate voltage=1.0 should be merged (average of 2.0 and 1.5 = 1.75)
        let idx = rv.iter().position(|&x| (x - 1.0).abs() < 1e-10).unwrap();
        assert!((ri[idx] - 1.75).abs() < 1e-10);
    }

    #[test]
    fn test_rectify_iv_curve_empty() {
        let (v, i) = rectify_iv_curve(&[], &[]);
        assert!(v.is_empty());
        assert!(i.is_empty());
    }

    #[test]
    fn test_fit_sandia_simple_roundtrip() {
        // Known parameters
        let il = 9.0;
        let i0 = 1e-10;
        let rs = 0.3;
        let rsh = 500.0;
        let a = 1.6; // nNsVth

        let (voltage, current) = generate_iv_curve(il, i0, rs, rsh, a, 200);
        let v_oc = *voltage.last().unwrap();
        let i_sc = current[0];

        let result = fit_sandia_simple(&voltage, &current, Some(v_oc), Some(i_sc), None, 0.2, 0.1);
        assert!(result.is_ok(), "fit_sandia_simple failed: {:?}", result.err());
        let r = result.unwrap();

        // Check recovered parameters are within reasonable tolerance
        let il_err = (r.photocurrent - il).abs() / il;
        let rsh_err = (r.resistance_shunt - rsh).abs() / rsh;
        let a_err = (r.n_ns_vth - a).abs() / a;

        assert!(il_err < 0.05, "IL error too large: {:.4}", il_err);
        assert!(rsh_err < 0.5, "Rsh error too large: {:.4}", rsh_err);
        assert!(a_err < 0.15, "nNsVth error too large: {:.4}", a_err);
        assert!(r.saturation_current > 0.0, "I0 should be positive");
        assert!(r.resistance_series >= 0.0, "Rs should be non-negative");
    }

    #[test]
    fn test_fit_desoto_typical_module() {
        // Typical 60-cell silicon module specs (similar to CS5P-220M)
        let v_mp = 29.0;
        let i_mp = 7.6;
        let v_oc = 36.3;
        let i_sc = 8.1;
        let alpha_sc = 0.003; // A/K
        let beta_voc = -0.125; // V/K
        let cells_in_series = 60;
        let eg_ref = 1.121;
        let d_eg_dt = -0.0002677;

        let result = fit_desoto(
            v_mp,
            i_mp,
            v_oc,
            i_sc,
            alpha_sc,
            beta_voc,
            cells_in_series,
            eg_ref,
            d_eg_dt,
        );
        assert!(result.is_ok(), "fit_desoto failed: {:?}", result.err());
        let r = result.unwrap();

        // Sanity checks on fitted parameters
        assert!(r.photocurrent > 0.0, "IL should be positive: {}", r.photocurrent);
        assert!(
            (r.photocurrent - i_sc).abs() < 1.0,
            "IL should be close to Isc: {}",
            r.photocurrent
        );
        assert!(r.saturation_current > 0.0, "I0 should be positive: {}", r.saturation_current);
        assert!(
            r.saturation_current < 1e-5,
            "I0 should be very small: {}",
            r.saturation_current
        );
        assert!(r.resistance_series > 0.0, "Rs should be positive: {}", r.resistance_series);
        assert!(r.resistance_series < 5.0, "Rs should be reasonable: {}", r.resistance_series);
        assert!(r.resistance_shunt > 10.0, "Rsh should be large: {}", r.resistance_shunt);
        assert!(r.n_ns_vth > 0.0, "a_ref should be positive: {}", r.n_ns_vth);
    }

    #[test]
    fn test_fit_sandia_simple_too_few_points() {
        let v = vec![0.0, 1.0, 2.0];
        let i = vec![5.0, 4.0, 0.0];
        let result = fit_sandia_simple(&v, &i, None, None, None, 0.2, 0.1);
        assert!(result.is_err());
    }
}
