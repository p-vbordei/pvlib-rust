/// NREL's PVWatts inverter model.
///
/// Calculates inverter efficiency as a function of input DC power using the
/// PVWatts efficiency curve.
///
/// # Arguments
/// * `pdc` - DC power input to the inverter (W).
/// * `pdc0` - DC input limit of the inverter (W).
/// * `eta_inv_nom` - Nominal inverter efficiency (default 0.96).
/// * `eta_inv_ref` - Reference inverter efficiency (default 0.9637).
///
/// # Returns
/// AC power output (W).
pub fn pvwatts_ac(pdc: f64, pdc0: f64, eta_inv_nom: f64, eta_inv_ref: f64) -> f64 {
    if pdc <= 0.0 {
        return 0.0;
    }

    let pac0 = eta_inv_nom * pdc0;
    let zeta = pdc / pdc0;
    let eta = (eta_inv_nom / eta_inv_ref) * (-0.0162 * zeta - 0.0059 / zeta + 0.9858);
    let power_ac = eta * pdc;

    power_ac.clamp(0.0, pac0)
}

/// Extend PVWatts inverter model for multiple MPPT inputs.
///
/// Sums DC power from all MPPT inputs, then applies the PVWatts model.
///
/// # Arguments
/// * `pdc` - DC power on each MPPT input (W).
/// * `pdc0` - Total DC power limit of the inverter (W).
/// * `eta_inv_nom` - Nominal inverter efficiency.
/// * `eta_inv_ref` - Reference inverter efficiency.
///
/// # Returns
/// AC power output (W).
pub fn pvwatts_multi(pdc: &[f64], pdc0: f64, eta_inv_nom: f64, eta_inv_ref: f64) -> f64 {
    let total_pdc: f64 = pdc.iter().sum();
    pvwatts_ac(total_pdc, pdc0, eta_inv_nom, eta_inv_ref)
}

/// Internal Sandia efficiency calculation (without clipping/limits).
///
/// Matches pvlib-python `_sandia_eff`.
fn sandia_eff(
    v_dc: f64, p_dc: f64,
    p_aco: f64, p_dco: f64, v_dco: f64, p_so: f64,
    c0: f64, c1: f64, c2: f64, c3: f64,
) -> f64 {
    let a = p_dco * (1.0 + c1 * (v_dc - v_dco));
    let b = p_so * (1.0 + c2 * (v_dc - v_dco));
    let c = c0 * (1.0 + c3 * (v_dc - v_dco));

    (p_aco / (a - b) - c * (a - b)) * (p_dc - b) + c * (p_dc - b).powi(2)
}

/// Sandia National Laboratories Grid-Connected Inverter Model.
///
/// # References
/// King, D. et al, 2007, "Performance Model for Grid-Connected Photovoltaic Inverters".
///
/// # Arguments
/// * `v_dc` - DC voltage input to the inverter (V).
/// * `p_dc` - DC power input to the inverter (W).
/// * `p_aco` - Maximum AC power rating of inverter (W).
/// * `p_dco` - DC power level at which AC power rating is achieved (W).
/// * `v_dco` - DC voltage at which AC power rating is achieved (V).
/// * `p_so` - DC power required to start the inversion process (W).
/// * `c0` - Parameter defining the curvature of the relationship (1/W).
/// * `c1` - Empirical coefficient (1/V).
/// * `c2` - Empirical coefficient (1/V).
/// * `c3` - Empirical coefficient (1/V).
/// * `p_nt` - AC power consumed by inverter at night (W).
///
/// # Returns
/// AC power output in Watts.
pub fn sandia(
    v_dc: f64, p_dc: f64,
    p_aco: f64, p_dco: f64, v_dco: f64, p_so: f64,
    c0: f64, c1: f64, c2: f64, c3: f64, p_nt: f64,
) -> f64 {
    let power_ac = sandia_eff(v_dc, p_dc, p_aco, p_dco, v_dco, p_so, c0, c1, c2, c3);
    sandia_limits(power_ac, p_dc, p_aco, p_nt, p_so)
}

/// Apply minimum and maximum power limits to Sandia AC power output.
fn sandia_limits(power_ac: f64, p_dc: f64, p_aco: f64, p_nt: f64, p_so: f64) -> f64 {
    if p_dc < p_so {
        return -p_nt.abs();
    }
    power_ac.min(p_aco)
}

/// Sandia inverter model for multiple MPPT inputs.
///
/// Extension of the Sandia model to inverters with multiple, unbalanced
/// MPPT inputs as described in Hansen et al., 2022.
///
/// # Arguments
/// * `v_dc` - DC voltage on each MPPT input (V).
/// * `p_dc` - DC power on each MPPT input (W).
/// * `p_aco` - Maximum AC power rating of inverter (W).
/// * `p_dco` - DC power level at which AC power rating is achieved (W).
/// * `v_dco` - DC voltage at which AC power rating is achieved (V).
/// * `p_so` - DC power required to start the inversion process (W).
/// * `c0` - Parameter defining the curvature (1/W).
/// * `c1` - Empirical coefficient (1/V).
/// * `c2` - Empirical coefficient (1/V).
/// * `c3` - Empirical coefficient (1/V).
/// * `p_nt` - AC power consumed by inverter at night (W).
///
/// # Returns
/// AC power output in Watts.
///
/// # Panics
/// Panics if `v_dc` and `p_dc` have different lengths.
pub fn sandia_multi(
    v_dc: &[f64], p_dc: &[f64],
    p_aco: f64, p_dco: f64, v_dco: f64, p_so: f64,
    c0: f64, c1: f64, c2: f64, c3: f64, p_nt: f64,
) -> f64 {
    assert_eq!(v_dc.len(), p_dc.len(), "v_dc and p_dc must have the same length");

    let power_dc_total: f64 = p_dc.iter().sum();

    if power_dc_total == 0.0 {
        return sandia_limits(0.0, power_dc_total, p_aco, p_nt, p_so);
    }

    let mut power_ac = 0.0;
    for (vdc_i, pdc_i) in v_dc.iter().zip(p_dc.iter()) {
        let weight = pdc_i / power_dc_total;
        power_ac += weight * sandia_eff(*vdc_i, power_dc_total, p_aco, p_dco, v_dco, p_so, c0, c1, c2, c3);
    }

    sandia_limits(power_ac, power_dc_total, p_aco, p_nt, p_so)
}

/// Anton Driesse (ADR) Inverter Model.
///
/// A highly robust parameter model used often in place of Sandia.
///
/// # References
/// Driesse, A. et al., 2020. "Beyond PVWatts: A streamlined PV grid-connected inverter model."
pub fn adr(v_dc: f64, p_dc: f64, p_nom: f64, v_nom: f64, adr_coeffs: [f64; 6]) -> f64 {
    let v_ratio = v_dc / v_nom.max(0.1);
    let p_ratio = p_dc / p_nom.max(0.1);

    // Empirical ADR efficiency polynomial
    let eff = adr_coeffs[0]
            + adr_coeffs[1] * v_ratio
            + adr_coeffs[2] * p_ratio
            + adr_coeffs[3] * v_ratio.powi(2)
            + adr_coeffs[4] * p_ratio.powi(2)
            + adr_coeffs[5] * v_ratio * p_ratio;

    let p_ac = p_dc * eff.clamp(0.0, 1.0);
    p_ac.clamp(0.0, p_nom)
}
