/// Calculate PV module efficiency using the ADR model.
///
/// The efficiency varies with irradiance and operating temperature
/// and is determined by 5 model parameters (Driesse et al., 2020).
///
/// # Arguments
/// * `effective_irradiance` - Effective irradiance on PV module in W/m^2.
/// * `temp_cell` - PV module operating temperature in degrees C.
/// * `k_a` - Absolute scaling factor (efficiency at reference conditions).
/// * `k_d` - Dark irradiance / diode coefficient (negative).
/// * `tc_d` - Temperature coefficient of diode coefficient.
/// * `k_rs` - Series resistance loss coefficient.
/// * `k_rsh` - Shunt resistance loss coefficient.
///
/// # Returns
/// Module efficiency (same units/scale as k_a).
pub fn pvefficiency_adr(
    effective_irradiance: f64,
    temp_cell: f64,
    k_a: f64,
    k_d: f64,
    tc_d: f64,
    k_rs: f64,
    k_rsh: f64,
) -> f64 {
    let g_ref = 1000.0;
    let t_ref = 25.0;

    let s = effective_irradiance / g_ref;
    let dt = temp_cell - t_ref;

    // Eq 29: s_o and s_o_ref using 10^x
    let s_o = 10.0_f64.powf(k_d + dt * tc_d);
    let s_o_ref = 10.0_f64.powf(k_d);

    // Eq 28, 30: normalized voltage
    let v = (s / s_o + 1.0).ln() / (1.0 / s_o_ref + 1.0).ln();

    // Eq 25: efficiency
    k_a * ((1.0 + k_rs + k_rsh) * v - k_rs * s - k_rsh * v * v)
}

/// DC power using the Huld model (used by PVGIS).
///
/// P_dc = G' * (pdc0 + k1*ln(G') + k2*ln(G')^2 + k3*T' + k4*T'*ln(G') + k5*T'*ln(G')^2 + k6*T'^2)
///
/// where G' = irradiance/1000, T' = temp_module - 25.
///
/// # Arguments
/// * `effective_irradiance` - Irradiance converted to photocurrent in W/m^2.
/// * `temp_module` - Module back-surface temperature in degrees C.
/// * `pdc0` - Power at reference conditions (1000 W/m^2, 25C) in W.
/// * `k` - Six empirical coefficients [k1, k2, k3, k4, k5, k6].
///
/// # Returns
/// DC power in W.
pub fn huld(
    effective_irradiance: f64,
    temp_module: f64,
    pdc0: f64,
    k: [f64; 6],
) -> f64 {
    let gprime = effective_irradiance / 1000.0;
    let tprime = temp_module - 25.0;

    let log_g = if gprime > 0.0 { gprime.ln() } else { 0.0 };

    gprime
        * (pdc0
            + k[0] * log_g
            + k[1] * log_g * log_g
            + k[2] * tprime
            + k[3] * tprime * log_g
            + k[4] * tprime * log_g * log_g
            + k[5] * tprime * tprime)
}
