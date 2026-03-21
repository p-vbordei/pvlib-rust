/// Spectral Mismatch Modifier (SMM)
/// 
/// Applies a polynomial modifier based on Air Mass to account for spectral variations
/// in module efficiency compared to STC (AM 1.5).
/// 
/// # Arguments
/// * `airmass_absolute` - The absolute air mass.
/// * `coefficients` - 6 empirical polynomial coefficients [C0, C1, C2, C3, C4, C5].
/// 
/// # Returns
/// A unitless modifier to multiply against DC power or current.
pub fn spectral_mismatch_modifier(airmass_absolute: f64, coefficients: [f64; 6]) -> f64 {
    let x = airmass_absolute.max(1.0);
    
    let ans = coefficients[0]
        + coefficients[1] * x
        + coefficients[2] * x.powi(2)
        + coefficients[3] * x.powi(3)
        + coefficients[4] * x.powi(4)
        + coefficients[5] * x.powi(5);
        
    ans.max(0.0)
}

/// First Solar empirical spectral mismatch modifier.
/// 
/// Evaluates the impact of precipitable water and airmass on CdTe and generic silicon spectrums.
/// 
/// # References
/// Lee, M. et al., 2018. "A new empirical model for spectral correction of PV performance."
pub fn first_solar_spectral_correction(precipitable_water: f64, airmass_absolute: f64, coefficients: [f64; 6]) -> f64 {
    let pw = precipitable_water;
    let am = airmass_absolute;
    let m = coefficients[0] 
          + coefficients[1]*am 
          + coefficients[2]*pw 
          + coefficients[3]*am.powi(2) 
          + coefficients[4]*pw.powi(2) 
          + coefficients[5]*am*pw;
          
    m.max(0.0)
}

/// SAPM spectral loss factor.
///
/// Fourth-order polynomial of absolute airmass: f1 = A0 + A1*AM + A2*AM^2 + A3*AM^3 + A4*AM^4.
///
/// # Arguments
/// * `airmass_absolute` - Absolute airmass. NaN values result in 0.
/// * `coefficients` - Array of 5 coefficients [A0, A1, A2, A3, A4].
///
/// # Returns
/// Spectral mismatch factor (unitless, non-negative).
pub fn spectral_factor_sapm(airmass_absolute: f64, coefficients: [f64; 5]) -> f64 {
    if airmass_absolute.is_nan() {
        return 0.0;
    }
    let am = airmass_absolute;
    let f1 = coefficients[0]
        + coefficients[1] * am
        + coefficients[2] * am.powi(2)
        + coefficients[3] * am.powi(3)
        + coefficients[4] * am.powi(4);
    f1.max(0.0)
}

/// Caballero (2018) spectral mismatch factor.
///
/// Estimates spectral modifier from airmass, AOD at 500nm, and precipitable water.
///
/// # Arguments
/// * `precipitable_water` - Atmospheric precipitable water in cm.
/// * `airmass_absolute` - Absolute airmass.
/// * `aod500` - Aerosol optical depth at 500 nm.
/// * `coefficients` - Array of 12 coefficients from Caballero (2018).
///
/// # Returns
/// Spectral mismatch modifier (unitless).
pub fn spectral_factor_caballero(
    precipitable_water: f64,
    airmass_absolute: f64,
    aod500: f64,
    coefficients: [f64; 12],
) -> f64 {
    let ama = airmass_absolute;
    let aod500_ref = 0.084;
    let pw_ref = 1.4164;

    let c = coefficients;

    // f_AM: 4th order polynomial in airmass
    let f_am = c[0] + c[1] * ama + c[2] * ama.powi(2) + c[3] * ama.powi(3) + c[4] * ama.powi(4);

    // f_AOD: Eq 6 with Table 1 selectors (c[10], c[11])
    let f_aod = (aod500 - aod500_ref)
        * (c[5] + c[10] * c[6] * ama + c[11] * c[6] * ama.ln() + c[7] * ama.powi(2));

    // f_PW: Eq 7 with Table 1
    let f_pw = (precipitable_water - pw_ref) * (c[8] + c[9] * ama.ln());

    f_am + f_aod + f_pw
}

/// Built-in Caballero coefficients for common module types.
pub fn caballero_coefficients(module_type: &str) -> Option<[f64; 12]> {
    match module_type {
        "cdte" => Some([
            1.0044, 0.0095, -0.0037, 0.0002, 0.0000, -0.0046, -0.0182, 0.0, 0.0095, 0.0068,
            0.0, 1.0,
        ]),
        "monosi" => Some([
            0.9706, 0.0377, -0.0123, 0.0025, -0.0002, 0.0159, -0.0165, 0.0, -0.0016, -0.0027,
            1.0, 0.0,
        ]),
        "multisi" => Some([
            0.9836, 0.0254, -0.0085, 0.0016, -0.0001, 0.0094, -0.0132, 0.0, -0.0002, -0.0011,
            1.0, 0.0,
        ]),
        "cigs" => Some([
            0.9801, 0.0283, -0.0092, 0.0019, -0.0001, 0.0117, -0.0126, 0.0, -0.0011, -0.0019,
            1.0, 0.0,
        ]),
        "asi" => Some([
            1.1060, -0.0848, 0.0302, -0.0076, 0.0006, -0.1283, 0.0986, -0.0254, 0.0156, 0.0146,
            1.0, 0.0,
        ]),
        "perovskite" => Some([
            1.0637, -0.0491, 0.0180, -0.0047, 0.0004, -0.0773, 0.0583, -0.0159, 0.01251,
            0.0109, 1.0, 0.0,
        ]),
        _ => None,
    }
}

