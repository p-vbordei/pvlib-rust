use pvlib::spectrum::*;

// --- spectral_factor_sapm ---

#[test]
fn test_spectral_factor_sapm_at_reference() {
    // Coefficients that give 1.0 at AM=1.5 (typical crystalline Si)
    // Using simple coefficients: [1, 0, 0, 0, 0] -> always 1.0
    let f = spectral_factor_sapm(1.5, [1.0, 0.0, 0.0, 0.0, 0.0]);
    assert!((f - 1.0).abs() < 1e-10, "Should be 1.0, got {}", f);
}

#[test]
fn test_spectral_factor_sapm_nan_input() {
    let f = spectral_factor_sapm(f64::NAN, [1.0, 0.0, 0.0, 0.0, 0.0]);
    assert!((f - 0.0).abs() < 1e-10, "NaN input should return 0");
}

#[test]
fn test_spectral_factor_sapm_non_negative() {
    // Even with coefficients that could go negative, result is clamped to 0
    let f = spectral_factor_sapm(5.0, [-10.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(f >= 0.0, "Should be non-negative, got {}", f);
}

// --- spectral_mismatch_modifier ---

#[test]
fn test_spectral_mismatch_modifier_constant() {
    let m = spectral_mismatch_modifier(1.5, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!((m - 1.0).abs() < 1e-10);
}

#[test]
fn test_spectral_mismatch_modifier_clamps_airmass() {
    // AM < 1 should be clamped to 1
    let m1 = spectral_mismatch_modifier(0.5, [0.5, 0.5, 0.0, 0.0, 0.0, 0.0]);
    let m2 = spectral_mismatch_modifier(1.0, [0.5, 0.5, 0.0, 0.0, 0.0, 0.0]);
    assert!((m1 - m2).abs() < 1e-10, "AM < 1 should clamp to AM=1");
}

// --- first_solar_spectral_correction ---

#[test]
fn test_first_solar_constant() {
    let m = first_solar_spectral_correction(1.0, 1.5, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!((m - 1.0).abs() < 1e-10);
}

#[test]
fn test_first_solar_non_negative() {
    let m = first_solar_spectral_correction(1.0, 1.5, [-100.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(m >= 0.0, "Should clamp to 0");
}

// --- spectral_factor_caballero ---

#[test]
fn test_caballero_with_builtin_coeffs() {
    let coeffs = caballero_coefficients("monosi").unwrap();
    let f = spectral_factor_caballero(1.4164, 1.5, 0.084, coeffs);
    // At reference conditions (pw_ref, aod_ref), f_aod and f_pw should be ~0
    // so result should be ~f_am at AM=1.5
    assert!(f > 0.8 && f < 1.2, "Should be near 1.0 at reference conditions, got {}", f);
}

#[test]
fn test_caballero_coefficients_known_types() {
    assert!(caballero_coefficients("cdte").is_some());
    assert!(caballero_coefficients("monosi").is_some());
    assert!(caballero_coefficients("multisi").is_some());
    assert!(caballero_coefficients("cigs").is_some());
    assert!(caballero_coefficients("asi").is_some());
    assert!(caballero_coefficients("perovskite").is_some());
}

#[test]
fn test_caballero_coefficients_unknown() {
    assert!(caballero_coefficients("unknown_type").is_none());
}
