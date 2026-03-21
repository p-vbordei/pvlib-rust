use pvlib::albedo::*;
use pvlib::pvarray::*;
use pvlib::transformer::*;
use pvlib::spectrum::*;

// --- albedo ---

#[test]
fn test_inland_water_dvoracek_high_sun() {
    // High sun elevation -> lower albedo (more absorption)
    let albedo_high = inland_water_dvoracek(60.0, "clear_water_no_waves");
    let albedo_low = inland_water_dvoracek(10.0, "clear_water_no_waves");
    assert!(albedo_high < albedo_low,
        "Higher sun should give lower albedo: high={}, low={}", albedo_high, albedo_low);
}

#[test]
fn test_inland_water_dvoracek_range() {
    let albedo = inland_water_dvoracek(45.0, "clear_water_no_waves");
    assert!(albedo > 0.0 && albedo < 1.0, "Albedo should be in (0,1), got {}", albedo);
}

#[test]
fn test_inland_water_dvoracek_negative_elevation() {
    // Negative elevation should be clamped to 0
    let albedo_neg = inland_water_dvoracek(-10.0, "clear_water_no_waves");
    let albedo_zero = inland_water_dvoracek(0.0, "clear_water_no_waves");
    assert!((albedo_neg - albedo_zero).abs() < 1e-10,
        "Negative elevation should equal zero elevation");
}

#[test]
fn test_inland_water_dvoracek_custom() {
    let albedo = inland_water_dvoracek_custom(45.0, 0.13, 0.29);
    let albedo_named = inland_water_dvoracek(45.0, "clear_water_no_waves");
    assert!((albedo - albedo_named).abs() < 1e-10,
        "Custom coefficients should match named condition");
}

#[test]
fn test_surface_albedo() {
    assert_eq!(surface_albedo("grass"), Some(0.20));
    assert_eq!(surface_albedo("snow"), Some(0.65));
    assert_eq!(surface_albedo("unknown"), None);
}

// --- pvarray ---

#[test]
fn test_pvefficiency_adr_at_stc() {
    // At STC (1000 W/m2, 25C), efficiency should equal k_a
    let eta = pvefficiency_adr(1000.0, 25.0, 100.0, -6.0, 0.02, 0.05, 0.10);
    assert!((eta - 100.0).abs() < 0.01, "ADR at STC should equal k_a=100, got {}", eta);
}

#[test]
fn test_pvefficiency_adr_low_irradiance() {
    // At lower irradiance, efficiency should decrease
    let eta_1000 = pvefficiency_adr(1000.0, 25.0, 100.0, -6.0, 0.02, 0.05, 0.10);
    let eta_200 = pvefficiency_adr(200.0, 25.0, 100.0, -6.0, 0.02, 0.05, 0.10);
    assert!(eta_200 < eta_1000, "Efficiency should decrease at lower irradiance");
}

#[test]
fn test_pvefficiency_adr_reference_value() {
    // Python: pvefficiency_adr(200, 25, k_a=100, k_d=-6.0, tc_d=0.02, k_rs=0.05, k_rsh=0.10)
    // = 92.79729308
    let eta = pvefficiency_adr(200.0, 25.0, 100.0, -6.0, 0.02, 0.05, 0.10);
    assert!((eta - 92.797).abs() < 0.01,
        "ADR at 200 W/m2 should be ~92.797, got {}", eta);
}

#[test]
fn test_huld_at_stc() {
    // At STC with zero k coefficients, power should be pdc0
    let pdc = huld(1000.0, 25.0, 250.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!((pdc - 250.0).abs() < 0.01, "Huld at STC with k=0 should give pdc0, got {}", pdc);
}

#[test]
fn test_huld_scales_with_irradiance() {
    let k = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let pdc_1000 = huld(1000.0, 25.0, 250.0, k);
    let pdc_500 = huld(500.0, 25.0, 250.0, k);
    assert!((pdc_500 / pdc_1000 - 0.5).abs() < 0.01,
        "With k=0, power should scale linearly with irradiance");
}

#[test]
fn test_huld_zero_irradiance() {
    let pdc = huld(0.0, 25.0, 250.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!((pdc - 0.0).abs() < 1e-10, "Zero irradiance should give zero power");
}

// --- transformer ---

#[test]
fn test_simple_efficiency_no_losses() {
    // With zero losses, output should equal input
    let output = simple_efficiency(1000.0, 5000.0, 0.0, 0.0);
    assert!((output - 1000.0).abs() < 0.01,
        "With no losses, output should equal input, got {}", output);
}

#[test]
fn test_simple_efficiency_losses() {
    // With losses, output should be less than input
    let output = simple_efficiency(1000.0, 5000.0, 0.005, 0.01);
    assert!(output < 1000.0, "Output should be less than input with losses");
    assert!(output > 0.0, "Output should still be positive");
}

#[test]
fn test_simple_efficiency_zero_input() {
    // At zero input, output should be negative (no-load losses)
    let output = simple_efficiency(0.0, 5000.0, 0.005, 0.01);
    assert!(output < 0.0, "Zero input with no-load losses should give negative output, got {}", output);
}

// --- spectrum ---

#[test]
fn test_spectral_factor_sapm_basic() {
    // Typical coefficients near 1 at AM=1.5
    let coeffs = [1.0, 0.0, 0.0, 0.0, 0.0]; // constant 1.0
    let f1 = spectral_factor_sapm(1.5, coeffs);
    assert!((f1 - 1.0).abs() < 1e-10, "SAPM with constant coeffs should be 1.0");
}

#[test]
fn test_spectral_factor_sapm_nan() {
    let f1 = spectral_factor_sapm(f64::NAN, [1.0, 0.0, 0.0, 0.0, 0.0]);
    assert_eq!(f1, 0.0, "NaN airmass should return 0");
}

#[test]
fn test_spectral_factor_sapm_nonnegative() {
    // Even with coefficients that could make it negative, result is clamped
    let f1 = spectral_factor_sapm(10.0, [-100.0, 0.0, 0.0, 0.0, 0.0]);
    assert_eq!(f1, 0.0, "SAPM should never return negative");
}

#[test]
fn test_spectral_factor_caballero_monosi() {
    let coeffs = caballero_coefficients("monosi").unwrap();
    let modifier = spectral_factor_caballero(1.4164, 1.5, 0.084, coeffs);
    // At reference conditions (pw=pw_ref, aod=aod_ref), modifier should be ~f_AM only
    assert!((modifier - 1.0).abs() < 0.05,
        "Caballero at reference conditions should be near 1.0, got {}", modifier);
}

#[test]
fn test_caballero_coefficients_known_types() {
    assert!(caballero_coefficients("cdte").is_some());
    assert!(caballero_coefficients("monosi").is_some());
    assert!(caballero_coefficients("multisi").is_some());
    assert!(caballero_coefficients("cigs").is_some());
    assert!(caballero_coefficients("asi").is_some());
    assert!(caballero_coefficients("perovskite").is_some());
    assert!(caballero_coefficients("unknown").is_none());
}
