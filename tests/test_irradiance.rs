use pvlib::irradiance::*;

// --- King diffuse model tests ---

#[test]
fn test_king_horizontal_surface() {
    // For a horizontal surface (tilt=0), cos(0)=1, so (1-cos)/2=0
    // and the GHI term vanishes. Result = dhi * (1+1)/2 = dhi.
    let result = king(0.0, 100.0, 500.0, 30.0);
    assert!((result - 100.0).abs() < 0.01, "Horizontal tilt should return dhi, got {}", result);
}

#[test]
fn test_king_tilted_surface() {
    // For a 30 degree tilt, zenith=45:
    // dhi*(1+cos(30))/2 + ghi*(0.012*45-0.04)*(1-cos(30))/2
    let result = king(30.0, 100.0, 500.0, 45.0);
    assert!(result > 0.0, "King should return positive diffuse for valid inputs, got {}", result);
    // dhi term: 100 * (1 + 0.8660) / 2 = 93.30
    // ghi term: 500 * (0.54 - 0.04) * (1 - 0.8660) / 2 = 500 * 0.50 * 0.0670 = 16.75
    // total ~ 110.05
    assert!((result - 110.05).abs() < 1.0, "Expected ~110, got {}", result);
}

#[test]
fn test_king_high_zenith() {
    // At high zenith, the ghi correction term is larger
    let result_low = king(30.0, 100.0, 500.0, 20.0);
    let result_high = king(30.0, 100.0, 500.0, 70.0);
    assert!(result_high > result_low,
        "Higher zenith should increase king diffuse: low={}, high={}", result_low, result_high);
}

#[test]
fn test_king_negative_clamp() {
    let result = king(10.0, 0.0, 0.0, 0.0);
    assert!(result >= 0.0, "King should never return negative, got {}", result);
}

#[test]
fn test_king_zero_irradiance() {
    let result = king(30.0, 0.0, 0.0, 45.0);
    assert!((result - 0.0).abs() < 1e-10, "Zero irradiance should give zero diffuse, got {}", result);
}

// --- DIRINDEX model tests ---

#[test]
fn test_dirindex_clear_sky_conditions() {
    // Under clear sky conditions (ghi == ghi_clearsky), DIRINDEX should
    // return approximately dni_clearsky since the ratio is ~1.
    let ghi = 500.0;
    let ghi_cs = 500.0;
    let dni_cs = 700.0;
    let zenith = 30.0;
    let doy = 172;

    let result = dirindex(ghi, ghi_cs, dni_cs, zenith, doy, Some(101325.0));
    assert!((result - dni_cs).abs() < 1.0,
        "Under clear sky, dirindex should return ~dni_clearsky ({}), got {}", dni_cs, result);
}

#[test]
fn test_dirindex_cloudy_conditions() {
    // When measured GHI is less than clearsky GHI, DNI should be reduced
    let ghi = 300.0;
    let ghi_cs = 500.0;
    let dni_cs = 700.0;
    let zenith = 30.0;
    let doy = 172;

    let result = dirindex(ghi, ghi_cs, dni_cs, zenith, doy, Some(101325.0));
    assert!(result < dni_cs,
        "Cloudy conditions should reduce DNI below clearsky: got {} vs {}", result, dni_cs);
    assert!(result >= 0.0, "DIRINDEX should not return negative values, got {}", result);
}

#[test]
fn test_dirindex_zero_ghi() {
    let result = dirindex(0.0, 500.0, 700.0, 30.0, 172, Some(101325.0));
    assert!((result - 0.0).abs() < 1e-6, "Zero GHI should give zero DNI, got {}", result);
}

#[test]
fn test_dirindex_high_zenith() {
    let result = dirindex(100.0, 200.0, 300.0, 88.0, 172, Some(101325.0));
    assert!(result >= 0.0, "High zenith should return non-negative, got {}", result);
}

#[test]
fn test_dirindex_none_pressure() {
    // Using None for pressure should default to 101325 Pa
    let result = dirindex(500.0, 500.0, 700.0, 30.0, 172, None);
    assert!(result > 0.0, "With default pressure, should return positive DNI, got {}", result);
}
