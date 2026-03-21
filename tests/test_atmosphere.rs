use pvlib::atmosphere::*;

#[test]
fn test_pres2alt_sea_level() {
    let alt = pres2alt(101325.0);
    assert!((alt - 0.0).abs() < 1.0, "Sea level pressure -> ~0 m altitude, got {}", alt);
}

#[test]
fn test_pres2alt_inverse_of_alt2pres() {
    let altitude = 1500.0;
    let pressure = alt2pres(altitude);
    let recovered = pres2alt(pressure);
    assert!((recovered - altitude).abs() < 1.0,
        "pres2alt should be inverse of alt2pres: expected {}, got {}", altitude, recovered);
}

#[test]
fn test_gueymard94_pw() {
    // At 20C, 50% RH, precipitable water should be a reasonable value
    let pw = gueymard94_pw(20.0, 50.0);
    assert!(pw > 0.1 && pw < 10.0, "PW should be reasonable, got {}", pw);

    // Higher temp and humidity should give more precipitable water
    let pw_high = gueymard94_pw(30.0, 80.0);
    assert!(pw_high > pw, "Higher T and RH should give more PW");
}

#[test]
fn test_gueymard94_pw_minimum() {
    // Very cold and dry should clamp to 0.1
    let pw = gueymard94_pw(-20.0, 5.0);
    assert!((pw - 0.1).abs() < 1e-6 || pw > 0.1, "PW minimum should be 0.1");
}

#[test]
fn test_tdew_from_rh() {
    // At 100% RH, dew point equals air temperature
    let tdew = tdew_from_rh(25.0, 100.0);
    assert!((tdew - 25.0).abs() < 0.01, "At 100% RH, Tdew should equal Tair, got {}", tdew);

    // At lower RH, dew point should be below air temp
    let tdew2 = tdew_from_rh(25.0, 50.0);
    assert!(tdew2 < 25.0, "Dew point should be below air temp at 50% RH");
}

#[test]
fn test_rh_from_tdew() {
    // When dew point equals air temp, RH should be 100%
    let rh = rh_from_tdew(25.0, 25.0);
    assert!((rh - 100.0).abs() < 0.01, "RH should be 100% when Tdew == Tair, got {}", rh);

    // When dew point is below air temp, RH should be < 100%
    let rh2 = rh_from_tdew(25.0, 15.0);
    assert!(rh2 < 100.0 && rh2 > 0.0, "RH should be between 0 and 100, got {}", rh2);
}

#[test]
fn test_tdew_rh_roundtrip() {
    let temp = 25.0;
    let rh_original = 60.0;
    let tdew = tdew_from_rh(temp, rh_original);
    let rh_recovered = rh_from_tdew(temp, tdew);
    assert!((rh_recovered - rh_original).abs() < 0.01,
        "Roundtrip RH should match: expected {}, got {}", rh_original, rh_recovered);
}

#[test]
fn test_bird_hulstrom80_aod_bb() {
    let aod_bb = bird_hulstrom80_aod_bb(0.1, 0.1);
    let expected = 0.27583 * 0.1 + 0.35 * 0.1;
    assert!((aod_bb - expected).abs() < 1e-10, "AOD_bb mismatch: expected {}, got {}", expected, aod_bb);
}

#[test]
fn test_kasten96_lt() {
    let lt = kasten96_lt(1.0, 1.0, 0.1);
    assert!(lt > 0.0, "Linke turbidity should be positive, got {}", lt);

    // Higher AOD should give higher turbidity
    let lt2 = kasten96_lt(1.0, 1.0, 0.5);
    assert!(lt2 > lt, "Higher AOD should give higher turbidity");
}

#[test]
fn test_angstrom_aod_at_lambda() {
    // AOD at same wavelength should be unchanged
    let aod = angstrom_aod_at_lambda(0.1, 500.0, 1.14, 500.0);
    assert!((aod - 0.1).abs() < 1e-10, "AOD at same wavelength should be unchanged");

    // AOD at longer wavelength should be smaller (for positive alpha)
    let aod700 = angstrom_aod_at_lambda(0.1, 500.0, 1.14, 700.0);
    assert!(aod700 < 0.1, "AOD at longer wavelength should be smaller");
}

#[test]
fn test_angstrom_alpha() {
    // Compute alpha from two known AOD values, then verify roundtrip
    let aod1 = 0.2;
    let lambda1 = 380.0;
    let aod2 = 0.1;
    let lambda2 = 500.0;
    let alpha = angstrom_alpha(aod1, lambda1, aod2, lambda2);
    assert!(alpha > 0.0, "Alpha should be positive for decreasing AOD with wavelength");

    // Use computed alpha to recover aod2 from aod1
    let aod2_recovered = angstrom_aod_at_lambda(aod1, lambda1, alpha, lambda2);
    assert!((aod2_recovered - aod2).abs() < 1e-10,
        "Roundtrip AOD should match: expected {}, got {}", aod2, aod2_recovered);
}

#[test]
fn test_windspeed_powerlaw() {
    // Same height should give same speed
    let ws = windspeed_powerlaw(10.0, 10.0, 10.0, 1.0 / 7.0);
    assert!((ws - 10.0).abs() < 1e-10, "Same height should give same speed");

    // Higher height should give higher speed
    let ws_high = windspeed_powerlaw(10.0, 10.0, 50.0, 1.0 / 7.0);
    assert!(ws_high > 10.0, "Wind speed should increase with height");

    // Lower height should give lower speed
    let ws_low = windspeed_powerlaw(10.0, 10.0, 2.0, 1.0 / 7.0);
    assert!(ws_low < 10.0, "Wind speed should decrease with lower height");
}

#[test]
fn test_windspeed_powerlaw_negative_inputs() {
    let ws = windspeed_powerlaw(-5.0, 10.0, 20.0, 1.0 / 7.0);
    assert!(ws.is_nan(), "Negative wind speed should return NaN");

    let ws2 = windspeed_powerlaw(10.0, -10.0, 20.0, 1.0 / 7.0);
    assert!(ws2.is_nan(), "Negative reference height should return NaN");
}
