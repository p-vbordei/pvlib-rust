use pvlib::shading::*;

// --- ground_angle ---

#[test]
fn test_ground_angle_zero_height() {
    // At slant_height=0, ground angle should be 0 (looking from bottom of row)
    let psi = ground_angle(30.0, 0.5, 0.0);
    assert!((psi - 0.0).abs() < 1e-6, "Ground angle at slant_height=0 should be 0, got {}", psi);
}

#[test]
fn test_ground_angle_increases_with_height() {
    let psi_low = ground_angle(30.0, 0.5, 0.2);
    let psi_high = ground_angle(30.0, 0.5, 0.8);
    assert!(psi_high > psi_low, "Ground angle should increase with slant height");
}

#[test]
fn test_ground_angle_increases_with_gcr() {
    let psi_low = ground_angle(30.0, 0.3, 0.5);
    let psi_high = ground_angle(30.0, 0.7, 0.5);
    assert!(psi_high > psi_low, "Ground angle should increase with GCR");
}

// --- masking_angle ---

#[test]
fn test_masking_angle_bottom() {
    // At slant_height=0, masking angle should be largest (worst case)
    let ma_bottom = masking_angle(30.0, 0.5, 0.0);
    let ma_top = masking_angle(30.0, 0.5, 0.9);
    assert!(ma_bottom > ma_top, "Masking angle should be larger at bottom of row");
}

#[test]
fn test_masking_angle_top() {
    // At slant_height=1, masking angle should be 0 (no shading from top)
    let ma = masking_angle(30.0, 0.5, 1.0);
    assert!((ma - 0.0).abs() < 1e-6, "Masking angle at top of row should be 0, got {}", ma);
}

#[test]
fn test_masking_angle_zero_gcr() {
    let ma = masking_angle(30.0, 0.0, 0.0);
    assert_eq!(ma, 0.0, "Masking angle with zero GCR should be 0");
}

// --- masking_angle_passias ---

#[test]
fn test_masking_angle_passias_zero_tilt() {
    let ma = masking_angle_passias(0.0, 0.5);
    assert!((ma - 0.0).abs() < 0.01, "Passias masking angle at tilt=0 should be ~0, got {}", ma);
}

#[test]
fn test_masking_angle_passias_positive() {
    let ma = masking_angle_passias(30.0, 0.5);
    assert!(ma > 0.0, "Passias masking angle should be positive for tilted rows, got {}", ma);
}

#[test]
fn test_masking_angle_passias_increases_with_gcr() {
    let ma_low = masking_angle_passias(30.0, 0.3);
    let ma_high = masking_angle_passias(30.0, 0.7);
    assert!(ma_high > ma_low, "Passias masking angle should increase with GCR");
}

// --- projected_solar_zenith_angle ---

#[test]
fn test_psza_sun_overhead() {
    // Sun directly overhead (zenith=0), PSZA should be 0
    let psza = projected_solar_zenith_angle(0.0, 180.0, 0.0, 180.0);
    assert!((psza - 0.0).abs() < 1e-6, "PSZA should be 0 when sun is overhead, got {}", psza);
}

#[test]
fn test_psza_sun_along_axis() {
    // Sun along the axis direction: projection onto perpendicular plane should be small
    // Sun azimuth = axis azimuth = 180 (south), zenith = 30
    let psza = projected_solar_zenith_angle(30.0, 180.0, 0.0, 180.0);
    // When sun is along axis, sx_prime ~ 0, so PSZA ~ 0
    assert!(psza.abs() < 1.0, "PSZA should be near 0 when sun is along axis, got {}", psza);
}

#[test]
fn test_psza_sun_perpendicular_to_axis() {
    // Sun perpendicular to axis: N-S axis (azimuth=180), sun from east (azimuth=90)
    let psza = projected_solar_zenith_angle(30.0, 90.0, 0.0, 180.0);
    assert!(psza.abs() > 10.0, "PSZA should be significant when sun is perpendicular to axis, got {}", psza);
}

// --- shaded_fraction1d ---

#[test]
fn test_shaded_fraction_no_shade() {
    // Sun high overhead, wide pitch -> no shading
    let sf = shaded_fraction1d(10.0, 180.0, 90.0, 30.0, 2.0, 10.0, 0.0, 0.0, 0.0);
    assert!((sf - 0.0).abs() < 0.1, "Should have minimal shading with wide pitch, got {}", sf);
}

#[test]
fn test_shaded_fraction_clamped() {
    // Result should always be between 0 and 1
    let sf = shaded_fraction1d(80.0, 135.0, 90.0, 30.0, 2.0, 3.0, 0.0, 0.0, 0.0);
    assert!(sf >= 0.0 && sf <= 1.0, "Shaded fraction should be in [0,1], got {}", sf);
}

#[test]
fn test_shaded_fraction_reference_value() {
    // Compare with Python reference: fixed-tilt south-facing, morning
    // Python: shaded_fraction1d(solar_zenith=80, solar_azimuth=135,
    //     axis_azimuth=90, shaded_row_rotation=30, shading_row_rotation=30,
    //     collector_width=2, pitch=3, axis_tilt=0,
    //     surface_to_axis_offset=0.05, cross_axis_slope=0)
    // = 0.47755694708090535
    let sf = shaded_fraction1d(80.0, 135.0, 90.0, 30.0, 2.0, 3.0, 0.0, 0.05, 0.0);
    assert!((sf - 0.4776).abs() < 0.01,
        "Should match Python reference ~0.4776, got {}", sf);
}

// --- sky_diffuse_passias ---

#[test]
fn test_sky_diffuse_passias_zero_angle() {
    let loss = sky_diffuse_passias(0.0);
    assert!((loss - 0.0).abs() < 1e-6, "No masking angle should mean no loss, got {}", loss);
}

#[test]
fn test_sky_diffuse_passias_increases_with_angle() {
    let loss_10 = sky_diffuse_passias(10.0);
    let loss_30 = sky_diffuse_passias(30.0);
    assert!(loss_30 > loss_10, "Larger masking angle should mean more loss");
}

// --- sky_diffuse_pass_equation ---

#[test]
fn test_sky_diffuse_pass_zero_angle() {
    let pass = sky_diffuse_pass_equation(0.0);
    assert!((pass - 1.0).abs() < 1e-6, "No masking should pass all sky diffuse");
}
