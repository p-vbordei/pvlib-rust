use pvlib::inverter::*;
use pvlib::tracking::calc_surface_orientation;

// --- pvwatts_ac ---

#[test]
fn test_pvwatts_ac_zero_input() {
    let pac = pvwatts_ac(0.0, 5000.0, 0.96, 0.9637);
    assert_eq!(pac, 0.0);
}

#[test]
fn test_pvwatts_ac_negative_input() {
    let pac = pvwatts_ac(-100.0, 5000.0, 0.96, 0.9637);
    assert_eq!(pac, 0.0);
}

#[test]
fn test_pvwatts_ac_clipping() {
    // pac0 = 0.96 * 5000 = 4800
    let pac = pvwatts_ac(5000.0, 5000.0, 0.96, 0.9637);
    assert!(pac <= 4800.0, "AC power should be clipped at pac0=4800, got {}", pac);
    assert!(pac > 0.0);
}

#[test]
fn test_pvwatts_ac_mid_power() {
    // At pdc = pdc0/2, zeta = 0.5
    // eta = (0.96/0.9637) * (-0.0162*0.5 - 0.0059/0.5 + 0.9858)
    //     = 0.99627 * (-0.0081 - 0.0118 + 0.9858)
    //     = 0.99627 * 0.9659
    //     = 0.96229
    let pdc = 2500.0;
    let pac = pvwatts_ac(pdc, 5000.0, 0.96, 0.9637);
    let expected = 0.96229 * 2500.0; // ~2405.7
    assert!((pac - expected).abs() < 5.0,
        "PVWatts at 50% load should be ~{}, got {}", expected, pac);
}

// --- pvwatts_multi ---

#[test]
fn test_pvwatts_multi_single_input() {
    let pdc = [2500.0];
    let pac_multi = pvwatts_multi(&pdc, 5000.0, 0.96, 0.9637);
    let pac_single = pvwatts_ac(2500.0, 5000.0, 0.96, 0.9637);
    assert!((pac_multi - pac_single).abs() < 1e-10,
        "pvwatts_multi with one input should match pvwatts_ac");
}

#[test]
fn test_pvwatts_multi_two_inputs() {
    let pdc = [1500.0, 1000.0];
    let pac_multi = pvwatts_multi(&pdc, 5000.0, 0.96, 0.9637);
    let pac_single = pvwatts_ac(2500.0, 5000.0, 0.96, 0.9637);
    assert!((pac_multi - pac_single).abs() < 1e-10,
        "pvwatts_multi should sum DC then apply model");
}

#[test]
fn test_pvwatts_multi_zero_inputs() {
    let pdc = [0.0, 0.0];
    let pac = pvwatts_multi(&pdc, 5000.0, 0.96, 0.9637);
    assert_eq!(pac, 0.0);
}

// --- sandia ---

#[test]
fn test_sandia_below_startup() {
    // p_dc < p_so should return -|p_nt|
    let pac = sandia(240.0, 5.0, 4000.0, 4200.0, 240.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.5);
    assert!((pac - (-0.5)).abs() < 1e-10,
        "Below startup power should return -p_nt, got {}", pac);
}

#[test]
fn test_sandia_clipping() {
    // With high DC power, AC should be capped at p_aco
    let pac = sandia(240.0, 5000.0, 4000.0, 4200.0, 240.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.5);
    assert!(pac <= 4000.0, "AC power should be clipped at p_aco=4000, got {}", pac);
}

#[test]
fn test_sandia_at_rated() {
    // At rated conditions (p_dc = p_dco, v_dc = v_dco), with c0=0, output should be near p_aco
    let pac = sandia(240.0, 4200.0, 4000.0, 4200.0, 240.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.5);
    // With c0=0: pac = (Paco/(A-B)) * (Pdc-B) = (4000/(4200-10)) * (4200-10) = 4000
    assert!((pac - 4000.0).abs() < 1.0,
        "At rated conditions with c0=0, AC power should be ~p_aco, got {}", pac);
}

// --- sandia_multi ---

#[test]
fn test_sandia_multi_single_input() {
    let v_dc = [240.0];
    let p_dc = [2000.0];
    let pac_multi = sandia_multi(&v_dc, &p_dc, 4000.0, 4200.0, 240.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.5);
    let pac_single = sandia(240.0, 2000.0, 4000.0, 4200.0, 240.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.5);
    assert!((pac_multi - pac_single).abs() < 1e-6,
        "sandia_multi with one input should match sandia: {} vs {}", pac_multi, pac_single);
}

#[test]
fn test_sandia_multi_balanced_inputs() {
    // Two balanced inputs: each contributes 50% weight
    let v_dc = [240.0, 240.0];
    let p_dc = [1000.0, 1000.0];
    let pac_multi = sandia_multi(&v_dc, &p_dc, 4000.0, 4200.0, 240.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.5);
    let pac_single = sandia(240.0, 2000.0, 4000.0, 4200.0, 240.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.5);
    assert!((pac_multi - pac_single).abs() < 1e-6,
        "Balanced sandia_multi should match single sandia: {} vs {}", pac_multi, pac_single);
}

#[test]
fn test_sandia_multi_below_startup() {
    let v_dc = [240.0, 240.0];
    let p_dc = [3.0, 2.0]; // total = 5.0 < p_so = 10.0
    let pac = sandia_multi(&v_dc, &p_dc, 4000.0, 4200.0, 240.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.5);
    assert!((pac - (-0.5)).abs() < 1e-10,
        "Below startup should return -p_nt, got {}", pac);
}

#[test]
#[should_panic(expected = "v_dc and p_dc must have the same length")]
fn test_sandia_multi_mismatched_lengths() {
    let v_dc = [240.0, 240.0];
    let p_dc = [1000.0];
    sandia_multi(&v_dc, &p_dc, 4000.0, 4200.0, 240.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.5);
}

// --- calc_surface_orientation ---

#[test]
fn test_calc_surface_orientation_horizontal() {
    // tracker_theta=0, axis_tilt=0 => surface_tilt=0
    let (tilt, _azimuth) = calc_surface_orientation(0.0, 0.0, 180.0);
    assert!((tilt - 0.0).abs() < 1e-6,
        "Horizontal tracker should have tilt=0, got {}", tilt);
}

#[test]
fn test_calc_surface_orientation_tilted_tracker() {
    // tracker_theta=30, axis_tilt=0, axis_azimuth=180
    // surface_tilt = acos(cos(30)*cos(0)) = acos(cos(30)) = 30
    let (tilt, azimuth) = calc_surface_orientation(30.0, 0.0, 180.0);
    assert!((tilt - 30.0).abs() < 1e-4,
        "Surface tilt should be 30 deg, got {}", tilt);
    // With positive tracker_theta and axis_azimuth=180:
    // azimuth_delta = asin(sin(30)/sin(30)) = asin(1) = 90
    // but actually: azimuth_delta = asin(sin(30)/sin(30)) = 90
    // surface_azimuth = (180 + 90) % 360 = 270
    // This means the panel faces west, which is correct for a south-facing
    // axis rotated westward.
    assert!((azimuth - 270.0).abs() < 1e-4,
        "Surface azimuth should be 270 (west), got {}", azimuth);
}

#[test]
fn test_calc_surface_orientation_negative_theta() {
    // tracker_theta=-30, axis_tilt=0, axis_azimuth=180
    let (tilt, azimuth) = calc_surface_orientation(-30.0, 0.0, 180.0);
    assert!((tilt - 30.0).abs() < 1e-4,
        "Surface tilt should be 30 deg for theta=-30, got {}", tilt);
    // azimuth_delta = asin(sin(-30)/sin(30)) = asin(-1) = -90
    // surface_azimuth = (180 - 90) % 360 = 90 (east)
    assert!((azimuth - 90.0).abs() < 1e-4,
        "Surface azimuth should be 90 (east), got {}", azimuth);
}

#[test]
fn test_calc_surface_orientation_with_axis_tilt() {
    // tracker_theta=0, axis_tilt=10, axis_azimuth=180
    // surface_tilt = acos(cos(0)*cos(10)) = acos(cos(10)) = 10
    let (tilt, azimuth) = calc_surface_orientation(0.0, 10.0, 180.0);
    assert!((tilt - 10.0).abs() < 1e-4,
        "Surface tilt should be 10 when axis_tilt=10 and theta=0, got {}", tilt);
    // When tracker_theta=0, azimuth_delta = asin(0/sin(10)) = 0
    // But the function handles sin_st != 0 case, so azimuth_delta = asin(0) = 0
    // surface_azimuth = (180 + 0) % 360 = 180
    assert!((azimuth - 180.0).abs() < 1e-4,
        "Surface azimuth should be 180 when theta=0 and axis_azimuth=180, got {}", azimuth);
}

#[test]
fn test_calc_surface_orientation_symmetry() {
    // Positive and negative theta should give same tilt
    let (tilt_pos, _) = calc_surface_orientation(45.0, 0.0, 180.0);
    let (tilt_neg, _) = calc_surface_orientation(-45.0, 0.0, 180.0);
    assert!((tilt_pos - tilt_neg).abs() < 1e-10,
        "Tilt should be symmetric: {} vs {}", tilt_pos, tilt_neg);
}
