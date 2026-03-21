use pvlib::temperature::{
    sapm_module, sapm_module_default, sapm_cell_from_module,
    sapm_cell_from_module_default, noct_sam, noct_sam_default,
    generic_linear, generic_linear_default,
};

#[test]
fn test_sapm_module_known_values() {
    // poa=1000, temp_air=25, wind_speed=2, a=-3.56, b=-0.075
    // Expected: 1000*exp(-3.56 + -0.075*2) + 25 = 49.4775
    let tm = sapm_module(1000.0, 25.0, 2.0, -3.56, -0.075);
    assert!(
        (tm - 49.4775).abs() < 0.01,
        "SAPM module temp expected ~49.48, got {}",
        tm
    );
}

#[test]
fn test_sapm_module_default_matches() {
    let tm = sapm_module_default(1000.0, 25.0, 2.0);
    let tm_explicit = sapm_module(1000.0, 25.0, 2.0, -3.56, -0.0750);
    assert!(
        (tm - tm_explicit).abs() < 1e-10,
        "Default should match explicit parameters"
    );
}

#[test]
fn test_sapm_module_zero_irradiance() {
    // With zero irradiance, module temp should equal air temp
    let tm = sapm_module(0.0, 20.0, 5.0, -3.56, -0.075);
    assert!(
        (tm - 20.0).abs() < 0.01,
        "At zero irradiance, module temp should be ~air temp, got {}",
        tm
    );
}

#[test]
fn test_sapm_cell_from_module_known() {
    // module_temp=49.4775, poa=1000, delta_t=3, irrad_ref=1000
    // Expected: 49.4775 + 1000/1000 * 3 = 52.4775
    let tc = sapm_cell_from_module(49.4775, 1000.0, 3.0, 1000.0);
    assert!(
        (tc - 52.4775).abs() < 0.01,
        "SAPM cell from module expected ~52.48, got {}",
        tc
    );
}

#[test]
fn test_sapm_cell_from_module_default() {
    let tc = sapm_cell_from_module_default(49.4775, 1000.0);
    let tc_explicit = sapm_cell_from_module(49.4775, 1000.0, 3.0, 1000.0);
    assert!((tc - tc_explicit).abs() < 1e-10);
}

#[test]
fn test_noct_sam_known_values() {
    // poa=1000, temp_air=25, wind_speed=2, noct=45, eff=0.20,
    // tau_alpha=0.9, height=1, standoff=4 (>3.5 so adj=0)
    // Expected: ~49.11
    let tc = noct_sam(1000.0, 25.0, 2.0, 45.0, 0.20, 0.9, 1, 4.0);
    assert!(
        (tc - 49.1127).abs() < 0.01,
        "NOCT SAM expected ~49.11, got {}",
        tc
    );
}

#[test]
fn test_noct_sam_default() {
    let tc = noct_sam_default(1000.0, 25.0, 2.0, 45.0, 0.20);
    let tc_explicit = noct_sam(1000.0, 25.0, 2.0, 45.0, 0.20, 0.9, 1, 4.0);
    assert!((tc - tc_explicit).abs() < 1e-10);
}

#[test]
fn test_noct_sam_height_2() {
    // Height 2 uses wind_adj = 0.61 * ws instead of 0.51
    let tc_h1 = noct_sam(1000.0, 25.0, 5.0, 45.0, 0.15, 0.9, 1, 4.0);
    let tc_h2 = noct_sam(1000.0, 25.0, 5.0, 45.0, 0.15, 0.9, 2, 4.0);
    // Higher array -> more wind cooling -> lower temp
    assert!(
        tc_h2 < tc_h1,
        "Height 2 should be cooler than height 1: h1={}, h2={}",
        tc_h1,
        tc_h2
    );
}

#[test]
fn test_generic_linear_known_values() {
    // poa=1000, temp_air=25, ws=2, u_c=29, u_v=0, eff=0.20, abs=0.9
    // heat_input = 1000*(0.9-0.20) = 700, total_loss = 29
    // T = 25 + 700/29 = 49.1379
    let tc = generic_linear(1000.0, 25.0, 2.0, 29.0, 0.0, 0.20, 0.9);
    assert!(
        (tc - 49.1379).abs() < 0.01,
        "Generic linear expected ~49.14, got {}",
        tc
    );
}

#[test]
fn test_generic_linear_with_wind() {
    // With u_v > 0, wind should reduce temperature
    let tc_no_wind = generic_linear(1000.0, 25.0, 0.0, 29.0, 6.0, 0.15, 0.9);
    let tc_with_wind = generic_linear(1000.0, 25.0, 5.0, 29.0, 6.0, 0.15, 0.9);
    assert!(
        tc_with_wind < tc_no_wind,
        "Wind should reduce cell temp: no_wind={}, with_wind={}",
        tc_no_wind,
        tc_with_wind
    );
}

#[test]
fn test_generic_linear_default() {
    let tc = generic_linear_default(1000.0, 25.0, 2.0);
    let tc_explicit = generic_linear(1000.0, 25.0, 2.0, 29.0, 0.0, 0.15, 0.9);
    assert!((tc - tc_explicit).abs() < 1e-10);
}
