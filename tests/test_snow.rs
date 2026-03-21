use pvlib::snow::*;

// --- fully_covered_nrel ---

#[test]
fn test_fully_covered_heavy_snowfall() {
    assert!(fully_covered_nrel(2.0, 1.0, 1.0));
}

#[test]
fn test_fully_covered_light_snowfall() {
    assert!(!fully_covered_nrel(0.5, 1.0, 1.0));
}

#[test]
fn test_fully_covered_zero_timestep() {
    assert!(!fully_covered_nrel(2.0, 0.0, 1.0));
}

#[test]
fn test_fully_covered_at_threshold() {
    assert!(fully_covered_nrel(1.0, 1.0, 1.0));
}

// --- coverage_nrel ---

#[test]
fn test_coverage_nrel_new_snowfall() {
    let cov = coverage_nrel(2.0, 100.0, -5.0, 30.0, 0.5, 1.0, 1.0, -80.0, 0.197);
    assert!((cov - 1.0).abs() < 1e-10, "Should be fully covered, got {}", cov);
}

#[test]
fn test_coverage_nrel_sliding() {
    let cov = coverage_nrel(0.0, 200.0, 5.0, 30.0, 0.8, 1.0, 1.0, -80.0, 0.197);
    assert!(cov < 0.8, "Coverage should decrease from sliding, got {}", cov);
    assert!(cov >= 0.0, "Coverage should not be negative, got {}", cov);
}

#[test]
fn test_coverage_nrel_no_sliding_cold() {
    let cov = coverage_nrel(0.0, 800.0, -10.0, 30.0, 0.8, 1.0, 1.0, -80.0, 0.197);
    assert!(
        (cov - 0.8).abs() < 1e-10,
        "Coverage should stay at 0.8 when can't slide, got {}",
        cov
    );
}

#[test]
fn test_coverage_nrel_clamps_to_zero() {
    let cov = coverage_nrel(0.0, 200.0, 5.0, 90.0, 0.01, 1.0, 1.0, -80.0, 10.0);
    assert!((cov - 0.0).abs() < 1e-10, "Coverage should clamp to 0, got {}", cov);
}

// --- dc_loss_nrel ---

#[test]
fn test_dc_loss_no_coverage() {
    assert!((dc_loss_nrel(0.0, 4) - 0.0).abs() < 1e-10);
}

#[test]
fn test_dc_loss_full_coverage() {
    assert!((dc_loss_nrel(1.0, 4) - 1.0).abs() < 1e-10);
}

#[test]
fn test_dc_loss_partial_coverage() {
    let loss = dc_loss_nrel(0.3, 4);
    assert!((loss - 0.5).abs() < 1e-10, "Expected 0.5, got {}", loss);
}

#[test]
fn test_dc_loss_single_string() {
    let loss = dc_loss_nrel(0.01, 1);
    assert!((loss - 1.0).abs() < 1e-10, "Single string: any coverage -> full loss");
}

#[test]
fn test_dc_loss_zero_strings() {
    assert!((dc_loss_nrel(0.5, 0) - 0.0).abs() < 1e-10);
}

// --- snow_slides ---

#[test]
fn test_snow_slides_steep_warm() {
    assert!(snow_slides(30.0, 10.0, 500.0));
}

#[test]
fn test_snow_slides_flat() {
    assert!(!snow_slides(5.0, 10.0, 500.0));
}

#[test]
fn test_snow_slides_cold() {
    assert!(!snow_slides(30.0, -5.0, 100.0));
}
