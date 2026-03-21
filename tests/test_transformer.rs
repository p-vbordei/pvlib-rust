use pvlib::transformer::*;

#[test]
fn test_simple_efficiency_no_losses() {
    // With zero losses, output should equal input
    let out = simple_efficiency(500.0, 1000.0, 0.0, 0.0);
    assert!((out - 500.0).abs() < 1e-6, "Zero losses: output should equal input, got {}", out);
}

#[test]
fn test_simple_efficiency_with_losses() {
    let out = simple_efficiency(1000.0, 1000.0, 0.01, 0.01);
    // Output should be less than input
    assert!(out < 1000.0, "Should have losses, got {}", out);
    assert!(out > 900.0, "Losses shouldn't be extreme, got {}", out);
}

#[test]
fn test_simple_efficiency_scales() {
    // Higher input -> higher output
    let out1 = simple_efficiency(500.0, 1000.0, 0.01, 0.01);
    let out2 = simple_efficiency(800.0, 1000.0, 0.01, 0.01);
    assert!(out2 > out1, "Higher input should give higher output");
}

#[test]
fn test_simple_efficiency_no_load_loss_constant() {
    // No-load loss is constant regardless of power level
    let out_low = simple_efficiency(100.0, 1000.0, 0.02, 0.0);
    let out_high = simple_efficiency(900.0, 1000.0, 0.02, 0.0);
    // With only no_load loss and no load-dependent loss:
    // loss is a constant fraction
    let loss_low = 100.0 - out_low;
    let loss_high = 900.0 - out_high;
    // Both should have similar absolute loss
    assert!(
        (loss_low - loss_high).abs() < 50.0,
        "No-load losses should be roughly constant: {} vs {}",
        loss_low, loss_high
    );
}
