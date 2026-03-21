use pvlib::pvarray::*;

// --- pvefficiency_adr ---

#[test]
fn test_adr_at_stc() {
    // At STC (1000 W/m2, 25C), efficiency should be close to k_a
    let eff = pvefficiency_adr(1000.0, 25.0, 0.20, -6.0, 0.02, 0.05, 0.05);
    // At reference: s=1, dt=0, s_o_ref = 10^(-6), v = ln(1/s_o_ref+1)/ln(1/s_o_ref+1) = 1
    // eff = k_a * ((1 + k_rs + k_rsh) * 1 - k_rs * 1 - k_rsh * 1) = k_a * 1 = 0.20
    assert!(
        (eff - 0.20).abs() < 0.01,
        "At STC, efficiency should be ~k_a=0.20, got {}",
        eff
    );
}

#[test]
fn test_adr_low_irradiance() {
    let eff_low = pvefficiency_adr(100.0, 25.0, 0.20, -6.0, 0.02, 0.05, 0.05);
    let eff_stc = pvefficiency_adr(1000.0, 25.0, 0.20, -6.0, 0.02, 0.05, 0.05);
    // At low irradiance, efficiency typically drops
    assert!(eff_low < eff_stc, "Low irradiance should have lower efficiency");
}

#[test]
fn test_adr_high_temp() {
    let eff_25 = pvefficiency_adr(1000.0, 25.0, 0.20, -6.0, 0.02, 0.05, 0.05);
    let eff_60 = pvefficiency_adr(1000.0, 60.0, 0.20, -6.0, 0.02, 0.05, 0.05);
    // Higher temp with positive tc_d increases dark current -> lower efficiency
    assert!(eff_60 < eff_25, "Higher temp should reduce efficiency");
}

// --- huld ---

#[test]
fn test_huld_at_stc() {
    // At STC with zero k coefficients, P = G' * pdc0 = 1 * 300 = 300
    let p = huld(1000.0, 25.0, 300.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!((p - 300.0).abs() < 1e-6, "At STC with zero k, power should be pdc0, got {}", p);
}

#[test]
fn test_huld_scales_with_irradiance() {
    let p500 = huld(500.0, 25.0, 300.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let p1000 = huld(1000.0, 25.0, 300.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(
        (p500 / p1000 - 0.5).abs() < 1e-6,
        "Power should scale linearly with G when k=0"
    );
}

#[test]
fn test_huld_zero_irradiance() {
    let p = huld(0.0, 25.0, 300.0, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    assert!((p - 0.0).abs() < 1e-10, "Zero irradiance should give zero power");
}
