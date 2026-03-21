use pvlib::singlediode::*;

// Standard test parameters from pvlib-python docstrings
fn test_params() -> Bishop88Params {
    Bishop88Params::new(1.0, 9e-10, 4.0, 5000.0, 4.0)
}

// --- estimate_voc ---

#[test]
fn test_estimate_voc_positive() {
    let voc = estimate_voc(1.0, 9e-10, 4.0);
    // Voc_est = 4 * ln(1/9e-10 + 1) ~ 4 * ln(1.11e9) ~ 4 * 20.83 ~ 83.3
    assert!(voc > 50.0 && voc < 100.0, "Voc estimate should be reasonable, got {}", voc);
}

#[test]
fn test_estimate_voc_proportional_to_nnsvth() {
    let voc1 = estimate_voc(1.0, 9e-10, 4.0);
    let voc2 = estimate_voc(1.0, 9e-10, 8.0);
    assert!((voc2 / voc1 - 2.0).abs() < 0.01, "Voc should scale linearly with nNsVth");
}

// --- bishop88 basic ---

#[test]
fn test_bishop88_short_circuit() {
    // At diode_voltage = 0, current should be close to photocurrent
    let p = test_params();
    let out = bishop88(0.0, &p);
    assert!(
        (out.current - 1.0).abs() < 0.01,
        "Isc should be ~photocurrent at Vd=0, got {}",
        out.current
    );
    // voltage should be negative (V = 0 - I*Rs)
    assert!(out.voltage < 0.0, "V should be negative at Vd=0 with Rs>0");
}

#[test]
fn test_bishop88_at_voc() {
    // At open circuit, current should be ~0
    let p = test_params();
    let voc_est = estimate_voc(p.photocurrent, p.saturation_current, p.n_ns_vth);
    let out = bishop88(voc_est, &p);
    // Current should be very close to zero (since Rsh is large)
    // With Rs=4 and Rsh=5000, some leakage is expected
    assert!(
        out.current.abs() < 0.02,
        "Current at Voc should be ~0, got {}",
        out.current
    );
}

#[test]
fn test_bishop88_power_positive_midpoint() {
    // At a midpoint diode voltage, power should be positive
    let p = test_params();
    let voc_est = estimate_voc(p.photocurrent, p.saturation_current, p.n_ns_vth);
    let out = bishop88(voc_est * 0.5, &p);
    assert!(out.power > 0.0, "Power should be positive at midpoint, got {}", out.power);
}

// --- bishop88_i_from_v ---

#[test]
fn test_bishop88_i_from_v_at_zero() {
    let p = test_params();
    let i = bishop88_i_from_v(0.0, &p);
    // At V=0, current should be close to Isc ~ photocurrent
    assert!(
        (i - 1.0).abs() < 0.01,
        "I at V=0 should be ~Isc=1A, got {}",
        i
    );
}

#[test]
fn test_bishop88_i_from_v_at_voc() {
    let p = test_params();
    // Find Voc first
    let voc = bishop88_v_from_i(0.0, &p);
    let i = bishop88_i_from_v(voc, &p);
    assert!(
        i.abs() < 0.01,
        "I at Voc should be ~0, got {}",
        i
    );
}

// --- bishop88_v_from_i ---

#[test]
fn test_bishop88_v_from_i_at_zero() {
    let p = test_params();
    let v = bishop88_v_from_i(0.0, &p);
    // Voc should be positive and reasonable
    assert!(v > 0.0, "Voc should be positive, got {}", v);
    assert!(v < 100.0, "Voc should be < 100V, got {}", v);
}

#[test]
fn test_bishop88_v_from_i_at_isc() {
    let p = test_params();
    let v = bishop88_v_from_i(1.0, &p);
    // At Isc, voltage should be ~0 (or slightly negative due to Rs)
    assert!(
        v.abs() < 5.0,
        "V at Isc should be near 0, got {}",
        v
    );
}

// --- bishop88_mpp ---

#[test]
fn test_bishop88_mpp_positive_power() {
    let p = test_params();
    let mpp = bishop88_mpp(&p);
    assert!(mpp.power > 0.0, "MPP power should be positive, got {}", mpp.power);
    assert!(mpp.current > 0.0, "MPP current should be positive, got {}", mpp.current);
    assert!(mpp.voltage > 0.0, "MPP voltage should be positive, got {}", mpp.voltage);
}

#[test]
fn test_bishop88_mpp_less_than_voc() {
    let p = test_params();
    let mpp = bishop88_mpp(&p);
    let voc = bishop88_v_from_i(0.0, &p);
    assert!(
        mpp.voltage < voc,
        "MPP voltage ({}) should be less than Voc ({})",
        mpp.voltage, voc
    );
}

#[test]
fn test_bishop88_mpp_current_less_than_isc() {
    let p = test_params();
    let mpp = bishop88_mpp(&p);
    let isc = bishop88_i_from_v(0.0, &p);
    assert!(
        mpp.current < isc,
        "MPP current ({}) should be less than Isc ({})",
        mpp.current, isc
    );
}

#[test]
fn test_bishop88_mpp_is_maximum() {
    // Verify that the MPP power is indeed a maximum by checking nearby points
    let p = test_params();
    let mpp = bishop88_mpp(&p);
    let voc_est = estimate_voc(p.photocurrent, p.saturation_current, p.n_ns_vth);

    // Sample 100 points along the IV curve and verify none has higher power
    for i in 0..100 {
        let vd = voc_est * (i as f64) / 100.0;
        let out = bishop88(vd, &p);
        assert!(
            out.power <= mpp.power + 1e-6,
            "Found higher power ({}) than MPP ({}) at vd={}",
            out.power, mpp.power, vd
        );
    }
}

// --- Consistency tests ---

#[test]
fn test_bishop88_i_from_v_roundtrip() {
    let p = test_params();
    let v_target = 30.0;
    let i = bishop88_i_from_v(v_target, &p);
    let v_back = bishop88_v_from_i(i, &p);
    assert!(
        (v_back - v_target).abs() < 0.01,
        "Roundtrip V: expected {}, got {}",
        v_target, v_back
    );
}

// --- Legacy i_from_v / v_from_i ---

#[test]
fn test_legacy_i_from_v() {
    let i = i_from_v(0.0, 1.0, 9e-10, 4.0, 5000.0, 4.0);
    assert!((i - 1.0).abs() < 0.01, "Legacy i_from_v at V=0 should give ~Isc, got {}", i);
}

#[test]
fn test_legacy_v_from_i() {
    let v = v_from_i(0.0, 1.0, 9e-10, 4.0, 5000.0, 4.0);
    assert!(v > 0.0, "Legacy v_from_i at I=0 should give positive Voc, got {}", v);
}

#[test]
fn test_bishop88_gradients_computed() {
    let p = test_params();
    let voc_est = estimate_voc(p.photocurrent, p.saturation_current, p.n_ns_vth);
    let (out, grad) = bishop88_with_gradients(voc_est * 0.5, &p);
    assert!(out.power > 0.0);
    // dv/dvd should be > 0 (monotonic mapping)
    assert!(grad.dv_dvd > 0.0, "dv/dvd should be positive, got {}", grad.dv_dvd);
    // di/dvd should be negative (current decreases as voltage increases)
    assert!(grad.di_dvd < 0.0, "di/dvd should be negative, got {}", grad.di_dvd);
}
