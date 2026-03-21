use pvlib::soiling::*;

// --- HSU model ---

#[test]
fn test_hsu_clean_start() {
    let (ratio, mass) = hsu(0.0, 5.0, 20.0, 30e-6, 60e-6, 0.0009, 0.004, 3600.0, 0.0);
    assert!(ratio <= 1.0 && ratio > 0.9, "Ratio should be near 1 for small accumulation, got {}", ratio);
    assert!(mass > 0.0, "Mass should accumulate");
}

#[test]
fn test_hsu_rain_cleaning() {
    let (ratio, mass) = hsu(10.0, 5.0, 20.0, 30e-6, 60e-6, 0.0009, 0.004, 3600.0, 1.0);
    assert!(ratio > 0.99, "Should be nearly clean after rain, got {}", ratio);
    assert!(mass < 0.01, "Mass should be near zero after cleaning, got {}", mass);
}

#[test]
fn test_hsu_mass_accumulates() {
    let (ratio_low, _) = hsu(0.0, 5.0, 20.0, 30e-6, 60e-6, 0.0009, 0.004, 3600.0, 0.0);
    let (ratio_high, _) = hsu(0.0, 5.0, 20.0, 30e-6, 60e-6, 0.0009, 0.004, 3600.0, 5.0);
    assert!(
        ratio_high < ratio_low,
        "More accumulated mass should give lower ratio: {} vs {}",
        ratio_high, ratio_low
    );
}

#[test]
fn test_hsu_soiling_ratio_clean() {
    assert!((hsu_soiling_ratio(0.0) - 1.0).abs() < 1e-6);
}

#[test]
fn test_hsu_soiling_ratio_decreases_with_mass() {
    let r1 = hsu_soiling_ratio(0.1);
    let r2 = hsu_soiling_ratio(1.0);
    let r3 = hsu_soiling_ratio(10.0);
    assert!(r1 > r2, "More mass -> lower ratio");
    assert!(r2 > r3, "More mass -> lower ratio");
}

#[test]
fn test_hsu_mass_rate_horizontal() {
    let m0 = hsu_mass_rate(30e-6, 60e-6, 0.0009, 0.004, 0.0, 3600.0);
    let m45 = hsu_mass_rate(30e-6, 60e-6, 0.0009, 0.004, 45.0, 3600.0);
    assert!(m0 > m45, "Horizontal should collect more than tilted");
}

// --- Kimber model ---

#[test]
fn test_kimber_rain_cleans() {
    let soiling = kimber(10.0, 6.0, 0.0015, 0.3, 1.0, 0.1, false);
    assert!((soiling - 0.0).abs() < 1e-10, "Rain above threshold should clean");
}

#[test]
fn test_kimber_accumulates() {
    let soiling = kimber(0.0, 6.0, 0.0015, 0.3, 1.0, 0.01, false);
    let expected = 0.01 + 0.0015;
    assert!(
        (soiling - expected).abs() < 1e-10,
        "Should accumulate: expected {}, got {}",
        expected, soiling
    );
}

#[test]
fn test_kimber_grace_period() {
    let soiling = kimber(0.0, 6.0, 0.0015, 0.3, 1.0, 0.05, true);
    assert!((soiling - 0.0).abs() < 1e-10, "Grace period should return 0");
}

#[test]
fn test_kimber_max_soiling() {
    let soiling = kimber(0.0, 6.0, 0.0015, 0.3, 1.0, 0.299, false);
    let expected = (0.299 + 0.0015_f64).min(0.3);
    assert!(
        (soiling - expected).abs() < 1e-10,
        "Should cap at max_soiling: expected {}, got {}",
        expected, soiling
    );
}

#[test]
fn test_kimber_no_rain_no_grace() {
    let soiling = kimber(0.0, 6.0, 0.002, 0.5, 1.0, 0.0, false);
    assert!(
        (soiling - 0.002).abs() < 1e-10,
        "One day soiling from zero: expected 0.002, got {}",
        soiling
    );
}

// --- accumulation_model ---

#[test]
fn test_accumulation_model_rain_cleans() {
    let acc = accumulation_model(10.0, 5.0, 20.0, 30.0, 60.0);
    assert!((acc - 0.0).abs() < 1e-10, "Rain should clean");
}

#[test]
fn test_accumulation_model_positive() {
    let acc = accumulation_model(0.0, 5.0, 20.0, 30.0, 60.0);
    assert!(acc > 0.0, "Should accumulate without rain");
}
