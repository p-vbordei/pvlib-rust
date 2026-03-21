use pvlib::clearsky::{bird, bird_default, lookup_linke_turbidity};

#[test]
fn test_bird_known_values() {
    // Reference values computed from the Python pvlib Bird model formulas:
    // zenith=30, airmass_relative=1.5, aod380=0.15, aod500=0.1,
    // precipitable_water=1.0, ozone=0.3, pressure=101325, dni_extra=1364,
    // asymmetry=0.85, albedo=0.2
    let result = bird(30.0, 1.5, 0.15, 0.1, 1.0, 0.3, 101325.0, 1364.0, 0.85, 0.2);

    assert!(
        (result.dni - 872.1025).abs() < 1.0,
        "DNI expected ~872.1, got {}",
        result.dni
    );
    assert!(
        (result.direct_horizontal - 755.2629).abs() < 1.0,
        "Direct horizontal expected ~755.3, got {}",
        result.direct_horizontal
    );
    assert!(
        (result.ghi - 896.6138).abs() < 1.0,
        "GHI expected ~896.6, got {}",
        result.ghi
    );
    assert!(
        (result.dhi - 141.3509).abs() < 1.0,
        "DHI expected ~141.4, got {}",
        result.dhi
    );
}

#[test]
fn test_bird_nighttime() {
    let result = bird(95.0, 1.0, 0.1, 0.1, 1.0, 0.3, 101325.0, 1364.0, 0.85, 0.2);
    assert_eq!(result.ghi, 0.0);
    assert_eq!(result.dni, 0.0);
    assert_eq!(result.dhi, 0.0);
    assert_eq!(result.direct_horizontal, 0.0);
}

#[test]
fn test_bird_default_wrapper() {
    let result = bird_default(30.0, 1.5, 0.15, 0.1, 1.0);
    // Should produce the same values as calling bird() with defaults
    let result_full = bird(30.0, 1.5, 0.15, 0.1, 1.0, 0.3, 101325.0, 1364.0, 0.85, 0.2);
    assert_eq!(result, result_full);
}

#[test]
fn test_bird_zenith_zero() {
    // Sun directly overhead
    let result = bird(0.0, 1.0, 0.1, 0.05, 1.0, 0.3, 101325.0, 1364.0, 0.85, 0.2);
    assert!(result.ghi > 900.0, "GHI at zenith=0 should be high, got {}", result.ghi);
    assert!(result.dni > 900.0, "DNI at zenith=0 should be high, got {}", result.dni);
    // At zenith=0, direct_horizontal should equal DNI
    assert!(
        (result.direct_horizontal - result.dni).abs() < 0.1,
        "At zenith=0, direct_horizontal should equal DNI"
    );
}

#[test]
fn test_lookup_linke_turbidity_polar() {
    let lt = lookup_linke_turbidity(70.0, 25.0, 1); // winter in Arctic
    assert!(lt >= 1.5 && lt <= 3.0, "Polar winter turbidity expected 2.0, got {}", lt);
}

#[test]
fn test_lookup_linke_turbidity_tropical() {
    let lt = lookup_linke_turbidity(5.0, 30.0, 7); // tropical summer
    assert!(lt >= 3.5 && lt <= 5.5, "Tropical summer turbidity expected ~4.5, got {}", lt);
}

#[test]
fn test_lookup_linke_turbidity_desert() {
    let lt = lookup_linke_turbidity(28.0, 30.0, 7); // Sahara summer
    assert!(lt >= 3.0 && lt <= 5.0, "Desert turbidity expected ~3.5-4.0, got {}", lt);
}

#[test]
fn test_lookup_linke_turbidity_seasonal() {
    let winter = lookup_linke_turbidity(50.0, 10.0, 1);
    let summer = lookup_linke_turbidity(50.0, 10.0, 7);
    assert!(summer > winter, "Summer turbidity should exceed winter");
}
