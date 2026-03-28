use pvlib::albedo::*;

#[test]
fn test_surface_albedo_known() {
    assert!((surface_albedo("grass").unwrap() - 0.20).abs() < 1e-10);
    assert!((surface_albedo("snow").unwrap() - 0.65).abs() < 1e-10);
    assert!((surface_albedo("sand").unwrap() - 0.40).abs() < 1e-10);
}

#[test]
fn test_surface_albedo_unknown() {
    assert!(surface_albedo("mars_regolith").is_none());
}

#[test]
fn test_surface_albedo_range() {
    for surface in ["urban", "grass", "soil", "sand", "snow", "asphalt", "concrete", "sea"] {
        let a = surface_albedo(surface).unwrap();
        assert!((0.0..=1.0).contains(&a), "{} albedo out of range: {}", surface, a);
    }
}

#[test]
fn test_inland_water_dvoracek_high_elevation() {
    // At high sun elevation, albedo should be low (less reflection)
    let a = inland_water_dvoracek(60.0, "clear_water_no_waves");
    assert!(a > 0.0 && a < 0.5, "High elevation albedo should be low, got {}", a);
}

#[test]
fn test_inland_water_dvoracek_low_elevation() {
    // At low sun elevation, more reflection
    let a_low = inland_water_dvoracek(5.0, "clear_water_no_waves");
    let a_high = inland_water_dvoracek(60.0, "clear_water_no_waves");
    assert!(a_low > a_high, "Lower elevation should have higher albedo: {} vs {}", a_low, a_high);
}

#[test]
fn test_inland_water_dvoracek_negative_elevation() {
    // Negative elevation gets clamped to 0
    let a0 = inland_water_dvoracek(0.0, "clear_water_no_waves");
    let a_neg = inland_water_dvoracek(-10.0, "clear_water_no_waves");
    assert!((a0 - a_neg).abs() < 1e-10, "Negative elevation should clamp to 0");
}

#[test]
fn test_inland_water_dvoracek_custom() {
    let a = inland_water_dvoracek_custom(45.0, 0.13, 0.29);
    assert!(a > 0.0 && a < 1.0, "Custom albedo should be in range, got {}", a);
}
