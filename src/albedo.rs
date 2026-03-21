/// Surface albedo values for common surfaces.
pub fn surface_albedo(surface: &str) -> Option<f64> {
    match surface {
        "urban" => Some(0.18),
        "grass" => Some(0.20),
        "fresh grass" => Some(0.26),
        "soil" => Some(0.17),
        "sand" => Some(0.40),
        "snow" => Some(0.65),
        "fresh snow" => Some(0.75),
        "asphalt" => Some(0.12),
        "concrete" => Some(0.30),
        "aluminum" => Some(0.85),
        "copper" => Some(0.74),
        "fresh steel" => Some(0.35),
        "dirty steel" => Some(0.08),
        "sea" => Some(0.06),
        _ => None,
    }
}

/// Water color and roughness coefficients for the Dvoracek inland water albedo model.
fn water_coefficients(surface_condition: &str) -> Option<(f64, f64)> {
    match surface_condition {
        "clear_water_no_waves" => Some((0.13, 0.29)),
        "clear_water_ripples_up_to_2.5cm" => Some((0.16, 0.70)),
        "clear_water_ripples_larger_than_2.5cm_occasional_whitecaps" => Some((0.23, 1.25)),
        "clear_water_frequent_whitecaps" => Some((0.30, 2.0)),
        "green_water_ripples_up_to_2.5cm" => Some((0.22, 0.70)),
        "muddy_water_no_waves" => Some((0.19, 0.29)),
        _ => None,
    }
}

/// Estimation of albedo for inland water bodies using the Dvoracek model.
///
/// albedo = c^(r * sin(elevation) + 1)
///
/// # Arguments
/// * `solar_elevation` - Sun elevation angle in degrees.
/// * `surface_condition` - Water surface condition string.
///
/// # Returns
/// Albedo value (dimensionless).
pub fn inland_water_dvoracek(solar_elevation: f64, surface_condition: &str) -> f64 {
    let (c, r) = water_coefficients(surface_condition)
        .unwrap_or_else(|| panic!("Unknown surface condition: {}", surface_condition));

    let elev = solar_elevation.max(0.0);
    c.powf(r * elev.to_radians().sin() + 1.0)
}

/// Estimation of albedo for inland water bodies using custom coefficients.
///
/// albedo = color_coeff^(wave_roughness_coeff * sin(elevation) + 1)
///
/// # Arguments
/// * `solar_elevation` - Sun elevation angle in degrees.
/// * `color_coeff` - Water color coefficient.
/// * `wave_roughness_coeff` - Water wave roughness coefficient.
///
/// # Returns
/// Albedo value (dimensionless).
pub fn inland_water_dvoracek_custom(
    solar_elevation: f64,
    color_coeff: f64,
    wave_roughness_coeff: f64,
) -> f64 {
    let elev = solar_elevation.max(0.0);
    color_coeff.powf(wave_roughness_coeff * elev.to_radians().sin() + 1.0)
}
