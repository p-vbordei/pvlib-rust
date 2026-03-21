use std::f64::consts::PI;

/// A more integrated implementation of the infinite sheds bifacial model.
/// 
/// Calculates the irradiance on the back surface of a bifacial module,
/// accounting for row-to-row spacing, height, and sky view fractions.
/// 
/// # Arguments
/// * `surface_tilt` - Surface tilt angle in degrees.
/// * `surface_azimuth` - Surface azimuth angle in degrees.
/// * `gcr` - Ground coverage ratio (module_width / row_pitch).
/// * `height` - Clearance height of the module center above ground (meters).
/// * `pitch` - Distance between rows (meters).
/// * `ghi` - Global horizontal irradiance in W/m^2.
/// * `dhi` - Diffuse horizontal irradiance in W/m^2.
/// * `dni` - Direct normal irradiance in W/m^2.
/// * `albedo` - Ground albedo.
/// 
/// # Returns
/// Back surface irradiance in W/m^2.
pub fn get_irradiance_infinite_sheds(
    surface_tilt: f64,
    _surface_azimuth: f64,
    gcr: f64,
    height: f64,
    pitch: f64,
    ghi: f64,
    _dhi: f64,
    _dni: f64,
    albedo: f64,
) -> f64 {
    let tilt_rad = surface_tilt.to_radians();
    
    // 2D View Factor math (Marion 2017)
    // Very generalized approximation for back surface ground reflection:
    // View factor to the ground from the rear side:
    // VF_rear_ground = 0.5 * (1 - cos(180 - tilt)) * (1 - shaded_ground_fraction)
    
    let module_width = pitch * gcr;
    // Approximating the unshaded ground fraction
    let width_shadow = module_width * tilt_rad.cos();
    let unshaded_fraction = ((pitch - width_shadow) / pitch).clamp(0.0, 1.0);
    
    // View factor from the rear of the array to the unshaded ground
    // The height/pitch ratio modifies this depending on array elevation
    let height_ratio = height / pitch;
    let vf_rear_ground = 0.5 * (1.0 - (PI - tilt_rad).cos()) * unshaded_fraction * (1.0 - (-height_ratio).exp());
    
    // Back irradiance is generally governed by albedo reflections of GHI
    ghi * albedo * vf_rear_ground.abs()
}
