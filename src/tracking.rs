/// Calculate single-axis tracker positions with backtracking support.
/// 
/// # Arguments
/// * `solar_zenith` - Apparent solar zenith angle in degrees.
/// * `solar_azimuth` - Apparent solar azimuth angle in degrees.
/// * `axis_tilt` - Tilt of the tracker axis from horizontal in degrees.
/// * `axis_azimuth` - Azimuth of the tracker axis in degrees (North=0, East=90).
/// * `max_angle` - Maximum rotation angle of the tracker from horizontal (e.g. 45.0 or 60.0 degrees).
/// * `backtrack` - Enable backtracking (True).
/// * `gcr` - Ground Coverage Ratio (module width / row pitch).
/// 
/// # Returns
/// A tuple containing `(surface_tilt, surface_azimuth, aoi)`. All in degrees.
pub fn singleaxis(
    solar_zenith: f64,
    solar_azimuth: f64,
    axis_tilt: f64,
    axis_azimuth: f64,
    max_angle: f64,
    backtrack: bool,
    gcr: f64,
) -> (f64, f64, f64) {
    if solar_zenith >= 90.0 { return (0.0, axis_azimuth, 90.0); }

    let sz_rad = solar_zenith.to_radians();
    let sa_rad = solar_azimuth.to_radians();
    let aa_rad = axis_azimuth.to_radians();
    let _at_rad = axis_tilt.to_radians(); // Assuming HSAT for MVP

    // Ideal Rotational angle (un-limited)
    let rot_angle_rad = (sz_rad.tan() * (sa_rad - aa_rad).sin()).atan();
    let mut rot_deg = rot_angle_rad.to_degrees();

    // Backtracking logic (standard simplified geometric implementation)
    if backtrack && gcr > 0.0 {
        let temp = (rot_angle_rad.cos() * gcr).clamp(-1.0, 1.0);
        let shade_angle = temp.acos();
        
        if rot_angle_rad.abs() > shade_angle {
            let bt_angle_rad = rot_angle_rad.signum() * (rot_angle_rad.abs() - shade_angle);
            rot_deg -= bt_angle_rad.to_degrees();
        }
    }

    // Apply hardware limits
    rot_deg = rot_deg.clamp(-max_angle, max_angle);

    let surface_tilt = rot_deg.abs();
    let surface_azimuth = if rot_deg >= 0.0 {
        (axis_azimuth - 90.0).rem_euclid(360.0)
    } else {
        (axis_azimuth + 90.0).rem_euclid(360.0)
    };

    let cos_aoi = sz_rad.cos() * surface_tilt.to_radians().cos()
        + sz_rad.sin() * surface_tilt.to_radians().sin() * (sa_rad - surface_azimuth.to_radians()).cos();
    
    let aoi = cos_aoi.clamp(-1.0, 1.0).acos().to_degrees();

    (surface_tilt, surface_azimuth, aoi)
}

/// Calculate tracking axis tilt from slope tilt and azimuths.
pub fn calc_axis_tilt(slope_azimuth: f64, slope_tilt: f64, axis_azimuth: f64) -> f64 {
    let sa_rad = slope_azimuth.to_radians();
    let st_rad = slope_tilt.to_radians();
    let aa_rad = axis_azimuth.to_radians();
    
    let axis_tilt_rad = (st_rad.tan() * (aa_rad - sa_rad).cos()).atan();
    axis_tilt_rad.to_degrees()
}

/// Calculate the surface tilt and azimuth angles for a given tracker rotation.
///
/// # Arguments
/// * `tracker_theta` - Tracker rotation angle (degrees). Right-handed rotation
///   around the axis defined by `axis_tilt` and `axis_azimuth`.
/// * `axis_tilt` - Tilt of the axis of rotation with respect to horizontal (degrees).
/// * `axis_azimuth` - Compass direction along which the axis of rotation lies (degrees).
///
/// # Returns
/// A tuple `(surface_tilt, surface_azimuth)` in degrees.
///
/// # References
/// Marion, W.F. and Dobos, A.P., 2013, "Rotation Angle for the Optimum Tracking
/// of One-Axis Trackers", NREL/TP-6A20-58891.
pub fn calc_surface_orientation(tracker_theta: f64, axis_tilt: f64, axis_azimuth: f64) -> (f64, f64) {
    let tt_rad = tracker_theta.to_radians();
    let at_rad = axis_tilt.to_radians();

    // Surface tilt: acos(cos(tracker_theta) * cos(axis_tilt))
    let surface_tilt_rad = (tt_rad.cos() * at_rad.cos()).clamp(-1.0, 1.0).acos();
    let surface_tilt = surface_tilt_rad.to_degrees();

    // Surface azimuth: axis_azimuth + azimuth_delta
    let sin_st = surface_tilt_rad.sin();

    let azimuth_delta = if sin_st.abs() < 1e-10 {
        // surface_tilt ~= 0, azimuth is arbitrary; use 90 per pvlib convention
        90.0
    } else {
        // azimuth_delta = asin(sin(tracker_theta) / sin(surface_tilt))
        let raw = (tt_rad.sin() / sin_st).clamp(-1.0, 1.0).asin().to_degrees();

        if tracker_theta.abs() < 90.0 {
            raw
        } else {
            -raw + tracker_theta.signum() * 180.0
        }
    };

    let surface_azimuth = (axis_azimuth + azimuth_delta).rem_euclid(360.0);

    (surface_tilt, surface_azimuth)
}

/// Calculate cross-axis tilt.
pub fn calc_cross_axis_tilt(slope_azimuth: f64, slope_tilt: f64, axis_azimuth: f64, axis_tilt: f64) -> f64 {
    let sa_rad = slope_azimuth.to_radians();
    let st_rad = slope_tilt.to_radians();
    let aa_rad = axis_azimuth.to_radians();
    let at_rad = axis_tilt.to_radians();
    
    let cross_axis_tilt_rad = (st_rad.tan() * (aa_rad - sa_rad).sin() * at_rad.cos()).atan();
    cross_axis_tilt_rad.to_degrees()
}

