/// Calculate single-axis tracker positions with backtracking support.
///
/// Determines the rotation angle of a single-axis tracker using the
/// projected solar zenith angle approach from Anderson & Mikofski (2020).
/// Backtracking uses the slope-aware method (Eq. 14 from the same reference).
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
///
/// # References
/// Anderson, K., and Mikofski, M., 2020, "Slope-Aware Backtracking for
/// Single-Axis Trackers", Technical Report NREL/TP-5K00-76626.
pub fn singleaxis(
    solar_zenith: f64,
    solar_azimuth: f64,
    axis_tilt: f64,
    axis_azimuth: f64,
    max_angle: f64,
    backtrack: bool,
    gcr: f64,
) -> (f64, f64, f64) {
    use crate::shading::projected_solar_zenith_angle;
    use crate::irradiance::aoi as irr_aoi;

    // Sun below horizon — return stow position
    if solar_zenith >= 90.0 {
        let (st, sa) = calc_surface_orientation(0.0, axis_tilt, axis_azimuth);
        let a = irr_aoi(st, sa, solar_zenith, solar_azimuth);
        return (st, sa, a);
    }

    // Calculate cross-axis tilt from the axis geometry.
    // For a flat site (slope_tilt=0) cross_axis_tilt is 0.
    // Users with sloped terrain should use calc_cross_axis_tilt externally.
    let cross_axis_tilt: f64 = 0.0;

    // Ideal rotation angle via projected solar zenith angle
    // (Anderson & Mikofski 2020, handles arbitrary axis_tilt)
    let mut tracker_theta = projected_solar_zenith_angle(
        solar_zenith,
        solar_azimuth,
        axis_tilt,
        axis_azimuth,
    );

    // Backtracking — slope-aware method (Anderson & Mikofski 2020, Eq. 14)
    if backtrack && gcr > 0.0 {
        let axes_distance = 1.0 / (gcr * cross_axis_tilt.to_radians().cos());

        // temp = |axes_distance * cos(tracker_theta - cross_axis_tilt)|
        let temp = (axes_distance * (tracker_theta - cross_axis_tilt).to_radians().cos()).abs();

        if temp < 1.0 {
            // Backtracking correction needed
            let omega_correction = -tracker_theta.signum() * temp.acos().to_degrees();
            tracker_theta += omega_correction;
        }
        // else: no row-to-row shade, no correction needed (Eqs. 15-16)
    }

    // Apply hardware limits
    tracker_theta = tracker_theta.clamp(-max_angle, max_angle);

    // Calculate surface tilt & azimuth using full orientation model
    let (surface_tilt, surface_azimuth) =
        calc_surface_orientation(tracker_theta, axis_tilt, axis_azimuth);

    // Angle of incidence
    let aoi = irr_aoi(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth);

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

