/// Angle from horizontal of the line from a point on the row slant length
/// to the bottom of the facing row.
///
/// # Arguments
/// * `surface_tilt` - Surface tilt angle in degrees from horizontal.
/// * `gcr` - Ground coverage ratio (row slant length / row spacing).
/// * `slant_height` - Fraction [0-1] of the module slant height.
///
/// # Returns
/// Ground angle (psi) in degrees.
pub fn ground_angle(surface_tilt: f64, gcr: f64, slant_height: f64) -> f64 {
    let x1 = gcr * slant_height * surface_tilt.to_radians().sin();
    let x2 = gcr * slant_height * surface_tilt.to_radians().cos() + 1.0;
    x1.atan2(x2).to_degrees()
}

/// Calculates the masking angle for rows in a standard PV array.
///
/// The masking angle is the elevation angle below which diffuse irradiance
/// is blocked by the preceding row.
///
/// # Arguments
/// * `surface_tilt` - Tilt of the PV modules in degrees.
/// * `gcr` - Ground Coverage Ratio (module length / row pitch).
/// * `slant_height` - Fraction [0-1] of the module slant height to evaluate.
///
/// # Returns
/// Masking angle in degrees.
pub fn masking_angle(surface_tilt: f64, gcr: f64, slant_height: f64) -> f64 {
    if gcr <= 0.0 {
        return 0.0;
    }
    let numerator = gcr * (1.0 - slant_height) * surface_tilt.to_radians().sin();
    let denominator = 1.0 - gcr * (1.0 - slant_height) * surface_tilt.to_radians().cos();
    (numerator / denominator).atan().to_degrees()
}

/// Average masking angle over the slant height of a row (Passias 1984).
///
/// # Arguments
/// * `surface_tilt` - Panel tilt from horizontal in degrees.
/// * `gcr` - Ground coverage ratio.
///
/// # Returns
/// Average masking angle in degrees.
pub fn masking_angle_passias(surface_tilt: f64, gcr: f64) -> f64 {
    let beta = surface_tilt.to_radians();
    let sin_b = beta.sin();
    let cos_b = beta.cos();

    if sin_b.abs() < 1e-10 || gcr <= 0.0 {
        return 0.0;
    }

    let x = 1.0 / gcr;

    let term1 = -x * sin_b * (2.0 * x * cos_b - (x * x + 1.0)).abs().ln() / 2.0;
    let term2 = (x * cos_b - 1.0) * ((x * cos_b - 1.0) / (x * sin_b)).atan();
    let term3 = (1.0 - x * cos_b) * (cos_b / sin_b).atan();
    let term4 = x * x.ln() * sin_b;

    let psi_avg = term1 + term2 + term3 + term4;

    if psi_avg.is_finite() {
        psi_avg.to_degrees()
    } else {
        0.0
    }
}

/// Projected solar zenith angle onto the plane perpendicular to a tracker axis.
///
/// # Arguments
/// * `solar_zenith` - Sun's apparent zenith in degrees.
/// * `solar_azimuth` - Sun's azimuth in degrees.
/// * `axis_tilt` - Axis tilt angle from horizontal in degrees.
/// * `axis_azimuth` - Axis azimuth angle in degrees (N=0, E=90, S=180, W=270).
///
/// # Returns
/// Projected solar zenith angle in degrees.
pub fn projected_solar_zenith_angle(
    solar_zenith: f64,
    solar_azimuth: f64,
    axis_tilt: f64,
    axis_azimuth: f64,
) -> f64 {
    let sin_sz = solar_zenith.to_radians().sin();
    let cos_aa = axis_azimuth.to_radians().cos();
    let sin_aa = axis_azimuth.to_radians().sin();
    let sin_at = axis_tilt.to_radians().sin();

    // Sun's x, y, z coordinates
    let sx = sin_sz * solar_azimuth.to_radians().sin();
    let sy = sin_sz * solar_azimuth.to_radians().cos();
    let sz = solar_zenith.to_radians().cos();

    // Project onto surface: Eq. (4) from Anderson & Mikofski (2020)
    let sx_prime = sx * cos_aa - sy * sin_aa;
    let sz_prime = sx * sin_aa * sin_at + sy * sin_at * cos_aa + sz * axis_tilt.to_radians().cos();

    // Eq. (5)
    sx_prime.atan2(sz_prime).to_degrees()
}

/// 1D shaded fraction for rows (e.g. single-axis trackers).
///
/// Based on Anderson & Jensen (2024).
///
/// # Arguments
/// * `solar_zenith` - Solar zenith angle in degrees.
/// * `solar_azimuth` - Solar azimuth angle in degrees.
/// * `axis_azimuth` - Axis azimuth in degrees.
/// * `shaded_row_rotation` - Rotation of the shaded row in degrees.
/// * `collector_width` - Vertical length of a tilted row.
/// * `pitch` - Axis-to-axis horizontal spacing.
/// * `axis_tilt` - Tilt of the rows axis from horizontal in degrees.
/// * `surface_to_axis_offset` - Distance between rotating axis and collector surface.
/// * `cross_axis_slope` - Angle of the plane containing row axes from horizontal in degrees.
///
/// # Returns
/// Shaded fraction [0-1].
pub fn shaded_fraction1d(
    solar_zenith: f64,
    solar_azimuth: f64,
    axis_azimuth: f64,
    shaded_row_rotation: f64,
    collector_width: f64,
    pitch: f64,
    axis_tilt: f64,
    surface_to_axis_offset: f64,
    cross_axis_slope: f64,
) -> f64 {
    // Use same rotation for shading row
    let shading_row_rotation = shaded_row_rotation;

    let psza = projected_solar_zenith_angle(solar_zenith, solar_azimuth, axis_tilt, axis_azimuth);

    let thetas_1_s_diff = shading_row_rotation - psza;
    let thetas_2_s_diff = shaded_row_rotation - psza;
    let theta_s_rotation_diff = psza - cross_axis_slope;

    let cos_theta_2_s_diff_abs = thetas_2_s_diff.to_radians().cos().abs().max(1e-6);
    let collector_width_safe = collector_width.max(1e-6);
    let cross_axis_cos = cross_axis_slope.to_radians().cos().max(1e-6);

    // Eq. (12) from Anderson & Jensen (2024)
    let t_asterisk = 0.5
        + thetas_1_s_diff.to_radians().cos().abs() / cos_theta_2_s_diff_abs / 2.0
        + (psza.signum()
            * surface_to_axis_offset
            / collector_width_safe
            / cos_theta_2_s_diff_abs
            * (thetas_2_s_diff.to_radians().sin() - thetas_1_s_diff.to_radians().sin()))
        - (pitch / collector_width_safe
            * theta_s_rotation_diff.to_radians().cos()
            / cos_theta_2_s_diff_abs
            / cross_axis_cos);

    t_asterisk.clamp(0.0, 1.0)
}

/// Computes the fraction of sky diffuse irradiance that passes the masking angle.
pub fn sky_diffuse_pass_equation(masking_angle: f64) -> f64 {
    (1.0 + masking_angle.to_radians().cos()) / 2.0
}

/// The diffuse irradiance loss caused by row-to-row sky diffuse shading.
///
/// Uses the Passias model assuming isotropic sky diffuse irradiance.
///
/// # Arguments
/// * `masking_angle` - Elevation angle below which diffuse irradiance is blocked in degrees.
///
/// # Returns
/// Fraction [0-1] of blocked sky diffuse irradiance.
pub fn sky_diffuse_passias(masking_angle: f64) -> f64 {
    1.0 - (masking_angle / 2.0).to_radians().cos().powi(2)
}
