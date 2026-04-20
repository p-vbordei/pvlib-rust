/// ASHRAE Incidence Angle Modifier (IAM) model.
/// 
/// `iam = 1 - b0 * (1/cos(aoi) - 1)`
/// 
/// # Arguments
/// * `aoi` - Angle of incidence in degrees.
/// * `b0` - IAM parameter (typically 0.05).
/// 
/// # Returns
/// IAM multiplier (dimensionless).
#[inline]
pub fn ashrae(aoi: f64, b0: f64) -> f64 {
    let aoi_rad = aoi.to_radians();
    if aoi_rad >= std::f64::consts::PI / 2.0 {
        return 0.0;
    }
    let cos_aoi = aoi_rad.cos();
    if cos_aoi <= 0.0 {
        return 0.0;
    }
    
    let iam = 1.0 - b0 * (1.0 / cos_aoi - 1.0);
    iam.clamp(0.0, 1.0)
}

/// Martin and Ruiz Incidence Angle Modifier (IAM) model.
/// 
/// # Arguments
/// * `aoi` - Angle of incidence in degrees.
/// * `a_r` - Angular losses coefficient.
/// 
/// # Returns
/// IAM multiplier (dimensionless).
#[inline]
pub fn martin_ruiz(aoi: f64, a_r: f64) -> f64 {
    if aoi >= 90.0 {
        return 0.0;
    }
    let aoi_rad = aoi.to_radians();
    // IAM = (1 - exp(-cos(aoi) / a_r)) / (1 - exp(-1 / a_r))
    let numerator = 1.0 - (-aoi_rad.cos() / a_r).exp();
    let denominator = 1.0 - (-1.0 / a_r).exp();
    if denominator == 0.0 {
        return 1.0;
    }
    (numerator / denominator).clamp(0.0, 1.0)
}

/// Physical IAM model using Fresnel/Snell's law.
///
/// Calculates transmission through a glass cover accounting for
/// reflection and absorption losses.
///
/// # Arguments
/// * `aoi` - Angle of incidence in degrees.
/// * `n` - Effective index of refraction (default 1.526 for glass).
/// * `k` - Glazing extinction coefficient in 1/meters (default 4.0).
/// * `l` - Glazing thickness in meters (default 0.002).
///
/// # Returns
/// IAM multiplier (dimensionless).
#[inline]
pub fn physical(aoi: f64, n: f64, k: f64, l: f64) -> f64 {
    if aoi.is_nan() {
        return f64::NAN;
    }
    if aoi.abs() >= 90.0 {
        return 0.0;
    }

    let aoi_rad = aoi.to_radians();
    let cos_aoi = aoi_rad.cos().max(0.0);
    let sin_aoi = (1.0 - cos_aoi * cos_aoi).sqrt();

    // Snell's law: refraction angle
    let sin_refr = sin_aoi / n;
    let cos_refr = (1.0 - sin_refr * sin_refr).sqrt();

    // Fresnel reflectance for s and p polarization
    let n1_cos1 = cos_aoi; // n1=1
    let n2_cos1 = n * cos_aoi;
    let n1_cos2 = cos_refr; // n1=1
    let n2_cos2 = n * cos_refr;

    let rho_s = ((n1_cos1 - n2_cos2) / (n1_cos1 + n2_cos2)).powi(2);
    let rho_p = ((n1_cos2 - n2_cos1) / (n1_cos2 + n2_cos1)).powi(2);
    let rho_0 = ((1.0 - n) / (1.0 + n)).powi(2);

    // Transmittance through the interface
    let tau_s = 1.0 - rho_s;
    let tau_p = 1.0 - rho_p;
    let tau_0 = 1.0 - rho_0;

    // Absorption through the glass
    let tau_s = tau_s * (-k * l / cos_refr).exp();
    let tau_p = tau_p * (-k * l / cos_refr).exp();
    let tau_0 = tau_0 * (-k * l).exp();

    // IAM = average of s and p, normalized by normal incidence. Clamp
    // to [0, 1] so tiny floating-point drift at grazing angles cannot
    // produce an IAM > 1 or < 0 (matches `pvlib.iam.physical`).
    ((tau_s + tau_p) / (2.0 * tau_0)).clamp(0.0, 1.0)
}

/// Schlick approximation of Fresnel reflection as IAM.
///
/// `iam = 1 - (1 - cos(aoi))^5`
///
/// # Arguments
/// * `aoi` - Angle of incidence in degrees.
///
/// # Returns
/// IAM multiplier (dimensionless).
#[inline]
pub fn schlick(aoi: f64) -> f64 {
    if aoi.abs() >= 90.0 {
        return 0.0;
    }
    let cos_aoi = aoi.to_radians().cos();
    1.0 - (1.0 - cos_aoi).powi(5)
}

/// Diffuse IAM using the Schlick model (analytically integrated).
///
/// Returns (iam_sky, iam_ground) for isotropic diffuse irradiance
/// on a tilted surface.
///
/// # Arguments
/// * `surface_tilt` - Surface tilt angle in degrees from horizontal.
///
/// # Returns
/// Tuple of (iam_sky, iam_ground).
#[inline]
pub fn schlick_diffuse(surface_tilt: f64) -> (f64, f64) {
    let cos_b = surface_tilt.to_radians().cos();
    let sin_b = surface_tilt.to_radians().sin();
    let beta_rad = surface_tilt.to_radians();
    let pi = std::f64::consts::PI;

    // Eq 4 from Xie et al. (2022): relative transmittance of sky diffuse
    let cuk = (2.0 / (pi * (1.0 + cos_b)))
        * ((30.0 / 7.0) * pi - (160.0 / 21.0) * beta_rad - (10.0 / 3.0) * pi * cos_b
            + (160.0 / 21.0) * cos_b * sin_b
            - (5.0 / 3.0) * pi * cos_b * sin_b.powi(2)
            + (20.0 / 7.0) * cos_b * sin_b.powi(3)
            - (5.0 / 16.0) * pi * cos_b * sin_b.powi(4)
            + (16.0 / 105.0) * cos_b * sin_b.powi(5));

    // Eq 6: relative transmittance of ground-reflected radiation
    let cug = if surface_tilt < 1e-6 {
        0.0
    } else {
        40.0 / (21.0 * (1.0 - cos_b)) - (1.0 + cos_b) / (1.0 - cos_b) * cuk
    };

    (cuk, cug)
}

/// Diffuse IAM using the Martin and Ruiz model.
///
/// Returns (iam_sky, iam_ground) for diffuse irradiance.
///
/// # Arguments
/// * `surface_tilt` - Surface tilt angle in degrees from horizontal.
/// * `a_r` - Angular losses coefficient (default 0.16).
///
/// # Returns
/// Tuple of (iam_sky, iam_ground).
#[inline]
pub fn martin_ruiz_diffuse(surface_tilt: f64, a_r: f64) -> (f64, f64) {
    let pi = std::f64::consts::PI;

    // Avoid undefined results for horizontal or upside-down surfaces
    let tilt = if surface_tilt == 0.0 {
        1e-6
    } else if surface_tilt == 180.0 {
        180.0 - 1e-6
    } else {
        surface_tilt
    };

    let c1 = 4.0 / 3.0 / pi; // 0.4244
    let c2 = 0.5 * a_r - 0.154;

    let beta = tilt.to_radians();
    let cos_beta = beta.cos();
    let sin_beta = if tilt < 90.0 {
        beta.sin()
    } else {
        (pi - beta).sin()
    };

    let trig_sky = sin_beta + (pi - beta - sin_beta) / (1.0 + cos_beta);
    let trig_gnd = sin_beta + (beta - sin_beta) / (1.0 - cos_beta);

    let iam_sky = 1.0 - (-(c1 + c2 * trig_sky) * trig_sky / a_r).exp();
    let iam_gnd = 1.0 - (-(c1 + c2 * trig_gnd) * trig_gnd / a_r).exp();

    (iam_sky, iam_gnd)
}

/// Sandia Array Performance Model (SAPM) Incidence Angle Modifier.
/// 
/// # References
/// King, D. et al., 2004, "Sandia Photovoltaic Array Performance Model".
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn sapm(aoi: f64, b0: f64, b1: f64, b2: f64, b3: f64, b4: f64, b5: f64) -> f64 {
    if !(0.0..=90.0).contains(&aoi) { return 0.0; }
    
    let iam = b0 
            + b1 * aoi 
            + b2 * aoi.powi(2) 
            + b3 * aoi.powi(3) 
            + b4 * aoi.powi(4) 
            + b5 * aoi.powi(5);
            
    iam.clamp(0.0, 1.0)
}
