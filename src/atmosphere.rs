/// Calculate the relative airmass using the Kasten-Young 1989 model.
/// 
/// # Arguments
/// * `zenith` - True or apparent zenith angle in degrees.
pub fn get_relative_airmass(zenith: f64) -> f64 {
    // Avoid calculating for zenith severely > 90 to prevent math errors.
    // However, the formula can technically handle slightly above 90 before going negative.
    let z = zenith;
    let cos_z = z.to_radians().cos();
    
    // Kasten-Young 1989 formula
    let c = 96.07995 - z;
    if c <= 0.0 {
        return f64::NAN; // Airmass is undefined for zenith angles >= 96.07995
    }
    
    let term = 0.50572 * c.powf(-1.6364);
    let am_rel = 1.0 / (cos_z + term);

    // Limit AM to reasonable numbers for horizon
    if !(0.0..=100.0).contains(&am_rel) {
        f64::NAN
    } else {
        am_rel
    }
}

/// Calculate atmospheric pressure (Pascals) from altitude (meters)
/// using the standard atmosphere model.
pub fn alt2pres(altitude: f64) -> f64 {
    101325.0 * (1.0 - 2.25577e-5 * altitude).powf(5.25588)
}

/// Calculate absolute airmass given relative airmass and pressure (Pascals).
pub fn get_absolute_airmass(airmass_relative: f64, pressure: f64) -> f64 {
    airmass_relative * pressure / 101325.0
}

/// Determine altitude from site pressure.
///
/// # Arguments
/// * `pressure` - Atmospheric pressure in Pascals.
///
/// # Returns
/// Altitude above sea level in meters.
pub fn pres2alt(pressure: f64) -> f64 {
    44331.5 - 4946.62 * pressure.powf(0.190263)
}

/// Calculate precipitable water (cm) from ambient air temperature (C)
/// and relative humidity (%) using Gueymard (1994) model.
///
/// # Arguments
/// * `temp_air` - Ambient air temperature in degrees C.
/// * `relative_humidity` - Relative humidity in percent (0-100).
///
/// # Returns
/// Precipitable water in cm (minimum 0.1).
pub fn gueymard94_pw(temp_air: f64, relative_humidity: f64) -> f64 {
    let t = temp_air + 273.15; // Kelvin
    let rh = relative_humidity;
    let theta = t / 273.15;

    let hv = 0.4976 + 1.5265 * theta
        + (13.6897 * theta - 14.9188 * theta.powi(3)).exp();

    let es = (22.330 - 49.140 * (100.0 / t) - 10.922 * (100.0 / t).powi(2)
        - 0.39015 * t / 100.0)
        .exp();

    let rho_v = 216.7 * rh / (100.0 * t) * es;

    let pw = 0.1 * hv * rho_v;
    pw.max(0.1)
}

/// Calculate dew point temperature from air temperature and relative humidity
/// using the Magnus formula.
///
/// # Arguments
/// * `temperature` - Air temperature in degrees C.
/// * `rh` - Relative humidity in percent (0-100).
///
/// # Returns
/// Dew point temperature in degrees C.
pub fn tdew_from_rh(temperature: f64, rh: f64) -> f64 {
    const B: f64 = 17.62;
    const C: f64 = 243.12;

    let ln_term = (B * temperature) / (C + temperature) + (rh / 100.0).ln();
    C * ln_term / (B - ln_term)
}

/// Calculate relative humidity from air temperature and dew point temperature
/// using the Magnus formula.
///
/// # Arguments
/// * `temperature` - Air temperature in degrees C.
/// * `temp_dew` - Dew point temperature in degrees C.
///
/// # Returns
/// Relative humidity in percent (0-100).
pub fn rh_from_tdew(temperature: f64, temp_dew: f64) -> f64 {
    const B: f64 = 17.62;
    const C: f64 = 243.12;

    let es = (B * temp_dew) / (C + temp_dew);
    let e = (B * temperature) / (C + temperature);
    100.0 * (es - e).exp()
}

/// Approximate broadband aerosol optical depth using Bird and Hulstrom (1980).
///
/// # Arguments
/// * `aod380` - AOD measured at 380 nm.
/// * `aod500` - AOD measured at 500 nm.
///
/// # Returns
/// Broadband AOD.
pub fn bird_hulstrom80_aod_bb(aod380: f64, aod500: f64) -> f64 {
    0.27583 * aod380 + 0.35 * aod500
}

/// Calculate Linke turbidity using Kasten pyrheliometric formula (1996).
///
/// # Arguments
/// * `airmass` - Pressure-adjusted (absolute) airmass.
/// * `precipitable_water` - Precipitable water in cm.
/// * `aod_bb` - Broadband aerosol optical depth.
///
/// # Returns
/// Linke turbidity factor.
pub fn kasten96_lt(airmass: f64, precipitable_water: f64, aod_bb: f64) -> f64 {
    let delta_cda = -0.101 + 0.235 * airmass.powf(-0.16);
    let delta_w = 0.112 * airmass.powf(-0.55) * precipitable_water.powf(0.34);
    let delta_a = aod_bb;

    -(9.4 + 0.9 * airmass)
        * (-airmass * (delta_cda + delta_w + delta_a)).exp().ln()
        / airmass
}

/// Get AOD at specified wavelength using Angstrom turbidity model.
///
/// # Arguments
/// * `aod0` - AOD measured at wavelength `lambda0`.
/// * `lambda0` - Wavelength corresponding to `aod0` in nm.
/// * `alpha` - Angstrom alpha exponent.
/// * `lambda1` - Desired wavelength in nm.
///
/// # Returns
/// AOD at desired wavelength.
pub fn angstrom_aod_at_lambda(aod0: f64, lambda0: f64, alpha: f64, lambda1: f64) -> f64 {
    aod0 * (lambda1 / lambda0).powf(-alpha)
}

/// Calculate Angstrom alpha exponent from two wavelength-AOD pairs.
///
/// # Arguments
/// * `aod1` - AOD at wavelength `lambda1`.
/// * `lambda1` - First wavelength in nm.
/// * `aod2` - AOD at wavelength `lambda2`.
/// * `lambda2` - Second wavelength in nm.
///
/// # Returns
/// Angstrom alpha exponent.
pub fn angstrom_alpha(aod1: f64, lambda1: f64, aod2: f64, lambda2: f64) -> f64 {
    -(aod1 / aod2).ln() / (lambda1 / lambda2).ln()
}

/// Estimate wind speed at a different height using the power law (Hellmann).
///
/// # Arguments
/// * `wind_speed` - Measured wind speed in m/s.
/// * `height_ref` - Reference height in meters.
/// * `height_desired` - Desired height in meters.
/// * `alpha` - Power law exponent (default 1/7 ≈ 0.143).
///
/// # Returns
/// Wind speed at desired height in m/s.
pub fn windspeed_powerlaw(wind_speed: f64, height_ref: f64, height_desired: f64, alpha: f64) -> f64 {
    if wind_speed < 0.0 || height_ref <= 0.0 || height_desired <= 0.0 {
        return f64::NAN;
    }
    wind_speed * (height_desired / height_ref).powf(alpha)
}

/// Atmospheric refraction calculation.
///
/// Converts true elevation to apparent elevation.
/// 
/// # Arguments
/// * `true_elevation` - True solar elevation in degrees.
/// 
/// # Returns
/// Atmospheric refraction in degrees to add to true elevation.
pub fn get_refraction(true_elevation: f64) -> f64 {
    let e = true_elevation;
    if e >= 85.0 {
        0.0
    } else if e > -0.575 {
        // Simple approximate refraction formula (Saemundsson, 1986)
        let r = 1.02 / (e + 10.3 / (e + 5.11)).to_radians().tan();
        r / 60.0 // Convert arcminutes to degrees
    } else {
        // Below horizon
        34.5 / 60.0 
    }
}

