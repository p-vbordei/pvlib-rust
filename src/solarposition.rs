use chrono::{DateTime, Utc};
use chrono_tz::Tz;
use spa::{solar_position, SpaError, StdFloatOps};
use crate::location::Location;

use std::f64::consts::PI;

/// Output of the solar position calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct SolarPosition {
    /// True zenith angle in degrees
    pub zenith: f64,
    /// True azimuth angle in degrees (North = 0, East = 90)
    pub azimuth: f64,
    /// True elevation angle in degrees
    pub elevation: f64,
}

/// Calculate the solar position for a given location and time.
/// This uses the NREL SPA (Solar Position Algorithm).
#[inline]
pub fn get_solarposition(location: &Location, time: DateTime<Tz>) -> Result<SolarPosition, SpaError> {
    let utc_time: DateTime<Utc> = time.with_timezone(&Utc);

    let result = solar_position::<StdFloatOps>(utc_time, location.latitude, location.longitude)?;

    Ok(SolarPosition {
        zenith: result.zenith_angle,
        azimuth: result.azimuth,
        elevation: 90.0 - result.zenith_angle,
    })
}

// --- Helper ---

/// Calculates the day angle for the Earth's orbit around the Sun.
/// For the Spencer method, offset=1.
fn calculate_simple_day_angle(day_of_year: f64, offset: f64) -> f64 {
    (2.0 * PI / 365.0) * (day_of_year - offset)
}

// --- Equation of Time ---

/// Equation of time in minutes using Spencer (1971) / Iqbal (1983) Fourier series.
///
/// Returns the difference between solar time and mean solar time in minutes.
#[inline]
pub fn equation_of_time_spencer71(day_of_year: f64) -> f64 {
    let day_angle = calculate_simple_day_angle(day_of_year, 1.0);
    (1440.0 / (2.0 * PI))
        * (0.0000075
            + 0.001868 * day_angle.cos()
            - 0.032077 * day_angle.sin()
            - 0.014615 * (2.0 * day_angle).cos()
            - 0.040849 * (2.0 * day_angle).sin())
}

/// Equation of time in minutes from PVCDROM / PV Education.
///
/// Returns the difference between solar time and mean solar time in minutes.
#[inline]
pub fn equation_of_time_pvcdrom(day_of_year: f64) -> f64 {
    let bday = calculate_simple_day_angle(day_of_year, 1.0) - (2.0 * PI / 365.0) * 80.0;
    9.87 * (2.0 * bday).sin() - 7.53 * bday.cos() - 1.5 * bday.sin()
}

// --- Declination ---

/// Solar declination in radians using Spencer (1971) Fourier series.
///
/// Returns the angular position of the sun at solar noon relative to
/// the plane of the equator, approximately +/-23.45 degrees.
#[inline]
pub fn declination_spencer71(day_of_year: f64) -> f64 {
    let day_angle = calculate_simple_day_angle(day_of_year, 1.0);
    0.006918
        - 0.399912 * day_angle.cos()
        + 0.070257 * day_angle.sin()
        - 0.006758 * (2.0 * day_angle).cos()
        + 0.000907 * (2.0 * day_angle).sin()
        - 0.002697 * (3.0 * day_angle).cos()
        + 0.00148 * (3.0 * day_angle).sin()
}

/// Solar declination in radians using Cooper (1969).
///
/// delta = 23.45 * sin(2*pi*(284+doy)/365) converted to radians.
#[inline]
pub fn declination_cooper69(day_of_year: f64) -> f64 {
    let day_angle = calculate_simple_day_angle(day_of_year, 1.0);
    (23.45_f64).to_radians() * (day_angle + (2.0 * PI / 365.0) * 285.0).sin()
}

// --- Hour Angle ---

/// Hour angle in degrees.
///
/// `hours_from_midnight_utc` is fractional hours since midnight UTC.
/// `longitude` is in degrees. `equation_of_time` is in minutes.
///
/// HA = 15 * (hours - 12) + longitude + eot/4
#[inline]
pub fn hour_angle(hours_from_midnight_utc: f64, longitude: f64, equation_of_time: f64) -> f64 {
    15.0 * (hours_from_midnight_utc - 12.0) + longitude + equation_of_time / 4.0
}

// --- Zenith and Azimuth (analytical, all inputs in radians) ---

/// Solar zenith angle in radians from analytical spherical trigonometry.
///
/// All inputs are in radians.
#[inline]
pub fn solar_zenith_analytical(latitude: f64, hour_angle: f64, declination: f64) -> f64 {
    let cos_zenith = declination.cos() * latitude.cos() * hour_angle.cos()
        + declination.sin() * latitude.sin();
    // clamp to [-1, 1] before acos to avoid NaN from floating point errors
    cos_zenith.clamp(-1.0, 1.0).acos()
}

/// Solar azimuth angle in radians from analytical spherical trigonometry.
///
/// All inputs are in radians. Returns azimuth measured from north (0) clockwise.
#[inline]
pub fn solar_azimuth_analytical(
    latitude: f64,
    hour_angle: f64,
    declination: f64,
    zenith: f64,
) -> f64 {
    let numer = zenith.cos() * latitude.sin() - declination.sin();
    let denom = zenith.sin() * latitude.cos();

    let mut cos_azi = if denom.abs() < 1e-8 {
        1.0
    } else {
        numer / denom
    };

    // clamp to [-1, 1]
    if (cos_azi - 1.0).abs() < 1e-8 {
        cos_azi = 1.0;
    }
    if (cos_azi + 1.0).abs() < 1e-8 {
        cos_azi = -1.0;
    }
    cos_azi = cos_azi.clamp(-1.0, 1.0);

    let sign_ha = if hour_angle > 0.0 {
        1.0
    } else if hour_angle < 0.0 {
        -1.0
    } else {
        0.0
    };

    sign_ha * cos_azi.acos() + PI
}

// --- Sunrise, Sunset, Transit ---

/// Result of geometric sunrise/sunset/transit calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct SunRiseSetTransit {
    /// Sunrise as fractional hours from midnight (local time)
    pub sunrise: f64,
    /// Sunset as fractional hours from midnight (local time)
    pub sunset: f64,
    /// Solar transit (noon) as fractional hours from midnight (local time)
    pub transit: f64,
}

/// Geometric calculation of sunrise, sunset, and transit times.
///
/// * `latitude` - degrees, positive north
/// * `longitude` - degrees, positive east
/// * `declination` - radians
/// * `equation_of_time` - minutes
/// * `utc_offset_hours` - timezone offset from UTC in hours (e.g. -5 for EST)
///
/// Returns sunrise, sunset, and transit as fractional hours from local midnight.
/// Returns `None` if the sun never rises or never sets (polar day/night).
#[inline]
pub fn sun_rise_set_transit_geometric(
    latitude: f64,
    longitude: f64,
    declination: f64,
    equation_of_time: f64,
    utc_offset_hours: f64,
) -> Option<SunRiseSetTransit> {
    let latitude_rad = latitude.to_radians();
    let tan_product = -declination.tan() * latitude_rad.tan();

    // If |tan_product| > 1, the sun never rises or never sets
    if tan_product.abs() > 1.0 {
        return None;
    }

    let sunset_angle = tan_product.acos().to_degrees(); // degrees
    let sunrise_angle = -sunset_angle;

    // Convert hour angles to hours from local midnight
    let sunrise_hour =
        (sunrise_angle - longitude - equation_of_time / 4.0) / 15.0 + 12.0 + utc_offset_hours;
    let sunset_hour =
        (sunset_angle - longitude - equation_of_time / 4.0) / 15.0 + 12.0 + utc_offset_hours;
    let transit_hour =
        (0.0 - longitude - equation_of_time / 4.0) / 15.0 + 12.0 + utc_offset_hours;

    Some(SunRiseSetTransit {
        sunrise: sunrise_hour,
        sunset: sunset_hour,
        transit: transit_hour,
    })
}

// --- Earth-Sun Distance ---

/// Earth-sun distance in AU using Spencer (1971).
///
/// Uses the Fourier series for the reciprocal squared of the earth-sun distance.
#[inline]
pub fn nrel_earthsun_distance(day_of_year: f64) -> f64 {
    let day_angle = calculate_simple_day_angle(day_of_year, 1.0);
    let r_over_r0_squared = 1.000110
        + 0.034221 * day_angle.cos()
        + 0.001280 * day_angle.sin()
        + 0.000719 * (2.0 * day_angle).cos()
        + 0.000077 * (2.0 * day_angle).sin();
    // This gives 1/r^2 in AU, so distance = 1/sqrt(value)
    1.0 / r_over_r0_squared.sqrt()
}
