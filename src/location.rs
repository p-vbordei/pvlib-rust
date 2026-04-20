use chrono::DateTime;
use chrono_tz::Tz;

use crate::atmosphere;
use crate::clearsky;
use crate::solarposition::{self, SolarPosition};

/// Represents a physical location on Earth.
///
/// `pvlib-python` equivalently uses `Location` class with attributes:
/// latitude, longitude, tz, altitude, and name.
#[derive(Debug, Clone, PartialEq)]
pub struct Location {
    pub latitude: f64,
    pub longitude: f64,
    pub tz: Tz,
    pub altitude: f64,
    pub name: String,
}

/// Error produced by [`Location::try_new`] when a coordinate is out of range.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum LocationError {
    #[error("latitude {0} is outside the valid range [-90, 90]")]
    Latitude(f64),
    #[error("longitude {0} is outside the valid range [-180, 180]")]
    Longitude(f64),
    #[error("altitude {0} is not a finite number")]
    Altitude(f64),
}

impl Location {
    /// Create a new Location instance.
    ///
    /// Inputs are **not** validated — pass garbage latitude or longitude and
    /// you will get garbage solar-position output. Use [`Location::try_new`]
    /// at trust boundaries (weather API input, user-supplied config) for a
    /// validated constructor.
    ///
    /// # Arguments
    ///
    /// * `latitude` - Latitude in decimal degrees. Positive north of equator, negative to south.
    /// * `longitude` - Longitude in decimal degrees. Positive east of prime meridian, negative to west.
    /// * `tz` - Timezone as a `chrono_tz::Tz` enum variant.
    /// * `altitude` - Altitude from sea level in meters.
    /// * `name` - Name of the location.
    pub fn new(latitude: f64, longitude: f64, tz: Tz, altitude: f64, name: &str) -> Self {
        Self {
            latitude,
            longitude,
            tz,
            altitude,
            name: name.to_string(),
        }
    }

    /// Create a new Location with validated coordinates.
    ///
    /// Returns `Err(LocationError)` if latitude is outside `[-90, 90]`,
    /// longitude is outside `[-180, 180]`, or altitude is non-finite.
    pub fn try_new(
        latitude: f64,
        longitude: f64,
        tz: Tz,
        altitude: f64,
        name: &str,
    ) -> Result<Self, LocationError> {
        if !(-90.0..=90.0).contains(&latitude) || !latitude.is_finite() {
            return Err(LocationError::Latitude(latitude));
        }
        if !(-180.0..=180.0).contains(&longitude) || !longitude.is_finite() {
            return Err(LocationError::Longitude(longitude));
        }
        if !altitude.is_finite() {
            return Err(LocationError::Altitude(altitude));
        }
        Ok(Self::new(latitude, longitude, tz, altitude, name))
    }

    /// Calculate the solar position for this location at the given time.
    ///
    /// Convenience wrapper around `solarposition::get_solarposition`.
    pub fn get_solarposition(&self, time: DateTime<Tz>) -> Result<SolarPosition, spa::SpaError> {
        solarposition::get_solarposition(self, time)
    }

    /// Calculate clear sky irradiance for this location at the given time.
    ///
    /// # Arguments
    /// * `time` - Date and time with timezone.
    /// * `model` - Clear sky model: "ineichen", "haurwitz", or "simplified_solis".
    ///
    /// # Returns
    /// GHI, DNI, DHI in W/m^2. Returns zeros if the sun is below the horizon.
    pub fn get_clearsky(&self, time: DateTime<Tz>, model: &str) -> clearsky::ClearSkyIrradiance {
        let solar_pos = match self.get_solarposition(time) {
            Ok(sp) => sp,
            Err(_) => return clearsky::ClearSkyIrradiance { ghi: 0.0, dni: 0.0, dhi: 0.0 },
        };

        match model {
            "haurwitz" => {
                let ghi = clearsky::haurwitz(solar_pos.zenith);
                clearsky::ClearSkyIrradiance { ghi, dni: 0.0, dhi: 0.0 }
            }
            "simplified_solis" => {
                let apparent_elevation = 90.0 - solar_pos.zenith;
                clearsky::simplified_solis(apparent_elevation, 0.1, 1.0, atmosphere::alt2pres(self.altitude))
            }
            _ => {
                // Default to ineichen
                let (_am_rel, am_abs) = self.get_airmass(time);
                if am_abs.is_nan() || am_abs <= 0.0 {
                    return clearsky::ClearSkyIrradiance { ghi: 0.0, dni: 0.0, dhi: 0.0 };
                }
                let month = {
                    use chrono::Datelike;
                    time.month()
                };
                let linke_turbidity = clearsky::lookup_linke_turbidity(self.latitude, self.longitude, month);
                clearsky::ineichen(solar_pos.zenith, am_abs, linke_turbidity, self.altitude, 1364.0)
            }
        }
    }

    /// Calculate relative and absolute airmass for this location at the given time.
    ///
    /// Uses Kasten-Young model for relative airmass and site pressure derived
    /// from altitude for absolute airmass.
    ///
    /// # Returns
    /// `(airmass_relative, airmass_absolute)`. Values may be NaN if the sun is
    /// below the horizon.
    pub fn get_airmass(&self, time: DateTime<Tz>) -> (f64, f64) {
        let solar_pos = match self.get_solarposition(time) {
            Ok(sp) => sp,
            Err(_) => return (f64::NAN, f64::NAN),
        };

        let am_rel = atmosphere::get_relative_airmass(solar_pos.zenith);
        let pressure = atmosphere::alt2pres(self.altitude);
        let am_abs = atmosphere::get_absolute_airmass(am_rel, pressure);
        (am_rel, am_abs)
    }
}

/// Lookup altitude for a given latitude and longitude.
///
/// This is a simplified approximation. Most populated areas are near sea level,
/// so this returns 0.0 as a default. For accurate altitude data, use SRTM or
/// similar elevation datasets.
///
/// # Arguments
/// * `_latitude` - Latitude in decimal degrees.
/// * `_longitude` - Longitude in decimal degrees.
///
/// # Returns
/// Estimated altitude in meters above sea level.
pub fn lookup_altitude(_latitude: f64, _longitude: f64) -> f64 {
    // A proper implementation would query SRTM or similar elevation data.
    // For now, return 0.0 (sea level) as a safe default.
    0.0
}
