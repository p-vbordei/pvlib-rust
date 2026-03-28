use rayon::prelude::*;
use chrono::TimeZone;
use crate::{solarposition, atmosphere, clearsky, irradiance, temperature, iam, inverter};

// ---------------------------------------------------------------------------
// Solar Position Batch
// ---------------------------------------------------------------------------

/// Batch solar position calculation for multiple timestamps.
/// Returns (zenith_vec, azimuth_vec, elevation_vec).
pub fn solar_position_batch(
    location: &crate::location::Location,
    times: &[chrono::DateTime<chrono_tz::Tz>],
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), spa::SpaError> {
    let results: Result<Vec<_>, _> = times.par_iter()
        .map(|t| solarposition::get_solarposition(location, *t))
        .collect();
    let results = results?;
    let zenith = results.iter().map(|r| r.zenith).collect();
    let azimuth = results.iter().map(|r| r.azimuth).collect();
    let elevation = results.iter().map(|r| r.elevation).collect();
    Ok((zenith, azimuth, elevation))
}

/// Convenience batch solar position for UTC `NaiveDateTime` timestamps.
///
/// Internally creates a `Location` with the UTC timezone, converts each
/// `NaiveDateTime` to `DateTime<Tz>`, and delegates to [`solar_position_batch`].
pub fn solar_position_batch_utc(
    latitude: f64,
    longitude: f64,
    altitude: f64,
    times: &[chrono::NaiveDateTime],
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), spa::SpaError> {
    let location = crate::location::Location::new(latitude, longitude, chrono_tz::UTC, altitude, "UTC");
    let datetimes: Vec<chrono::DateTime<chrono_tz::Tz>> = times
        .iter()
        .map(|ndt| chrono::Utc.from_utc_datetime(ndt).with_timezone(&chrono_tz::UTC))
        .collect();
    solar_position_batch(&location, &datetimes)
}

// ---------------------------------------------------------------------------
// Atmosphere Batch
// ---------------------------------------------------------------------------

/// Batch relative airmass for an array of zenith angles.
pub fn airmass_relative_batch(zenith: &[f64]) -> Vec<f64> {
    zenith.par_iter()
        .map(|z| atmosphere::get_relative_airmass(*z))
        .collect()
}

/// Batch absolute airmass.
pub fn airmass_absolute_batch(airmass_relative: &[f64], pressure: f64) -> Vec<f64> {
    airmass_relative.par_iter()
        .map(|am| atmosphere::get_absolute_airmass(*am, pressure))
        .collect()
}

// ---------------------------------------------------------------------------
// Clear Sky Batch
// ---------------------------------------------------------------------------

/// Batch Ineichen clear sky model.
/// Returns (ghi_vec, dni_vec, dhi_vec).
pub fn ineichen_batch(
    zenith: &[f64],
    airmass_absolute: &[f64],
    linke_turbidity: f64,
    altitude: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    assert_eq!(zenith.len(), airmass_absolute.len(), "zenith and airmass_absolute must have the same length");
    let results: Vec<_> = zenith.par_iter()
        .zip(airmass_absolute.par_iter())
        .map(|(z, am)| clearsky::ineichen(*z, *am, linke_turbidity, altitude, 1364.0))
        .collect();
    let ghi = results.iter().map(|r| r.ghi).collect();
    let dni = results.iter().map(|r| r.dni).collect();
    let dhi = results.iter().map(|r| r.dhi).collect();
    (ghi, dni, dhi)
}

/// Batch Bird clear sky model.
pub fn bird_batch(
    zenith: &[f64],
    airmass_relative: &[f64],
    aod380: f64,
    aod500: f64,
    precipitable_water: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    assert_eq!(zenith.len(), airmass_relative.len(), "zenith and airmass_relative must have the same length");
    let results: Vec<_> = zenith.par_iter()
        .zip(airmass_relative.par_iter())
        .map(|(z, am)| clearsky::bird_default(*z, *am, aod380, aod500, precipitable_water))
        .collect();
    let ghi = results.iter().map(|r| r.ghi).collect();
    let dni = results.iter().map(|r| r.dni).collect();
    let dhi = results.iter().map(|r| r.dhi).collect();
    (ghi, dni, dhi)
}

// ---------------------------------------------------------------------------
// Irradiance Batch
// ---------------------------------------------------------------------------

/// Batch AOI calculation.
pub fn aoi_batch(
    surface_tilt: f64,
    surface_azimuth: f64,
    solar_zenith: &[f64],
    solar_azimuth: &[f64],
) -> Vec<f64> {
    assert_eq!(solar_zenith.len(), solar_azimuth.len(), "solar_zenith and solar_azimuth must have the same length");
    solar_zenith.par_iter()
        .zip(solar_azimuth.par_iter())
        .map(|(z, a)| irradiance::aoi(surface_tilt, surface_azimuth, *z, *a))
        .collect()
}

/// Batch extraterrestrial radiation.
pub fn extra_radiation_batch(day_of_year: &[i32]) -> Vec<f64> {
    day_of_year.par_iter()
        .map(|d| irradiance::get_extra_radiation(*d))
        .collect()
}

/// Batch Erbs decomposition. Returns (dni_vec, dhi_vec).
pub fn erbs_batch(
    ghi: &[f64],
    zenith: &[f64],
    day_of_year: &[u32],
    dni_extra: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let n = ghi.len();
    assert_eq!(zenith.len(), n, "zenith len mismatch");
    assert_eq!(day_of_year.len(), n, "day_of_year len mismatch");
    assert_eq!(dni_extra.len(), n, "dni_extra len mismatch");
    let results: Vec<_> = ghi.par_iter()
        .zip(zenith.par_iter())
        .zip(day_of_year.par_iter())
        .zip(dni_extra.par_iter())
        .map(|(((g, z), d), e)| irradiance::erbs(*g, *z, *d, *e))
        .collect();
    let dni = results.iter().map(|r| r.0).collect();
    let dhi = results.iter().map(|r| r.1).collect();
    (dni, dhi)
}

/// Batch DISC decomposition. Returns (dni_vec, kt_vec, airmass_vec).
pub fn disc_batch(
    ghi: &[f64],
    solar_zenith: &[f64],
    day_of_year: &[i32],
    pressure: Option<f64>,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = ghi.len();
    assert_eq!(solar_zenith.len(), n, "solar_zenith len mismatch");
    assert_eq!(day_of_year.len(), n, "day_of_year len mismatch");
    let results: Vec<_> = ghi.par_iter()
        .zip(solar_zenith.par_iter())
        .zip(day_of_year.par_iter())
        .map(|((g, z), d)| irradiance::disc(*g, *z, *d, pressure))
        .collect();
    let dni = results.iter().map(|r| r.dni).collect();
    let kt = results.iter().map(|r| r.kt).collect();
    let am = results.iter().map(|r| r.airmass).collect();
    (dni, kt, am)
}

/// Batch Perez transposition model.
#[allow(clippy::too_many_arguments)]
pub fn perez_batch(
    surface_tilt: f64,
    surface_azimuth: f64,
    dhi: &[f64],
    dni: &[f64],
    dni_extra: &[f64],
    solar_zenith: &[f64],
    solar_azimuth: &[f64],
    airmass: &[f64],
    aoi_vals: &[f64],
) -> Vec<f64> {
    let n = dhi.len();
    assert_eq!(dni.len(), n, "dni len mismatch");
    assert_eq!(dni_extra.len(), n, "dni_extra len mismatch");
    assert_eq!(solar_zenith.len(), n, "solar_zenith len mismatch");
    assert_eq!(solar_azimuth.len(), n, "solar_azimuth len mismatch");
    assert_eq!(airmass.len(), n, "airmass len mismatch");
    assert_eq!(aoi_vals.len(), n, "aoi_vals len mismatch");
    (0..n).into_par_iter()
        .map(|i| irradiance::perez(
            surface_tilt, surface_azimuth,
            dhi[i], dni[i], dni_extra[i],
            solar_zenith[i], solar_azimuth[i],
            airmass[i], aoi_vals[i],
        ))
        .collect()
}

/// Batch get_total_irradiance. Returns PoaComponents for each timestep.
#[allow(clippy::too_many_arguments)]
pub fn total_irradiance_batch(
    surface_tilt: f64,
    surface_azimuth: f64,
    solar_zenith: &[f64],
    solar_azimuth: &[f64],
    dni: &[f64],
    ghi: &[f64],
    dhi: &[f64],
    albedo: f64,
    model: irradiance::DiffuseModel,
    dni_extra: &[f64],
    airmass: &[f64],
) -> Vec<irradiance::PoaComponents> {
    let n = solar_zenith.len();
    assert_eq!(solar_azimuth.len(), n, "solar_azimuth len mismatch");
    assert_eq!(dni.len(), n, "dni len mismatch");
    assert_eq!(ghi.len(), n, "ghi len mismatch");
    assert_eq!(dhi.len(), n, "dhi len mismatch");
    assert_eq!(dni_extra.len(), n, "dni_extra len mismatch");
    assert_eq!(airmass.len(), n, "airmass len mismatch");
    (0..n).into_par_iter()
        .map(|i| irradiance::get_total_irradiance(
            surface_tilt, surface_azimuth,
            solar_zenith[i], solar_azimuth[i],
            dni[i], ghi[i], dhi[i],
            albedo, model,
            Some(dni_extra[i]),
            Some(airmass[i]),
        ))
        .collect()
}

// ---------------------------------------------------------------------------
// Temperature Batch
// ---------------------------------------------------------------------------

/// Batch SAPM cell temperature. Returns (cell_temp_vec, module_temp_vec).
#[allow(clippy::too_many_arguments)]
pub fn sapm_cell_temperature_batch(
    poa_global: &[f64],
    temp_air: &[f64],
    wind_speed: &[f64],
    a: f64, b: f64, delta_t: f64, irrad_ref: f64,
) -> (Vec<f64>, Vec<f64>) {
    let n = poa_global.len();
    assert_eq!(temp_air.len(), n, "temp_air len mismatch");
    assert_eq!(wind_speed.len(), n, "wind_speed len mismatch");
    let results: Vec<_> = (0..n).into_par_iter()
        .map(|i| temperature::sapm_cell_temperature(
            poa_global[i], temp_air[i], wind_speed[i], a, b, delta_t, irrad_ref,
        ))
        .collect();
    let cell = results.iter().map(|r| r.0).collect();
    let module = results.iter().map(|r| r.1).collect();
    (cell, module)
}

/// Batch Faiman cell temperature.
pub fn faiman_batch(
    poa_global: &[f64],
    temp_air: &[f64],
    wind_speed: &[f64],
    u0: f64,
    u1: f64,
) -> Vec<f64> {
    let n = poa_global.len();
    assert_eq!(temp_air.len(), n, "temp_air len mismatch");
    assert_eq!(wind_speed.len(), n, "wind_speed len mismatch");
    (0..n).into_par_iter()
        .map(|i| temperature::faiman(poa_global[i], temp_air[i], wind_speed[i], u0, u1))
        .collect()
}

// ---------------------------------------------------------------------------
// IAM Batch
// ---------------------------------------------------------------------------

/// Batch Physical IAM.
pub fn iam_physical_batch(aoi: &[f64], n: f64, k: f64, l: f64) -> Vec<f64> {
    aoi.par_iter()
        .map(|a| iam::physical(*a, n, k, l))
        .collect()
}

/// Batch ASHRAE IAM.
pub fn iam_ashrae_batch(aoi: &[f64], b0: f64) -> Vec<f64> {
    aoi.par_iter()
        .map(|a| iam::ashrae(*a, b0))
        .collect()
}

// ---------------------------------------------------------------------------
// Inverter Batch
// ---------------------------------------------------------------------------

/// Batch PVWatts AC power.
pub fn pvwatts_ac_batch(
    pdc: &[f64],
    pdc0: f64,
    eta_inv_nom: f64,
    eta_inv_ref: f64,
) -> Vec<f64> {
    pdc.par_iter()
        .map(|p| inverter::pvwatts_ac(*p, pdc0, eta_inv_nom, eta_inv_ref))
        .collect()
}

// ---------------------------------------------------------------------------
// Full Pipeline Batch
// ---------------------------------------------------------------------------

/// Input data for batch simulation -- one value per timestep.
#[derive(Debug, Clone)]
pub struct WeatherSeries {
    pub times: Vec<chrono::DateTime<chrono_tz::Tz>>,
    pub ghi: Vec<f64>,
    pub dni: Vec<f64>,
    pub dhi: Vec<f64>,
    pub temp_air: Vec<f64>,
    pub wind_speed: Vec<f64>,
    pub albedo: Option<Vec<f64>>,
}

impl WeatherSeries {
    /// Construct a `WeatherSeries` from UTC `NaiveDateTime` timestamps and a
    /// timezone name (e.g. `"US/Eastern"`, `"UTC"`, `"Europe/Berlin"`).
    ///
    /// Each `NaiveDateTime` is interpreted as UTC and converted to the target
    /// timezone.  Returns `Err` if `tz_name` cannot be parsed.
    pub fn from_utc(
        timestamps: &[chrono::NaiveDateTime],
        tz_name: &str,
        ghi: Vec<f64>,
        dni: Vec<f64>,
        dhi: Vec<f64>,
        temp_air: Vec<f64>,
        wind_speed: Vec<f64>,
    ) -> Result<Self, String> {
        let tz: chrono_tz::Tz = tz_name
            .parse()
            .map_err(|_| format!("Unknown timezone: {}", tz_name))?;
        let times: Vec<chrono::DateTime<chrono_tz::Tz>> = timestamps
            .iter()
            .map(|ndt| chrono::Utc.from_utc_datetime(ndt).with_timezone(&tz))
            .collect();
        Ok(Self {
            times,
            ghi,
            dni,
            dhi,
            temp_air,
            wind_speed,
            albedo: None,
        })
    }
}

/// Output from batch simulation -- one value per timestep.
#[derive(Debug, Clone)]
pub struct SimulationSeries {
    pub solar_zenith: Vec<f64>,
    pub solar_elevation: Vec<f64>,
    pub solar_azimuth: Vec<f64>,
    pub airmass: Vec<f64>,
    pub aoi: Vec<f64>,
    pub poa_global: Vec<f64>,
    pub poa_direct: Vec<f64>,
    pub poa_diffuse: Vec<f64>,
    pub cell_temperature: Vec<f64>,
    pub effective_irradiance: Vec<f64>,
    pub dc_power: Vec<f64>,
    pub ac_power: Vec<f64>,
}

impl SimulationSeries {
    /// Total energy produced in Wh (assuming 1-hour timesteps).
    pub fn total_energy_wh(&self) -> f64 {
        self.ac_power.iter().filter(|p| **p > 0.0).sum()
    }

    /// Peak AC power in W.
    pub fn peak_power(&self) -> f64 {
        self.ac_power.iter().cloned().fold(0.0_f64, f64::max)
    }

    /// Capacity factor (ratio of actual energy to theoretical maximum).
    pub fn capacity_factor(&self, system_capacity_w: f64) -> f64 {
        let hours = self.ac_power.len() as f64;
        if hours == 0.0 || system_capacity_w == 0.0 { return 0.0; }
        self.total_energy_wh() / (system_capacity_w * hours)
    }
}

/// Batch ModelChain -- runs the full PV simulation pipeline on a time series
/// using rayon for parallel processing.
///
/// This is the main entry point for production batch simulations.
/// A typical TMY year (8760 hourly timesteps) completes in milliseconds.
pub struct BatchModelChain {
    pub location: crate::location::Location,
    pub surface_tilt: f64,
    pub surface_azimuth: f64,
    pub system_capacity_dc: f64,
    /// Temperature coefficient, e.g. -0.004
    pub gamma_pdc: f64,
    pub inverter_capacity: f64,
    pub inverter_efficiency: f64,
    pub albedo: f64,
    pub transposition_model: irradiance::DiffuseModel,
    /// When true, automatically decompose GHI into DNI/DHI using the Erbs model
    /// if DNI and DHI are both near zero but GHI is positive.
    pub auto_decomposition: bool,
    /// Bifaciality factor: 0.0 = monofacial, 0.65-0.85 typical for bifacial modules.
    pub bifaciality_factor: f64,
    /// Ground albedo used for rear-side irradiance calculation.
    pub bifacial_ground_albedo: f64,
}

impl BatchModelChain {
    /// Create a new BatchModelChain with PVWatts-style defaults.
    pub fn pvwatts(
        location: crate::location::Location,
        surface_tilt: f64,
        surface_azimuth: f64,
        system_capacity_dc: f64,
    ) -> Self {
        Self {
            location,
            surface_tilt,
            surface_azimuth,
            system_capacity_dc,
            gamma_pdc: -0.004,
            inverter_capacity: system_capacity_dc,
            inverter_efficiency: 0.96,
            albedo: 0.2,
            transposition_model: irradiance::DiffuseModel::Perez,
            auto_decomposition: false,
            bifaciality_factor: 0.0,
            bifacial_ground_albedo: 0.2,
        }
    }

    /// Builder: set temperature coefficient.
    pub fn with_gamma_pdc(mut self, gamma_pdc: f64) -> Self {
        self.gamma_pdc = gamma_pdc;
        self
    }

    /// Builder: set inverter parameters.
    pub fn with_inverter(mut self, capacity: f64, efficiency: f64) -> Self {
        self.inverter_capacity = capacity;
        self.inverter_efficiency = efficiency;
        self
    }

    /// Builder: set ground albedo.
    pub fn with_albedo(mut self, albedo: f64) -> Self {
        self.albedo = albedo;
        self
    }

    /// Builder: set transposition model.
    pub fn with_transposition(mut self, model: irradiance::DiffuseModel) -> Self {
        self.transposition_model = model;
        self
    }

    /// Builder: enable/disable automatic GHI to DNI/DHI decomposition via Erbs model.
    pub fn with_auto_decomposition(mut self, enabled: bool) -> Self {
        self.auto_decomposition = enabled;
        self
    }

    /// Builder: set bifacial parameters.
    ///
    /// `bifaciality_factor` is typically 0.65-0.85 for bifacial modules (0.0 = monofacial).
    /// `ground_albedo` is the albedo used for rear-side irradiance calculation.
    pub fn with_bifacial(mut self, bifaciality_factor: f64, ground_albedo: f64) -> Self {
        self.bifaciality_factor = bifaciality_factor;
        self.bifacial_ground_albedo = ground_albedo;
        self
    }

    /// Run the full simulation on a weather time series.
    ///
    /// Uses rayon for automatic parallelization across CPU cores.
    pub fn run(&self, weather: &WeatherSeries) -> Result<SimulationSeries, spa::SpaError> {
        let n = weather.times.len();
        assert_eq!(weather.ghi.len(), n);
        assert_eq!(weather.dni.len(), n, "dni len mismatch");
        assert_eq!(weather.dhi.len(), n, "dhi len mismatch");
        assert_eq!(weather.temp_air.len(), n, "temp_air len mismatch");
        assert_eq!(weather.wind_speed.len(), n, "wind_speed len mismatch");
        if let Some(albedos) = &weather.albedo {
            assert_eq!(albedos.len(), n, "albedo len mismatch");
        }

        let pressure = atmosphere::alt2pres(self.location.altitude);
        let albedo_default = self.albedo;

        // Run entire pipeline in parallel for each timestep
        let results: Result<Vec<_>, _> = (0..n).into_par_iter().map(|i| {
            // 1. Solar position
            let solpos = solarposition::get_solarposition(&self.location, weather.times[i])?;

            // 2. Airmass
            let am_rel = atmosphere::get_relative_airmass(solpos.zenith);
            let am_abs = if am_rel.is_nan() || am_rel <= 0.0 {
                0.0
            } else {
                atmosphere::get_absolute_airmass(am_rel, pressure)
            };

            // 3. AOI
            let aoi_val = irradiance::aoi(
                self.surface_tilt, self.surface_azimuth,
                solpos.zenith, solpos.azimuth,
            );

            // 4. Extraterrestrial irradiance
            let doy = weather.times[i].format("%j").to_string().parse::<i32>().unwrap_or(1);
            let dni_extra = irradiance::get_extra_radiation(doy);

            // 4b. Auto-decompose GHI → DNI/DHI if enabled and needed
            let (dni_i, dhi_i) = if self.auto_decomposition
                && weather.dni[i].abs() < 1.0
                && weather.dhi[i].abs() < 1.0
                && weather.ghi[i] > 0.0
            {
                let dni_extra_val = dni_extra;
                irradiance::erbs(weather.ghi[i], solpos.zenith, doy as u32, dni_extra_val)
            } else {
                (weather.dni[i], weather.dhi[i])
            };

            // 5. Transposition
            let poa = irradiance::get_total_irradiance(
                self.surface_tilt, self.surface_azimuth,
                solpos.zenith, solpos.azimuth,
                dni_i, weather.ghi[i], dhi_i,
                weather.albedo.as_ref().map_or(albedo_default, |a| a[i]),
                self.transposition_model,
                Some(dni_extra),
                if am_rel.is_nan() { None } else { Some(am_rel) },
            );

            // 6. IAM
            // Note: Uses Physical IAM model (Fresnel/Snell's law) matching the PVWatts configuration.
            let iam_val = iam::physical(aoi_val, 1.526, 4.0, 0.002);

            // 7. Effective irradiance (spectral modifier = 1.0 for NoLoss)
            let spectral_modifier = 1.0;
            let eff_irrad = ((poa.poa_direct * iam_val + poa.poa_diffuse) * spectral_modifier).max(0.0);

            // 8. Cell temperature
            // Note: Uses PVWatts temperature model (T_cell = T_air + POA*(NOCT-20)/800).
            let t_cell = weather.temp_air[i] + poa.poa_global * (45.0 - 20.0) / 800.0;

            // 9. DC power
            let pdc = self.system_capacity_dc * (eff_irrad / 1000.0)
                * (1.0 + self.gamma_pdc * (t_cell - 25.0));
            let pdc = pdc.max(0.0);

            // 10. AC power
            let pac = inverter::pvwatts_ac(
                pdc, self.system_capacity_dc,
                self.inverter_efficiency, 0.9637,
            );

            // 11. Bifacial rear-side gain
            let pac = if self.bifaciality_factor > 0.0 && poa.poa_global > 10.0 {
                let rear_gain = (self.bifaciality_factor * self.bifacial_ground_albedo
                    * weather.ghi[i] / poa.poa_global).min(0.25);
                pac * (1.0 + rear_gain)
            } else {
                pac
            };

            Ok((solpos.zenith, solpos.azimuth, am_abs, aoi_val,
                poa.poa_global, poa.poa_direct, poa.poa_diffuse,
                t_cell, eff_irrad, pdc, pac))
        }).collect();

        let results = results?;

        Ok(SimulationSeries {
            solar_zenith: results.iter().map(|r| r.0).collect(),
            solar_elevation: results.iter().map(|r| 90.0 - r.0).collect(),
            solar_azimuth: results.iter().map(|r| r.1).collect(),
            airmass: results.iter().map(|r| r.2).collect(),
            aoi: results.iter().map(|r| r.3).collect(),
            poa_global: results.iter().map(|r| r.4).collect(),
            poa_direct: results.iter().map(|r| r.5).collect(),
            poa_diffuse: results.iter().map(|r| r.6).collect(),
            cell_temperature: results.iter().map(|r| r.7).collect(),
            effective_irradiance: results.iter().map(|r| r.8).collect(),
            dc_power: results.iter().map(|r| r.9).collect(),
            ac_power: results.iter().map(|r| r.10).collect(),
        })
    }
}
