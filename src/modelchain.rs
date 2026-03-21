use chrono::DateTime;
use chrono_tz::Tz;
use crate::location::Location;
use crate::pvsystem::PVSystem;
use crate::solarposition::get_solarposition;
use crate::irradiance::{
    aoi, get_total_irradiance, get_extra_radiation, poa_direct, erbs,
    DiffuseModel, PoaComponents,
};
use crate::atmosphere::{get_relative_airmass, get_absolute_airmass, alt2pres};
use crate::iam;
use crate::temperature;
use crate::inverter;

// ---------------------------------------------------------------------------
// Model selection enums
// ---------------------------------------------------------------------------

/// DC power model selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DCModel {
    /// PVWatts simple linear model (nameplate * POA/1000 * temp correction)
    PVWatts,
}

/// AC power model selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ACModel {
    /// PVWatts simple inverter model
    PVWatts,
    /// Sandia inverter model
    Sandia,
    /// ADR inverter model
    ADR,
}

/// Angle-of-incidence (IAM) model selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AOIModel {
    /// Fresnel/Snell physical model
    Physical,
    /// ASHRAE model
    ASHRAE,
    /// SAPM polynomial model
    SAPM,
    /// Martin-Ruiz model
    MartinRuiz,
    /// No IAM losses
    NoLoss,
}

/// Spectral modifier model selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpectralModel {
    /// No spectral losses
    NoLoss,
}

/// Cell temperature model selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum TemperatureModel {
    /// SAPM cell temperature model
    SAPM,
    /// PVsyst cell temperature model
    PVSyst,
    /// Faiman cell temperature model
    Faiman,
    /// Fuentes cell temperature model
    Fuentes,
    /// NOCT SAM cell temperature model
    NOCT_SAM,
    /// Simple PVWatts-style (T_air + POA*(NOCT-20)/800)
    PVWatts,
}

/// Transposition model selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TranspositionModel {
    /// Isotropic sky model
    Isotropic,
    /// Hay-Davies model
    HayDavies,
    /// Perez 1990 model
    Perez,
    /// Klucher model
    Klucher,
    /// Reindl model
    Reindl,
}

/// Losses model selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LossesModel {
    /// PVWatts default losses
    PVWatts,
    /// No additional losses
    NoLoss,
}

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for all model selections in a ModelChain run.
#[derive(Debug, Clone)]
pub struct ModelChainConfig {
    pub dc_model: DCModel,
    pub ac_model: ACModel,
    pub aoi_model: AOIModel,
    pub spectral_model: SpectralModel,
    pub temperature_model: TemperatureModel,
    pub transposition_model: TranspositionModel,
    pub losses_model: LossesModel,
}

// ---------------------------------------------------------------------------
// Weather input
// ---------------------------------------------------------------------------

/// Weather data for a single timestep.
#[derive(Debug, Clone)]
pub struct WeatherInput {
    /// Timestamp with timezone
    pub time: DateTime<Tz>,
    /// Global horizontal irradiance [W/m2] (optional if DNI+DHI provided)
    pub ghi: Option<f64>,
    /// Direct normal irradiance [W/m2]
    pub dni: Option<f64>,
    /// Diffuse horizontal irradiance [W/m2]
    pub dhi: Option<f64>,
    /// Ambient air temperature [C]
    pub temp_air: f64,
    /// Wind speed [m/s]
    pub wind_speed: f64,
    /// Surface albedo (optional, default 0.25)
    pub albedo: Option<f64>,
}

/// POA (plane of array) input data for run_model_from_poa.
#[derive(Debug, Clone)]
pub struct POAInput {
    /// Timestamp with timezone
    pub time: DateTime<Tz>,
    /// Direct POA irradiance [W/m2]
    pub poa_direct: f64,
    /// Diffuse POA irradiance [W/m2]
    pub poa_diffuse: f64,
    /// Global POA irradiance [W/m2]
    pub poa_global: f64,
    /// Ambient air temperature [C]
    pub temp_air: f64,
    /// Wind speed [m/s]
    pub wind_speed: f64,
    /// Angle of incidence [degrees]
    pub aoi: f64,
}

/// Effective irradiance input data for run_model_from_effective_irradiance.
#[derive(Debug, Clone)]
pub struct EffectiveIrradianceInput {
    /// Timestamp with timezone
    pub time: DateTime<Tz>,
    /// Effective irradiance reaching cells [W/m2]
    pub effective_irradiance: f64,
    /// Global POA irradiance [W/m2] (for cell temperature calculation)
    pub poa_global: f64,
    /// Ambient air temperature [C]
    pub temp_air: f64,
    /// Wind speed [m/s]
    pub wind_speed: f64,
}

// ---------------------------------------------------------------------------
// Results
// ---------------------------------------------------------------------------

/// Full simulation result from ModelChain.
#[derive(Debug, Clone, PartialEq)]
pub struct ModelChainResult {
    /// Solar zenith angle [degrees]
    pub solar_zenith: f64,
    /// Solar azimuth angle [degrees]
    pub solar_azimuth: f64,
    /// Absolute airmass
    pub airmass: f64,
    /// Angle of incidence [degrees]
    pub aoi: f64,
    /// POA irradiance components
    pub poa_global: f64,
    pub poa_direct: f64,
    pub poa_diffuse: f64,
    /// AOI modifier (IAM)
    pub aoi_modifier: f64,
    /// Spectral modifier
    pub spectral_modifier: f64,
    /// Effective irradiance reaching cells [W/m2]
    pub effective_irradiance: f64,
    /// Cell temperature [C]
    pub cell_temperature: f64,
    /// DC power output [W]
    pub dc_power: f64,
    /// AC power output [W]
    pub ac_power: f64,
}

/// Legacy simulation result for backward compatibility.
#[derive(Debug, Clone, PartialEq)]
pub struct SimulationResult {
    pub poa_global: f64,
    pub temp_cell: f64,
    pub dc_power: f64,
    pub ac_power: f64,
}

// ---------------------------------------------------------------------------
// ModelChain
// ---------------------------------------------------------------------------

/// Orchestrates the full PV system performance simulation pipeline.
pub struct ModelChain {
    pub system: PVSystem,
    pub location: Location,
    pub surface_tilt: f64,
    pub surface_azimuth: f64,
    pub inverter_pac0: f64,
    pub inverter_eta: f64,
    pub config: ModelChainConfig,
}

impl ModelChain {
    /// Create a new ModelChain with explicit parameters (legacy constructor).
    pub fn new(
        system: PVSystem,
        location: Location,
        surface_tilt: f64,
        surface_azimuth: f64,
        inverter_pac0: f64,
        inverter_eta: f64,
    ) -> Self {
        Self {
            system,
            location,
            surface_tilt,
            surface_azimuth,
            inverter_pac0,
            inverter_eta,
            config: ModelChainConfig {
                dc_model: DCModel::PVWatts,
                ac_model: ACModel::PVWatts,
                aoi_model: AOIModel::ASHRAE,
                spectral_model: SpectralModel::NoLoss,
                temperature_model: TemperatureModel::PVWatts,
                transposition_model: TranspositionModel::Isotropic,
                losses_model: LossesModel::NoLoss,
            },
        }
    }

    /// Create a ModelChain with a custom configuration.
    pub fn with_config(
        system: PVSystem,
        location: Location,
        surface_tilt: f64,
        surface_azimuth: f64,
        inverter_pac0: f64,
        inverter_eta: f64,
        config: ModelChainConfig,
    ) -> Self {
        Self {
            system,
            location,
            surface_tilt,
            surface_azimuth,
            inverter_pac0,
            inverter_eta,
            config,
        }
    }

    /// Factory: PVWatts-style configuration.
    ///
    /// Uses PVWatts DC/AC, Physical AOI, PVWatts temperature, Perez transposition.
    pub fn with_pvwatts(
        system: PVSystem,
        location: Location,
        surface_tilt: f64,
        surface_azimuth: f64,
        inverter_pac0: f64,
        inverter_eta: f64,
    ) -> Self {
        Self {
            system,
            location,
            surface_tilt,
            surface_azimuth,
            inverter_pac0,
            inverter_eta,
            config: ModelChainConfig {
                dc_model: DCModel::PVWatts,
                ac_model: ACModel::PVWatts,
                aoi_model: AOIModel::Physical,
                spectral_model: SpectralModel::NoLoss,
                temperature_model: TemperatureModel::PVWatts,
                transposition_model: TranspositionModel::Perez,
                losses_model: LossesModel::PVWatts,
            },
        }
    }

    /// Factory: SAPM-style configuration.
    ///
    /// Uses PVWatts DC, PVWatts AC, ASHRAE AOI, SAPM temperature, HayDavies transposition.
    pub fn with_sapm(
        system: PVSystem,
        location: Location,
        surface_tilt: f64,
        surface_azimuth: f64,
        inverter_pac0: f64,
        inverter_eta: f64,
    ) -> Self {
        Self {
            system,
            location,
            surface_tilt,
            surface_azimuth,
            inverter_pac0,
            inverter_eta,
            config: ModelChainConfig {
                dc_model: DCModel::PVWatts,
                ac_model: ACModel::PVWatts,
                aoi_model: AOIModel::ASHRAE,
                spectral_model: SpectralModel::NoLoss,
                temperature_model: TemperatureModel::SAPM,
                transposition_model: TranspositionModel::HayDavies,
                losses_model: LossesModel::NoLoss,
            },
        }
    }

    // -----------------------------------------------------------------------
    // Legacy run_model (backward compatible)
    // -----------------------------------------------------------------------

    /// Run the model for a single timestep with given weather data (legacy API).
    pub fn run_model(
        &self,
        time: DateTime<Tz>,
        _ghi: f64,
        dni: f64,
        dhi: f64,
        temp_air: f64,
        _wind_speed: f64,
    ) -> Result<SimulationResult, spa::SpaError> {
        let solpos = get_solarposition(&self.location, time)?;
        let incidence = aoi(self.surface_tilt, self.surface_azimuth, solpos.zenith, solpos.azimuth);
        let iam_mult = iam::ashrae(incidence, 0.05);
        let poa_diffuse_val = crate::irradiance::isotropic(self.surface_tilt, dhi);
        let poa_dir = poa_direct(incidence, dni);
        let poa_global = poa_dir * iam_mult + poa_diffuse_val;
        let temp_cell = temp_air + poa_global * (45.0 - 20.0) / 800.0;
        let pdc = self.system.get_dc_power_total(poa_global, temp_cell);
        let pdc0 = self.system.get_nameplate_dc_total();
        let eta_inv_nom = self.inverter_eta;
        let eta_inv_ref = 0.9637;
        let pac = inverter::pvwatts_ac(pdc, pdc0, eta_inv_nom, eta_inv_ref);

        Ok(SimulationResult {
            poa_global,
            temp_cell,
            dc_power: pdc,
            ac_power: pac,
        })
    }

    // -----------------------------------------------------------------------
    // Enhanced run methods
    // -----------------------------------------------------------------------

    /// Run the full simulation pipeline from weather data.
    ///
    /// Steps: solar position -> airmass -> AOI -> transposition -> AOI modifier
    /// -> spectral modifier -> effective irradiance -> cell temperature
    /// -> DC power -> AC power.
    pub fn run_model_from_weather(
        &self,
        weather: &WeatherInput,
    ) -> Result<ModelChainResult, spa::SpaError> {
        // Complete irradiance if needed
        let (ghi, dni, dhi) = self.resolve_irradiance(weather)?;
        let albedo = weather.albedo.unwrap_or(0.25);

        // 1. Solar position
        let solpos = get_solarposition(&self.location, weather.time)?;

        // 2. Airmass
        let am_rel = get_relative_airmass(solpos.zenith);
        let pressure = alt2pres(self.location.altitude);
        let am_abs = if am_rel.is_nan() {
            0.0
        } else {
            get_absolute_airmass(am_rel, pressure)
        };

        // 3. AOI
        let aoi_val = aoi(self.surface_tilt, self.surface_azimuth, solpos.zenith, solpos.azimuth);

        // 4. Day of year for extraterrestrial irradiance
        let day_of_year = weather.time.format("%j").to_string().parse::<i32>().unwrap_or(1);
        let dni_extra = get_extra_radiation(day_of_year);

        // 5. Transposition -> POA components
        let diffuse_model = match self.config.transposition_model {
            TranspositionModel::Isotropic => DiffuseModel::Isotropic,
            TranspositionModel::HayDavies => DiffuseModel::HayDavies,
            TranspositionModel::Perez => DiffuseModel::Perez,
            TranspositionModel::Klucher => DiffuseModel::Klucher,
            TranspositionModel::Reindl => DiffuseModel::Reindl,
        };

        let poa = get_total_irradiance(
            self.surface_tilt,
            self.surface_azimuth,
            solpos.zenith,
            solpos.azimuth,
            dni,
            ghi,
            dhi,
            albedo,
            diffuse_model,
            Some(dni_extra),
            if am_rel.is_nan() { None } else { Some(am_rel) },
        );

        // Continue from POA
        self.compute_from_poa(
            solpos.zenith,
            solpos.azimuth,
            am_abs,
            aoi_val,
            &poa,
            weather.temp_air,
            weather.wind_speed,
        )
    }

    /// Run the simulation starting from plane-of-array irradiance data.
    ///
    /// Skips transposition step but still computes solar position for airmass.
    pub fn run_model_from_poa(
        &self,
        input: &POAInput,
    ) -> Result<ModelChainResult, spa::SpaError> {
        let solpos = get_solarposition(&self.location, input.time)?;
        let am_rel = get_relative_airmass(solpos.zenith);
        let pressure = alt2pres(self.location.altitude);
        let am_abs = if am_rel.is_nan() { 0.0 } else { get_absolute_airmass(am_rel, pressure) };

        let poa = PoaComponents {
            poa_global: input.poa_global,
            poa_direct: input.poa_direct,
            poa_diffuse: input.poa_diffuse,
            poa_sky_diffuse: input.poa_diffuse,
            poa_ground_diffuse: 0.0,
        };

        self.compute_from_poa(
            solpos.zenith,
            solpos.azimuth,
            am_abs,
            input.aoi,
            &poa,
            input.temp_air,
            input.wind_speed,
        )
    }

    /// Run the simulation starting from effective irradiance.
    ///
    /// Skips solar position, transposition, AOI, and spectral steps.
    pub fn run_model_from_effective_irradiance(
        &self,
        input: &EffectiveIrradianceInput,
    ) -> Result<ModelChainResult, spa::SpaError> {
        let solpos = get_solarposition(&self.location, input.time)?;
        let am_rel = get_relative_airmass(solpos.zenith);
        let pressure = alt2pres(self.location.altitude);
        let am_abs = if am_rel.is_nan() { 0.0 } else { get_absolute_airmass(am_rel, pressure) };

        let temp_cell = self.calc_cell_temperature(
            input.poa_global,
            input.temp_air,
            input.wind_speed,
        );
        let pdc = self.calc_dc_power(input.effective_irradiance, temp_cell);
        let pac = self.calc_ac_power(pdc);

        Ok(ModelChainResult {
            solar_zenith: solpos.zenith,
            solar_azimuth: solpos.azimuth,
            airmass: am_abs,
            aoi: 0.0,
            poa_global: input.poa_global,
            poa_direct: 0.0,
            poa_diffuse: 0.0,
            aoi_modifier: 1.0,
            spectral_modifier: 1.0,
            effective_irradiance: input.effective_irradiance,
            cell_temperature: temp_cell,
            dc_power: pdc,
            ac_power: pac,
        })
    }

    /// Fill missing GHI, DNI, or DHI using the Erbs decomposition model.
    ///
    /// If all three are provided, returns them as-is. If GHI is provided but
    /// DNI or DHI are missing, uses Erbs to decompose. If GHI is missing but
    /// DNI and DHI are provided, computes GHI = DNI * cos(zenith) + DHI.
    pub fn complete_irradiance(
        &self,
        weather: &WeatherInput,
    ) -> Result<(f64, f64, f64), spa::SpaError> {
        self.resolve_irradiance(weather)
    }

    // -----------------------------------------------------------------------
    // Internal pipeline steps
    // -----------------------------------------------------------------------

    fn resolve_irradiance(
        &self,
        weather: &WeatherInput,
    ) -> Result<(f64, f64, f64), spa::SpaError> {
        match (weather.ghi, weather.dni, weather.dhi) {
            (Some(ghi), Some(dni), Some(dhi)) => Ok((ghi, dni, dhi)),
            (Some(ghi), _, _) => {
                // Decompose GHI -> DNI + DHI using Erbs
                let solpos = get_solarposition(&self.location, weather.time)?;
                let day_of_year = weather.time.format("%j").to_string().parse::<i32>().unwrap_or(1);
                let dni_extra = get_extra_radiation(day_of_year);
                let (dni, dhi) = erbs(ghi, solpos.zenith, day_of_year as u32, dni_extra);
                Ok((ghi, dni, dhi))
            }
            (None, Some(dni), Some(dhi)) => {
                // Compute GHI from DNI + DHI
                let solpos = get_solarposition(&self.location, weather.time)?;
                let cos_z = solpos.zenith.to_radians().cos().max(0.0);
                let ghi = dni * cos_z + dhi;
                Ok((ghi, dni, dhi))
            }
            _ => {
                // Not enough data, return zeros
                Ok((0.0, 0.0, 0.0))
            }
        }
    }

    fn compute_from_poa(
        &self,
        solar_zenith: f64,
        solar_azimuth: f64,
        airmass: f64,
        aoi_val: f64,
        poa: &PoaComponents,
        temp_air: f64,
        wind_speed: f64,
    ) -> Result<ModelChainResult, spa::SpaError> {
        // AOI modifier (IAM)
        let aoi_modifier = self.calc_aoi_modifier(aoi_val);

        // Spectral modifier
        let spectral_modifier = self.calc_spectral_modifier();

        // Effective irradiance
        let effective_irradiance = (poa.poa_direct * aoi_modifier + poa.poa_diffuse)
            * spectral_modifier;

        // Cell temperature
        let temp_cell = self.calc_cell_temperature(poa.poa_global, temp_air, wind_speed);

        // DC power
        let pdc = self.calc_dc_power(effective_irradiance, temp_cell);

        // AC power
        let pac = self.calc_ac_power(pdc);

        Ok(ModelChainResult {
            solar_zenith,
            solar_azimuth,
            airmass,
            aoi: aoi_val,
            poa_global: poa.poa_global,
            poa_direct: poa.poa_direct,
            poa_diffuse: poa.poa_diffuse,
            aoi_modifier,
            spectral_modifier,
            effective_irradiance,
            cell_temperature: temp_cell,
            dc_power: pdc,
            ac_power: pac,
        })
    }

    fn calc_aoi_modifier(&self, aoi_val: f64) -> f64 {
        match self.config.aoi_model {
            AOIModel::Physical => iam::physical(aoi_val, 1.526, 4.0, 0.002),
            AOIModel::ASHRAE => iam::ashrae(aoi_val, 0.05),
            AOIModel::MartinRuiz => iam::martin_ruiz(aoi_val, 0.16),
            AOIModel::SAPM => {
                iam::sapm(aoi_val, 1.0, -0.002438, 3.103e-4, -1.246e-5, 2.112e-7, -1.359e-9)
            }
            AOIModel::NoLoss => 1.0,
        }
    }

    fn calc_spectral_modifier(&self) -> f64 {
        match self.config.spectral_model {
            SpectralModel::NoLoss => 1.0,
        }
    }

    fn calc_cell_temperature(&self, poa_global: f64, temp_air: f64, wind_speed: f64) -> f64 {
        match self.config.temperature_model {
            TemperatureModel::SAPM => {
                let (temp_cell, _) = temperature::sapm_cell_temperature(
                    poa_global, temp_air, wind_speed, -3.56, -0.075, 3.0, 1000.0,
                );
                temp_cell
            }
            TemperatureModel::PVSyst => {
                temperature::pvsyst_cell_temperature(
                    poa_global, temp_air, wind_speed, 29.0, 0.0, 0.15, 0.9,
                )
            }
            TemperatureModel::Faiman => {
                temperature::faiman(poa_global, temp_air, wind_speed, 25.0, 6.84)
            }
            TemperatureModel::Fuentes => {
                temperature::fuentes(poa_global, temp_air, wind_speed, 45.0)
            }
            TemperatureModel::NOCT_SAM => {
                temperature::noct_sam_default(poa_global, temp_air, wind_speed, 45.0, 0.15)
            }
            TemperatureModel::PVWatts => {
                temp_air + poa_global * (45.0 - 20.0) / 800.0
            }
        }
    }

    fn calc_dc_power(&self, effective_irradiance: f64, temp_cell: f64) -> f64 {
        match self.config.dc_model {
            DCModel::PVWatts => self.system.get_dc_power_total(effective_irradiance, temp_cell),
        }
    }

    fn calc_ac_power(&self, pdc: f64) -> f64 {
        match self.config.ac_model {
            ACModel::PVWatts | ACModel::Sandia | ACModel::ADR => {
                let pdc0 = self.system.get_nameplate_dc_total();
                let eta_inv_nom = self.inverter_eta;
                let eta_inv_ref = 0.9637;
                inverter::pvwatts_ac(pdc, pdc0, eta_inv_nom, eta_inv_ref)
            }
        }
    }
}
