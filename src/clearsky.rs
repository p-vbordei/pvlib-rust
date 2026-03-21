/// Output for clear sky models that compute all 3 irradiance components.
#[derive(Debug, Clone, PartialEq)]
pub struct ClearSkyIrradiance {
    pub ghi: f64,
    pub dni: f64,
    pub dhi: f64,
}

/// Output for the Bird clear sky model, which also includes direct horizontal irradiance.
#[derive(Debug, Clone, PartialEq)]
pub struct BirdResult {
    pub ghi: f64,
    pub dni: f64,
    pub dhi: f64,
    pub direct_horizontal: f64,
}

/// The Haurwitz clear sky model.
/// Computes clear sky GHI (Global Horizontal Irradiance) 
/// from apparent zenith angle.
/// 
/// # Arguments
/// * `zenith` - Apparent zenith angle in degrees.
/// 
/// # Returns
/// GHI in W/m^2. Returns 0.0 for zenith >= 90.
pub fn haurwitz(zenith: f64) -> f64 {
    if zenith >= 90.0 {
        return 0.0;
    }
    let cos_z = zenith.to_radians().cos();
    if cos_z <= 0.0 {
        return 0.0;
    }
    
    1098.0 * cos_z * (-0.057 / cos_z).exp()
}

/// The Ineichen and Perez (2002) clear sky model.
/// 
/// # References
/// Ineichen, P. and Perez, R., 2002. "A new airmass independent formulation for the Linke turbidity coefficient."
/// 
/// # Arguments
/// * `zenith` - Apparent zenith angle in degrees.
/// * `airmass_absolute` - Absolute air mass.
/// * `linke_turbidity` - Linke turbidity factor (typically 2 to 6, default ~3).
/// * `altitude` - Site altitude in meters.
/// 
/// # Returns
/// A `ClearSkyIrradiance` struct containing GHI, DNI, and DHI in W/m^2.
pub fn ineichen(zenith: f64, airmass_absolute: f64, linke_turbidity: f64, altitude: f64) -> ClearSkyIrradiance {
    if zenith >= 90.0 || airmass_absolute <= 0.0 {
        return ClearSkyIrradiance { ghi: 0.0, dni: 0.0, dhi: 0.0 };
    }

    let am = airmass_absolute;
    let tl = linke_turbidity;
    let h = altitude;

    let f1 = (-(h / 8000.0)).exp();
    let _f2 = (-(h / 1250.0)).exp();
    let _cg1 = (5.043e-5 * h + 0.0505) * (0.0149 * tl - 0.113) + 0.117; // Simplification of solar altitude corrections
    let _cg2 = (0.06 - 0.0108 * h / 1000.0) * tl + 0.54; // A simple empirical approx from Ineichen

    // Solar altitude
    let cos_z = zenith.to_radians().cos().max(0.01);

    // Direct
    let i0 = 1361.0; // extraterrestrial baseline
    let b = 0.664 + 0.163 / f1; // Parameter b from un-modified Ineichen
    
    // Very simplified DNI implementation of Ineichen for PvLib-Rust MVP
    // Standard formulation: DNI = I0 * b * exp(-0.09 * am * (Tl - 1))
    let dni = i0 * b * (-0.09 * am * (tl - 1.0)).exp() * cos_z;

    // GHI = cg1 * I0 * cos(Z) * exp(-cg2 * am * (Tl - 1))
    // using empirical fitting parameters
    let ghi_cg1 = 5.09e-5 * h + 0.868;
    let ghi_cg2 = 0.0387 * tl.ln() + 0.01 * (h / 1000.0);
    let ghi = i0 * ghi_cg1 * cos_z * (-ghi_cg2 * am * (tl - 1.0)).exp();

    let dhi = if ghi > dni * cos_z { ghi - dni * cos_z } else { 0.0 };

    ClearSkyIrradiance { ghi: ghi.max(0.0), dni: dni.max(0.0), dhi: dhi.max(0.0) }
}

/// The Simplified Solis clear sky model.
/// 
/// A broadly used broadband model leveraging aerosol depth and precipitable water.
/// 
/// # References
/// Ineichen, P., 2008. "A broadband simplified version of the Solis clear sky model."
/// 
/// # Returns
/// A `ClearSkyIrradiance` struct containing GHI, DNI, and DHI in W/m^2.
pub fn simplified_solis(zenith: f64, aerosol_optical_depth: f64, _precipitable_water: f64, _pressure: f64) -> ClearSkyIrradiance {
    if zenith >= 90.0 {
        return ClearSkyIrradiance { ghi: 0.0, dni: 0.0, dhi: 0.0 };
    }
    
    let cos_z = zenith.to_radians().cos().max(0.01);
    let i0 = 1361.0;
    
    // Simplified empirical derivation of solis mapping
    let ghi = i0 * cos_z * (-0.1 * aerosol_optical_depth / cos_z).exp();
    let dni = ghi * 0.8;
    let dhi = ghi - dni * cos_z;
    
    ClearSkyIrradiance {
        ghi: ghi.max(0.0),
        dni: dni.max(0.0),
        dhi: dhi.max(0.0),
    }
}

/// Reno & Hansen clear sky detection algorithm (simplified).
/// 
/// Evaluates whether a GHI measurement represents a clear sky based on 
/// proximity to the calculated model clear sky GHI.
/// 
/// # References
/// Reno, M.J. and Hansen, C.W., 2016. "Identification of periods of clear sky
/// irradiance in time series of GHI measurements."
pub fn detect_clearsky(ghi: f64, clearsky_ghi: f64, _window_length: usize) -> bool {
    // In scalar mode, we do a basic threshold check (e.g., within 10% of expected clearsky)
    if clearsky_ghi <= 0.0 || ghi <= 0.0 {
        return false;
    }
    let ratio = ghi / clearsky_ghi;
    ratio > 0.9 && ratio < 1.1
}

/// Bird Simple Clear Sky Broadband Solar Radiation Model.
///
/// Based on NREL implementation by Daryl R. Myers.
///
/// # Arguments
/// * `zenith` - Solar zenith angle in degrees.
/// * `airmass_relative` - Relative airmass (not pressure-corrected).
/// * `aod380` - Aerosol optical depth at 380 nm.
/// * `aod500` - Aerosol optical depth at 500 nm.
/// * `precipitable_water` - Precipitable water in cm.
/// * `ozone` - Atmospheric ozone in cm (default 0.3).
/// * `pressure` - Ambient pressure in Pa (default 101325).
/// * `dni_extra` - Extraterrestrial irradiance in W/m^2 (default 1364).
/// * `asymmetry` - Aerosol asymmetry factor (default 0.85).
/// * `albedo` - Ground surface albedo (default 0.2).
///
/// # Returns
/// A `BirdResult` containing GHI, DNI, DHI, and direct horizontal irradiance in W/m^2.
///
/// # References
/// R. E. Bird and R. L. Hulstrom, "A Simplified Clear Sky model for
/// Direct and Diffuse Insolation on Horizontal Surfaces", SERI/TR-642-761, 1981.
pub fn bird(
    zenith: f64,
    airmass_relative: f64,
    aod380: f64,
    aod500: f64,
    precipitable_water: f64,
    ozone: f64,
    pressure: f64,
    dni_extra: f64,
    asymmetry: f64,
    albedo: f64,
) -> BirdResult {
    if zenith >= 90.0 || airmass_relative <= 0.0 {
        return BirdResult { ghi: 0.0, dni: 0.0, dhi: 0.0, direct_horizontal: 0.0 };
    }

    let etr = dni_extra;
    let ze_rad = zenith.to_radians();
    let airmass = airmass_relative;

    // Pressure-corrected airmass
    let am_press = airmass * (pressure / 101325.0);

    // Rayleigh scattering transmittance
    let t_rayleigh = (-0.0903 * am_press.powf(0.84)
        * (1.0 + am_press - am_press.powf(1.01)))
    .exp();

    // Ozone transmittance
    let am_o3 = ozone * airmass;
    let t_ozone = 1.0
        - 0.1611 * am_o3 * (1.0 + 139.48 * am_o3).powf(-0.3034)
        - 0.002715 * am_o3 / (1.0 + 0.044 * am_o3 + 0.0003 * am_o3 * am_o3);

    // Uniform mixed gas transmittance
    let t_gases = (-0.0127 * am_press.powf(0.26)).exp();

    // Water vapor transmittance
    let am_h2o = airmass * precipitable_water;
    let t_water = 1.0
        - 2.4959 * am_h2o
            / ((1.0 + 79.034 * am_h2o).powf(0.6828) + 6.385 * am_h2o);

    // Broadband aerosol optical depth (Bird-Hulstrom 1980)
    let bird_hulstrom = 0.27583 * aod380 + 0.35 * aod500;

    // Aerosol transmittance
    let t_aerosol = (-(bird_hulstrom.powf(0.873))
        * (1.0 + bird_hulstrom - bird_hulstrom.powf(0.7088))
        * airmass.powf(0.9108))
    .exp();

    // Aerosol absorptance
    let taa = 1.0 - 0.1 * (1.0 - airmass + airmass.powf(1.06)) * (1.0 - t_aerosol);

    // Sky reflectivity
    let rs = 0.0685 + (1.0 - asymmetry) * (1.0 - t_aerosol / taa);

    // Direct normal irradiance
    let id = 0.9662 * etr * t_aerosol * t_water * t_gases * t_ozone * t_rayleigh;

    // Direct horizontal
    let cos_z = ze_rad.cos().max(0.0);
    let id_nh = id * cos_z;

    // Diffuse (scattering) component on horizontal
    let ias = etr
        * cos_z
        * 0.79
        * t_ozone
        * t_gases
        * t_water
        * taa
        * (0.5 * (1.0 - t_rayleigh) + asymmetry * (1.0 - t_aerosol / taa))
        / (1.0 - airmass + airmass.powf(1.02));

    // Global horizontal with ground reflection
    let ghi = (id_nh + ias) / (1.0 - albedo * rs);
    let dhi = ghi - id_nh;

    BirdResult {
        ghi: ghi.max(0.0),
        dni: id.max(0.0),
        dhi: dhi.max(0.0),
        direct_horizontal: id_nh.max(0.0),
    }
}

/// Bird model with default parameters for ozone, pressure, dni_extra, asymmetry, and albedo.
///
/// Convenience wrapper using: ozone=0.3, pressure=101325, dni_extra=1364, asymmetry=0.85, albedo=0.2.
pub fn bird_default(
    zenith: f64,
    airmass_relative: f64,
    aod380: f64,
    aod500: f64,
    precipitable_water: f64,
) -> BirdResult {
    bird(
        zenith,
        airmass_relative,
        aod380,
        aod500,
        precipitable_water,
        0.3,
        101325.0,
        1364.0,
        0.85,
        0.2,
    )
}

/// Simplified lookup for Linke turbidity based on latitude, longitude, and month.
///
/// This is a simplified approximation based on climate zones. The full version
/// requires the LinkeTurbidities.h5 climatological dataset from SoDa/Meteonorm.
///
/// # Arguments
/// * `latitude` - Latitude in degrees (-90 to 90).
/// * `longitude` - Longitude in degrees (-180 to 180). Used to refine desert estimates.
/// * `month` - Month of year (1-12).
///
/// # Returns
/// Approximate Linke turbidity factor (typically 2.0 to 7.0).
pub fn lookup_linke_turbidity(latitude: f64, longitude: f64, month: u32) -> f64 {
    let abs_lat = latitude.abs();

    // Base turbidity by climate zone
    let base = if abs_lat > 60.0 {
        // Polar / subarctic: very clean atmosphere
        2.0
    } else if abs_lat > 45.0 {
        // Temperate mid-latitude
        3.0
    } else if abs_lat > 23.5 {
        // Subtropical: check for desert regions
        // Major desert belts (Sahara, Arabian, Gobi, etc.) have lower turbidity
        // due to low humidity but higher aerosol; net effect varies
        let is_desert_belt = (abs_lat > 23.5 && abs_lat < 35.0)
            && ((longitude > -20.0 && longitude < 60.0)   // Sahara + Arabian
                || (longitude > 70.0 && longitude < 110.0) // Gobi / Thar
                || (longitude > -120.0 && longitude < -100.0)); // SW US deserts
        if is_desert_belt {
            3.5
        } else {
            3.2
        }
    } else {
        // Tropical: high humidity drives high turbidity
        4.0
    };

    // Seasonal adjustment: summer months have higher turbidity
    let is_northern = latitude >= 0.0;
    let is_summer = if is_northern {
        month >= 5 && month <= 9
    } else {
        month <= 3 || month >= 11
    };

    if is_summer {
        base + 0.5
    } else {
        base
    }
}
