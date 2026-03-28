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
#[inline]
pub fn haurwitz(zenith: f64) -> f64 {
    if zenith >= 90.0 {
        return 0.0;
    }
    let cos_z = zenith.to_radians().cos();
    if cos_z <= 0.0 {
        return 0.0;
    }
    
    1098.0 * cos_z * (-0.059 / cos_z).exp()
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
/// * `dni_extra` - Extraterrestrial normal irradiance in W/m^2.
///
/// # Returns
/// A `ClearSkyIrradiance` struct containing GHI, DNI, and DHI in W/m^2.
#[inline]
pub fn ineichen(zenith: f64, airmass_absolute: f64, linke_turbidity: f64, altitude: f64, dni_extra: f64) -> ClearSkyIrradiance {
    if zenith >= 90.0 || airmass_absolute <= 0.0 {
        return ClearSkyIrradiance { ghi: 0.0, dni: 0.0, dhi: 0.0 };
    }

    let am = airmass_absolute;
    let tl = linke_turbidity;
    let h = altitude;

    let fh1 = (-h / 8000.0).exp();
    let fh2 = (-h / 1250.0).exp();

    let cos_z = zenith.to_radians().cos().max(0.01);

    let i0 = dni_extra;

    let cg1 = 5.09e-05 * h + 0.868;
    let cg2 = 3.92e-05 * h + 0.0387;

    // GHI
    let ghi = cg1 * i0 * cos_z * (-cg2 * am * (fh1 + fh2 * (tl - 1.0))).exp().max(0.0);

    // DNI
    let b = 0.664 + 0.163 / fh1;
    let bnci = i0 * (b * (-0.09 * am * (tl - 1.0)).exp()).max(0.0);
    let denom = cos_z;
    let ratio = if denom > 0.0 {
        ((1.0 - (0.1 - 0.2 * (-tl).exp()) / (0.1 + 0.882 / fh1)) / denom).clamp(0.0, 1e20)
    } else {
        0.0
    };
    let bnci_2 = ghi * ratio;
    let dni = bnci.min(bnci_2);

    // DHI
    let dhi = ghi - dni * cos_z;

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
/// Simplified Solis clear sky model.
///
/// Note: this function takes `apparent_elevation` (degrees above horizon),
/// NOT zenith angle, matching the pvlib-python API.
///
/// # Arguments
/// * `apparent_elevation` - Apparent solar elevation in degrees above the horizon.
/// * `aod700` - Aerosol optical depth at 700 nm (default 0.1).
/// * `precipitable_water` - Precipitable water in cm (default 1.0, minimum 0.2).
/// * `pressure` - Atmospheric pressure in Pascals (default 101325).
///
/// # Returns
/// A `ClearSkyIrradiance` struct containing GHI, DNI, and DHI in W/m^2.
#[inline]
pub fn simplified_solis(apparent_elevation: f64, aod700: f64, precipitable_water: f64, pressure: f64) -> ClearSkyIrradiance {
    if apparent_elevation <= 0.0 {
        return ClearSkyIrradiance { ghi: 0.0, dni: 0.0, dhi: 0.0 };
    }

    let dni_extra = 1364.0;
    let p = pressure;
    let p0 = 101325.0_f64;
    let w = precipitable_water.max(0.2);
    let aod = aod700;
    let ln_w = w.ln();
    let ln_p = (p / p0).ln();

    // i0p: enhanced extraterrestrial irradiance
    let io0 = 1.08 * w.powf(0.0051);
    let i01 = 0.97 * w.powf(0.032);
    let i02 = 0.12 * w.powf(0.56);
    let i0p = dni_extra * (i02 * aod * aod + i01 * aod + io0 + 0.071 * ln_p);

    // taub, b coefficients (DNI)
    let tb1 = 1.82 + 0.056 * ln_w + 0.0071 * ln_w * ln_w;
    let tb0 = 0.33 + 0.045 * ln_w + 0.0096 * ln_w * ln_w;
    let tbp = 0.0089 * w + 0.13;
    let taub = tb1 * aod + tb0 + tbp * ln_p;

    let b1 = 0.00925 * aod * aod + 0.0148 * aod - 0.0172;
    let b0 = -0.7565 * aod * aod + 0.5057 * aod + 0.4557;
    let b = b1 * ln_w + b0;

    // taug, g coefficients (GHI)
    let tg1 = 1.24 + 0.047 * ln_w + 0.0061 * ln_w * ln_w;
    let tg0 = 0.27 + 0.043 * ln_w + 0.0090 * ln_w * ln_w;
    let tgp = 0.0079 * w + 0.1;
    let taug = tg1 * aod + tg0 + tgp * ln_p;

    let g = -0.0147 * ln_w - 0.3079 * aod * aod + 0.2846 * aod + 0.3798;

    // taud, d coefficients (DHI)
    let (td4, td3, td2, td1, td0, tdp) = if aod < 0.05 {
        (
            86.0 * w - 13800.0,
            -3.11 * w + 79.4,
            -0.23 * w + 74.8,
            0.092 * w - 8.86,
            0.0042 * w + 3.12,
            -0.83 * (1.0 + aod).powf(-17.2),
        )
    } else {
        (
            -0.21 * w + 11.6,
            0.27 * w - 20.7,
            -0.134 * w + 15.5,
            0.0554 * w - 5.71,
            0.0057 * w + 2.94,
            -0.71 * (1.0 + aod).powf(-15.0),
        )
    };
    let taud = td4 * aod.powi(4) + td3 * aod.powi(3) + td2 * aod.powi(2) + td1 * aod + td0 + tdp * ln_p;

    let dp = 1.0 / (18.0 + 152.0 * aod);
    let d = -0.337 * aod * aod + 0.63 * aod + 0.116 + dp * ln_p;

    // Compute irradiances
    let sin_elev = apparent_elevation.to_radians().sin().max(1e-30);

    let dni = i0p * (-taub / sin_elev.powf(b)).exp();
    let ghi = i0p * (-taug / sin_elev.powf(g)).exp() * sin_elev;
    let dhi = i0p * (-taud / sin_elev.powf(d)).exp();

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
#[inline]
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
#[allow(clippy::too_many_arguments)]
#[inline]
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
#[inline]
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
#[inline]
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
        (5..=9).contains(&month)
    } else {
        month <= 3 || month >= 11
    };

    if is_summer {
        base + 0.5
    } else {
        base
    }
}
