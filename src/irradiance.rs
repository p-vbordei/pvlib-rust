use std::f64::consts::PI;
use crate::atmosphere;

/// Perez et al. (1990) Table 1 — brightness coefficients
/// Columns: f11, f12, f13, f21, f22, f23
/// Rows: 8 sky clearness bins (epsilon ranges)
const PEREZ_COEFFICIENTS: [[f64; 6]; 8] = [
    [-0.0083117, 0.5877277, -0.0620636, -0.0596012, 0.0721249, -0.0220216],
    [0.1299457, 0.6825954, -0.1513752, -0.0189325, 0.0659650, -0.0288748],
    [0.3296958, 0.4868735, -0.2210958, 0.0554140, -0.0639588, -0.0260542],
    [0.5682053, 0.1874990, -0.2951290, 0.1088631, -0.1519229, -0.0139754],
    [0.8730280, -0.3920403, -0.3616149, 0.2255647, -0.4620442,  0.0012448],
    [1.1326077, -1.2367284, -0.4118494, 0.2877813, -0.8230357,  0.0558225],
    [1.0601591, -1.5999137, -0.3589221, 0.2642124, -1.1272340,  0.1310694],
    [0.6777470, -0.3272588, -0.2504286, 0.1561313, -1.3765031,  0.2506212],
];

/// Calculate the angle of incidence (AOI) of the solar vector on a surface.
pub fn aoi(surface_tilt: f64, surface_azimuth: f64, solar_zenith: f64, solar_azimuth: f64) -> f64 {
    let tilt_rad = surface_tilt.to_radians();
    let surf_az_rad = surface_azimuth.to_radians();
    let zen_rad = solar_zenith.to_radians();
    let sol_az_rad = solar_azimuth.to_radians();

    let cos_aoi = zen_rad.cos() * tilt_rad.cos()
        + zen_rad.sin() * tilt_rad.sin() * (sol_az_rad - surf_az_rad).cos();
    
    let cos_aoi = cos_aoi.clamp(-1.0, 1.0);
    cos_aoi.acos().to_degrees()
}

/// Calculate extraterrestrial solar irradiance for a day of year (Spencer 1971).
pub fn get_extra_radiation(dayofyear: i32) -> f64 {
    let b = 2.0 * PI * ((dayofyear - 1) as f64) / 365.0;
    let rover_r0_sqrd = 1.00011
        + 0.034221 * b.cos()
        + 0.00128 * b.sin()
        + 0.000719 * (2.0 * b).cos()
        + 0.000077 * (2.0 * b).sin();
    1366.1 * rover_r0_sqrd
}

/// Isotropic diffuse model.
pub fn isotropic(surface_tilt: f64, dhi: f64) -> f64 {
    dhi * (1.0 + surface_tilt.to_radians().cos()) / 2.0
}

/// Hay-Davies diffuse sky model.
/// 
/// # References
/// Hay, J.E. and Davies, J.A., 1980, "Calculations of the solar radiation incident on an inclined surface", 
/// in Proceedings of the First Canadian Solar Radiation Data Workshop.
pub fn haydavies(surface_tilt: f64, _surface_azimuth: f64, dhi: f64, dni: f64, dni_extra: f64, solar_zenith: f64, _solar_azimuth: f64, aoi_in: f64) -> f64 {
    let mut a = 0.0;
    if dni_extra > 0.0 {
        a = dni / dni_extra;
    }
    let a = a.clamp(0.0, 1.0);
    let mut cos_z = solar_zenith.to_radians().cos();
    if cos_z < 0.0436 { cos_z = 0.0436; }
    
    let cos_aoi = aoi_in.to_radians().cos().max(0.0);
    let r_b = cos_aoi / cos_z;
    
    dhi * ((1.0 - a) * (1.0 + surface_tilt.to_radians().cos()) / 2.0 + a * r_b)
}

/// Klucher diffuse sky model.
/// 
/// # References
/// Klucher, T.M., 1979, "Evaluation of models to predict insolation on tilted surfaces," 
/// Solar Energy, 23(2), pp. 111-114.
pub fn klucher(surface_tilt: f64, _surface_azimuth: f64, dhi: f64, ghi: f64, solar_zenith: f64, _solar_azimuth: f64, aoi_in: f64) -> f64 {
    let mut f = 0.0;
    if ghi > 0.0 {
        let frac = dhi / ghi;
        f = 1.0 - frac * frac;
    }
    let f = f.clamp(0.0, 1.0);
    
    let _cos_z = solar_zenith.to_radians().cos();
    let cos_aoi = aoi_in.to_radians().cos().max(0.0);
    let tilt_rad = surface_tilt.to_radians();
    
    let term1 = 1.0 + f * (tilt_rad / 2.0).sin().powi(3);
    let term2 = 1.0 + f * cos_aoi.powi(2) * (solar_zenith.to_radians().sin()).powi(3);
    
    dhi * ((1.0 + tilt_rad.cos()) / 2.0) * term1 * term2
}

/// Perez diffuse sky model.
/// 
/// # References
/// Perez, R., Ineichen, P., Seals, R., Michalsky, J. and Stewart, R., 1990, 
/// "Modeling daylight availability and irradiance components from direct and global irradiance," 
/// Solar Energy, 44(5), pp. 271-289.
pub fn perez(surface_tilt: f64, _surface_azimuth: f64, dhi: f64, dni: f64, dni_extra: f64, solar_zenith: f64, _solar_azimuth: f64, airmass: f64, aoi_in: f64) -> f64 {
    let mut cos_z = solar_zenith.to_radians().cos();
    if cos_z < 0.0436 { cos_z = 0.0436; } // cap to ~87.5 deg
    let cos_aoi = aoi_in.to_radians().cos().max(0.0); // beam parallel to surface if >90

    // sky clearness epsilon
    let _a = (dni_extra * 1e-6).max(1.0); // essentially 1.0 for bounds, simplified for delta
    let delta = dhi * airmass / dni_extra;
    
    let mut epsilon = 1.0;
    if dhi > 0.0 {
        epsilon = ((dhi + dni) / dhi + 1.041 * solar_zenith.to_radians().powi(3)) / 
                  (1.0 + 1.041 * solar_zenith.to_radians().powi(3));
    }

    let bin = if epsilon < 1.065 { 0 }
    else if epsilon < 1.230 { 1 }
    else if epsilon < 1.500 { 2 }
    else if epsilon < 1.950 { 3 }
    else if epsilon < 2.800 { 4 }
    else if epsilon < 4.500 { 5 }
    else if epsilon < 6.200 { 6 }
    else { 7 };

    let coeffs = PEREZ_COEFFICIENTS[bin];
    let mut f1 = coeffs[0] + coeffs[1] * delta + coeffs[2] * solar_zenith.to_radians();
    f1 = f1.max(0.0);
    let f2 = coeffs[3] + coeffs[4] * delta + coeffs[5] * solar_zenith.to_radians();

    let a_perez = cos_aoi;
    let b_perez = cos_z;

    dhi * ((1.0 - f1) * (1.0 + surface_tilt.to_radians().cos()) / 2.0 + f1 * a_perez / b_perez + f2 * surface_tilt.to_radians().sin())
}

/// Erbs decomposition model.
/// 
/// # References
/// Erbs, D.G., Klein, S.A. and Duffie, J.A., 1982, 
/// "Estimation of the diffuse radiation fraction for hourly, daily and monthly-average global radiation," 
/// Solar Energy, 28(4), pp. 293-302.
pub fn erbs(ghi: f64, zenith: f64, _day_of_year: u32, dni_extra: f64) -> (f64, f64) {
    if ghi <= 0.0 || zenith >= 90.0 { return (0.0, 0.0); }
    let mut cos_z = zenith.to_radians().cos();
    if cos_z < 0.0436 { cos_z = 0.0436; }
    
    let kt = ghi / (dni_extra * cos_z);
    
    let kd = if kt <= 0.22 {
        1.0 - 0.09 * kt
    } else if kt <= 0.80 {
        0.9511 - 0.1604 * kt + 4.388 * kt.powi(2) - 16.638 * kt.powi(3) + 12.336 * kt.powi(4)
    } else {
        0.165
    };
    
    let dhi = ghi * kd.clamp(0.0, 1.0);
    let dni = ((ghi - dhi) / cos_z).max(0.0);
    
    (dni, dhi)
}

/// Boland (2008) decomposition model.
/// Logistic regression model for continuous diffuse fraction estimation.
/// 
/// # References
/// Boland, J., Scott, L. and Luther, M., 2008. 
/// "Modelling the diffuse fraction of global solar radiation on a horizontal surface."
pub fn boland(ghi: f64, zenith: f64, dni_extra: f64) -> (f64, f64) {
    if ghi <= 0.0 || zenith >= 90.0 { return (0.0, 0.0); }
    let cos_z = zenith.to_radians().cos().max(0.0436);
    
    let kt = ghi / (dni_extra * cos_z);
    
    // Boland logistic equation: DF = 1 / (1 + exp(a*(kt - b)))
    // Default coefficients: a=8.645, b=0.613 (15-minute data, Boland et al.)
    let a_coeff = 8.645;
    let b_coeff = 0.613;
    let kd = 1.0 / (1.0 + (a_coeff * (kt - b_coeff)).exp());
    let dhi = ghi * kd.clamp(0.0, 1.0);
    let dni = ((ghi - dhi) / cos_z).max(0.0);
    
    (dni, dhi)
}

/// DIRINT (Perez 1992) decomposition model.
/// 
/// Note: This is a highly simplified representation of DIRINT for estimating DNI 
/// from GHI without full climatic parameter timeseries tracking.
/// 
/// # References
/// Perez, R., Ineichen, P., Maxwell, E., Seals, R. and Zelenka, A., 1992. 
/// "Dynamic global-to-direct irradiance conversion models."
pub fn dirint(ghi: f64, zenith: f64, _dew_point: f64, _pressure: f64, dni_extra: f64) -> (f64, f64) {
    // In a full time-series context, DIRINT uses persistence bins. 
    // Here we approximate it by defaulting to a slightly more aggressive Erbs.
    if ghi <= 0.0 || zenith >= 90.0 { return (0.0, 0.0); }
    let cos_z = zenith.to_radians().cos().max(0.0436);
    
    let kt = ghi / (dni_extra * cos_z);
    
    // Approximate diffuse fraction
    let kd = if kt <= 0.2 {
        0.99
    } else if kt <= 0.8 {
        0.95 - 0.9 * (kt - 0.2)
    } else {
        0.15
    };
    
    let dhi = ghi * kd.clamp(0.0, 1.0);
    let dni = ((ghi - dhi) / cos_z).max(0.0);
    (dni, dhi)
}

/// POA direct beam.
pub fn poa_direct(aoi_in: f64, dni: f64) -> f64 {
    let aoi_rad = aoi_in.to_radians();
    if aoi_rad.abs() > std::f64::consts::PI / 2.0 {
        0.0
    } else {
        (dni * aoi_rad.cos()).max(0.0)
    }
}

/// Reindl transposition model (anisotropic sky).
/// 
/// A highly cited mathematical model bridging the gap between Hay-Davies and Perez models.
/// 
/// # References
/// Reindl, D.T., Beckman, W.A. and Duffie, J.A., 1990. "Evaluation of hourly tilt data models".
#[allow(clippy::too_many_arguments)]
pub fn reindl(surface_tilt: f64, dhi: f64, ghi: f64, dni: f64, dni_extra: f64, solar_zenith: f64, aoi_in: f64) -> f64 {
    let mut a = 0.0;
    if dni_extra > 0.0 { a = dni / dni_extra; }
    let a = a.clamp(0.0, 1.0);
    
    let cos_z = solar_zenith.to_radians().cos().max(0.0436);
    let cos_aoi = aoi_in.to_radians().cos().max(0.0);
    let r_b = cos_aoi / cos_z;
    
    let f = if ghi > 0.0 { (dni / ghi).powi(2) } else { 0.0 };
    
    let tilt_rad = surface_tilt.to_radians();
    let term1 = dhi * (1.0 - a) * (1.0 + tilt_rad.cos()) / 2.0 * (1.0 + f * (tilt_rad / 2.0).sin().powi(3));
    let term2 = dhi * a * r_b;
    
    term1 + term2
}

/// Clearness Index (Kt).
/// 
/// The ratio of global horizontal irradiance to extraterrestrial horizontal irradiance.
pub fn clearness_index(ghi: f64, solar_zenith: f64, dni_extra: f64) -> f64 {
    let cos_z = solar_zenith.to_radians().cos().max(0.01);
    let ghi_extra = dni_extra * cos_z;
    if ghi_extra <= 0.0 { 0.0 } else { (ghi / ghi_extra).clamp(0.0, 1.0) }
}

/// Zenith-independent clearness index (Kt*).
///
/// # References
/// Perez, R. et al., 1990. "Making full use of the clearness index for parameterizing hourly insolation conditions."
pub fn clearness_index_zenith_independent(clearness_idx: f64, _solar_zenith: f64, airmass_absolute: f64) -> f64 {
    let am = airmass_absolute.max(1.0);
    // Approximation of the geometric zenith independence formula
    let denominator = 1.031 * (-1.4 / (0.9 + 9.4 / am)).exp() + 0.1;
    (clearness_idx / denominator).max(0.0)
}

/// Cosine of the angle of incidence (AOI projection).
///
/// Calculates the dot product of the sun position unit vector and the surface
/// normal unit vector. When the sun is behind the surface, the returned value
/// is negative. Input all angles in degrees.
///
/// # References
/// Same geometry as [`aoi`], but returns cos(AOI) without taking arccos.
pub fn aoi_projection(surface_tilt: f64, surface_azimuth: f64, solar_zenith: f64, solar_azimuth: f64) -> f64 {
    let tilt_rad = surface_tilt.to_radians();
    let surf_az_rad = surface_azimuth.to_radians();
    let zen_rad = solar_zenith.to_radians();
    let sol_az_rad = solar_azimuth.to_radians();

    let projection = zen_rad.cos() * tilt_rad.cos()
        + zen_rad.sin() * tilt_rad.sin() * (sol_az_rad - surf_az_rad).cos();

    projection.clamp(-1.0, 1.0)
}

/// Beam component of plane-of-array irradiance.
///
/// Calculates `DNI * max(cos(AOI), 0)`.
///
/// # Parameters
/// - `surface_tilt`: Panel tilt from horizontal [degrees]
/// - `surface_azimuth`: Panel azimuth [degrees]
/// - `solar_zenith`: Solar zenith angle [degrees]
/// - `solar_azimuth`: Solar azimuth angle [degrees]
/// - `dni`: Direct normal irradiance [W/m²]
pub fn beam_component(surface_tilt: f64, surface_azimuth: f64, solar_zenith: f64, solar_azimuth: f64, dni: f64) -> f64 {
    let proj = aoi_projection(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth);
    (dni * proj).max(0.0)
}

/// Ground-reflected diffuse irradiance on a tilted surface.
///
/// Calculated as `GHI * albedo * (1 - cos(tilt)) / 2`.
///
/// # Parameters
/// - `surface_tilt`: Panel tilt from horizontal [degrees]
/// - `ghi`: Global horizontal irradiance [W/m²]
/// - `albedo`: Ground surface albedo (typically 0.1–0.4) [unitless]
///
/// # References
/// Loutzenhiser P.G. et al., 2007, "Empirical validation of models to compute
/// solar irradiance on inclined surfaces for building energy simulation",
/// Solar Energy vol. 81, pp. 254-267.
pub fn get_ground_diffuse(surface_tilt: f64, ghi: f64, albedo: f64) -> f64 {
    ghi * albedo * (1.0 - surface_tilt.to_radians().cos()) * 0.5
}

/// Components of plane-of-array irradiance.
#[derive(Debug, Clone, Copy)]
pub struct PoaComponents {
    /// Total in-plane irradiance [W/m²]
    pub poa_global: f64,
    /// Total in-plane beam irradiance [W/m²]
    pub poa_direct: f64,
    /// Total in-plane diffuse irradiance [W/m²]
    pub poa_diffuse: f64,
    /// In-plane diffuse irradiance from sky [W/m²]
    pub poa_sky_diffuse: f64,
    /// In-plane diffuse irradiance from ground [W/m²]
    pub poa_ground_diffuse: f64,
}

/// Determine in-plane irradiance components.
///
/// Combines DNI with sky diffuse and ground-reflected irradiance to calculate
/// total, direct, and diffuse irradiance components in the plane of array.
/// Negative beam irradiation due to AOI > 90° is set to zero.
///
/// # Parameters
/// - `aoi_val`: Angle of incidence [degrees]
/// - `dni`: Direct normal irradiance [W/m²]
/// - `poa_sky_diffuse`: Sky diffuse irradiance in the plane of array [W/m²]
/// - `poa_ground_diffuse`: Ground-reflected irradiance in the plane of array [W/m²]
pub fn poa_components(aoi_val: f64, dni: f64, poa_sky_diffuse: f64, poa_ground_diffuse: f64) -> PoaComponents {
    let poa_direct = (dni * aoi_val.to_radians().cos()).max(0.0);
    let poa_diffuse = poa_sky_diffuse + poa_ground_diffuse;
    let poa_global = poa_direct + poa_diffuse;

    PoaComponents {
        poa_global,
        poa_direct,
        poa_diffuse,
        poa_sky_diffuse,
        poa_ground_diffuse,
    }
}

/// Result of total irradiance calculation.
pub type TotalIrradiance = PoaComponents;

/// Sky diffuse irradiance model selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiffuseModel {
    Isotropic,
    Klucher,
    HayDavies,
    Reindl,
    Perez,
}

/// Determine in-plane sky diffuse irradiance using the specified model.
///
/// Dispatches to the appropriate diffuse sky model: isotropic, klucher,
/// haydavies, reindl, or perez.
///
/// # Parameters
/// - `surface_tilt`: Panel tilt from horizontal [degrees]
/// - `surface_azimuth`: Panel azimuth [degrees]
/// - `solar_zenith`: Solar zenith angle [degrees]
/// - `solar_azimuth`: Solar azimuth angle [degrees]
/// - `dni`: Direct normal irradiance [W/m²]
/// - `ghi`: Global horizontal irradiance [W/m²]
/// - `dhi`: Diffuse horizontal irradiance [W/m²]
/// - `model`: Sky diffuse irradiance model
/// - `dni_extra`: Extraterrestrial DNI [W/m²] (required for HayDavies, Reindl, Perez)
/// - `airmass`: Relative airmass (required for Perez)
#[allow(clippy::too_many_arguments)]
pub fn get_sky_diffuse(
    surface_tilt: f64,
    surface_azimuth: f64,
    solar_zenith: f64,
    solar_azimuth: f64,
    dni: f64,
    ghi: f64,
    dhi: f64,
    model: DiffuseModel,
    dni_extra: Option<f64>,
    airmass: Option<f64>,
) -> f64 {
    let aoi_val = aoi(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth);

    match model {
        DiffuseModel::Isotropic => isotropic(surface_tilt, dhi),
        DiffuseModel::Klucher => klucher(surface_tilt, surface_azimuth, dhi, ghi, solar_zenith, solar_azimuth, aoi_val),
        DiffuseModel::HayDavies => {
            let extra = dni_extra.unwrap_or(0.0);
            haydavies(surface_tilt, surface_azimuth, dhi, dni, extra, solar_zenith, solar_azimuth, aoi_val)
        }
        DiffuseModel::Reindl => {
            let extra = dni_extra.unwrap_or(0.0);
            reindl(surface_tilt, dhi, ghi, dni, extra, solar_zenith, aoi_val)
        }
        DiffuseModel::Perez => {
            let extra = dni_extra.unwrap_or(0.0);
            let am = airmass.unwrap_or_else(|| atmosphere::get_relative_airmass(solar_zenith));
            perez(surface_tilt, surface_azimuth, dhi, dni, extra, solar_zenith, solar_azimuth, am, aoi_val)
        }
    }
}

/// Determine total in-plane irradiance and its beam, sky diffuse, and ground
/// reflected components using the specified sky diffuse irradiance model.
///
/// # Parameters
/// - `surface_tilt`: Panel tilt from horizontal [degrees]
/// - `surface_azimuth`: Panel azimuth [degrees]
/// - `solar_zenith`: Solar zenith angle [degrees]
/// - `solar_azimuth`: Solar azimuth angle [degrees]
/// - `dni`: Direct normal irradiance [W/m²]
/// - `ghi`: Global horizontal irradiance [W/m²]
/// - `dhi`: Diffuse horizontal irradiance [W/m²]
/// - `albedo`: Ground surface albedo [unitless]
/// - `model`: Sky diffuse irradiance model
/// - `dni_extra`: Extraterrestrial DNI [W/m²] (required for HayDavies, Reindl, Perez)
/// - `airmass`: Relative airmass (required for Perez)
#[allow(clippy::too_many_arguments)]
pub fn get_total_irradiance(
    surface_tilt: f64,
    surface_azimuth: f64,
    solar_zenith: f64,
    solar_azimuth: f64,
    dni: f64,
    ghi: f64,
    dhi: f64,
    albedo: f64,
    model: DiffuseModel,
    dni_extra: Option<f64>,
    airmass: Option<f64>,
) -> TotalIrradiance {
    let aoi_val = aoi(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth);

    let sky_diffuse = get_sky_diffuse(
        surface_tilt, surface_azimuth, solar_zenith, solar_azimuth,
        dni, ghi, dhi, model, dni_extra, airmass,
    );

    let ground_diffuse = get_ground_diffuse(surface_tilt, ghi, albedo);

    poa_components(aoi_val, dni, sky_diffuse, ground_diffuse)
}

/// Output of the DISC decomposition model.
#[derive(Debug, Clone, Copy)]
pub struct DiscOutput {
    /// Direct normal irradiance [W/m²]
    pub dni: f64,
    /// Clearness index [unitless]
    pub kt: f64,
    /// Airmass used in the calculation [unitless]
    pub airmass: f64,
}

/// DISC model helper: calculate Kn from clearness index and airmass.
fn disc_kn(kt: f64, am: f64) -> (f64, f64) {
    let am = am.min(12.0);

    let (a, b, c) = if kt <= 0.6 {
        (
            0.512 + kt * (-1.56 + kt * (2.286 - 2.222 * kt)),
            0.37 + 0.962 * kt,
            -0.28 + kt * (0.932 - 2.048 * kt),
        )
    } else {
        (
            -5.743 + kt * (21.77 + kt * (-27.49 + 11.56 * kt)),
            41.4 + kt * (-118.5 + kt * (66.05 + 31.9 * kt)),
            -47.01 + kt * (184.2 + kt * (-222.0 + 73.81 * kt)),
        )
    };

    let delta_kn = a + b * (c * am).exp();
    let knc = 0.866 + am * (-0.122 + am * (0.0121 + am * (-0.000653 + 1.4e-05 * am)));
    let kn = knc - delta_kn;

    (kn, am)
}

/// Estimate Direct Normal Irradiance from Global Horizontal Irradiance
/// using the DISC model.
///
/// The DISC algorithm converts GHI to DNI through empirical relationships
/// between the global and direct clearness indices.
///
/// # Parameters
/// - `ghi`: Global horizontal irradiance [W/m²]
/// - `solar_zenith`: True (not refraction-corrected) solar zenith angle [degrees]
/// - `day_of_year`: Day of year (1–365)
/// - `pressure`: Site pressure [Pa]. Use `None` for relative airmass only.
///
/// # References
/// Maxwell, E. L., 1987, "A Quasi-Physical Model for Converting Hourly
/// Global Horizontal to Direct Normal Insolation", Technical Report
/// No. SERI/TR-215-3087, Golden, CO: Solar Energy Research Institute.
pub fn disc(ghi: f64, solar_zenith: f64, day_of_year: i32, pressure: Option<f64>) -> DiscOutput {
    let max_zenith = 87.0;
    let min_cos_zenith = 0.065;

    // DISC uses solar constant = 1370 with Spencer 1971 full Fourier series
    let b = 2.0 * PI * ((day_of_year - 1) as f64) / 365.0;
    let rover = 1.00011 + 0.034221 * b.cos() + 0.00128 * b.sin()
        + 0.000719 * (2.0 * b).cos() + 0.000077 * (2.0 * b).sin();
    let i0 = 1370.0 * rover;

    // Clearness index
    let cos_z = solar_zenith.to_radians().cos().max(min_cos_zenith);
    let ghi_extra = i0 * cos_z;
    let kt = if ghi_extra > 0.0 { (ghi / ghi_extra).clamp(0.0, 1.0) } else { 0.0 };

    // Airmass — DISC was calibrated against Kasten 1966, not Kasten-Young 1989
    // Kasten 1966: AM = 1 / (cos(z) + 0.15 * (93.885 - z)^(-1.253))
    let mut am = {
        let z = solar_zenith;
        let cos_z = z.to_radians().cos();
        let c = 93.885 - z;
        if c <= 0.0 {
            f64::NAN
        } else {
            1.0 / (cos_z + 0.15 * c.powf(-1.253))
        }
    };
    if let Some(p) = pressure {
        am = atmosphere::get_absolute_airmass(am, p);
    }

    let (kn, am) = disc_kn(kt, am);
    let mut dni = kn * i0;

    if solar_zenith > max_zenith || ghi < 0.0 || dni < 0.0 {
        dni = 0.0;
    }

    DiscOutput { dni, kt, airmass: am }
}

/// Output of the Erbs-Driesse decomposition model.
#[derive(Debug, Clone, Copy)]
pub struct ErbsDriesseOutput {
    /// Direct normal irradiance [W/m²]
    pub dni: f64,
    /// Diffuse horizontal irradiance [W/m²]
    pub dhi: f64,
    /// Clearness index [unitless]
    pub kt: f64,
}

/// Estimate DNI and DHI from GHI using the continuous Erbs-Driesse model.
///
/// The Erbs-Driesse model is a reformulation of the original Erbs model
/// that provides continuity of the function and its first derivative at
/// the two transition points.
///
/// # Parameters
/// - `ghi`: Global horizontal irradiance [W/m²]
/// - `solar_zenith`: True (not refraction-corrected) zenith angle [degrees]
/// - `day_of_year`: Day of year (1–365)
///
/// # References
/// Driesse, A., Jensen, A., Perez, R., 2024. A Continuous form of the
/// Perez diffuse sky model for forward and reverse transposition.
/// Solar Energy vol. 267. doi:10.1016/j.solener.2023.112093
pub fn erbs_driesse(ghi: f64, solar_zenith: f64, day_of_year: i32) -> ErbsDriesseOutput {
    let max_zenith = 87.0;
    let min_cos_zenith = 0.065;

    let ghi = ghi.max(0.0);

    let dni_extra = get_extra_radiation(day_of_year);

    // Clearness index
    let cos_z = solar_zenith.to_radians().cos().max(min_cos_zenith);
    let ghi_extra = dni_extra * cos_z;
    let kt = if ghi_extra > 0.0 { (ghi / ghi_extra).clamp(0.0, 1.0) } else { 0.0 };

    // Central polynomial coefficients
    let p = [12.26911439571261, -16.4705084246973, 4.24692671521831700,
             -0.11390583806313881, 0.946296633571001];

    // Diffuse fraction
    let df = if kt <= 0.216 {
        1.0 - 0.09 * kt
    } else if kt <= 0.792 {
        // np.polyval evaluates p[0]*x^4 + p[1]*x^3 + ...
        p[0] * kt.powi(4) + p[1] * kt.powi(3) + p[2] * kt.powi(2) + p[3] * kt + p[4]
    } else {
        0.165
    };

    let dhi = df * ghi;
    let mut dni = (ghi - dhi) / solar_zenith.to_radians().cos();

    let bad = solar_zenith > max_zenith || ghi < 0.0 || dni < 0.0;
    let dhi = if bad { ghi } else { dhi };
    if bad {
        dni = 0.0;
    }

    ErbsDriesseOutput { dni, dhi, kt }
}

/// King diffuse sky model.
///
/// Determines the diffuse irradiance from the sky on a tilted surface using
/// the King model. Ground-reflected irradiance is not included.
///
/// # Parameters
/// - `surface_tilt`: Panel tilt from horizontal [degrees]
/// - `dhi`: Diffuse horizontal irradiance [W/m²]
/// - `ghi`: Global horizontal irradiance [W/m²]
/// - `solar_zenith`: Apparent (refraction-corrected) solar zenith angle [degrees]
pub fn king(surface_tilt: f64, dhi: f64, ghi: f64, solar_zenith: f64) -> f64 {
    let cos_tilt = surface_tilt.to_radians().cos();
    let sky_diffuse = dhi * (1.0 + cos_tilt) / 2.0
        + ghi * (0.012 * solar_zenith - 0.04) * (1.0 - cos_tilt) / 2.0;
    sky_diffuse.max(0.0)
}

/// DIRINDEX model for estimating DNI from GHI using clearsky information.
///
/// The DIRINDEX model modifies the DIRINT model by incorporating information
/// from a clear sky model. It computes:
/// `DNI = DNI_clear * DIRINT(GHI) / DIRINT(GHI_clear)`
///
/// # Parameters
/// - `ghi`: Global horizontal irradiance [W/m²]
/// - `ghi_clearsky`: Clear-sky global horizontal irradiance [W/m²]
/// - `dni_clearsky`: Clear-sky direct normal irradiance [W/m²]
/// - `zenith`: True (not refraction-corrected) zenith angle [degrees]
/// - `day_of_year`: Day of year (1–365)
/// - `pressure`: Site pressure [Pa]. Use `None` for standard pressure (101325 Pa).
///
/// # References
/// Perez, R., Ineichen, P., Moore, K., Kmiecik, M., Chain, C., George, R.,
/// & Vignola, F. (2002). A new operational model for satellite-derived
/// irradiances: description and validation. Solar Energy, 73(5), 307-317.
pub fn dirindex(
    ghi: f64,
    ghi_clearsky: f64,
    dni_clearsky: f64,
    zenith: f64,
    day_of_year: i32,
    pressure: Option<f64>,
) -> f64 {
    let dni_extra = get_extra_radiation(day_of_year);
    let p = pressure.unwrap_or(101325.0);

    let (dni_dirint, _) = dirint(ghi, zenith, 0.0, p, dni_extra);
    let (dni_dirint_clear, _) = dirint(ghi_clearsky, zenith, 0.0, p, dni_extra);

    if dni_dirint_clear <= 0.0 {
        return 0.0;
    }

    let dni = dni_clearsky * dni_dirint / dni_dirint_clear;
    dni.max(0.0)
}

// ---------------------------------------------------------------------------
// Perez-Driesse transposition model
// ---------------------------------------------------------------------------

/// Knot vector for the Perez-Driesse quadratic B-splines.
const PD_KNOTS: [f64; 13] = [
    0.000, 0.000, 0.000,
    0.061, 0.187, 0.333, 0.487, 0.643, 0.778, 0.839,
    1.000, 1.000, 1.000,
];

/// Coefficient table for the Perez-Driesse splines.
/// Original layout: 13 rows x 6 columns (f11,f12,f13, f21,f22,f23).
/// After transpose+reshape to (2,3,13), index as COEFS[i-1][j-1].
const PD_COEFS: [[[f64; 13]; 3]; 2] = [
    // i=1 (F1 coefficients)
    [
        // j=1: f11
        [-0.053, -0.008,  0.131,  0.328,  0.557,  0.861,  1.212,  1.099,  0.544,  0.544,  0.000,  0.000,  0.000],
        // j=2: f12
        [ 0.529,  0.588,  0.770,  0.471,  0.241, -0.323, -1.239, -1.847,  0.157,  0.157,  0.000,  0.000,  0.000],
        // j=3: f13
        [-0.028, -0.062, -0.167, -0.216, -0.300, -0.355, -0.444, -0.365, -0.213, -0.213,  0.000,  0.000,  0.000],
    ],
    // i=2 (F2 coefficients)
    [
        // j=1: f21
        [-0.071, -0.060, -0.026,  0.069,  0.086,  0.240,  0.305,  0.275,  0.118,  0.118,  0.000,  0.000,  0.000],
        // j=2: f22
        [ 0.061,  0.072,  0.106, -0.105, -0.085, -0.467, -0.797, -1.132, -1.455, -1.455,  0.000,  0.000,  0.000],
        // j=3: f23
        [-0.019, -0.022, -0.032, -0.028, -0.012, -0.008,  0.047,  0.124,  0.292,  0.292,  0.000,  0.000,  0.000],
    ],
];

/// Evaluate a quadratic B-spline defined by the Perez-Driesse knots and coefficients.
///
/// This is equivalent to `scipy.interpolate.splev(x, (knots, coefs, 2))`.
fn pd_splev(x: f64, coefs: &[f64; 13]) -> f64 {
    let t = &PD_KNOTS;
    let k = 2_usize; // quadratic
    let n = t.len() - k - 1; // 10 basis functions

    // Clamp x to knot domain [t[k], t[n]]
    let x = x.clamp(t[k], t[n]);

    // De Boor's algorithm for evaluating B-spline at x
    // Find knot span: largest i such that t[i] <= x < t[i+1], with i in [k, n-1]
    let mut span = k;
    for i in k..n {
        if t[i + 1] > x {
            span = i;
            break;
        }
        span = i;
    }

    // Initialize: d[j] = coefs[span - k + j] for j = 0..=k
    let mut d = [0.0_f64; 3]; // k+1 = 3
    for j in 0..=k {
        let idx = span - k + j;
        if idx < 13 {
            d[j] = coefs[idx];
        }
    }

    // Triangular computation
    for r in 1..=k {
        for j in (r..=k).rev() {
            let left = span + j - k;
            let right = span + 1 + j - r;
            let denom = t[right] - t[left];
            if denom.abs() < 1e-15 {
                d[j] = 0.0;
            } else {
                let alpha = (x - t[left]) / denom;
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
            }
        }
    }

    d[k]
}

/// Compute the delta parameter (sky brightness) for Perez-Driesse.
fn pd_calc_delta(dhi: f64, dni_extra: f64, solar_zenith: f64, airmass: Option<f64>) -> f64 {
    let am = match airmass {
        Some(a) => {
            if solar_zenith >= 90.0 {
                // Use max airmass at horizon
                atmosphere::get_relative_airmass(89.999)
            } else {
                a
            }
        }
        None => {
            if solar_zenith >= 90.0 {
                atmosphere::get_relative_airmass(89.999)
            } else {
                atmosphere::get_relative_airmass(solar_zenith)
            }
        }
    };

    let am = if am.is_nan() { atmosphere::get_relative_airmass(89.999) } else { am };

    if dni_extra <= 0.0 || am <= 0.0 {
        return 0.0;
    }

    dhi / (dni_extra / am)
}

/// Compute the zeta parameter (sky clearness) for Perez-Driesse.
fn pd_calc_zeta(dhi: f64, dni: f64, zenith: f64) -> f64 {
    if dhi <= 0.0 && dni <= 0.0 {
        return 0.0;
    }

    let sum = dhi + dni;
    let mut zeta = if sum > 0.0 { dni / sum } else { 0.0 };

    if dhi == 0.0 {
        zeta = 0.0;
    }

    // Apply kappa correction (analogous to eq. 7)
    let kappa = 1.041;
    let kterm = kappa * zenith.to_radians().powi(3);
    let denom = 1.0 - kterm * (zeta - 1.0);
    if denom.abs() > 1e-15 {
        zeta /= denom;
    }

    zeta
}

/// Evaluate the Perez-Driesse spline function f(i,j,zeta).
fn pd_f(i: usize, j: usize, zeta: f64) -> f64 {
    pd_splev(zeta, &PD_COEFS[i - 1][j - 1])
}

/// Continuous Perez-Driesse diffuse sky model.
///
/// The Perez-Driesse model is a reformulation of the 1990 Perez model
/// that provides continuity of the function and of its first derivatives
/// by replacing the look-up table of coefficients with quadratic splines.
///
/// # Parameters
/// - `surface_tilt`: Panel tilt from horizontal [degrees]
/// - `surface_azimuth`: Panel azimuth [degrees]
/// - `dhi`: Diffuse horizontal irradiance [W/m²]
/// - `dni`: Direct normal irradiance [W/m²]
/// - `dni_extra`: Extraterrestrial normal irradiance [W/m²]
/// - `solar_zenith`: Apparent (refraction-corrected) zenith angle [degrees]
/// - `solar_azimuth`: Solar azimuth angle [degrees]
/// - `airmass`: Relative (not pressure-corrected) airmass [unitless].
///   If `None`, calculated internally using Kasten-Young 1989.
///
/// # References
/// Driesse, A., Jensen, A., Perez, R., 2024. A Continuous form of the
/// Perez diffuse sky model for forward and reverse transposition.
/// Solar Energy vol. 267. doi:10.1016/j.solener.2023.112093
#[allow(clippy::too_many_arguments)]
pub fn perez_driesse(
    surface_tilt: f64,
    surface_azimuth: f64,
    dhi: f64,
    dni: f64,
    dni_extra: f64,
    solar_zenith: f64,
    solar_azimuth: f64,
    airmass: Option<f64>,
) -> f64 {
    let delta = pd_calc_delta(dhi, dni_extra, solar_zenith, airmass);
    let zeta = pd_calc_zeta(dhi, dni, solar_zenith);

    let z = solar_zenith.to_radians();

    let f1 = pd_f(1, 1, zeta) + pd_f(1, 2, zeta) * delta + pd_f(1, 3, zeta) * z;
    let f2 = pd_f(2, 1, zeta) + pd_f(2, 2, zeta) * delta + pd_f(2, 3, zeta) * z;

    // Clip F1 to [0, 0.9] as recommended
    let f1 = f1.clamp(0.0, 0.9);

    // A = max(cos(AOI), 0)
    let a = aoi_projection(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth).max(0.0);

    // B = max(cos(zenith), cos(85))
    let b = solar_zenith.to_radians().cos().max(85.0_f64.to_radians().cos());

    let term1 = 0.5 * (1.0 - f1) * (1.0 + surface_tilt.to_radians().cos());
    let term2 = f1 * a / b;
    let term3 = f2 * surface_tilt.to_radians().sin();

    (dhi * (term1 + term2 + term3)).max(0.0)
}

