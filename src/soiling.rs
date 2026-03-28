/// Simplified soiling loss mass accumulation model (inspired by Hsu et al).
///
/// Estimates daily soiling accumulation based on particulate matter concentrations,
/// assuming natural rain washing.
///
/// # Arguments
/// * `rainfall` - Daily rainfall in mm.
/// * `cleaning_threshold` - Rainfall threshold to clean the panels (e.g., 5.0 mm).
/// * `tilt` - Surface tilt in degrees.
/// * `pm2_5` - PM 2.5 concentration in ug/m^3.
/// * `pm10` - PM 10 concentration in ug/m^3.
///
/// # Returns
/// Soiling mass accumulation fraction for the day.
pub fn accumulation_model(rainfall: f64, cleaning_threshold: f64, tilt: f64, pm2_5: f64, pm10: f64) -> f64 {
    if rainfall >= cleaning_threshold {
        return 0.0; // Cleaned by rain
    }

    // Rough empirical accumulation rate combined with tilt gravity slide-off
    let accumulation_rate = (0.001 * pm2_5 + 0.002 * pm10) / 100.0;

    // Closer to horizontal accumulates more
    let tilt_factor = tilt.clamp(0.0, 90.0).to_radians().cos();

    accumulation_rate * tilt_factor
}

/// HSU soiling model - single timestep mass accumulation update.
///
/// Computes the soiling ratio using the Fixed Velocity model from Humboldt State
/// University (HSU). The soiling ratio ranges from 0 to 1, where 1 means no soiling.
///
/// This is a single-timestep function. To simulate over time, accumulate mass across
/// timesteps, resetting to 0 when rainfall exceeds the cleaning threshold, then
/// convert accumulated mass to soiling ratio with `hsu_soiling_ratio`.
///
/// # Arguments
/// * `pm2_5` - PM2.5 concentration [g/m^3].
/// * `pm10` - PM10 concentration [g/m^3].
/// * `depo_veloc_2_5` - Deposition velocity for PM2.5 [m/s], default 0.0009.
/// * `depo_veloc_10` - Deposition velocity for PM10 [m/s], default 0.004.
/// * `surface_tilt` - Module tilt from horizontal [degrees].
/// * `dt_sec` - Timestep duration [seconds].
///
/// # Returns
/// Mass deposited on the tilted surface during this timestep [g/m^2].
///
/// # References
/// Coello, M. and Boyle, L. (2019). "Simple Model For Predicting Time Series
/// Soiling of Photovoltaic Panels." IEEE Journal of Photovoltaics.
pub fn hsu_mass_rate(
    pm2_5: f64,
    pm10: f64,
    depo_veloc_2_5: f64,
    depo_veloc_10: f64,
    surface_tilt: f64,
    dt_sec: f64,
) -> f64 {
    let horiz_mass = (pm2_5 * depo_veloc_2_5 + (pm10 - pm2_5).max(0.0) * depo_veloc_10) * dt_sec;
    horiz_mass * surface_tilt.to_radians().cos()
}

/// Converts accumulated soiling mass to a soiling ratio using the HSU model.
///
/// soiling_ratio = 1 - 0.3437 * erf(0.17 * accum_mass^0.8473)
///
/// # Arguments
/// * `accumulated_mass` - Total accumulated soiling mass on the surface [g/m^2].
///
/// # Returns
/// Soiling ratio (0 to 1), where 1 means perfectly clean.
pub fn hsu_soiling_ratio(accumulated_mass: f64) -> f64 {
    // Approximate erf using Abramowitz and Stegun formula 7.1.26
    let x = 0.17 * accumulated_mass.powf(0.8473);
    1.0 - 0.3437 * erf_approx(x)
}

/// HSU soiling model convenience function for a single timestep.
///
/// Returns the soiling ratio given current accumulated mass and whether
/// rainfall has cleaned the panels.
///
/// # Arguments
/// * `rainfall` - Rainfall accumulated in the current period [mm].
/// * `cleaning_threshold` - Rainfall needed to clean panels [mm].
/// * `surface_tilt` - Module tilt [degrees].
/// * `pm2_5` - PM2.5 concentration [g/m^3].
/// * `pm10` - PM10 concentration [g/m^3].
/// * `depo_veloc_2_5` - PM2.5 deposition velocity [m/s], default 0.0009.
/// * `depo_veloc_10` - PM10 deposition velocity [m/s], default 0.004.
/// * `dt_sec` - Timestep duration [seconds].
/// * `previous_mass` - Accumulated mass from previous timestep [g/m^2].
///
/// # Returns
/// A tuple of (soiling_ratio, new_accumulated_mass).
#[allow(clippy::too_many_arguments)]
pub fn hsu(
    rainfall: f64,
    cleaning_threshold: f64,
    surface_tilt: f64,
    pm2_5: f64,
    pm10: f64,
    depo_veloc_2_5: f64,
    depo_veloc_10: f64,
    dt_sec: f64,
    previous_mass: f64,
) -> (f64, f64) {
    let is_cleaned = rainfall >= cleaning_threshold;
    let mass_this_step = hsu_mass_rate(pm2_5, pm10, depo_veloc_2_5, depo_veloc_10, surface_tilt, dt_sec);

    let accum_mass = if is_cleaned {
        mass_this_step // reset after cleaning
    } else {
        previous_mass + mass_this_step
    };

    (hsu_soiling_ratio(accum_mass), accum_mass)
}

/// Kimber soiling model - single timestep update.
///
/// Linear soiling accumulation with rain cleaning resets. Soiling builds up
/// at a daily rate unless rainfall exceeds the cleaning threshold.
///
/// # Arguments
/// * `rainfall` - Rainfall accumulated in the current period [mm].
/// * `cleaning_threshold` - Daily rainfall needed to clean panels [mm], default 6.0.
/// * `soiling_loss_rate` - Fraction of energy lost per day of soiling [unitless], default 0.0015.
/// * `max_soiling` - Maximum soiling loss fraction [unitless], default 0.3.
/// * `timestep_days` - Duration of the current timestep as a fraction of a day.
/// * `previous_soiling` - Soiling loss from previous timestep (0 to max_soiling).
/// * `is_grace_period` - True if within the grace period after a rain event (ground too damp for soiling).
///
/// # Returns
/// Updated soiling loss fraction (0.0 to max_soiling).
///
/// # References
/// Kimber, A. et al. (2006). "The Effect of Soiling on Large Grid-Connected
/// Photovoltaic Systems in California and the Southwest Region of the United States."
/// IEEE 4th World Conference on Photovoltaic Energy Conference.
pub fn kimber(
    rainfall: f64,
    cleaning_threshold: f64,
    soiling_loss_rate: f64,
    max_soiling: f64,
    timestep_days: f64,
    previous_soiling: f64,
    is_grace_period: bool,
) -> f64 {
    // Rain event cleans the panels
    if rainfall >= cleaning_threshold {
        return 0.0;
    }

    // During grace period (ground damp), no soiling accumulation
    if is_grace_period {
        return 0.0;
    }

    // Accumulate soiling
    let new_soiling = previous_soiling + soiling_loss_rate * timestep_days;
    new_soiling.min(max_soiling)
}

/// Approximate error function using Abramowitz and Stegun formula 7.1.26.
/// Maximum error ~1.5e-7.
fn erf_approx(x: f64) -> f64 {
    let sign = if x >= 0.0 { 1.0 } else { -1.0 };
    let x = x.abs();

    let p = 0.3275911;
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;

    let t = 1.0 / (1.0 + p * x);
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    let result = 1.0 - (a1 * t + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5) * (-x * x).exp();
    sign * result
}
