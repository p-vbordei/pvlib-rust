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

// ---------------------------------------------------------------------------
// Reno–Hansen clear-sky detection (ported from pvlib-python)
// ---------------------------------------------------------------------------

/// Thresholds for [`detect_clearsky`] (Reno–Hansen 2016).
///
/// Default values match pvlib-python defaults for 10-minute windows of
/// 1-minute GHI data. Use [`ClearSkyThresholds::from_sample_interval`] to
/// get the Jordan–Hansen (2023) values interpolated for a specific
/// sample interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClearSkyThresholds {
    /// Centered sliding-window length, in minutes.
    pub window_length_minutes: f64,
    /// Max |mean(meas) - α·mean(clear)| per window [W/m²].
    pub mean_diff: f64,
    /// Max |max(meas) - α·max(clear)| per window [W/m²].
    pub max_diff: f64,
    /// Lower bound on `line_length(meas) - line_length(α·clear)`.
    pub lower_line_length: f64,
    /// Upper bound on `line_length(meas) - line_length(α·clear)`.
    pub upper_line_length: f64,
    /// Upper bound on the normalised std-dev of slopes (Hz⁻¹).
    pub var_diff: f64,
    /// Upper bound on the max |Δ(meas − α·clear)| per window.
    pub slope_dev: f64,
    /// Maximum iterations for the α rescaling fixed-point.
    pub max_iterations: usize,
}

impl Default for ClearSkyThresholds {
    fn default() -> Self {
        // Reno-Hansen 2016 defaults for 10-min window / 1-min data.
        Self {
            window_length_minutes: 10.0,
            mean_diff: 75.0,
            max_diff: 75.0,
            lower_line_length: -5.0,
            upper_line_length: 10.0,
            var_diff: 0.005,
            slope_dev: 8.0,
            max_iterations: 20,
        }
    }
}

impl ClearSkyThresholds {
    /// Interpolate thresholds for a given sample interval (Jordan &
    /// Hansen 2023, Table 1). Valid for `sample_interval_minutes ∈ [1, 30]`.
    pub fn from_sample_interval(sample_interval_minutes: f64) -> Self {
        let si = sample_interval_minutes.clamp(1.0, 30.0);
        let interp = |xs: &[f64], ys: &[f64]| -> f64 {
            // linear interpolation over 4 breakpoints [1, 5, 15, 30]
            debug_assert_eq!(xs.len(), ys.len());
            if si <= xs[0] {
                return ys[0];
            }
            for i in 0..xs.len() - 1 {
                if si <= xs[i + 1] {
                    let t = (si - xs[i]) / (xs[i + 1] - xs[i]);
                    return ys[i] + t * (ys[i + 1] - ys[i]);
                }
            }
            *ys.last().unwrap()
        };
        let breakpoints = [1.0, 5.0, 15.0, 30.0];
        Self {
            window_length_minutes: interp(&breakpoints, &[50.0, 60.0, 90.0, 120.0]),
            mean_diff: 75.0,
            max_diff: interp(&breakpoints, &[60.0, 65.0, 75.0, 90.0]),
            lower_line_length: -45.0,
            upper_line_length: 80.0,
            var_diff: interp(&breakpoints, &[0.005, 0.01, 0.032, 0.07]),
            slope_dev: interp(&breakpoints, &[50.0, 60.0, 75.0, 96.0]),
            max_iterations: 20,
        }
    }
}

/// Full result from [`detect_clearsky_detail`]: per-sample clear flag plus
/// the fitted α rescaling factor.
#[derive(Debug, Clone)]
pub struct ClearSkyDetectionResult {
    pub clear_samples: Vec<bool>,
    pub alpha: f64,
    pub iterations: usize,
}

/// Detect clear-sky periods using the **Reno–Hansen (2016) windowed
/// 5-criterion algorithm**.
///
/// Faithful port of `pvlib.clearsky.detect_clearsky`. Samples that fall
/// inside any centered window passing the five criteria (mean diff, max
/// diff, line-length diff, slope-n-std, slope-dev) are flagged as clear.
/// The algorithm iterates an α rescaling of the clear-sky series so
/// that the detected clear samples best match the clear-sky model.
///
/// # Parameters
/// - `measured`, `clearsky` — parallel time-series samples of measured
///   GHI and expected clearsky GHI [W/m²]. Must have equal length and be
///   equally spaced in time at `sample_interval_minutes`.
/// - `sample_interval_minutes` — spacing between consecutive samples
///   (e.g. 1.0 for 1-minute data, 60.0 for hourly).
/// - `thresholds` — detection thresholds; [`ClearSkyThresholds::default`]
///   matches pvlib-python defaults. For non-default sample intervals,
///   use [`ClearSkyThresholds::from_sample_interval`] to get
///   Jordan–Hansen 2023 values.
///
/// # Panics
///
/// Panics if `measured.len() != clearsky.len()` or if the window would
/// contain fewer than 3 samples.
///
/// # References
/// - Reno, M.J. and Hansen, C.W., 2016. "Identification of periods of
///   clear sky irradiance in time series of GHI measurements."
///   Renewable Energy, v90, p. 520-531.
/// - Jordan, D.C. and Hansen, C., 2023. "Clear-sky detection for PV
///   degradation analysis using multiple regression." Renewable Energy,
///   v209, p. 393-400.
pub fn detect_clearsky(
    measured: &[f64],
    clearsky: &[f64],
    sample_interval_minutes: f64,
    thresholds: ClearSkyThresholds,
) -> Vec<bool> {
    detect_clearsky_detail(measured, clearsky, sample_interval_minutes, thresholds).clear_samples
}

/// Same as [`detect_clearsky`] but also returns the fitted α and the
/// iteration count. Useful for diagnostics or for chaining into a second
/// pass with the fitted scaling.
pub fn detect_clearsky_detail(
    measured: &[f64],
    clearsky: &[f64],
    sample_interval_minutes: f64,
    thresholds: ClearSkyThresholds,
) -> ClearSkyDetectionResult {
    let n = measured.len();
    assert_eq!(clearsky.len(), n, "detect_clearsky: length mismatch");

    let samples_per_window =
        (thresholds.window_length_minutes / sample_interval_minutes).round() as usize;
    assert!(
        samples_per_window >= 3,
        "detect_clearsky: samples_per_window ({samples_per_window}) must be ≥ 3; \
         increase window_length_minutes or reduce sample_interval_minutes"
    );

    if n < samples_per_window {
        return ClearSkyDetectionResult {
            clear_samples: vec![false; n],
            alpha: 1.0,
            iterations: 0,
        };
    }

    let num_windows = n + 1 - samples_per_window;

    // Per-window statistics for the measured series: mean, max, slope-n-std,
    // line length, and (for slope-dev) the raw sample-to-sample forward diffs.
    let meas_mean = window_mean(measured, samples_per_window, num_windows);
    let meas_max = window_max(measured, samples_per_window, num_windows);
    let meas_slope_nstd = window_slope_nstd(
        measured,
        sample_interval_minutes,
        samples_per_window,
        num_windows,
        &meas_mean,
    );
    let meas_line_length = window_line_length(
        measured,
        sample_interval_minutes,
        samples_per_window,
        num_windows,
    );

    // Clearsky window stats (α applied inside the loop).
    let clear_mean = window_mean(clearsky, samples_per_window, num_windows);
    let clear_max = window_max(clearsky, samples_per_window, num_windows);

    let mut alpha: f64 = 1.0;
    let mut clear_windows_flags = vec![false; num_windows];
    let mut iterations = 0;

    for it in 0..thresholds.max_iterations {
        iterations = it + 1;

        // Scaled clearsky stats; line length of α·clear equals α·line-length
        // only in the limit of negligible sample_interval — recompute fully.
        let scaled_clear: Vec<f64> = clearsky.iter().map(|v| alpha * v).collect();
        let clear_line_length = window_line_length(
            &scaled_clear,
            sample_interval_minutes,
            samples_per_window,
            num_windows,
        );
        let residual: Vec<f64> = measured
            .iter()
            .zip(&scaled_clear)
            .map(|(m, c)| m - c)
            .collect();
        let slope_max_diff = window_max_abs_diff(
            &residual,
            samples_per_window,
            num_windows,
        );

        for w in 0..num_windows {
            let line_diff = meas_line_length[w] - clear_line_length[w];
            let c1 = (meas_mean[w] - alpha * clear_mean[w]).abs() < thresholds.mean_diff;
            let c2 = (meas_max[w] - alpha * clear_max[w]).abs() < thresholds.max_diff;
            let c3 = line_diff > thresholds.lower_line_length
                && line_diff < thresholds.upper_line_length;
            let c4 = meas_slope_nstd[w] < thresholds.var_diff;
            let c5 = slope_max_diff[w] < thresholds.slope_dev;
            let c6 = clear_mean[w] != 0.0 && !clear_mean[w].is_nan();
            clear_windows_flags[w] = c1 && c2 && c3 && c4 && c5 && c6;
        }

        // Propagate window flags to samples — any sample inside a clear
        // window is clear.
        let mut sample_clear = vec![false; n];
        for w in 0..num_windows {
            if clear_windows_flags[w] {
                for j in 0..samples_per_window {
                    sample_clear[w + j] = true;
                }
            }
        }

        // Refit α over the identified clear samples.
        let mut num = 0.0;
        let mut den = 0.0;
        for i in 0..n {
            if sample_clear[i] && clearsky[i].is_finite() && measured[i].is_finite() {
                num += measured[i] * clearsky[i];
                den += clearsky[i] * clearsky[i];
            }
        }
        let previous_alpha = alpha;
        if den.is_finite() && den > 0.0 {
            alpha = num / den;
        }
        if (alpha * 10_000.0).round() == (previous_alpha * 10_000.0).round() {
            break;
        }
    }

    // Final sample-clear propagation after convergence.
    let mut clear_samples = vec![false; n];
    for w in 0..num_windows {
        if clear_windows_flags[w] {
            for j in 0..samples_per_window {
                clear_samples[w + j] = true;
            }
        }
    }

    ClearSkyDetectionResult {
        clear_samples,
        alpha,
        iterations,
    }
}

/// Window mean over a sliding window of `w` consecutive samples. Returns
/// `num_windows` values (one per starting position).
fn window_mean(data: &[f64], w: usize, num_windows: usize) -> Vec<f64> {
    // Prefix-sum for O(n) mean; NaN samples contribute NaN to the mean.
    let mut out = Vec::with_capacity(num_windows);
    for start in 0..num_windows {
        let mut s = 0.0_f64;
        for j in 0..w {
            s += data[start + j];
        }
        out.push(s / w as f64);
    }
    out
}

fn window_max(data: &[f64], w: usize, num_windows: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(num_windows);
    for start in 0..num_windows {
        let mut m = f64::NEG_INFINITY;
        for j in 0..w {
            let v = data[start + j];
            if v > m {
                m = v;
            }
        }
        out.push(m);
    }
    out
}

/// Σ √(Δ² + dt²) over each window (i.e. arc length on a regular-interval
/// time series). `w` samples per window → `w - 1` differences per window.
fn window_line_length(
    data: &[f64],
    sample_interval: f64,
    w: usize,
    num_windows: usize,
) -> Vec<f64> {
    let mut out = Vec::with_capacity(num_windows);
    let dt2 = sample_interval * sample_interval;
    for start in 0..num_windows {
        let mut s = 0.0_f64;
        for j in 0..w - 1 {
            let d = data[start + j + 1] - data[start + j];
            s += (d * d + dt2).sqrt();
        }
        out.push(s);
    }
    out
}

/// Standard deviation (sample, ddof=1) of `diff(data)/sample_interval`
/// within each window, normalised by the window mean of `data`
/// (matching pvlib-python's `_slope_nstd_windowed`).
fn window_slope_nstd(
    data: &[f64],
    sample_interval: f64,
    w: usize,
    num_windows: usize,
    window_means: &[f64],
) -> Vec<f64> {
    let mut out = Vec::with_capacity(num_windows);
    let n_slopes = w - 1; // per window
    let denom_ddof = (n_slopes - 1) as f64; // sample std
    for start in 0..num_windows {
        let mut sum_s = 0.0_f64;
        let mut sum_sq = 0.0_f64;
        for j in 0..n_slopes {
            let slope = (data[start + j + 1] - data[start + j]) / sample_interval;
            sum_s += slope;
            sum_sq += slope * slope;
        }
        let mean_s = sum_s / n_slopes as f64;
        let var = (sum_sq - n_slopes as f64 * mean_s * mean_s) / denom_ddof;
        let std = if var > 0.0 { var.sqrt() } else { 0.0 };
        let wm = window_means[start];
        out.push(if wm.abs() > 0.0 { std / wm } else { f64::NAN });
    }
    out
}

/// Max |diff(data)| within each window (= `w - 1` successive-sample deltas).
fn window_max_abs_diff(data: &[f64], w: usize, num_windows: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(num_windows);
    for start in 0..num_windows {
        let mut m = 0.0_f64;
        for j in 0..w - 1 {
            let d = (data[start + j + 1] - data[start + j]).abs();
            if d > m {
                m = d;
            }
        }
        out.push(m);
    }
    out
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
