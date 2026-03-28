# Batch API Improvements Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Improve the pvlib-rust library's batch API based on real-world friction points discovered in the rust-solar-worker production consumer, eliminating common workarounds and code duplication.

**Architecture:** All changes target `src/batch.rs` (the batch processing module) and `tests/test_batch.rs`. The `SimulationSeries` struct gains a `solar_elevation` field. `WeatherSeries` gains a `from_utc` constructor and supports `Option<Vec<f64>>` for `dni`/`dhi` to enable auto-decomposition. `BatchModelChain` gains builder methods for bifacial support, system losses, and auto-decomposition. A new `solar_position_batch_utc` convenience function is added.

**Tech Stack:** Rust, chrono, chrono_tz, rayon, spa crate

---

### Task 1: Add `solar_elevation` to `SimulationSeries`

**Files:**
- Modify: `src/batch.rs` (struct `SimulationSeries` at line 309, method `run` at line 487)
- Modify: `tests/test_batch.rs` (all tests constructing `SimulationSeries`)

- [ ] **Step 1: Write the failing test**

Add to `tests/test_batch.rs`:

```rust
#[test]
fn test_simulation_series_has_elevation() {
    let loc = tucson();
    let chain = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0);

    let weather = batch::WeatherSeries {
        times: vec![
            Eastern.with_ymd_and_hms(2020, 6, 15, 10, 0, 0).unwrap(),
            summer_noon(),
            Eastern.with_ymd_and_hms(2020, 6, 15, 14, 0, 0).unwrap(),
        ],
        ghi: vec![600.0, 900.0, 700.0],
        dni: vec![500.0, 800.0, 600.0],
        dhi: vec![100.0, 100.0, 100.0],
        temp_air: vec![28.0, 32.0, 30.0],
        wind_speed: vec![2.0, 1.5, 2.5],
        albedo: None,
    };

    let result = chain.run(&weather).unwrap();

    assert_eq!(result.solar_elevation.len(), 3);
    // Elevation + zenith = 90
    for i in 0..3 {
        let sum = result.solar_elevation[i] + result.solar_zenith[i];
        assert!((sum - 90.0).abs() < 1e-10,
            "elevation + zenith should be 90, got {}", sum);
    }
    // Noon elevation should be highest
    assert!(result.solar_elevation[1] > result.solar_elevation[0]);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test test_simulation_series_has_elevation -- --nocapture 2>&1 | head -20`
Expected: FAIL — `SimulationSeries` has no field `solar_elevation`

- [ ] **Step 3: Add `solar_elevation` field to `SimulationSeries` and populate it in `run()`**

In `src/batch.rs`, add the field to the struct:

```rust
pub struct SimulationSeries {
    pub solar_zenith: Vec<f64>,
    pub solar_azimuth: Vec<f64>,
    pub solar_elevation: Vec<f64>,  // NEW
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
```

In `BatchModelChain::run()`, populate it after collecting results:

```rust
Ok(SimulationSeries {
    solar_zenith: results.iter().map(|r| r.0).collect(),
    solar_azimuth: results.iter().map(|r| r.1).collect(),
    solar_elevation: results.iter().map(|r| 90.0 - r.0).collect(),  // NEW
    airmass: results.iter().map(|r| r.2).collect(),
    // ... rest unchanged
})
```

Update ALL existing tests in `tests/test_batch.rs` that construct `SimulationSeries` manually (the `test_simulation_series_helpers` and `test_simulation_series_edge_cases` tests) to include the new field:

```rust
// In test_simulation_series_helpers:
let series = batch::SimulationSeries {
    solar_zenith: vec![30.0, 40.0, 50.0],
    solar_azimuth: vec![180.0, 190.0, 200.0],
    solar_elevation: vec![60.0, 50.0, 40.0],  // NEW: 90 - zenith
    airmass: vec![1.1, 1.3, 1.5],
    // ... rest unchanged
};

// In test_simulation_series_edge_cases:
let empty = batch::SimulationSeries {
    solar_zenith: vec![],
    solar_azimuth: vec![],
    solar_elevation: vec![],  // NEW
    airmass: vec![],
    // ... rest unchanged
};
```

- [ ] **Step 4: Run all batch tests to verify they pass**

Run: `cargo test test_batch -- --nocapture 2>&1 | tail -5`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/batch.rs tests/test_batch.rs
git commit -m "feat: add solar_elevation field to SimulationSeries"
```

---

### Task 2: Add `solar_position_batch_utc` convenience function

**Files:**
- Modify: `src/batch.rs` (add new function after `solar_position_batch`)
- Modify: `tests/test_batch.rs` (add tests)

- [ ] **Step 1: Write the failing test**

Add to `tests/test_batch.rs`:

```rust
#[test]
fn test_solar_position_batch_utc() {
    use chrono::NaiveDate;

    let timestamps = vec![
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(14, 0, 0).unwrap(),
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(16, 0, 0).unwrap(),
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(18, 0, 0).unwrap(),
    ];

    let (zenith, azimuth, elevation) = batch::solar_position_batch_utc(
        32.2, -110.9, 700.0, &timestamps,
    ).unwrap();

    assert_eq!(zenith.len(), 3);
    assert_eq!(azimuth.len(), 3);
    assert_eq!(elevation.len(), 3);

    for i in 0..3 {
        let sum = zenith[i] + elevation[i];
        assert!((sum - 90.0).abs() < 1e-10,
            "elevation + zenith should be 90 at index {}", i);
    }

    // All should be daytime (UTC 14-18 = local morning/midday in Tucson)
    for z in &zenith {
        assert!(*z < 90.0, "zenith {} should be < 90 for daytime", z);
    }
}

#[test]
fn test_solar_position_batch_utc_matches_tz_version() {
    use chrono::{NaiveDate, TimeZone};

    let naive_ts = vec![
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(16, 0, 0).unwrap(),
    ];

    let (utc_zen, utc_az, utc_el) = batch::solar_position_batch_utc(
        32.2, -110.9, 700.0, &naive_ts,
    ).unwrap();

    // Same moment via Location + DateTime<Tz>
    let loc = tucson();
    let tz_time = chrono::Utc.from_utc_datetime(&naive_ts[0])
        .with_timezone(&Eastern);
    let (tz_zen, tz_az, tz_el) = batch::solar_position_batch(&loc, &[tz_time]).unwrap();

    assert!((utc_zen[0] - tz_zen[0]).abs() < 1e-6,
        "zenith mismatch: utc={} tz={}", utc_zen[0], tz_zen[0]);
    assert!((utc_az[0] - tz_az[0]).abs() < 1e-6,
        "azimuth mismatch: utc={} tz={}", utc_az[0], tz_az[0]);
    assert!((utc_el[0] - tz_el[0]).abs() < 1e-6,
        "elevation mismatch: utc={} tz={}", utc_el[0], tz_el[0]);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test test_solar_position_batch_utc -- --nocapture 2>&1 | head -20`
Expected: FAIL — `solar_position_batch_utc` not found

- [ ] **Step 3: Implement `solar_position_batch_utc`**

Add to `src/batch.rs` right after `solar_position_batch`:

```rust
/// Batch solar position from UTC NaiveDateTime timestamps.
///
/// Convenience function that accepts NaiveDateTime (assumed UTC) instead of
/// requiring callers to construct DateTime<Tz>. This is the common case for
/// weather API data which is typically timestamped in UTC.
///
/// Returns (zenith_vec, azimuth_vec, elevation_vec).
pub fn solar_position_batch_utc(
    latitude: f64,
    longitude: f64,
    altitude: f64,
    timestamps: &[chrono::NaiveDateTime],
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), spa::SpaError> {
    let location = crate::location::Location::new(
        latitude, longitude, chrono_tz::UTC, altitude, "",
    );
    let times: Vec<chrono::DateTime<chrono_tz::Tz>> = timestamps
        .iter()
        .map(|ndt| {
            chrono::Utc
                .from_utc_datetime(ndt)
                .with_timezone(&chrono_tz::UTC)
        })
        .collect();
    solar_position_batch(&location, &times)
}
```

Note: requires adding `use chrono::TimeZone;` to the top of `src/batch.rs` if not already present.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test test_solar_position_batch_utc -- --nocapture 2>&1 | tail -5`
Expected: Both tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/batch.rs tests/test_batch.rs
git commit -m "feat: add solar_position_batch_utc for NaiveDateTime input"
```

---

### Task 3: Add `WeatherSeries::from_utc` constructor

**Files:**
- Modify: `src/batch.rs` (add `impl WeatherSeries` block)
- Modify: `tests/test_batch.rs` (add test)

- [ ] **Step 1: Write the failing test**

Add to `tests/test_batch.rs`:

```rust
#[test]
fn test_weather_series_from_utc() {
    use chrono::NaiveDate;

    let timestamps = vec![
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(14, 0, 0).unwrap(),
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(15, 0, 0).unwrap(),
    ];

    let ws = batch::WeatherSeries::from_utc(
        &timestamps,
        "US/Eastern",
        vec![600.0, 700.0],
        vec![500.0, 600.0],
        vec![100.0, 100.0],
        vec![28.0, 30.0],
        vec![2.0, 2.5],
    ).unwrap();

    assert_eq!(ws.times.len(), 2);
    assert_eq!(ws.ghi.len(), 2);
    assert_eq!(ws.dni.len(), 2);
    assert_eq!(ws.dhi.len(), 2);
    assert_eq!(ws.temp_air.len(), 2);
    assert_eq!(ws.wind_speed.len(), 2);
    assert!(ws.albedo.is_none());
}

#[test]
fn test_weather_series_from_utc_invalid_tz() {
    use chrono::NaiveDate;

    let timestamps = vec![
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(14, 0, 0).unwrap(),
    ];

    let result = batch::WeatherSeries::from_utc(
        &timestamps,
        "Invalid/Timezone",
        vec![600.0],
        vec![500.0],
        vec![100.0],
        vec![28.0],
        vec![2.0],
    );

    assert!(result.is_err());
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test test_weather_series_from_utc -- --nocapture 2>&1 | head -20`
Expected: FAIL — `from_utc` method not found

- [ ] **Step 3: Implement `WeatherSeries::from_utc`**

Add to `src/batch.rs` after the `WeatherSeries` struct definition:

```rust
impl WeatherSeries {
    /// Create a WeatherSeries from UTC NaiveDateTime timestamps and a timezone name.
    ///
    /// Converts each NaiveDateTime (assumed UTC) to DateTime<Tz> in the specified timezone.
    /// This is the common case when working with weather API data timestamped in UTC.
    ///
    /// # Arguments
    /// * `timestamps` - NaiveDateTime values assumed to be UTC
    /// * `tz_name` - IANA timezone name (e.g., "US/Eastern", "Europe/Bucharest")
    /// * `ghi` - Global Horizontal Irradiance [W/m²]
    /// * `dni` - Direct Normal Irradiance [W/m²]
    /// * `dhi` - Diffuse Horizontal Irradiance [W/m²]
    /// * `temp_air` - Ambient temperature [°C]
    /// * `wind_speed` - Wind speed [m/s]
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
            .map_err(|_| format!("Invalid timezone: '{}'", tz_name))?;

        let times: Vec<chrono::DateTime<chrono_tz::Tz>> = timestamps
            .iter()
            .map(|ndt| {
                chrono::Utc
                    .from_utc_datetime(ndt)
                    .with_timezone(&tz)
            })
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test test_weather_series_from_utc -- --nocapture 2>&1 | tail -5`
Expected: Both tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/batch.rs tests/test_batch.rs
git commit -m "feat: add WeatherSeries::from_utc constructor for UTC timestamps"
```

---

### Task 4: Add auto GHI-to-DNI/DHI decomposition in `BatchModelChain`

**Files:**
- Modify: `src/batch.rs` (`BatchModelChain` struct + `run` method)
- Modify: `tests/test_batch.rs`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_batch.rs`:

```rust
#[test]
fn test_batch_modelchain_auto_decomposition() {
    let loc = tucson();
    let chain = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0)
        .with_auto_decomposition(true);

    // Provide GHI only (DNI/DHI all zeros — simulating weather API that only gives GHI)
    let weather = batch::WeatherSeries {
        times: vec![
            Eastern.with_ymd_and_hms(2020, 6, 15, 10, 0, 0).unwrap(),
            summer_noon(),
            Eastern.with_ymd_and_hms(2020, 6, 15, 14, 0, 0).unwrap(),
        ],
        ghi: vec![600.0, 900.0, 700.0],
        dni: vec![0.0, 0.0, 0.0],
        dhi: vec![0.0, 0.0, 0.0],
        temp_air: vec![28.0, 32.0, 30.0],
        wind_speed: vec![2.0, 1.5, 2.5],
        albedo: None,
    };

    let result = chain.run(&weather).unwrap();

    // With auto-decomposition, should still produce positive power
    for i in 0..3 {
        assert!(result.ac_power[i] > 0.0,
            "ac_power[{}] should be > 0 with auto-decomposition, got {}", i, result.ac_power[i]);
        assert!(result.poa_global[i] > 0.0,
            "poa_global[{}] should be > 0 with auto-decomposition, got {}", i, result.poa_global[i]);
    }
}

#[test]
fn test_batch_modelchain_no_decomposition_zero_power() {
    // Without auto-decomposition, zero DNI/DHI should produce near-zero POA direct
    let loc = tucson();
    let chain = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0);

    let weather = batch::WeatherSeries {
        times: vec![summer_noon()],
        ghi: vec![900.0],
        dni: vec![0.0],
        dhi: vec![0.0],
        temp_air: vec![32.0],
        wind_speed: vec![1.5],
        albedo: None,
    };

    let result = chain.run(&weather).unwrap();

    // Without decomposition, direct component is zero
    assert!(result.poa_direct[0] < 1.0,
        "poa_direct should be ~0 without decomposition, got {}", result.poa_direct[0]);
}

#[test]
fn test_auto_decomposition_matches_manual_erbs() {
    let loc = tucson();

    // Run with auto-decomposition
    let chain_auto = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0)
        .with_auto_decomposition(true);

    let time = Eastern.with_ymd_and_hms(2020, 6, 15, 12, 0, 0).unwrap();
    let ghi_val = 900.0;

    let weather_auto = batch::WeatherSeries {
        times: vec![time],
        ghi: vec![ghi_val],
        dni: vec![0.0],
        dhi: vec![0.0],
        temp_air: vec![32.0],
        wind_speed: vec![1.5],
        albedo: None,
    };

    let result_auto = chain_auto.run(&weather_auto).unwrap();

    // Manually decompose and run without auto
    let solpos = pvlib::solarposition::get_solarposition(&loc, time).unwrap();
    let doy: i32 = time.format("%j").to_string().parse().unwrap();
    let dni_extra = irradiance::get_extra_radiation(doy);
    let (dni_erbs, dhi_erbs) = irradiance::erbs(ghi_val, solpos.zenith, doy as u32, dni_extra);

    let chain_manual = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0);
    let weather_manual = batch::WeatherSeries {
        times: vec![time],
        ghi: vec![ghi_val],
        dni: vec![dni_erbs],
        dhi: vec![dhi_erbs],
        temp_air: vec![32.0],
        wind_speed: vec![1.5],
        albedo: None,
    };

    let result_manual = chain_manual.run(&weather_manual).unwrap();

    // Results should match
    assert!((result_auto.ac_power[0] - result_manual.ac_power[0]).abs() < 1.0,
        "auto={} manual={}", result_auto.ac_power[0], result_manual.ac_power[0]);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test test_batch_modelchain_auto_decomposition -- --nocapture 2>&1 | head -20`
Expected: FAIL — `with_auto_decomposition` not found

- [ ] **Step 3: Implement auto-decomposition**

In `src/batch.rs`, add the field to `BatchModelChain`:

```rust
pub struct BatchModelChain {
    pub location: crate::location::Location,
    pub surface_tilt: f64,
    pub surface_azimuth: f64,
    pub system_capacity_dc: f64,
    pub gamma_pdc: f64,
    pub inverter_capacity: f64,
    pub inverter_efficiency: f64,
    pub albedo: f64,
    pub transposition_model: irradiance::DiffuseModel,
    pub auto_decomposition: bool,  // NEW
}
```

Update the `pvwatts` constructor to set `auto_decomposition: false`.

Add the builder method:

```rust
/// Builder: enable automatic Erbs GHI→DNI/DHI decomposition.
///
/// When enabled, if both DNI and DHI are zero for a timestep but GHI is positive,
/// the Erbs (1982) model is applied to decompose GHI into DNI and DHI components.
/// This is essential when using weather APIs that only provide GHI.
pub fn with_auto_decomposition(mut self, enabled: bool) -> Self {
    self.auto_decomposition = enabled;
    self
}
```

In the `run()` method, add decomposition logic inside the per-timestep parallel closure, right after computing `solpos` and before computing `aoi_val`. Insert after the airmass calculation (step 2) and before step 3:

```rust
// 2b. Auto-decompose GHI → DNI/DHI if enabled and needed
let (dni_i, dhi_i) = if self.auto_decomposition
    && weather.dni[i].abs() < 1.0
    && weather.dhi[i].abs() < 1.0
    && weather.ghi[i] > 0.0
{
    let doy = weather.times[i].format("%j").to_string().parse::<i32>().unwrap_or(1);
    let dni_extra_val = irradiance::get_extra_radiation(doy);
    irradiance::erbs(weather.ghi[i], solpos.zenith, doy as u32, dni_extra_val)
} else {
    (weather.dni[i], weather.dhi[i])
};
```

Then replace `weather.dni[i]` and `weather.dhi[i]` with `dni_i` and `dhi_i` in the transposition call (step 5).

- [ ] **Step 4: Run all batch tests to verify they pass**

Run: `cargo test test_batch -- --nocapture 2>&1 | tail -5`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/batch.rs tests/test_batch.rs
git commit -m "feat: add auto GHI→DNI/DHI Erbs decomposition to BatchModelChain"
```

---

### Task 5: Add bifacial rear-side gain to `BatchModelChain`

**Files:**
- Modify: `src/batch.rs` (`BatchModelChain` struct + `run` method)
- Modify: `tests/test_batch.rs`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_batch.rs`:

```rust
#[test]
fn test_batch_modelchain_bifacial() {
    let loc = tucson();

    // Without bifacial
    let chain_mono = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0);

    // With bifacial
    let chain_bifi = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0)
        .with_bifacial(0.7, 0.25);

    let weather = batch::WeatherSeries {
        times: vec![summer_noon()],
        ghi: vec![900.0],
        dni: vec![800.0],
        dhi: vec![100.0],
        temp_air: vec![32.0],
        wind_speed: vec![1.5],
        albedo: None,
    };

    let result_mono = chain_mono.run(&weather).unwrap();
    let result_bifi = chain_bifi.run(&weather).unwrap();

    // Bifacial should produce more power than monofacial
    assert!(result_bifi.ac_power[0] > result_mono.ac_power[0],
        "bifacial {} should exceed monofacial {}", result_bifi.ac_power[0], result_mono.ac_power[0]);

    // But not more than 25% more (the cap)
    let ratio = result_bifi.ac_power[0] / result_mono.ac_power[0];
    assert!(ratio < 1.26,
        "bifacial gain ratio {} should be < 1.26", ratio);
}

#[test]
fn test_batch_modelchain_bifacial_zero_means_no_gain() {
    let loc = tucson();

    let chain_zero = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0)
        .with_bifacial(0.0, 0.2);
    let chain_none = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0);

    let weather = batch::WeatherSeries {
        times: vec![summer_noon()],
        ghi: vec![900.0],
        dni: vec![800.0],
        dhi: vec![100.0],
        temp_air: vec![32.0],
        wind_speed: vec![1.5],
        albedo: None,
    };

    let result_zero = chain_zero.run(&weather).unwrap();
    let result_none = chain_none.run(&weather).unwrap();

    assert!((result_zero.ac_power[0] - result_none.ac_power[0]).abs() < 1e-6,
        "zero bifaciality should match non-bifacial: {} vs {}",
        result_zero.ac_power[0], result_none.ac_power[0]);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test test_batch_modelchain_bifacial -- --nocapture 2>&1 | head -20`
Expected: FAIL — `with_bifacial` not found

- [ ] **Step 3: Implement bifacial support**

Add two fields to `BatchModelChain`:

```rust
pub struct BatchModelChain {
    // ... existing fields ...
    pub auto_decomposition: bool,
    pub bifaciality_factor: f64,  // NEW: 0.0 = monofacial, 0.65-0.85 typical for bifacial
    pub bifacial_ground_albedo: f64,  // NEW: albedo specifically for rear-side calculation
}
```

Update `pvwatts` constructor to set both to `0.0` and `0.2` respectively.

Add builder method:

```rust
/// Builder: enable bifacial rear-side gain.
///
/// Models additional energy from rear-side irradiance using a simplified
/// ground-reflected model. The rear-side gain is capped at 25% to prevent
/// unrealistic values.
///
/// # Arguments
/// * `bifaciality_factor` - Ratio of rear-to-front efficiency (0.0-1.0, typical 0.65-0.85)
/// * `ground_albedo` - Ground albedo for rear-side reflection (0.0-1.0, e.g. 0.25 for light gravel)
pub fn with_bifacial(mut self, bifaciality_factor: f64, ground_albedo: f64) -> Self {
    self.bifaciality_factor = bifaciality_factor;
    self.bifacial_ground_albedo = ground_albedo;
    self
}
```

In the `run()` method, after computing `pac` (step 10), apply bifacial gain:

```rust
// 11. Bifacial rear-side gain
let pac = if self.bifaciality_factor > 0.0 && poa.poa_global > 10.0 {
    let rear_gain = (self.bifaciality_factor * self.bifacial_ground_albedo
        * weather.ghi[i] / poa.poa_global).min(0.25);
    pac * (1.0 + rear_gain)
} else {
    pac
};
```

- [ ] **Step 4: Run all batch tests to verify they pass**

Run: `cargo test test_batch -- --nocapture 2>&1 | tail -5`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/batch.rs tests/test_batch.rs
git commit -m "feat: add bifacial rear-side gain to BatchModelChain"
```

---

### Task 6: Add system losses to `BatchModelChain`

**Files:**
- Modify: `src/batch.rs` (`BatchModelChain` struct + `run` method)
- Modify: `tests/test_batch.rs`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_batch.rs`:

```rust
#[test]
fn test_batch_modelchain_system_losses() {
    let loc = tucson();

    let chain_no_loss = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0);
    let chain_14pct = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0)
        .with_system_losses(0.14);

    let weather = batch::WeatherSeries {
        times: vec![summer_noon()],
        ghi: vec![900.0],
        dni: vec![800.0],
        dhi: vec![100.0],
        temp_air: vec![32.0],
        wind_speed: vec![1.5],
        albedo: None,
    };

    let result_no_loss = chain_no_loss.run(&weather).unwrap();
    let result_14pct = chain_14pct.run(&weather).unwrap();

    // 14% losses should reduce power by ~14%
    let ratio = result_14pct.ac_power[0] / result_no_loss.ac_power[0];
    assert!((ratio - 0.86).abs() < 0.01,
        "14% loss ratio should be ~0.86, got {}", ratio);
}

#[test]
fn test_batch_modelchain_zero_losses() {
    let loc = tucson();

    let chain_zero = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0)
        .with_system_losses(0.0);
    let chain_default = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0);

    let weather = batch::WeatherSeries {
        times: vec![summer_noon()],
        ghi: vec![900.0],
        dni: vec![800.0],
        dhi: vec![100.0],
        temp_air: vec![32.0],
        wind_speed: vec![1.5],
        albedo: None,
    };

    let r_zero = chain_zero.run(&weather).unwrap();
    let r_default = chain_default.run(&weather).unwrap();

    assert!((r_zero.ac_power[0] - r_default.ac_power[0]).abs() < 1e-6,
        "zero losses should match default: {} vs {}", r_zero.ac_power[0], r_default.ac_power[0]);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test test_batch_modelchain_system_losses -- --nocapture 2>&1 | head -20`
Expected: FAIL — `with_system_losses` not found

- [ ] **Step 3: Implement system losses**

Add field to `BatchModelChain`:

```rust
pub system_losses: f64,  // NEW: 0.0 = no losses, 0.14 = 14% losses
```

Set to `0.0` in `pvwatts` constructor.

Add builder:

```rust
/// Builder: set system losses factor.
///
/// Applied as a flat derating to DC power before AC conversion.
/// Accounts for wiring, connections, mismatch, soiling, etc.
///
/// # Arguments
/// * `losses` - Loss fraction (0.0-1.0, e.g. 0.14 for 14% total system losses)
pub fn with_system_losses(mut self, losses: f64) -> Self {
    self.system_losses = losses;
    self
}
```

In the `run()` method, apply losses to DC power after computing `pdc` (step 9):

```rust
// 9. DC power
let pdc = self.system_capacity_dc * (eff_irrad / 1000.0)
    * (1.0 + self.gamma_pdc * (t_cell - 25.0));
let pdc = pdc.max(0.0);

// 9b. System losses
let pdc = pdc * (1.0 - self.system_losses);
```

- [ ] **Step 4: Run all batch tests to verify they pass**

Run: `cargo test test_batch -- --nocapture 2>&1 | tail -5`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/batch.rs tests/test_batch.rs
git commit -m "feat: add system losses to BatchModelChain"
```

---

### Task 7: Integration test — full pipeline with all new features

**Files:**
- Modify: `tests/test_batch.rs`

- [ ] **Step 1: Write an integration test that exercises all new features together**

Add to `tests/test_batch.rs`:

```rust
#[test]
fn test_all_new_features_combined() {
    use chrono::NaiveDate;

    let loc = tucson();

    // Use all new builder options together
    let chain = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0)
        .with_gamma_pdc(-0.004)
        .with_inverter(5000.0, 0.96)
        .with_albedo(0.2)
        .with_auto_decomposition(true)
        .with_bifacial(0.7, 0.25)
        .with_system_losses(0.14);

    // Use from_utc constructor (GHI-only, DNI/DHI zero to trigger auto-decomposition)
    let timestamps = vec![
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(14, 0, 0).unwrap(),
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(16, 0, 0).unwrap(),
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(18, 0, 0).unwrap(),
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(20, 0, 0).unwrap(),
    ];

    let weather = batch::WeatherSeries::from_utc(
        &timestamps,
        "US/Eastern",
        vec![800.0, 600.0, 200.0, 0.0],   // GHI
        vec![0.0, 0.0, 0.0, 0.0],         // DNI (zeros — auto-decompose)
        vec![0.0, 0.0, 0.0, 0.0],         // DHI (zeros — auto-decompose)
        vec![32.0, 30.0, 26.0, 22.0],     // temp_air
        vec![2.0, 2.5, 3.0, 2.0],         // wind_speed
    ).unwrap();

    let result = chain.run(&weather).unwrap();

    // Check solar_elevation is populated
    assert_eq!(result.solar_elevation.len(), 4);
    for i in 0..4 {
        assert!((result.solar_elevation[i] + result.solar_zenith[i] - 90.0).abs() < 1e-10);
    }

    // Daytime timesteps should produce positive power
    assert!(result.ac_power[0] > 0.0, "midday should produce power");
    assert!(result.ac_power[1] > 0.0, "afternoon should produce power");

    // Night timestep should produce zero power
    assert!(result.ac_power[3] < 1.0, "night should produce ~0 power");

    // Power should be less than capacity (losses applied)
    for p in &result.ac_power {
        assert!(*p <= 5000.0 * 1.1, "power {} exceeds capacity", p);
    }

    // Verify solar_position_batch_utc works independently
    let (zen, _az, elev) = batch::solar_position_batch_utc(
        32.2, -110.9, 700.0, &timestamps,
    ).unwrap();
    assert_eq!(zen.len(), 4);
    assert_eq!(elev.len(), 4);
}
```

- [ ] **Step 2: Run integration test**

Run: `cargo test test_all_new_features_combined -- --nocapture 2>&1 | tail -10`
Expected: PASS

- [ ] **Step 3: Run the full test suite to verify no regressions**

Run: `cargo test 2>&1 | tail -10`
Expected: All tests PASS

- [ ] **Step 4: Commit**

```bash
git add tests/test_batch.rs
git commit -m "test: add integration test for all batch API improvements"
```
