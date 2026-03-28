use chrono::{NaiveDate, TimeZone};
use chrono_tz::US::Eastern;
use pvlib::location::Location;
use pvlib::batch;
use pvlib::{atmosphere, irradiance, iam, inverter};

/// Helper: create a test location (Tucson, AZ).
fn tucson() -> Location {
    Location::new(32.2, -110.9, Eastern, 700.0, "Tucson")
}

/// Helper: noon on a summer day.
fn summer_noon() -> chrono::DateTime<chrono_tz::Tz> {
    Eastern.with_ymd_and_hms(2020, 6, 15, 12, 0, 0).unwrap()
}

// ---------------------------------------------------------------------------
// Solar Position Batch
// ---------------------------------------------------------------------------

#[test]
fn test_batch_solar_position() {
    let loc = tucson();
    let times = vec![
        Eastern.with_ymd_and_hms(2020, 6, 15, 10, 0, 0).unwrap(),
        summer_noon(),
        Eastern.with_ymd_and_hms(2020, 6, 15, 14, 0, 0).unwrap(),
    ];

    let (zenith, azimuth, elevation) = batch::solar_position_batch(&loc, &times).unwrap();

    assert_eq!(zenith.len(), 3);
    assert_eq!(azimuth.len(), 3);
    assert_eq!(elevation.len(), 3);

    // Morning should have higher zenith than noon
    assert!(zenith[1] < zenith[0], "noon zenith {} should be less than morning {}", zenith[1], zenith[0]);

    // Zenith should be reasonable (0-90 for daytime)
    for z in &zenith {
        assert!(*z > 0.0 && *z < 90.0, "zenith {} out of range", z);
    }

    // Elevation + zenith = 90
    for (z, e) in zenith.iter().zip(elevation.iter()) {
        assert!((z + e - 90.0).abs() < 1e-10);
    }
}

// ---------------------------------------------------------------------------
// Atmosphere Batch
// ---------------------------------------------------------------------------

#[test]
fn test_batch_airmass_relative() {
    let zeniths = vec![0.0, 30.0, 60.0, 80.0];
    let am = batch::airmass_relative_batch(&zeniths);

    assert_eq!(am.len(), 4);
    // At zenith=0 airmass should be ~1.0
    assert!((am[0] - 1.0).abs() < 0.01);
    // Airmass increases with zenith
    assert!(am[1] > am[0]);
    assert!(am[2] > am[1]);
    assert!(am[3] > am[2]);
}

#[test]
fn test_batch_airmass_absolute() {
    let am_rel = vec![1.0, 2.0, 5.0];
    let pressure = 101325.0;
    let am_abs = batch::airmass_absolute_batch(&am_rel, pressure);

    // At standard pressure, absolute == relative
    for (r, a) in am_rel.iter().zip(am_abs.iter()) {
        assert!((r - a).abs() < 1e-10);
    }
}

// ---------------------------------------------------------------------------
// Clear Sky Batch
// ---------------------------------------------------------------------------

#[test]
fn test_batch_ineichen() {
    let zenith = vec![20.0, 40.0, 60.0];
    let am_abs = vec![1.1, 1.5, 2.0];
    let linke_turbidity = 3.0;
    let altitude = 700.0;

    let (ghi, dni, dhi) = batch::ineichen_batch(&zenith, &am_abs, linke_turbidity, altitude);

    assert_eq!(ghi.len(), 3);
    assert_eq!(dni.len(), 3);
    assert_eq!(dhi.len(), 3);

    // All should be positive for daytime zeniths
    for i in 0..3 {
        assert!(ghi[i] > 0.0, "ghi[{}] = {}", i, ghi[i]);
        assert!(dni[i] > 0.0, "dni[{}] = {}", i, dni[i]);
        assert!(dhi[i] >= 0.0, "dhi[{}] = {}", i, dhi[i]);
    }

    // GHI should decrease with zenith
    assert!(ghi[0] > ghi[1]);
    assert!(ghi[1] > ghi[2]);
}

#[test]
fn test_batch_bird() {
    let zenith = vec![20.0, 40.0, 60.0];
    let am_rel = vec![1.1, 1.5, 2.0];
    let aod380 = 0.15;
    let aod500 = 0.1;
    let pw = 1.0;

    let (ghi, dni, dhi) = batch::bird_batch(&zenith, &am_rel, aod380, aod500, pw);

    assert_eq!(ghi.len(), 3);
    for i in 0..3 {
        assert!(ghi[i] > 0.0);
        assert!(dni[i] > 0.0);
        assert!(dhi[i] >= 0.0);
    }
}

// ---------------------------------------------------------------------------
// Irradiance Batch
// ---------------------------------------------------------------------------

#[test]
fn test_batch_aoi() {
    let surface_tilt = 30.0;
    let surface_azimuth = 180.0;
    let solar_zenith = vec![20.0, 40.0, 60.0];
    let solar_azimuth = vec![180.0, 180.0, 180.0];

    let aoi = batch::aoi_batch(surface_tilt, surface_azimuth, &solar_zenith, &solar_azimuth);

    assert_eq!(aoi.len(), 3);
    for a in &aoi {
        assert!(*a >= 0.0 && *a <= 180.0);
    }
}

#[test]
fn test_batch_extra_radiation() {
    let doy = vec![1, 91, 182, 274];
    let extra = batch::extra_radiation_batch(&doy);

    assert_eq!(extra.len(), 4);
    for e in &extra {
        // Should be around 1366 W/m2 +/- ~3%
        assert!(*e > 1320.0 && *e < 1420.0, "extra radiation {} out of range", e);
    }
}

#[test]
fn test_batch_erbs() {
    let ghi = vec![500.0, 800.0, 300.0];
    let zenith = vec![30.0, 20.0, 50.0];
    let doy = vec![172, 172, 172];
    let dni_extra = vec![1366.0, 1366.0, 1366.0];

    let (dni, dhi) = batch::erbs_batch(&ghi, &zenith, &doy, &dni_extra);

    assert_eq!(dni.len(), 3);
    assert_eq!(dhi.len(), 3);

    for i in 0..3 {
        assert!(dni[i] >= 0.0);
        assert!(dhi[i] >= 0.0);
        // DHI should be less than GHI
        assert!(dhi[i] <= ghi[i] + 1.0);
    }
}

#[test]
fn test_batch_disc() {
    let ghi = vec![500.0, 800.0];
    let solar_zenith = vec![30.0, 20.0];
    let doy = vec![172, 172];

    let (dni, kt, am) = batch::disc_batch(&ghi, &solar_zenith, &doy, Some(101325.0));

    assert_eq!(dni.len(), 2);
    assert_eq!(kt.len(), 2);
    assert_eq!(am.len(), 2);

    for i in 0..2 {
        assert!(dni[i] >= 0.0);
        assert!(kt[i] >= 0.0 && kt[i] <= 1.0);
    }
}

// ---------------------------------------------------------------------------
// Batch results match scalar results
// ---------------------------------------------------------------------------

#[test]
fn test_batch_matches_scalar_airmass() {
    let zeniths = vec![10.0, 30.0, 50.0, 70.0, 85.0];
    let batch_result = batch::airmass_relative_batch(&zeniths);

    for (z, batch_am) in zeniths.iter().zip(batch_result.iter()) {
        let scalar_am = atmosphere::get_relative_airmass(*z);
        assert!(
            (batch_am - scalar_am).abs() < 1e-10,
            "Mismatch at zenith {}: batch={} scalar={}",
            z, batch_am, scalar_am
        );
    }
}

#[test]
fn test_batch_matches_scalar_iam() {
    let aoi_vals = vec![0.0, 15.0, 30.0, 45.0, 60.0, 75.0];
    let n = 1.526;
    let k = 4.0;
    let l = 0.002;

    let batch_result = batch::iam_physical_batch(&aoi_vals, n, k, l);

    for (a, batch_iam) in aoi_vals.iter().zip(batch_result.iter()) {
        let scalar_iam = iam::physical(*a, n, k, l);
        assert!(
            (batch_iam - scalar_iam).abs() < 1e-10,
            "Mismatch at aoi {}: batch={} scalar={}",
            a, batch_iam, scalar_iam
        );
    }
}

#[test]
fn test_batch_matches_scalar_inverter() {
    let pdc_vals = vec![0.0, 1000.0, 3000.0, 5000.0, 6000.0];
    let pdc0 = 5000.0;
    let eta_nom = 0.96;
    let eta_ref = 0.9637;

    let batch_result = batch::pvwatts_ac_batch(&pdc_vals, pdc0, eta_nom, eta_ref);

    for (p, batch_ac) in pdc_vals.iter().zip(batch_result.iter()) {
        let scalar_ac = inverter::pvwatts_ac(*p, pdc0, eta_nom, eta_ref);
        assert!(
            (batch_ac - scalar_ac).abs() < 1e-10,
            "Mismatch at pdc {}: batch={} scalar={}",
            p, batch_ac, scalar_ac
        );
    }
}

// ---------------------------------------------------------------------------
// BatchModelChain builder pattern
// ---------------------------------------------------------------------------

#[test]
fn test_batch_modelchain_builder() {
    let loc = tucson();
    let chain = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0)
        .with_gamma_pdc(-0.003)
        .with_inverter(5500.0, 0.97)
        .with_albedo(0.25)
        .with_transposition(irradiance::DiffuseModel::HayDavies);

    assert!((chain.gamma_pdc - (-0.003)).abs() < 1e-10);
    assert!((chain.inverter_capacity - 5500.0).abs() < 1e-10);
    assert!((chain.inverter_efficiency - 0.97).abs() < 1e-10);
    assert!((chain.albedo - 0.25).abs() < 1e-10);
    assert_eq!(chain.transposition_model, irradiance::DiffuseModel::HayDavies);
}

// ---------------------------------------------------------------------------
// SimulationSeries helper methods
// ---------------------------------------------------------------------------

#[test]
fn test_simulation_series_helpers() {
    let series = batch::SimulationSeries {
        solar_zenith: vec![30.0, 40.0, 50.0],
        solar_elevation: vec![60.0, 50.0, 40.0],
        solar_azimuth: vec![180.0, 190.0, 200.0],
        airmass: vec![1.1, 1.3, 1.5],
        aoi: vec![10.0, 20.0, 30.0],
        poa_global: vec![800.0, 700.0, 500.0],
        poa_direct: vec![600.0, 500.0, 300.0],
        poa_diffuse: vec![200.0, 200.0, 200.0],
        cell_temperature: vec![35.0, 33.0, 30.0],
        effective_irradiance: vec![780.0, 680.0, 480.0],
        dc_power: vec![3800.0, 3300.0, 2300.0],
        ac_power: vec![3600.0, 3100.0, 2100.0],
    };

    // Total energy = sum of positive ac_power
    let total = series.total_energy_wh();
    assert!((total - 8800.0).abs() < 1e-6);

    // Peak power = max ac_power
    let peak = series.peak_power();
    assert!((peak - 3600.0).abs() < 1e-6);

    // Capacity factor
    let cf = series.capacity_factor(5000.0);
    let expected_cf = 8800.0 / (5000.0 * 3.0);
    assert!((cf - expected_cf).abs() < 1e-6);
}

#[test]
fn test_simulation_series_edge_cases() {
    let empty = batch::SimulationSeries {
        solar_zenith: vec![],
        solar_elevation: vec![],
        solar_azimuth: vec![],
        airmass: vec![],
        aoi: vec![],
        poa_global: vec![],
        poa_direct: vec![],
        poa_diffuse: vec![],
        cell_temperature: vec![],
        effective_irradiance: vec![],
        dc_power: vec![],
        ac_power: vec![],
    };

    assert_eq!(empty.total_energy_wh(), 0.0);
    assert_eq!(empty.peak_power(), 0.0);
    assert_eq!(empty.capacity_factor(5000.0), 0.0);
    assert_eq!(empty.capacity_factor(0.0), 0.0);
}

// ---------------------------------------------------------------------------
// Full Pipeline (BatchModelChain)
// ---------------------------------------------------------------------------

#[test]
fn test_batch_modelchain_basic() {
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

    assert_eq!(result.ac_power.len(), 3);
    assert_eq!(result.dc_power.len(), 3);
    assert_eq!(result.solar_zenith.len(), 3);

    // All outputs should be positive for daytime inputs
    for i in 0..3 {
        assert!(result.ac_power[i] > 0.0, "ac_power[{}] = {}", i, result.ac_power[i]);
        assert!(result.dc_power[i] > 0.0, "dc_power[{}] = {}", i, result.dc_power[i]);
        assert!(result.poa_global[i] > 0.0);
        assert!(result.effective_irradiance[i] > 0.0);
    }

    // AC power should not exceed system capacity
    for p in &result.ac_power {
        assert!(*p <= 5000.0 * 1.1, "AC power {} exceeds capacity", p);
    }
}

// ---------------------------------------------------------------------------
// BatchModelChain vs scalar pipeline equivalence
// ---------------------------------------------------------------------------

/// Assert two f64 values are bitwise identical, treating NaN == NaN as true.
fn assert_f64_identical(left: f64, right: f64, label: &str, index: usize) {
    if left.is_nan() && right.is_nan() {
        return; // both NaN is considered equal here
    }
    assert_eq!(left, right, "{} mismatch at timestep {}", label, index);
}

#[test]
fn test_batch_modelchain_matches_scalar_pipeline() {
    let loc = tucson();
    let surface_tilt = 30.0;
    let surface_azimuth = 180.0;
    let system_capacity_dc = 5000.0;
    let gamma_pdc = -0.004;
    let inverter_efficiency = 0.96;
    let albedo = 0.2;
    let transposition_model = irradiance::DiffuseModel::Perez;

    let chain = batch::BatchModelChain::pvwatts(loc.clone(), surface_tilt, surface_azimuth, system_capacity_dc);

    // 5 distinct daytime timestamps with realistic weather.
    // All times are chosen so the sun is well above the horizon for
    // Tucson (lat 32.2) to avoid NaN from pre-dawn / post-sunset conditions.
    let times = vec![
        Eastern.with_ymd_and_hms(2020, 3, 21, 12, 0, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 6, 15, 10, 0, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 6, 15, 12, 0, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 9, 22, 13, 0, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 12, 21, 12, 0, 0).unwrap(),
    ];
    let ghi_vals   = vec![650.0, 700.0, 950.0, 700.0, 500.0];
    let dni_vals   = vec![550.0, 600.0, 850.0, 600.0, 400.0];
    let dhi_vals   = vec![100.0, 100.0, 100.0, 100.0,  100.0];
    let temp_vals  = vec![ 18.0,  28.0,  32.0,  26.0,  12.0];
    let wind_vals  = vec![  3.0,   2.0,   1.5,   2.5,   4.0];

    let weather = batch::WeatherSeries {
        times: times.clone(),
        ghi: ghi_vals.clone(),
        dni: dni_vals.clone(),
        dhi: dhi_vals.clone(),
        temp_air: temp_vals.clone(),
        wind_speed: wind_vals.clone(),
        albedo: None,
    };

    let batch_result = chain.run(&weather).unwrap();

    // Reproduce the scalar pipeline for each timestep identically to
    // BatchModelChain::run() (see src/batch.rs lines 423-479).
    let pressure = atmosphere::alt2pres(loc.altitude);
    for i in 0..times.len() {
        let solpos = pvlib::solarposition::get_solarposition(&loc, times[i]).unwrap();

        let am_rel = atmosphere::get_relative_airmass(solpos.zenith);
        let am_abs = if am_rel.is_nan() || am_rel <= 0.0 {
            0.0
        } else {
            atmosphere::get_absolute_airmass(am_rel, pressure)
        };

        let aoi_val = irradiance::aoi(surface_tilt, surface_azimuth, solpos.zenith, solpos.azimuth);

        let doy: i32 = times[i].format("%j").to_string().parse().unwrap_or(1);
        let dni_extra = irradiance::get_extra_radiation(doy);

        let poa = irradiance::get_total_irradiance(
            surface_tilt, surface_azimuth,
            solpos.zenith, solpos.azimuth,
            dni_vals[i], ghi_vals[i], dhi_vals[i],
            albedo,
            transposition_model,
            Some(dni_extra),
            if am_rel.is_nan() { None } else { Some(am_rel) },
        );

        let iam_val = iam::physical(aoi_val, 1.526, 4.0, 0.002);
        let eff_irrad = (poa.poa_direct * iam_val + poa.poa_diffuse).max(0.0);
        let t_cell = temp_vals[i] + poa.poa_global * (45.0 - 20.0) / 800.0;
        let pdc = (system_capacity_dc * (eff_irrad / 1000.0)
            * (1.0 + gamma_pdc * (t_cell - 25.0))).max(0.0);
        let pac = pvlib::inverter::pvwatts_ac(pdc, system_capacity_dc, inverter_efficiency, 0.9637);

        assert_f64_identical(batch_result.solar_zenith[i], solpos.zenith, "zenith", i);
        assert_f64_identical(batch_result.solar_azimuth[i], solpos.azimuth, "azimuth", i);
        assert_f64_identical(batch_result.airmass[i], am_abs, "airmass", i);
        assert_f64_identical(batch_result.aoi[i], aoi_val, "aoi", i);
        assert_f64_identical(batch_result.poa_global[i], poa.poa_global, "poa_global", i);
        assert_f64_identical(batch_result.poa_direct[i], poa.poa_direct, "poa_direct", i);
        assert_f64_identical(batch_result.poa_diffuse[i], poa.poa_diffuse, "poa_diffuse", i);
        assert_f64_identical(batch_result.cell_temperature[i], t_cell, "cell_temperature", i);
        assert_f64_identical(batch_result.effective_irradiance[i], eff_irrad, "effective_irradiance", i);
        assert_f64_identical(batch_result.dc_power[i], pdc, "dc_power", i);
        assert_f64_identical(batch_result.ac_power[i], pac, "ac_power", i);
    }
}

// ---------------------------------------------------------------------------
// TMY Year Performance Benchmark
// ---------------------------------------------------------------------------

#[test]
fn test_batch_modelchain_tmy_year_performance() {
    use std::time::Instant;
    use chrono::Duration;

    let loc = tucson();
    let chain = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0);

    // Generate 8760 hours of synthetic weather data for a full TMY year
    let start_time = Eastern.with_ymd_and_hms(2020, 1, 1, 0, 0, 0).unwrap();
    let n = 8760;

    let times: Vec<_> = (0..n)
        .map(|h| start_time + Duration::hours(h as i64))
        .collect();

    // Synthetic sinusoidal weather data
    let ghi: Vec<f64> = (0..n).map(|h| {
        let hour_of_day = (h % 24) as f64;
        if hour_of_day >= 6.0 && hour_of_day <= 18.0 {
            let solar_fraction = ((hour_of_day - 6.0) / 6.0 * std::f64::consts::PI).sin();
            800.0 * solar_fraction
        } else {
            0.0
        }
    }).collect();

    let dni: Vec<f64> = ghi.iter().map(|g| g * 0.7).collect();
    let dhi: Vec<f64> = ghi.iter().map(|g| g * 0.3).collect();
    let temp_air: Vec<f64> = (0..n).map(|h| {
        let hour_of_day = (h % 24) as f64;
        20.0 + 10.0 * ((hour_of_day - 14.0) / 12.0 * std::f64::consts::PI).sin()
    }).collect();
    let wind_speed: Vec<f64> = vec![2.0; n];

    let weather = batch::WeatherSeries {
        times,
        ghi,
        dni,
        dhi,
        temp_air,
        wind_speed,
        albedo: None,
    };

    let start = Instant::now();
    let result = chain.run(&weather).unwrap();
    let elapsed = start.elapsed();

    println!("BatchModelChain TMY year (8760 hours): {:?}", elapsed);

    // Verify results are physically reasonable
    assert_eq!(result.ac_power.len(), n);

    let total_energy = result.total_energy_wh();
    let peak = result.peak_power();
    let cf = result.capacity_factor(5000.0);

    println!("Total energy: {:.0} Wh ({:.1} kWh)", total_energy, total_energy / 1000.0);
    println!("Peak AC power: {:.0} W", peak);
    println!("Capacity factor: {:.3}", cf);

    // Sanity checks
    assert!(total_energy > 0.0, "Total energy should be positive");
    assert!(peak > 0.0, "Peak power should be positive");
    assert!(peak <= 5000.0 * 1.1, "Peak power should not greatly exceed capacity");
    assert!(cf > 0.0 && cf < 1.0, "Capacity factor {} out of range", cf);

    // Typical capacity factor for Tucson ~20-25%, with synthetic data maybe different
    // but should be in a reasonable range
    assert!(cf > 0.01, "Capacity factor too low: {}", cf);
    assert!(cf < 0.5, "Capacity factor too high: {}", cf);
}

// ---------------------------------------------------------------------------
// Solar Position Batch UTC
// ---------------------------------------------------------------------------

#[test]
fn test_solar_position_batch_utc() {
    // Summer daytime timestamps in UTC for Tucson (lat 32.2, lon -110.9)
    let times = vec![
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(18, 0, 0).unwrap(), // ~11 AM local
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(19, 0, 0).unwrap(), // ~noon local
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(20, 0, 0).unwrap(), // ~1 PM local
    ];

    let (zenith, azimuth, elevation) = batch::solar_position_batch_utc(32.2, -110.9, 700.0, &times).unwrap();

    assert_eq!(zenith.len(), 3);
    assert_eq!(azimuth.len(), 3);
    assert_eq!(elevation.len(), 3);

    // Daytime: zenith should be < 90
    for z in &zenith {
        assert!(*z > 0.0 && *z < 90.0, "zenith {} out of daytime range", z);
    }

    // Elevation + zenith = 90
    for (z, e) in zenith.iter().zip(elevation.iter()) {
        assert!((z + e - 90.0).abs() < 1e-10, "zenith {} + elevation {} != 90", z, e);
    }
}

#[test]
fn test_solar_position_batch_utc_matches_tz_version() {
    // Create the same moment in time via both paths and verify identical results.
    let ndt = NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(19, 0, 0).unwrap();

    // UTC convenience function
    let (zen_utc, az_utc, el_utc) = batch::solar_position_batch_utc(32.2, -110.9, 700.0, &[ndt]).unwrap();

    // Existing tz-aware path: construct the same instant manually
    let loc = Location::new(32.2, -110.9, chrono_tz::UTC, 700.0, "test");
    let dt = chrono::Utc.from_utc_datetime(&ndt).with_timezone(&chrono_tz::UTC);
    let (zen_tz, az_tz, el_tz) = batch::solar_position_batch(&loc, &[dt]).unwrap();

    assert!((zen_utc[0] - zen_tz[0]).abs() < 1e-10, "zenith mismatch: {} vs {}", zen_utc[0], zen_tz[0]);
    assert!((az_utc[0] - az_tz[0]).abs() < 1e-10, "azimuth mismatch: {} vs {}", az_utc[0], az_tz[0]);
    assert!((el_utc[0] - el_tz[0]).abs() < 1e-10, "elevation mismatch: {} vs {}", el_utc[0], el_tz[0]);
}

// ---------------------------------------------------------------------------
// WeatherSeries::from_utc
// ---------------------------------------------------------------------------

#[test]
fn test_weather_series_from_utc() {
    let timestamps = vec![
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(18, 0, 0).unwrap(),
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(19, 0, 0).unwrap(),
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(20, 0, 0).unwrap(),
    ];
    let ghi = vec![600.0, 900.0, 700.0];
    let dni = vec![500.0, 800.0, 600.0];
    let dhi = vec![100.0, 100.0, 100.0];
    let temp_air = vec![28.0, 32.0, 30.0];
    let wind_speed = vec![2.0, 1.5, 2.5];

    let ws = batch::WeatherSeries::from_utc(
        &timestamps, "US/Eastern",
        ghi.clone(), dni.clone(), dhi.clone(),
        temp_air.clone(), wind_speed.clone(),
    ).unwrap();

    assert_eq!(ws.times.len(), 3);
    assert_eq!(ws.ghi.len(), 3);
    assert_eq!(ws.dni.len(), 3);
    assert_eq!(ws.dhi.len(), 3);
    assert_eq!(ws.temp_air.len(), 3);
    assert_eq!(ws.wind_speed.len(), 3);
    assert!(ws.albedo.is_none());
}

#[test]
fn test_weather_series_from_utc_invalid_tz() {
    let timestamps = vec![
        NaiveDate::from_ymd_opt(2020, 6, 15).unwrap().and_hms_opt(18, 0, 0).unwrap(),
    ];
    let result = batch::WeatherSeries::from_utc(
        &timestamps, "Invalid/Timezone",
        vec![600.0], vec![500.0], vec![100.0],
        vec![28.0], vec![2.0],
    );
    assert!(result.is_err());
}

// ---------------------------------------------------------------------------
// Solar Elevation
// ---------------------------------------------------------------------------

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

    // elevation + zenith = 90 for every timestep
    for i in 0..3 {
        let sum = result.solar_elevation[i] + result.solar_zenith[i];
        assert!(
            (sum - 90.0).abs() < 1e-10,
            "elevation[{}] ({}) + zenith[{}] ({}) = {}, expected 90.0",
            i, result.solar_elevation[i], i, result.solar_zenith[i], sum
        );
    }
}

// ---------------------------------------------------------------------------
// Auto GHI→DNI/DHI Decomposition
// ---------------------------------------------------------------------------

/// Helper: build a daytime weather series with GHI-only (DNI/DHI = 0).
fn ghi_only_weather() -> batch::WeatherSeries {
    let times = vec![
        Eastern.with_ymd_and_hms(2020, 6, 15, 10, 0, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 6, 15, 12, 0, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 6, 15, 14, 0, 0).unwrap(),
    ];
    batch::WeatherSeries {
        times,
        ghi: vec![600.0, 900.0, 700.0],
        dni: vec![0.0, 0.0, 0.0],
        dhi: vec![0.0, 0.0, 0.0],
        temp_air: vec![28.0, 32.0, 30.0],
        wind_speed: vec![2.0, 1.5, 2.5],
        albedo: None,
    }
}

#[test]
fn test_batch_modelchain_auto_decomposition() {
    let loc = tucson();
    let chain = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0)
        .with_auto_decomposition(true);

    assert!(chain.auto_decomposition);

    let weather = ghi_only_weather();
    let result = chain.run(&weather).unwrap();

    assert_eq!(result.ac_power.len(), 3);
    // With auto decomposition enabled, GHI-only data should produce positive power
    for i in 0..3 {
        assert!(
            result.ac_power[i] > 0.0,
            "ac_power[{}] = {} should be positive with auto decomposition",
            i, result.ac_power[i]
        );
        assert!(
            result.poa_global[i] > 0.0,
            "poa_global[{}] = {} should be positive with auto decomposition",
            i, result.poa_global[i]
        );
    }
}

#[test]
fn test_batch_modelchain_no_decomposition_zero_power() {
    let loc = tucson();
    let chain = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0);

    assert!(!chain.auto_decomposition);

    let weather = ghi_only_weather();
    let result = chain.run(&weather).unwrap();

    // Without auto decomposition, zero DNI/DHI should produce near-zero poa_direct
    for i in 0..3 {
        assert!(
            result.poa_direct[i].abs() < 1.0,
            "poa_direct[{}] = {} should be near-zero without decomposition",
            i, result.poa_direct[i]
        );
    }
}

#[test]
fn test_auto_decomposition_matches_manual_erbs() {
    let loc = tucson();

    // Step 1: Run with auto decomposition on GHI-only data
    let chain_auto = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0)
        .with_auto_decomposition(true);
    let weather_ghi_only = ghi_only_weather();
    let result_auto = chain_auto.run(&weather_ghi_only).unwrap();

    // Step 2: Manually decompose GHI via erbs, then run without auto decomposition
    let times = &weather_ghi_only.times;
    let ghi = &weather_ghi_only.ghi;
    let mut manual_dni = Vec::new();
    let mut manual_dhi = Vec::new();

    for i in 0..times.len() {
        let solpos = pvlib::solarposition::get_solarposition(&loc, times[i]).unwrap();
        let doy: i32 = times[i].format("%j").to_string().parse().unwrap_or(1);
        let dni_extra = irradiance::get_extra_radiation(doy);
        let (dni, dhi) = irradiance::erbs(ghi[i], solpos.zenith, doy as u32, dni_extra);
        manual_dni.push(dni);
        manual_dhi.push(dhi);
    }

    let weather_manual = batch::WeatherSeries {
        times: times.clone(),
        ghi: ghi.clone(),
        dni: manual_dni,
        dhi: manual_dhi,
        temp_air: weather_ghi_only.temp_air.clone(),
        wind_speed: weather_ghi_only.wind_speed.clone(),
        albedo: None,
    };

    let chain_manual = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0);
    let result_manual = chain_manual.run(&weather_manual).unwrap();

    // Step 3: Compare results -- should be identical
    for i in 0..times.len() {
        assert!(
            (result_auto.poa_global[i] - result_manual.poa_global[i]).abs() < 1e-6,
            "poa_global mismatch at {}: auto={} manual={}",
            i, result_auto.poa_global[i], result_manual.poa_global[i]
        );
        assert!(
            (result_auto.ac_power[i] - result_manual.ac_power[i]).abs() < 1e-6,
            "ac_power mismatch at {}: auto={} manual={}",
            i, result_auto.ac_power[i], result_manual.ac_power[i]
        );
        assert!(
            (result_auto.dc_power[i] - result_manual.dc_power[i]).abs() < 1e-6,
            "dc_power mismatch at {}: auto={} manual={}",
            i, result_auto.dc_power[i], result_manual.dc_power[i]
        );
    }
}

// ---------------------------------------------------------------------------
// System Losses
// ---------------------------------------------------------------------------

#[test]
fn test_batch_modelchain_system_losses() {
    let loc = tucson();
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

    // No losses baseline
    let chain_no_loss = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0);
    let result_no_loss = chain_no_loss.run(&weather).unwrap();

    // 14% losses
    let chain_14 = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0)
        .with_system_losses(0.14);
    let result_14 = chain_14.run(&weather).unwrap();

    // DC power with 14% losses should be ~86% of no-loss DC power
    for i in 0..3 {
        let ratio = result_14.dc_power[i] / result_no_loss.dc_power[i];
        assert!(
            (ratio - 0.86).abs() < 0.001,
            "dc_power ratio at [{}] is {}, expected ~0.86",
            i, ratio
        );
    }

    // AC power with losses should be less than without
    for i in 0..3 {
        assert!(
            result_14.ac_power[i] < result_no_loss.ac_power[i],
            "ac_power[{}] with losses ({}) should be less than without ({})",
            i, result_14.ac_power[i], result_no_loss.ac_power[i]
        );
    }
}

#[test]
fn test_batch_modelchain_zero_losses() {
    let loc = tucson();
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

    // Default (system_losses = 0.0)
    let chain_default = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0);
    let result_default = chain_default.run(&weather).unwrap();

    // Explicitly set losses to 0.0
    let chain_zero = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0)
        .with_system_losses(0.0);
    let result_zero = chain_zero.run(&weather).unwrap();

    // Should be exactly identical
    for i in 0..3 {
        assert_eq!(
            result_default.ac_power[i], result_zero.ac_power[i],
            "ac_power[{}] mismatch: default={} zero_losses={}",
            i, result_default.ac_power[i], result_zero.ac_power[i]
        );
        assert_eq!(
            result_default.dc_power[i], result_zero.dc_power[i],
            "dc_power[{}] mismatch: default={} zero_losses={}",
            i, result_default.dc_power[i], result_zero.dc_power[i]
        );
    }
}

// ---------------------------------------------------------------------------
// Bifacial Rear-Side Gain
// ---------------------------------------------------------------------------

#[test]
fn test_batch_modelchain_bifacial() {
    let loc = tucson();
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

    // Monofacial baseline
    let mono = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0);
    let result_mono = mono.run(&weather).unwrap();

    // Bifacial with typical parameters
    let bifi = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0)
        .with_bifacial(0.75, 0.25);
    let result_bifi = bifi.run(&weather).unwrap();

    // Bifacial should produce more power than monofacial
    for i in 0..3 {
        assert!(
            result_bifi.ac_power[i] > result_mono.ac_power[i],
            "bifacial ac_power[{}] ({}) should exceed monofacial ({})",
            i, result_bifi.ac_power[i], result_mono.ac_power[i]
        );
        // Ratio should be reasonable: less than 1.30 (25% DC cap, but inverter
        // nonlinear efficiency can slightly amplify the ratio at AC level)
        let ratio = result_bifi.ac_power[i] / result_mono.ac_power[i];
        assert!(
            ratio < 1.30,
            "bifacial/monofacial ratio at [{}] is {}, expected < 1.30",
            i, ratio
        );
    }

    // Total energy should also be higher
    assert!(result_bifi.total_energy_wh() > result_mono.total_energy_wh());
}

#[test]
fn test_batch_modelchain_bifacial_zero_means_no_gain() {
    let loc = tucson();
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

    // Default (bifaciality_factor = 0.0)
    let chain_default = batch::BatchModelChain::pvwatts(loc.clone(), 30.0, 180.0, 5000.0);
    let result_default = chain_default.run(&weather).unwrap();

    // Explicitly set bifaciality_factor = 0.0 via builder
    let chain_zero = batch::BatchModelChain::pvwatts(loc, 30.0, 180.0, 5000.0)
        .with_bifacial(0.0, 0.3);
    let result_zero = chain_zero.run(&weather).unwrap();

    // Should be exactly identical
    for i in 0..3 {
        assert_eq!(
            result_default.ac_power[i], result_zero.ac_power[i],
            "ac_power[{}] mismatch: default={} zero_bifacial={}",
            i, result_default.ac_power[i], result_zero.ac_power[i]
        );
    }
}

// ---------------------------------------------------------------------------
// Integration: All New Features Combined
// ---------------------------------------------------------------------------

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
