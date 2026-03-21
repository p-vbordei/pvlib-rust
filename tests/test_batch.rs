use chrono::TimeZone;
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
