/// Tests verifying that every batch function produces results identical to
/// calling the corresponding scalar function on each element individually.
///
/// All assertions use exact f64 equality (`==`) because the batch wrappers
/// must call the exact same scalar function -- no approximation is involved.

use chrono::TimeZone;
use chrono_tz::US::Eastern;

use pvlib::batch;
use pvlib::location::Location;
use pvlib::{atmosphere, clearsky, iam, inverter, irradiance, solarposition, temperature};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn tucson() -> Location {
    Location::new(32.2, -110.9, Eastern, 700.0, "Tucson")
}

/// Seven distinct daytime timestamps spanning morning to late afternoon.
fn test_times() -> Vec<chrono::DateTime<chrono_tz::Tz>> {
    vec![
        Eastern.with_ymd_and_hms(2020, 6, 15, 7, 0, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 6, 15, 9, 0, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 6, 15, 10, 30, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 6, 15, 12, 0, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 6, 15, 14, 0, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 6, 15, 16, 0, 0).unwrap(),
        Eastern.with_ymd_and_hms(2020, 6, 15, 17, 30, 0).unwrap(),
    ]
}

// ---------------------------------------------------------------------------
// solar_position_batch vs solarposition::get_solarposition
// ---------------------------------------------------------------------------

#[test]
fn test_solar_position_batch_vs_scalar() {
    let loc = tucson();
    let times = test_times();

    let (batch_zen, batch_azi, batch_elev) =
        batch::solar_position_batch(&loc, &times).unwrap();

    for (i, t) in times.iter().enumerate() {
        let scalar = solarposition::get_solarposition(&loc, *t).unwrap();
        assert_eq!(
            batch_zen[i], scalar.zenith,
            "zenith mismatch at index {}: batch={} scalar={}",
            i, batch_zen[i], scalar.zenith
        );
        assert_eq!(
            batch_azi[i], scalar.azimuth,
            "azimuth mismatch at index {}: batch={} scalar={}",
            i, batch_azi[i], scalar.azimuth
        );
        assert_eq!(
            batch_elev[i], scalar.elevation,
            "elevation mismatch at index {}: batch={} scalar={}",
            i, batch_elev[i], scalar.elevation
        );
    }
}

// ---------------------------------------------------------------------------
// airmass_relative_batch vs atmosphere::get_relative_airmass
// ---------------------------------------------------------------------------

#[test]
fn test_airmass_relative_batch_vs_scalar() {
    let zeniths = vec![5.0, 15.0, 25.5, 38.0, 52.3, 67.0, 78.5, 85.0];

    let batch_result = batch::airmass_relative_batch(&zeniths);

    for (i, z) in zeniths.iter().enumerate() {
        let scalar = atmosphere::get_relative_airmass(*z);
        assert_eq!(
            batch_result[i], scalar,
            "airmass_relative mismatch at zenith {}: batch={} scalar={}",
            z, batch_result[i], scalar
        );
    }
}

// ---------------------------------------------------------------------------
// airmass_absolute_batch vs atmosphere::get_absolute_airmass
// ---------------------------------------------------------------------------

#[test]
fn test_airmass_absolute_batch_vs_scalar() {
    let am_rel = vec![1.0, 1.15, 1.55, 2.0, 3.5, 5.7, 10.0];
    let pressure = 95000.0; // typical high-altitude station

    let batch_result = batch::airmass_absolute_batch(&am_rel, pressure);

    for (i, am) in am_rel.iter().enumerate() {
        let scalar = atmosphere::get_absolute_airmass(*am, pressure);
        assert_eq!(
            batch_result[i], scalar,
            "airmass_absolute mismatch at am_rel {}: batch={} scalar={}",
            am, batch_result[i], scalar
        );
    }
}

// ---------------------------------------------------------------------------
// ineichen_batch vs clearsky::ineichen
// ---------------------------------------------------------------------------

#[test]
fn test_ineichen_batch_vs_scalar() {
    let zeniths = vec![10.0, 22.5, 35.0, 48.0, 60.0, 72.0, 80.0];
    let am_abs = vec![1.02, 1.10, 1.22, 1.50, 2.00, 3.10, 5.60];
    let linke_turbidity = 3.5;
    let altitude = 500.0;

    let (batch_ghi, batch_dni, batch_dhi) =
        batch::ineichen_batch(&zeniths, &am_abs, linke_turbidity, altitude);

    for i in 0..zeniths.len() {
        let scalar = clearsky::ineichen(zeniths[i], am_abs[i], linke_turbidity, altitude, 1364.0);
        assert_eq!(
            batch_ghi[i], scalar.ghi,
            "ineichen GHI mismatch at index {}: batch={} scalar={}",
            i, batch_ghi[i], scalar.ghi
        );
        assert_eq!(
            batch_dni[i], scalar.dni,
            "ineichen DNI mismatch at index {}: batch={} scalar={}",
            i, batch_dni[i], scalar.dni
        );
        assert_eq!(
            batch_dhi[i], scalar.dhi,
            "ineichen DHI mismatch at index {}: batch={} scalar={}",
            i, batch_dhi[i], scalar.dhi
        );
    }
}

// ---------------------------------------------------------------------------
// bird_batch vs clearsky::bird_default
// ---------------------------------------------------------------------------

#[test]
fn test_bird_batch_vs_scalar() {
    let zeniths = vec![12.0, 25.0, 38.0, 50.0, 63.0, 75.0, 82.0];
    let am_rel = vec![1.02, 1.10, 1.27, 1.56, 2.25, 3.85, 7.30];
    let aod380 = 0.15;
    let aod500 = 0.10;
    let precipitable_water = 1.4;

    let (batch_ghi, batch_dni, batch_dhi) =
        batch::bird_batch(&zeniths, &am_rel, aod380, aod500, precipitable_water);

    for i in 0..zeniths.len() {
        let scalar =
            clearsky::bird_default(zeniths[i], am_rel[i], aod380, aod500, precipitable_water);
        assert_eq!(
            batch_ghi[i], scalar.ghi,
            "bird GHI mismatch at index {}: batch={} scalar={}",
            i, batch_ghi[i], scalar.ghi
        );
        assert_eq!(
            batch_dni[i], scalar.dni,
            "bird DNI mismatch at index {}: batch={} scalar={}",
            i, batch_dni[i], scalar.dni
        );
        assert_eq!(
            batch_dhi[i], scalar.dhi,
            "bird DHI mismatch at index {}: batch={} scalar={}",
            i, batch_dhi[i], scalar.dhi
        );
    }
}

// ---------------------------------------------------------------------------
// aoi_batch vs irradiance::aoi
// ---------------------------------------------------------------------------

#[test]
fn test_aoi_batch_vs_scalar() {
    let surface_tilt = 30.0;
    let surface_azimuth = 180.0;
    let solar_zenith = vec![15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0];
    let solar_azimuth = vec![120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 250.0];

    let batch_result =
        batch::aoi_batch(surface_tilt, surface_azimuth, &solar_zenith, &solar_azimuth);

    for i in 0..solar_zenith.len() {
        let scalar = irradiance::aoi(
            surface_tilt,
            surface_azimuth,
            solar_zenith[i],
            solar_azimuth[i],
        );
        assert_eq!(
            batch_result[i], scalar,
            "AOI mismatch at index {}: batch={} scalar={}",
            i, batch_result[i], scalar
        );
    }
}

// ---------------------------------------------------------------------------
// extra_radiation_batch vs irradiance::get_extra_radiation
// ---------------------------------------------------------------------------

#[test]
fn test_extra_radiation_batch_vs_scalar() {
    let days = vec![1, 32, 60, 91, 121, 152, 182, 213, 244, 274];

    let batch_result = batch::extra_radiation_batch(&days);

    for (i, d) in days.iter().enumerate() {
        let scalar = irradiance::get_extra_radiation(*d);
        assert_eq!(
            batch_result[i], scalar,
            "extra_radiation mismatch at doy {}: batch={} scalar={}",
            d, batch_result[i], scalar
        );
    }
}

// ---------------------------------------------------------------------------
// erbs_batch vs irradiance::erbs
// ---------------------------------------------------------------------------

#[test]
fn test_erbs_batch_vs_scalar() {
    let ghi = vec![200.0, 400.0, 550.0, 700.0, 850.0, 300.0, 100.0];
    let zenith = vec![60.0, 45.0, 35.0, 25.0, 15.0, 50.0, 70.0];
    let day_of_year: Vec<u32> = vec![15, 80, 120, 172, 200, 260, 350];
    let dni_extra = vec![1415.0, 1390.0, 1370.0, 1366.0, 1360.0, 1370.0, 1412.0];

    let (batch_dni, batch_dhi) =
        batch::erbs_batch(&ghi, &zenith, &day_of_year, &dni_extra);

    for i in 0..ghi.len() {
        let (scalar_dni, scalar_dhi) =
            irradiance::erbs(ghi[i], zenith[i], day_of_year[i], dni_extra[i]);
        assert_eq!(
            batch_dni[i], scalar_dni,
            "erbs DNI mismatch at index {}: batch={} scalar={}",
            i, batch_dni[i], scalar_dni
        );
        assert_eq!(
            batch_dhi[i], scalar_dhi,
            "erbs DHI mismatch at index {}: batch={} scalar={}",
            i, batch_dhi[i], scalar_dhi
        );
    }
}

// ---------------------------------------------------------------------------
// disc_batch vs irradiance::disc
// ---------------------------------------------------------------------------

#[test]
fn test_disc_batch_vs_scalar() {
    let ghi = vec![150.0, 350.0, 500.0, 650.0, 800.0, 250.0];
    let solar_zenith = vec![70.0, 55.0, 40.0, 30.0, 20.0, 60.0];
    let day_of_year: Vec<i32> = vec![30, 90, 150, 180, 240, 310];
    let pressure = Some(101325.0);

    let (batch_dni, batch_kt, batch_am) =
        batch::disc_batch(&ghi, &solar_zenith, &day_of_year, pressure);

    for i in 0..ghi.len() {
        let scalar = irradiance::disc(ghi[i], solar_zenith[i], day_of_year[i], pressure);
        assert_eq!(
            batch_dni[i], scalar.dni,
            "disc DNI mismatch at index {}: batch={} scalar={}",
            i, batch_dni[i], scalar.dni
        );
        assert_eq!(
            batch_kt[i], scalar.kt,
            "disc kt mismatch at index {}: batch={} scalar={}",
            i, batch_kt[i], scalar.kt
        );
        assert_eq!(
            batch_am[i], scalar.airmass,
            "disc airmass mismatch at index {}: batch={} scalar={}",
            i, batch_am[i], scalar.airmass
        );
    }
}

// ---------------------------------------------------------------------------
// sapm_cell_temperature_batch vs temperature::sapm_cell_temperature
// ---------------------------------------------------------------------------

#[test]
fn test_sapm_cell_temperature_batch_vs_scalar() {
    let poa_global = vec![100.0, 300.0, 500.0, 700.0, 900.0, 1050.0, 400.0];
    let temp_air = vec![15.0, 20.0, 25.0, 28.0, 32.0, 35.0, 22.0];
    let wind_speed = vec![1.0, 2.0, 3.0, 1.5, 0.5, 4.0, 2.5];
    let a = -3.56;
    let b = -0.075;
    let delta_t = 3.0;
    let irrad_ref = 1000.0;

    let (batch_cell, batch_module) = batch::sapm_cell_temperature_batch(
        &poa_global,
        &temp_air,
        &wind_speed,
        a,
        b,
        delta_t,
        irrad_ref,
    );

    for i in 0..poa_global.len() {
        let (scalar_cell, scalar_module) = temperature::sapm_cell_temperature(
            poa_global[i],
            temp_air[i],
            wind_speed[i],
            a,
            b,
            delta_t,
            irrad_ref,
        );
        assert_eq!(
            batch_cell[i], scalar_cell,
            "sapm_cell_temperature cell mismatch at index {}: batch={} scalar={}",
            i, batch_cell[i], scalar_cell
        );
        assert_eq!(
            batch_module[i], scalar_module,
            "sapm_cell_temperature module mismatch at index {}: batch={} scalar={}",
            i, batch_module[i], scalar_module
        );
    }
}

// ---------------------------------------------------------------------------
// faiman_batch vs temperature::faiman
// ---------------------------------------------------------------------------

#[test]
fn test_faiman_batch_vs_scalar() {
    let poa_global = vec![50.0, 200.0, 450.0, 650.0, 800.0, 1000.0, 350.0, 120.0];
    let temp_air = vec![10.0, 15.0, 22.0, 27.0, 30.0, 34.0, 20.0, 12.0];
    let wind_speed = vec![0.5, 1.0, 2.0, 3.0, 1.5, 0.8, 4.0, 2.5];
    let u0 = 25.0;
    let u1 = 6.84;

    let batch_result =
        batch::faiman_batch(&poa_global, &temp_air, &wind_speed, u0, u1);

    for i in 0..poa_global.len() {
        let scalar =
            temperature::faiman(poa_global[i], temp_air[i], wind_speed[i], u0, u1);
        assert_eq!(
            batch_result[i], scalar,
            "faiman mismatch at index {}: batch={} scalar={}",
            i, batch_result[i], scalar
        );
    }
}

// ---------------------------------------------------------------------------
// iam_physical_batch vs iam::physical
// ---------------------------------------------------------------------------

#[test]
fn test_iam_physical_batch_vs_scalar() {
    let aoi_vals = vec![0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0];
    let n = 1.526;
    let k = 4.0;
    let l = 0.002;

    let batch_result = batch::iam_physical_batch(&aoi_vals, n, k, l);

    for (i, a) in aoi_vals.iter().enumerate() {
        let scalar = iam::physical(*a, n, k, l);
        assert_eq!(
            batch_result[i], scalar,
            "iam_physical mismatch at aoi {}: batch={} scalar={}",
            a, batch_result[i], scalar
        );
    }
}

// ---------------------------------------------------------------------------
// iam_ashrae_batch vs iam::ashrae
// ---------------------------------------------------------------------------

#[test]
fn test_iam_ashrae_batch_vs_scalar() {
    let aoi_vals = vec![0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 89.0];
    let b0 = 0.05;

    let batch_result = batch::iam_ashrae_batch(&aoi_vals, b0);

    for (i, a) in aoi_vals.iter().enumerate() {
        let scalar = iam::ashrae(*a, b0);
        assert_eq!(
            batch_result[i], scalar,
            "iam_ashrae mismatch at aoi {}: batch={} scalar={}",
            a, batch_result[i], scalar
        );
    }
}

// ---------------------------------------------------------------------------
// pvwatts_ac_batch vs inverter::pvwatts_ac
// ---------------------------------------------------------------------------

#[test]
fn test_pvwatts_ac_batch_vs_scalar() {
    let pdc_vals = vec![0.0, 500.0, 1500.0, 2500.0, 3500.0, 4500.0, 5000.0, 5500.0, 6000.0];
    let pdc0 = 5000.0;
    let eta_inv_nom = 0.96;
    let eta_inv_ref = 0.9637;

    let batch_result = batch::pvwatts_ac_batch(&pdc_vals, pdc0, eta_inv_nom, eta_inv_ref);

    for (i, p) in pdc_vals.iter().enumerate() {
        let scalar = inverter::pvwatts_ac(*p, pdc0, eta_inv_nom, eta_inv_ref);
        assert_eq!(
            batch_result[i], scalar,
            "pvwatts_ac mismatch at pdc {}: batch={} scalar={}",
            p, batch_result[i], scalar
        );
    }
}

// ---------------------------------------------------------------------------
// disc_batch with None pressure
// ---------------------------------------------------------------------------

#[test]
fn test_disc_batch_vs_scalar_no_pressure() {
    let ghi = vec![300.0, 600.0, 900.0, 450.0, 750.0];
    let solar_zenith = vec![55.0, 40.0, 20.0, 50.0, 30.0];
    let day_of_year: Vec<i32> = vec![45, 100, 172, 220, 300];

    let (batch_dni, batch_kt, batch_am) =
        batch::disc_batch(&ghi, &solar_zenith, &day_of_year, None);

    for i in 0..ghi.len() {
        let scalar = irradiance::disc(ghi[i], solar_zenith[i], day_of_year[i], None);
        assert_eq!(
            batch_dni[i], scalar.dni,
            "disc (no pressure) DNI mismatch at index {}: batch={} scalar={}",
            i, batch_dni[i], scalar.dni
        );
        assert_eq!(
            batch_kt[i], scalar.kt,
            "disc (no pressure) kt mismatch at index {}: batch={} scalar={}",
            i, batch_kt[i], scalar.kt
        );
        assert_eq!(
            batch_am[i], scalar.airmass,
            "disc (no pressure) airmass mismatch at index {}: batch={} scalar={}",
            i, batch_am[i], scalar.airmass
        );
    }
}
