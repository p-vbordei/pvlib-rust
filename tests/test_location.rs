use chrono::TimeZone;
use chrono_tz::US::Eastern;
use chrono_tz::Europe::Berlin;
use pvlib::location::{Location, lookup_altitude};

fn tucson() -> Location {
    Location::new(32.2, -110.9, Eastern, 700.0, "Tucson")
}

fn berlin() -> Location {
    Location::new(52.52, 13.405, Berlin, 34.0, "Berlin")
}

#[test]
fn test_get_solarposition_daytime() {
    let loc = tucson();
    // June 21, 2020 at noon local time
    let time = Eastern.with_ymd_and_hms(2020, 6, 21, 12, 0, 0).unwrap();
    let pos = loc.get_solarposition(time).unwrap();

    // At solar noon in Tucson in summer, zenith should be small (sun is high)
    assert!(pos.zenith > 0.0 && pos.zenith < 60.0,
        "Zenith at noon in Tucson summer should be 0-60, got {}", pos.zenith);
    assert!(pos.elevation > 30.0 && pos.elevation < 90.0,
        "Elevation at noon should be 30-90, got {}", pos.elevation);
    assert!(pos.azimuth > 0.0 && pos.azimuth < 360.0,
        "Azimuth should be 0-360, got {}", pos.azimuth);
}

#[test]
fn test_get_solarposition_nighttime() {
    let loc = tucson();
    // Midnight
    let time = Eastern.with_ymd_and_hms(2020, 6, 21, 2, 0, 0).unwrap();
    let pos = loc.get_solarposition(time).unwrap();

    assert!(pos.zenith > 90.0, "Zenith at night should be > 90, got {}", pos.zenith);
    assert!(pos.elevation < 0.0, "Elevation at night should be negative, got {}", pos.elevation);
}

#[test]
fn test_get_clearsky_daytime_ineichen() {
    let loc = tucson();
    let time = Eastern.with_ymd_and_hms(2020, 6, 21, 12, 0, 0).unwrap();
    let cs = loc.get_clearsky(time, "ineichen");

    assert!(cs.ghi > 100.0, "GHI at noon should be > 100 W/m2, got {}", cs.ghi);
    assert!(cs.dni >= 0.0, "DNI should be non-negative, got {}", cs.dni);
    assert!(cs.dhi >= 0.0, "DHI should be non-negative, got {}", cs.dhi);
}

#[test]
fn test_get_clearsky_daytime_haurwitz() {
    let loc = tucson();
    let time = Eastern.with_ymd_and_hms(2020, 6, 21, 12, 0, 0).unwrap();
    let cs = loc.get_clearsky(time, "haurwitz");

    assert!(cs.ghi > 100.0, "Haurwitz GHI at noon should be > 100 W/m2, got {}", cs.ghi);
}

#[test]
fn test_get_clearsky_nighttime() {
    let loc = tucson();
    let time = Eastern.with_ymd_and_hms(2020, 6, 21, 2, 0, 0).unwrap();
    let cs = loc.get_clearsky(time, "ineichen");

    assert!(cs.ghi == 0.0, "GHI at night should be 0, got {}", cs.ghi);
    assert!(cs.dni == 0.0, "DNI at night should be 0, got {}", cs.dni);
}

#[test]
fn test_get_clearsky_simplified_solis() {
    let loc = tucson();
    let time = Eastern.with_ymd_and_hms(2020, 6, 21, 12, 0, 0).unwrap();
    let cs = loc.get_clearsky(time, "simplified_solis");

    assert!(cs.ghi > 100.0, "Simplified Solis GHI at noon should be > 100 W/m2, got {}", cs.ghi);
}

#[test]
fn test_get_airmass_noon() {
    let loc = tucson();
    let time = Eastern.with_ymd_and_hms(2020, 6, 21, 12, 0, 0).unwrap();
    let (am_rel, am_abs) = loc.get_airmass(time);

    // At noon in summer, airmass should be relatively low (close to 1)
    assert!(am_rel > 0.9 && am_rel < 5.0,
        "Relative airmass at noon should be 0.9-5.0, got {}", am_rel);
    assert!(am_abs > 0.0 && am_abs < 5.0,
        "Absolute airmass at noon should be > 0 and < 5, got {}", am_abs);
    // Absolute should be less than relative for high altitude locations
    assert!(am_abs < am_rel,
        "Absolute AM ({}) should be < relative AM ({}) at altitude", am_abs, am_rel);
}

#[test]
fn test_get_airmass_low_sun() {
    let loc = berlin();
    // Winter morning - sun is low
    let time = Berlin.with_ymd_and_hms(2020, 12, 21, 9, 0, 0).unwrap();
    let (am_rel, _am_abs) = loc.get_airmass(time);

    if !am_rel.is_nan() {
        assert!(am_rel > 3.0,
            "Airmass for low sun should be > 3, got {}", am_rel);
    }
}

#[test]
fn test_lookup_altitude() {
    let alt = lookup_altitude(32.2, -110.9);
    assert_eq!(alt, 0.0, "Default altitude lookup should return 0.0");
}

#[test]
fn test_location_new() {
    let loc = Location::new(40.0, -74.0, Eastern, 10.0, "NYC");
    assert_eq!(loc.latitude, 40.0);
    assert_eq!(loc.longitude, -74.0);
    assert_eq!(loc.altitude, 10.0);
    assert_eq!(loc.name, "NYC");
}
