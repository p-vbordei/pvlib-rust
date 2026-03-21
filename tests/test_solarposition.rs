use pvlib::solarposition::*;

// --- Equation of Time ---

#[test]
fn test_eot_spencer71_vernal_equinox() {
    // Around day 80 (vernal equinox), EoT should be close to -7 to -8 minutes
    let eot = equation_of_time_spencer71(80.0);
    assert!(eot.abs() < 15.0, "EoT at day 80 should be within +/-15 min, got {}", eot);
}

#[test]
fn test_eot_spencer71_summer_solstice() {
    // Around day 172 (June 21), EoT is roughly -1 to +1 minutes
    let eot = equation_of_time_spencer71(172.0);
    assert!(eot.abs() < 5.0, "EoT at day 172 should be small, got {}", eot);
}

#[test]
fn test_eot_pvcdrom_vernal_equinox() {
    let eot = equation_of_time_pvcdrom(80.0);
    assert!(eot.abs() < 15.0, "EoT PVCDROM at day 80 should be within +/-15 min, got {}", eot);
}

#[test]
fn test_eot_spencer71_vs_pvcdrom_agreement() {
    // Both methods should roughly agree across the year
    for doy in [1.0, 80.0, 172.0, 266.0, 355.0] {
        let s = equation_of_time_spencer71(doy);
        let p = equation_of_time_pvcdrom(doy);
        assert!(
            (s - p).abs() < 3.0,
            "Spencer and PVCDROM should roughly agree at doy {}: spencer={}, pvcdrom={}",
            doy, s, p
        );
    }
}

// --- Declination ---

#[test]
fn test_declination_spencer71_summer_solstice() {
    // Day ~172 (June 21), declination should be near +23.45 degrees (~0.4093 rad)
    let dec = declination_spencer71(172.0);
    let dec_deg = dec.to_degrees();
    assert!(
        (dec_deg - 23.45).abs() < 2.0,
        "Declination at summer solstice should be ~23.45 deg, got {} deg",
        dec_deg
    );
}

#[test]
fn test_declination_spencer71_winter_solstice() {
    // Day ~355 (Dec 21), declination should be near -23.45 degrees
    let dec = declination_spencer71(355.0);
    let dec_deg = dec.to_degrees();
    assert!(
        (dec_deg + 23.45).abs() < 2.0,
        "Declination at winter solstice should be ~-23.45 deg, got {} deg",
        dec_deg
    );
}

#[test]
fn test_declination_cooper69_summer_solstice() {
    let dec = declination_cooper69(172.0);
    let dec_deg = dec.to_degrees();
    assert!(
        (dec_deg - 23.45).abs() < 2.0,
        "Cooper69 declination at summer solstice should be ~23.45 deg, got {} deg",
        dec_deg
    );
}

#[test]
fn test_declination_spencer71_vs_cooper69_agreement() {
    for doy in [1.0, 80.0, 172.0, 266.0, 355.0] {
        let s = declination_spencer71(doy).to_degrees();
        let c = declination_cooper69(doy).to_degrees();
        assert!(
            (s - c).abs() < 3.0,
            "Spencer and Cooper declinations should roughly agree at doy {}: {} vs {}",
            doy, s, c
        );
    }
}

// --- Hour Angle ---

#[test]
fn test_hour_angle_solar_noon() {
    // At solar noon (12h UTC) at prime meridian (lon=0) with zero EoT,
    // hour angle should be 0
    let ha = hour_angle(12.0, 0.0, 0.0);
    assert!((ha - 0.0).abs() < 1e-10, "Hour angle at solar noon should be 0, got {}", ha);
}

#[test]
fn test_hour_angle_morning() {
    // 6h UTC, lon=0, EoT=0 -> HA = 15*(6-12) = -90 degrees
    let ha = hour_angle(6.0, 0.0, 0.0);
    assert!((ha - (-90.0)).abs() < 1e-10, "Hour angle at 6h should be -90, got {}", ha);
}

#[test]
fn test_hour_angle_with_longitude() {
    // 12h UTC, lon=15 (1h east), EoT=0 -> HA = 15*(12-12) + 15 = 15 degrees
    let ha = hour_angle(12.0, 15.0, 0.0);
    assert!((ha - 15.0).abs() < 1e-10, "Hour angle with lon=15 should be 15, got {}", ha);
}

// --- Solar Zenith Analytical ---

#[test]
fn test_solar_zenith_at_noon_equinox_equator() {
    // At equator, solar noon, equinox: zenith should be ~0
    let lat = 0.0_f64.to_radians();
    let ha = 0.0_f64.to_radians();
    let dec = 0.0_f64.to_radians();
    let z = solar_zenith_analytical(lat, ha, dec);
    assert!(z.to_degrees().abs() < 1e-6, "Zenith should be 0 at equator noon equinox, got {} deg", z.to_degrees());
}

#[test]
fn test_solar_zenith_at_noon_summer_solstice_tropic() {
    // At 23.45N, solar noon, summer solstice (dec=23.45): zenith ~0
    let lat = 23.45_f64.to_radians();
    let ha = 0.0_f64.to_radians();
    let dec = 23.45_f64.to_radians();
    let z = solar_zenith_analytical(lat, ha, dec);
    assert!(z.to_degrees().abs() < 1.0, "Zenith should be ~0 at Tropic of Cancer on solstice, got {} deg", z.to_degrees());
}

// --- Solar Azimuth Analytical ---

#[test]
fn test_solar_azimuth_at_noon() {
    // At noon (ha=0) in northern hemisphere, sun is due south -> azimuth = pi (180 deg)
    let lat = 45.0_f64.to_radians();
    let ha = 0.0_f64.to_radians();
    let dec = 10.0_f64.to_radians();
    let z = solar_zenith_analytical(lat, ha, dec);
    let az = solar_azimuth_analytical(lat, ha, dec, z);
    assert!(
        (az.to_degrees() - 180.0).abs() < 1.0,
        "Azimuth at solar noon should be ~180 (south), got {} deg",
        az.to_degrees()
    );
}

#[test]
fn test_solar_azimuth_morning_east() {
    // In morning (ha < 0), sun should be east of south (azimuth < 180)
    let lat = 45.0_f64.to_radians();
    let ha = (-45.0_f64).to_radians();
    let dec = 10.0_f64.to_radians();
    let z = solar_zenith_analytical(lat, ha, dec);
    let az = solar_azimuth_analytical(lat, ha, dec, z);
    assert!(
        az.to_degrees() < 180.0,
        "Morning azimuth should be < 180 (east of south), got {} deg",
        az.to_degrees()
    );
}

// --- Sun Rise/Set/Transit ---

#[test]
fn test_sunrise_sunset_equinox_equator() {
    // At equinox (dec=0), equator: sunrise ~6h, sunset ~18h, transit ~12h (local)
    let result = sun_rise_set_transit_geometric(0.0, 0.0, 0.0, 0.0, 0.0);
    let srs = result.expect("Should return Some for equinox at equator");
    assert!((srs.sunrise - 6.0).abs() < 0.5, "Sunrise should be ~6h, got {}", srs.sunrise);
    assert!((srs.sunset - 18.0).abs() < 0.5, "Sunset should be ~18h, got {}", srs.sunset);
    assert!((srs.transit - 12.0).abs() < 0.5, "Transit should be ~12h, got {}", srs.transit);
}

#[test]
fn test_sunrise_sunset_polar_night() {
    // At 80N with winter solstice declination (-23.45 deg), sun never rises
    let dec = (-23.45_f64).to_radians();
    let result = sun_rise_set_transit_geometric(80.0, 0.0, dec, 0.0, 0.0);
    assert!(result.is_none(), "Should be None for polar night");
}

// --- Earth-Sun Distance ---

#[test]
fn test_earthsun_distance_perihelion() {
    // Near perihelion (day ~3, early January), distance should be ~0.983 AU
    let d = nrel_earthsun_distance(3.0);
    assert!(
        (d - 0.983).abs() < 0.01,
        "Distance near perihelion should be ~0.983 AU, got {}",
        d
    );
}

#[test]
fn test_earthsun_distance_aphelion() {
    // Near aphelion (day ~185, early July), distance should be ~1.017 AU
    let d = nrel_earthsun_distance(185.0);
    assert!(
        (d - 1.017).abs() < 0.01,
        "Distance near aphelion should be ~1.017 AU, got {}",
        d
    );
}

#[test]
fn test_earthsun_distance_range() {
    // Distance should always be between 0.98 and 1.02 AU
    for doy in 1..=365 {
        let d = nrel_earthsun_distance(doy as f64);
        assert!(
            d > 0.98 && d < 1.02,
            "Distance at doy {} out of range: {}",
            doy, d
        );
    }
}
