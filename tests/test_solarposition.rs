use pvlib::solarposition::*;

// --- Equation of Time ---

#[test]
fn test_eot_spencer71_vernal_equinox() {
    // Around day 80 (vernal equinox), EoT should be close to -7 to -8 minutes
    let eot = equation_of_time_spencer71(80.0);
    assert!(
        eot.abs() < 15.0,
        "EoT at day 80 should be within +/-15 min, got {}",
        eot
    );
}

#[test]
fn test_eot_spencer71_summer_solstice() {
    // Around day 172 (June 21), EoT is roughly -1 to +1 minutes
    let eot = equation_of_time_spencer71(172.0);
    assert!(
        eot.abs() < 5.0,
        "EoT at day 172 should be small, got {}",
        eot
    );
}

#[test]
fn test_eot_pvcdrom_vernal_equinox() {
    let eot = equation_of_time_pvcdrom(80.0);
    assert!(
        eot.abs() < 15.0,
        "EoT PVCDROM at day 80 should be within +/-15 min, got {}",
        eot
    );
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
            doy,
            s,
            p
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
            doy,
            s,
            c
        );
    }
}

// --- Hour Angle ---

#[test]
fn test_hour_angle_solar_noon() {
    // At solar noon (12h UTC) at prime meridian (lon=0) with zero EoT,
    // hour angle should be 0
    let ha = hour_angle(12.0, 0.0, 0.0);
    assert!(
        (ha - 0.0).abs() < 1e-10,
        "Hour angle at solar noon should be 0, got {}",
        ha
    );
}

#[test]
fn test_hour_angle_morning() {
    // 6h UTC, lon=0, EoT=0 -> HA = 15*(6-12) = -90 degrees
    let ha = hour_angle(6.0, 0.0, 0.0);
    assert!(
        (ha - (-90.0)).abs() < 1e-10,
        "Hour angle at 6h should be -90, got {}",
        ha
    );
}

#[test]
fn test_hour_angle_with_longitude() {
    // 12h UTC, lon=15 (1h east), EoT=0 -> HA = 15*(12-12) + 15 = 15 degrees
    let ha = hour_angle(12.0, 15.0, 0.0);
    assert!(
        (ha - 15.0).abs() < 1e-10,
        "Hour angle with lon=15 should be 15, got {}",
        ha
    );
}

// --- Solar Zenith Analytical ---

#[test]
fn test_solar_zenith_at_noon_equinox_equator() {
    // At equator, solar noon, equinox: zenith should be ~0
    let lat = 0.0_f64.to_radians();
    let ha = 0.0_f64.to_radians();
    let dec = 0.0_f64.to_radians();
    let z = solar_zenith_analytical(lat, ha, dec);
    assert!(
        z.to_degrees().abs() < 1e-6,
        "Zenith should be 0 at equator noon equinox, got {} deg",
        z.to_degrees()
    );
}

#[test]
fn test_solar_zenith_at_noon_summer_solstice_tropic() {
    // At 23.45N, solar noon, summer solstice (dec=23.45): zenith ~0
    let lat = 23.45_f64.to_radians();
    let ha = 0.0_f64.to_radians();
    let dec = 23.45_f64.to_radians();
    let z = solar_zenith_analytical(lat, ha, dec);
    assert!(
        z.to_degrees().abs() < 1.0,
        "Zenith should be ~0 at Tropic of Cancer on solstice, got {} deg",
        z.to_degrees()
    );
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
    assert!(
        (srs.sunrise - 6.0).abs() < 0.5,
        "Sunrise should be ~6h, got {}",
        srs.sunrise
    );
    assert!(
        (srs.sunset - 18.0).abs() < 0.5,
        "Sunset should be ~18h, got {}",
        srs.sunset
    );
    assert!(
        (srs.transit - 12.0).abs() < 0.5,
        "Transit should be ~12h, got {}",
        srs.transit
    );
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
            doy,
            d
        );
    }
}

// --- Apparent (refraction-corrected) solar position ---

/// Canonical NREL SPA reference point (Reda & Andreas 2004, Table A4.1 —
/// the same case pvlib-python's test_spa.py uses):
/// 2003-10-17 19:30:30 UTC, lat 39.742476, lon -105.1786, elev 1830.14 m,
/// pressure 820 mbar, temperature 11 C.
/// Reference values: topocentric zenith (apparent) theta = 50.11162 deg,
/// true zenith theta0 = 50.12800 deg, refraction delta_e = 0.016380 deg.
///
/// As of the `solar-positioning` crate swap (replacing the earlier `spa`
/// crate, which was actually PSA/Blanco-Muriel, not NREL SPA -- see
/// solar-engine's `docs/decisions/spa-solar-position.md`), this exercises the
/// full Reda & Andreas SPA (±0.0003 deg). `get_solarposition` doesn't accept
/// a caller-supplied delta_t, so this runs through the crate's own
/// Espenak-Meeus estimate for 2003-10 (~64.51 s) rather than NREL's example's
/// exact delta_t=67 s -- measured directly, that ~2.5 s delta_t difference
/// contributes at most ~5e-5 deg of zenith/azimuth error here, negligible
/// next to SPA's own ±0.0003 deg accuracy claim. Tolerance tightened from the
/// old PSA-era 0.01 deg to 0.001 deg (roughly 20x the largest error actually
/// observed), not loosened to match whatever the code outputs.
#[test]
fn test_apparent_zenith_nrel_spa_reference() {
    use chrono::TimeZone;
    let loc = pvlib::location::Location::new(39.742476, -105.1786, chrono_tz::UTC, 1830.14, "NREL");
    let t = chrono_tz::UTC
        .with_ymd_and_hms(2003, 10, 17, 19, 30, 30)
        .unwrap();
    // 820 mbar = 82000 Pa (the API takes Pascals, like pvlib spa_python).
    let sp = get_solarposition_apparent(&loc, t, Some(82000.0), Some(11.0)).unwrap();

    assert!(
        (sp.apparent_zenith - 50.11162).abs() < 0.001,
        "apparent zenith {} should match NREL SPA reference 50.11162 deg",
        sp.apparent_zenith
    );
    assert!(
        (sp.zenith - 50.12800).abs() < 0.001,
        "true zenith {} should match NREL SPA reference 50.12800 deg",
        sp.zenith
    );
    assert!(
        (sp.azimuth - 194.340241).abs() < 0.001,
        "azimuth {} should match NREL SPA reference 194.340241 deg",
        sp.azimuth
    );
    // Refraction correction itself matches NREL's 0.016380 deg closely.
    let delta_e = sp.apparent_elevation - sp.elevation;
    assert!(
        (delta_e - 0.016380).abs() < 0.001,
        "refraction correction {} should be ~0.01638 deg",
        delta_e
    );
    // Apparent zenith is always <= true zenith while the sun is up.
    assert!(sp.apparent_zenith < sp.zenith);
}

#[test]
fn test_refraction_negligible_at_high_sun() {
    // At 60 deg elevation the SPA refraction is under ~0.01 deg.
    let delta_e = atmospheric_refraction_correction(1013.25, 12.0, 60.0, DEFAULT_ATMOS_REFRACT);
    assert!(delta_e > 0.0);
    assert!(
        delta_e < 0.01,
        "refraction at e0=60 should be ~0.006 deg, got {}",
        delta_e
    );
}

#[test]
fn test_refraction_near_horizon() {
    // At the horizon (e0 = 0) the standard-atmosphere refraction is
    // ~0.48 deg; anywhere in 0.4-0.6 deg is physically right.
    let delta_e = atmospheric_refraction_correction(1013.25, 12.0, 0.0, DEFAULT_ATMOS_REFRACT);
    assert!(
        (0.4..=0.6).contains(&delta_e),
        "refraction at the horizon should be 0.4-0.6 deg, got {}",
        delta_e
    );
}

#[test]
fn test_refraction_zero_below_spa_horizon_threshold() {
    // pvlib/NREL SPA applies the correction only when
    // e0 >= -(SUN_RADIUS + atmos_refract) = -0.83337 deg with defaults.
    let threshold = -(SUN_RADIUS + DEFAULT_ATMOS_REFRACT);
    assert!((threshold - (-0.83337)).abs() < 1e-9);

    // Below the threshold: no correction.
    let below =
        atmospheric_refraction_correction(1013.25, 12.0, threshold - 0.001, DEFAULT_ATMOS_REFRACT);
    assert_eq!(below, 0.0);
    // At/above the threshold: correction applies.
    let at = atmospheric_refraction_correction(1013.25, 12.0, threshold, DEFAULT_ATMOS_REFRACT);
    assert!(at > 0.0);
}

#[test]
fn test_apparent_defaults_match_explicit_pvlib_defaults() {
    use chrono::TimeZone;
    let loc = pvlib::location::Location::new(46.52, 21.52, chrono_tz::UTC, 100.0, "Chisineu");
    let t = chrono_tz::UTC
        .with_ymd_and_hms(2025, 6, 21, 10, 0, 0)
        .unwrap();
    let defaulted = get_solarposition_apparent(&loc, t, None, None).unwrap();
    let explicit = get_solarposition_apparent(&loc, t, Some(101325.0), Some(12.0)).unwrap();
    assert_eq!(defaulted, explicit);

    // True fields are unchanged vs the non-refracted API.
    let plain = get_solarposition(&loc, t).unwrap();
    assert_eq!(defaulted.zenith, plain.zenith);
    assert_eq!(defaulted.azimuth, plain.azimuth);
    assert_eq!(defaulted.elevation, plain.elevation);
}
