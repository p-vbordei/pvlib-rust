use pvlib::iotools::{parse_pvgis_tmy_json, parse_pvgis_hourly_json, parse_pvgis_horizon_json};

/// Sample PVGIS TMY JSON response (minimal but structurally correct).
fn sample_tmy_json() -> &'static str {
    r#"{
        "inputs": {
            "location": {
                "latitude": 45.0,
                "longitude": 8.0,
                "elevation": 250.0
            }
        },
        "outputs": {
            "months_selected": [
                {"month": 1, "year": 2010},
                {"month": 2, "year": 2008},
                {"month": 3, "year": 2012},
                {"month": 4, "year": 2009},
                {"month": 5, "year": 2011},
                {"month": 6, "year": 2007},
                {"month": 7, "year": 2013},
                {"month": 8, "year": 2006},
                {"month": 9, "year": 2014},
                {"month": 10, "year": 2010},
                {"month": 11, "year": 2008},
                {"month": 12, "year": 2012}
            ],
            "tmy_hourly": [
                {
                    "time(UTC)": "20100101:0010",
                    "G(h)": 0.0,
                    "Gb(n)": 0.0,
                    "Gd(h)": 0.0,
                    "T2m": 2.5,
                    "WS10m": 1.3,
                    "SP": 101325.0,
                    "RH": 80.0,
                    "IR(h)": 250.0,
                    "WD10m": 180.0
                },
                {
                    "time(UTC)": "20100101:0110",
                    "G(h)": 0.0,
                    "Gb(n)": 0.0,
                    "Gd(h)": 0.0,
                    "T2m": 2.3,
                    "WS10m": 1.1,
                    "SP": 101320.0,
                    "RH": 82.0,
                    "IR(h)": 248.0,
                    "WD10m": 175.0
                },
                {
                    "time(UTC)": "20100615:1210",
                    "G(h)": 850.0,
                    "Gb(n)": 700.0,
                    "Gd(h)": 200.0,
                    "T2m": 28.0,
                    "WS10m": 3.5,
                    "SP": 100500.0,
                    "RH": 45.0,
                    "IR(h)": 350.0,
                    "WD10m": 220.0
                }
            ]
        },
        "meta": {}
    }"#
}

/// Sample PVGIS hourly radiation JSON response.
fn sample_hourly_json() -> &'static str {
    r#"{
        "inputs": {
            "location": {
                "latitude": 45.0,
                "longitude": 8.0,
                "elevation": 250.0
            }
        },
        "outputs": {
            "hourly": [
                {
                    "time": "20150101:0010",
                    "G(h)": 0.0,
                    "Gb(n)": 0.0,
                    "Gd(h)": 0.0,
                    "T2m": 1.0,
                    "WS10m": 2.0,
                    "SP": 101000.0,
                    "RH": 70.0
                },
                {
                    "time": "20150701:1210",
                    "G(h)": 900.0,
                    "Gb(n)": 750.0,
                    "Gd(h)": 180.0,
                    "T2m": 30.0,
                    "WS10m": 2.5,
                    "SP": 100200.0,
                    "RH": 40.0
                }
            ]
        },
        "meta": {}
    }"#
}

/// Sample PVGIS horizon JSON response.
fn sample_horizon_json() -> &'static str {
    r#"{
        "inputs": {
            "location": {
                "latitude": 45.0,
                "longitude": 8.0,
                "elevation": 250.0
            }
        },
        "outputs": {
            "horizon_profile": [
                {"A": -180.0, "H_hor": 5.0},
                {"A": -90.0,  "H_hor": 3.0},
                {"A": 0.0,    "H_hor": 10.0},
                {"A": 90.0,   "H_hor": 4.0},
                {"A": 180.0,  "H_hor": 5.0}
            ]
        },
        "meta": {}
    }"#
}

#[test]
fn test_parse_pvgis_tmy_json_metadata() {
    let data = parse_pvgis_tmy_json(sample_tmy_json()).unwrap();
    assert!((data.metadata.latitude - 45.0).abs() < 1e-6);
    assert!((data.metadata.longitude - 8.0).abs() < 1e-6);
    assert!((data.metadata.elevation.unwrap() - 250.0).abs() < 1e-6);
}

#[test]
fn test_parse_pvgis_tmy_json_months_selected() {
    let data = parse_pvgis_tmy_json(sample_tmy_json()).unwrap();
    let months = data.metadata.months_selected.as_ref().unwrap();
    assert_eq!(months.len(), 12);
    assert_eq!(months[0].month, 1);
    assert_eq!(months[0].year, 2010);
    assert_eq!(months[11].month, 12);
    assert_eq!(months[11].year, 2012);
}

#[test]
fn test_parse_pvgis_tmy_json_records() {
    let data = parse_pvgis_tmy_json(sample_tmy_json()).unwrap();
    assert_eq!(data.records.len(), 3);

    // First record: nighttime
    let r0 = &data.records[0];
    assert_eq!(r0.time, "20100101:0010");
    assert!((r0.ghi - 0.0).abs() < 1e-6);
    assert!((r0.temp_air - 2.5).abs() < 1e-6);
    assert!((r0.pressure - 101325.0).abs() < 1e-6);
    assert!((r0.relative_humidity - 80.0).abs() < 1e-6);

    // Third record: daytime with irradiance
    let r2 = &data.records[2];
    assert!((r2.ghi - 850.0).abs() < 1e-6);
    assert!((r2.dni - 700.0).abs() < 1e-6);
    assert!((r2.dhi - 200.0).abs() < 1e-6);
    assert!((r2.temp_air - 28.0).abs() < 1e-6);
    assert!((r2.wind_speed - 3.5).abs() < 1e-6);
}

#[test]
fn test_parse_pvgis_tmy_json_optional_fields() {
    let data = parse_pvgis_tmy_json(sample_tmy_json()).unwrap();
    let r0 = &data.records[0];
    assert!((r0.infrared.unwrap() - 250.0).abs() < 1e-6);
    assert!((r0.wind_direction.unwrap() - 180.0).abs() < 1e-6);
}

#[test]
fn test_parse_pvgis_hourly_json_records() {
    let data = parse_pvgis_hourly_json(sample_hourly_json()).unwrap();
    assert_eq!(data.records.len(), 2);

    let r1 = &data.records[1];
    assert_eq!(r1.time, "20150701:1210");
    assert!((r1.ghi - 900.0).abs() < 1e-6);
    assert!((r1.dni - 750.0).abs() < 1e-6);
    assert!((r1.dhi - 180.0).abs() < 1e-6);
    assert!((r1.temp_air - 30.0).abs() < 1e-6);
    assert!((r1.wind_speed - 2.5).abs() < 1e-6);
}

#[test]
fn test_parse_pvgis_hourly_json_no_months_selected() {
    let data = parse_pvgis_hourly_json(sample_hourly_json()).unwrap();
    assert!(data.metadata.months_selected.is_none());
}

#[test]
fn test_parse_pvgis_horizon_json() {
    let result = parse_pvgis_horizon_json(sample_horizon_json()).unwrap();
    // Original azimuths: -180, -90, 0, 90, 180
    // After +180: 0, 90, 180, 270, 360
    // 360 is filtered out (duplicate north).
    assert_eq!(result.len(), 4);

    // A=-180 -> pvlib az=0 (north), elevation=5
    assert!((result[0].0 - 0.0).abs() < 1e-6);
    assert!((result[0].1 - 5.0).abs() < 1e-6);

    // A=-90 -> pvlib az=90 (east), elevation=3
    assert!((result[1].0 - 90.0).abs() < 1e-6);
    assert!((result[1].1 - 3.0).abs() < 1e-6);

    // A=0 -> pvlib az=180 (south), elevation=10
    assert!((result[2].0 - 180.0).abs() < 1e-6);
    assert!((result[2].1 - 10.0).abs() < 1e-6);

    // A=90 -> pvlib az=270 (west), elevation=4
    assert!((result[3].0 - 270.0).abs() < 1e-6);
    assert!((result[3].1 - 4.0).abs() < 1e-6);
}

#[test]
fn test_parse_pvgis_tmy_json_invalid() {
    let result = parse_pvgis_tmy_json("not valid json");
    assert!(result.is_err());
}

#[test]
fn test_parse_pvgis_hourly_json_missing_hourly() {
    let json = r#"{"inputs": {"location": {}}, "outputs": {}, "meta": {}}"#;
    let result = parse_pvgis_hourly_json(json);
    assert!(result.is_err());
}

#[test]
fn test_parse_pvgis_horizon_json_missing_profile() {
    let json = r#"{"inputs": {}, "outputs": {}, "meta": {}}"#;
    let result = parse_pvgis_horizon_json(json);
    assert!(result.is_err());
}
