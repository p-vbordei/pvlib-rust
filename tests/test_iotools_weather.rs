use std::io::Write;
use tempfile::NamedTempFile;
use pvlib::iotools::{read_tmy3, read_epw};

#[test]
fn test_read_tmy3_basic() {
    let tmy3_data = "\
722780,Test Station,AZ,-7,33.45,-111.98,337
Date (MM/DD/YYYY),Time (HH:MM),GHI (W/m^2),DNI (W/m^2),DHI (W/m^2),Dry-bulb (C),Dew-point (C),RHum (%),Pressure (mbar),Wdir (degrees),Wspd (m/s),Alb (unitless),Pwat (cm)
01/01/1988,01:00,0,0,0,5.6,2.1,78,960,180,3.1,0.16,1.2
01/01/1988,02:00,0,0,0,4.8,1.5,75,961,200,2.8,0.16,1.1
01/01/1988,03:00,50,100,25,6.2,2.8,72,959,190,3.5,0.17,1.3
";
    let mut tmpfile = NamedTempFile::new().unwrap();
    write!(tmpfile, "{}", tmy3_data).unwrap();

    let wd = read_tmy3(tmpfile.path().to_str().unwrap()).unwrap();

    // Check metadata
    assert!((wd.metadata.latitude - 33.45).abs() < 1e-6);
    assert!((wd.metadata.longitude - (-111.98)).abs() < 1e-6);
    assert!((wd.metadata.elevation.unwrap() - 337.0).abs() < 1e-6);
    assert!((wd.metadata.tz_offset.unwrap() - (-7.0)).abs() < 1e-6);
    assert_eq!(wd.metadata.name.as_deref(), Some("Test Station"));
    assert_eq!(wd.metadata.state.as_deref(), Some("AZ"));

    // Check records
    assert_eq!(wd.records.len(), 3);

    let r0 = &wd.records[0];
    assert_eq!(r0.year, Some(1988));
    assert_eq!(r0.month, Some(1));
    assert_eq!(r0.day, Some(1));
    assert_eq!(r0.hour, Some(1));
    assert!((r0.ghi - 0.0).abs() < 1e-6);
    assert!((r0.temp_air - 5.6).abs() < 1e-6);
    assert!((r0.temp_dew.unwrap() - 2.1).abs() < 1e-6);
    assert!((r0.relative_humidity - 78.0).abs() < 1e-6);
    assert!((r0.pressure - 960.0).abs() < 1e-6);
    assert!((r0.wind_direction.unwrap() - 180.0).abs() < 1e-6);
    assert!((r0.wind_speed - 3.1).abs() < 1e-6);
    assert!((r0.albedo.unwrap() - 0.16).abs() < 1e-6);
    assert!((r0.precipitable_water.unwrap() - 1.2).abs() < 1e-6);

    let r2 = &wd.records[2];
    assert!((r2.ghi - 50.0).abs() < 1e-6);
    assert!((r2.dni - 100.0).abs() < 1e-6);
    assert!((r2.dhi - 25.0).abs() < 1e-6);
}

#[test]
fn test_read_tmy3_midnight_convention() {
    // TMY3 24:00 should become hour 0
    let tmy3_data = "\
722780,Test,AZ,-7,33.45,-111.98,337
Date (MM/DD/YYYY),Time (HH:MM),GHI (W/m^2),DNI (W/m^2),DHI (W/m^2),Dry-bulb (C),Dew-point (C),RHum (%),Pressure (mbar),Wdir (degrees),Wspd (m/s),Alb (unitless),Pwat (cm)
12/31/1988,24:00,0,0,0,3.0,1.0,80,1013,0,1.0,0.15,1.0
";
    let mut tmpfile = NamedTempFile::new().unwrap();
    write!(tmpfile, "{}", tmy3_data).unwrap();

    let wd = read_tmy3(tmpfile.path().to_str().unwrap()).unwrap();
    assert_eq!(wd.records[0].hour, Some(0));
}

#[test]
fn test_read_epw_basic() {
    // Build a minimal EPW file: 8 header lines + data rows
    // 35 columns per data row
    let mut epw = String::new();
    // Line 1: LOCATION
    epw.push_str("LOCATION,Phoenix,AZ,USA,TMY3,722780,33.45,-112.07,-7,337\n");
    // Lines 2-8: other header lines (content doesn't matter for parsing)
    epw.push_str("DESIGN CONDITIONS,0\n");
    epw.push_str("TYPICAL/EXTREME PERIODS,0\n");
    epw.push_str("GROUND TEMPERATURES,0\n");
    epw.push_str("HOLIDAYS/DAYLIGHT SAVINGS,No,0,0,0\n");
    epw.push_str("COMMENTS 1,test\n");
    epw.push_str("COMMENTS 2,test\n");
    epw.push_str("DATA PERIODS,1,1,Data,Sunday,1/1,12/31\n");

    // Data row: 35 comma-separated fields
    // year,month,day,hour,minute,source,temp_air,temp_dew,rh,pressure,
    // etr,etrn,ghi_ir,ghi,dni,dhi,gh_illum,dn_illum,dh_illum,zenith_lum,
    // wdir,wspd,sky_cover,opq_cover,vis,ceil,pwo,pwc,precip_water,aod,
    // snow,days_snow,albedo,liq_depth,liq_qty
    epw.push_str("1988,1,1,1,0,?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0?0,5.6,2.1,78,96000,0,0,300,0,0,0,0,0,0,0,180,3.1,0,0,20,77777,0,0,1.2,0.05,0,0,0.16,0,0\n");
    epw.push_str("1988,1,1,2,0,?0,4.8,1.5,75,96100,0,0,290,100,200,50,0,0,0,0,200,2.8,0,0,20,77777,0,0,1.1,0.05,0,0,0.15,0,0\n");

    let mut tmpfile = NamedTempFile::new().unwrap();
    write!(tmpfile, "{}", epw).unwrap();

    let wd = read_epw(tmpfile.path().to_str().unwrap()).unwrap();

    // Check metadata
    assert!((wd.metadata.latitude - 33.45).abs() < 1e-6);
    assert!((wd.metadata.longitude - (-112.07)).abs() < 1e-6);
    assert!((wd.metadata.elevation.unwrap() - 337.0).abs() < 1e-6);
    assert!((wd.metadata.tz_offset.unwrap() - (-7.0)).abs() < 1e-6);
    assert_eq!(wd.metadata.city.as_deref(), Some("Phoenix"));
    assert_eq!(wd.metadata.state.as_deref(), Some("AZ"));

    // Check records
    assert_eq!(wd.records.len(), 2);

    // EPW hour 1 -> internal hour 0
    let r0 = &wd.records[0];
    assert_eq!(r0.year, Some(1988));
    assert_eq!(r0.month, Some(1));
    assert_eq!(r0.day, Some(1));
    assert_eq!(r0.hour, Some(0));
    assert!((r0.temp_air - 5.6).abs() < 1e-6);
    assert!((r0.temp_dew.unwrap() - 2.1).abs() < 1e-6);
    assert!((r0.relative_humidity - 78.0).abs() < 1e-6);
    assert!((r0.pressure - 96000.0).abs() < 1e-6);
    assert!((r0.ghi - 0.0).abs() < 1e-6);
    assert!((r0.wind_direction.unwrap() - 180.0).abs() < 1e-6);
    assert!((r0.wind_speed - 3.1).abs() < 1e-6);
    assert!((r0.albedo.unwrap() - 0.16).abs() < 1e-6);
    assert!((r0.precipitable_water.unwrap() - 1.2).abs() < 1e-6);

    // Second record: EPW hour 2 -> internal hour 1
    let r1 = &wd.records[1];
    assert_eq!(r1.hour, Some(1));
    assert!((r1.ghi - 100.0).abs() < 1e-6);
    assert!((r1.dni - 200.0).abs() < 1e-6);
    assert!((r1.dhi - 50.0).abs() < 1e-6);
}

#[test]
fn test_read_tmy3_file_not_found() {
    let result = read_tmy3("/nonexistent/path/to/file.csv");
    assert!(result.is_err());
}

#[test]
fn test_read_epw_file_not_found() {
    let result = read_epw("/nonexistent/path/to/file.epw");
    assert!(result.is_err());
}
