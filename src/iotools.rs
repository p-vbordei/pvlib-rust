use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
#[cfg(feature = "pvgis")]
use std::sync::OnceLock;
#[cfg(feature = "pvgis")]
use std::time::Duration;
use serde_json::Value;

// ---------------------------------------------------------------------------
// Shared HTTP client
// ---------------------------------------------------------------------------
//
// A single `reqwest::blocking::Client` is reused across PVGIS / SAM calls
// to avoid rebuilding a TLS context on every request and — more importantly —
// to impose a bounded timeout. The previous `reqwest::blocking::get(url)`
// calls had no timeout and could hang a worker thread indefinitely on a
// stalled connection or a slow DNS resolver.

#[cfg(feature = "pvgis")]
fn http_client() -> &'static reqwest::blocking::Client {
    static CLIENT: OnceLock<reqwest::blocking::Client> = OnceLock::new();
    CLIENT.get_or_init(|| {
        reqwest::blocking::Client::builder()
            .timeout(Duration::from_secs(30))
            .user_agent(concat!("pvlib-rust/", env!("CARGO_PKG_VERSION")))
            .build()
            .expect("reqwest client build failed; rustls-tls is enabled by Cargo.toml")
    })
}

#[cfg(feature = "pvgis")]
fn http_get_text(url: &str) -> Result<String, Box<dyn Error>> {
    let resp = http_client().get(url).send()?.error_for_status()?;
    Ok(resp.text()?)
}

/// Retrieve a database from the SAM (System Advisor Model) library.
/// Modeled after `pvlib.pvsystem.retrieve_sam`.
/// 
/// Supported databases:
///  - "CEC Inverters"
///  - "CEC Modules"
/// 
/// Note: These files are downloaded from the NREL SAM GitHub repository.
/// 
/// # Arguments
/// * `name` - The name of the database to retrieve.
/// 
/// # Returns
/// A Vector of HashMaps, where each map corresponds to a row (usually an inverter or module)
/// keyed by the column headers.
#[cfg(feature = "pvgis")]
pub fn retrieve_sam(name: &str) -> Result<Vec<HashMap<String, String>>, Box<dyn Error>> {
    let url = match name {
        "CEC Inverters" | "cecinverter" => "https://raw.githubusercontent.com/NREL/SAM/patch/deploy/libraries/CEC%20Inverters.csv",
        "CEC Modules" | "cecmod" => "https://raw.githubusercontent.com/NREL/SAM/patch/deploy/libraries/CEC%20Modules.csv",
        _ => return Err(format!("Unknown SAM DB string. Please use 'CEC Inverters' or 'CEC Modules'. You provided: {}", name).into()),
    };

    let response = http_get_text(url)?;
    
    let mut reader = csv::ReaderBuilder::new()
        .flexible(true)
        .from_reader(response.as_bytes());
        
    let headers = reader.headers()?.clone();
    
    let mut records = Vec::new();
    for result in reader.records() {
        let record = result?;
        let mut map = HashMap::new();
        for (i, field) in record.iter().enumerate() {
            if let Some(header) = headers.get(i) {
                map.insert(header.to_string(), field.to_string());
            }
        }
        // SAM CSVs can have an initial row with units which we might want to skip.
        // If the 'Name' or similar field is blank or "Units", we can skip it.
        // For simplicity, we just return all rows. The user can filter.
        records.push(map);
    }
    
    Ok(records)
}

// ---------------------------------------------------------------------------
// Weather data types
// ---------------------------------------------------------------------------

/// A single hourly weather observation.
#[derive(Debug, Clone)]
pub struct WeatherRecord {
    /// Timestamp string as returned by PVGIS (e.g. "20050101:0010").
    pub time: String,
    /// Global horizontal irradiance (W/m²).
    pub ghi: f64,
    /// Direct normal irradiance (W/m²).
    pub dni: f64,
    /// Diffuse horizontal irradiance (W/m²).
    pub dhi: f64,
    /// Air temperature at 2 m (°C).
    pub temp_air: f64,
    /// Wind speed at 10 m (m/s).
    pub wind_speed: f64,
    /// Surface pressure (Pa or mbar depending on source).
    pub pressure: f64,
    /// Relative humidity (%).
    pub relative_humidity: f64,
    /// Infrared radiation downwards (W/m²), if available.
    pub infrared: Option<f64>,
    /// Wind direction at 10 m (°), if available.
    pub wind_direction: Option<f64>,
    /// Dew-point temperature (°C), if available.
    pub temp_dew: Option<f64>,
    /// Albedo (unitless), if available.
    pub albedo: Option<f64>,
    /// Precipitable water (cm), if available.
    pub precipitable_water: Option<f64>,
    /// Year from the data file.
    pub year: Option<i32>,
    /// Month from the data file.
    pub month: Option<u32>,
    /// Day from the data file.
    pub day: Option<u32>,
    /// Hour from the data file (0-23).
    pub hour: Option<u32>,
}

/// Metadata about the weather data source and location.
#[derive(Debug, Clone)]
pub struct WeatherMetadata {
    pub latitude: f64,
    pub longitude: f64,
    pub elevation: Option<f64>,
    /// Timezone offset from UTC (hours).
    pub tz_offset: Option<f64>,
    /// Site / station name.
    pub name: Option<String>,
    /// City name.
    pub city: Option<String>,
    /// State or province.
    pub state: Option<String>,
    /// Data source identifier.
    pub source: Option<String>,
    /// Months selected for TMY, if applicable.
    pub months_selected: Option<Vec<MonthYear>>,
    /// Additional key-value metadata from the API response.
    pub extra: HashMap<String, String>,
}

/// Month-year pair indicating which year was selected for a given month in TMY.
#[derive(Debug, Clone)]
pub struct MonthYear {
    pub month: i32,
    pub year: i32,
}

/// Container for weather time-series data plus metadata.
#[derive(Debug, Clone)]
pub struct WeatherData {
    pub records: Vec<WeatherRecord>,
    pub metadata: WeatherMetadata,
}

// ---------------------------------------------------------------------------
// PVGIS constants
// ---------------------------------------------------------------------------

#[cfg(feature = "pvgis")]
const PVGIS_BASE_URL: &str = "https://re.jrc.ec.europa.eu/api/v5_3/";

// ---------------------------------------------------------------------------
// PVGIS API functions
// ---------------------------------------------------------------------------

/// Retrieve Typical Meteorological Year (TMY) data from PVGIS.
///
/// # Arguments
/// * `latitude` – Latitude in decimal degrees (north positive).
/// * `longitude` – Longitude in decimal degrees (east positive).
/// * `outputformat` – Response format: `"json"`, `"csv"`, or `"epw"`.
/// * `startyear` – Optional first year for TMY calculation.
/// * `endyear` – Optional last year for TMY calculation.
#[cfg(feature = "pvgis")]
pub fn get_pvgis_tmy(
    latitude: f64,
    longitude: f64,
    outputformat: &str,
    startyear: Option<i32>,
    endyear: Option<i32>,
) -> Result<WeatherData, Box<dyn Error>> {
    let mut url = format!(
        "{}tmy?lat={}&lon={}&outputformat={}",
        PVGIS_BASE_URL, latitude, longitude, outputformat,
    );
    if let Some(sy) = startyear {
        url.push_str(&format!("&startyear={}", sy));
    }
    if let Some(ey) = endyear {
        url.push_str(&format!("&endyear={}", ey));
    }

    let response = http_get_text(&url)?;

    match outputformat {
        "json" => parse_pvgis_tmy_json(&response),
        _ => Err(format!("Unsupported PVGIS TMY outputformat: '{}'. Use 'json'.", outputformat).into()),
    }
}

/// Retrieve hourly solar radiation (and optionally PV power) data from PVGIS.
///
/// # Arguments
/// * `latitude` – Latitude in decimal degrees.
/// * `longitude` – Longitude in decimal degrees.
/// * `start` – First year of the time series.
/// * `end` – Last year of the time series.
/// * `pvcalculation` – If true, include estimated PV power output.
/// * `peakpower` – Nominal PV system power in kW (required if `pvcalculation` is true).
/// * `surface_tilt` – Tilt angle from horizontal (degrees).
/// * `surface_azimuth` – Orientation clockwise from north (degrees). Converted to PVGIS convention internally.
#[cfg(feature = "pvgis")]
#[allow(clippy::too_many_arguments)]
pub fn get_pvgis_hourly(
    latitude: f64,
    longitude: f64,
    start: i32,
    end: i32,
    pvcalculation: bool,
    peakpower: Option<f64>,
    surface_tilt: Option<f64>,
    surface_azimuth: Option<f64>,
) -> Result<WeatherData, Box<dyn Error>> {
    let tilt = surface_tilt.unwrap_or(0.0);
    // PVGIS uses south=0 convention; pvlib uses south=180, so subtract 180.
    let aspect = surface_azimuth.unwrap_or(180.0) - 180.0;
    let pvcalc_int = if pvcalculation { 1 } else { 0 };

    let mut url = format!(
        "{}seriescalc?lat={}&lon={}&startyear={}&endyear={}&pvcalculation={}&angle={}&aspect={}&outputformat=json",
        PVGIS_BASE_URL, latitude, longitude, start, end, pvcalc_int, tilt, aspect,
    );
    if let Some(pp) = peakpower {
        url.push_str(&format!("&peakpower={}", pp));
    }

    let response = http_get_text(&url)?;
    parse_pvgis_hourly_json(&response)
}

/// Retrieve horizon profile data from PVGIS.
///
/// Returns a vector of (azimuth, elevation) pairs where azimuth follows
/// the pvlib convention (north=0, clockwise, south=180).
#[cfg(feature = "pvgis")]
pub fn get_pvgis_horizon(
    latitude: f64,
    longitude: f64,
) -> Result<Vec<(f64, f64)>, Box<dyn Error>> {
    let url = format!(
        "{}printhorizon?lat={}&lon={}&outputformat=json",
        PVGIS_BASE_URL, latitude, longitude,
    );

    let response = http_get_text(&url)?;
    parse_pvgis_horizon_json(&response)
}

// ---------------------------------------------------------------------------
// JSON parsing helpers (public for testing)
// ---------------------------------------------------------------------------

/// Parse PVGIS TMY JSON response into `WeatherData`.
pub fn parse_pvgis_tmy_json(json_str: &str) -> Result<WeatherData, Box<dyn Error>> {
    let root: Value = serde_json::from_str(json_str)?;

    // Location metadata
    let inputs = &root["inputs"]["location"];
    let latitude = inputs["latitude"].as_f64().unwrap_or(0.0);
    let longitude = inputs["longitude"].as_f64().unwrap_or(0.0);
    let elevation = inputs["elevation"].as_f64();

    // Months selected
    let months_selected = root["outputs"]["months_selected"]
        .as_array()
        .map(|arr| {
            arr.iter()
                .map(|m| MonthYear {
                    month: m["month"].as_i64().unwrap_or(0) as i32,
                    year: m["year"].as_i64().unwrap_or(0) as i32,
                })
                .collect()
        });

    // Hourly records
    let hourly = root["outputs"]["tmy_hourly"]
        .as_array()
        .ok_or("Missing outputs.tmy_hourly in PVGIS TMY response")?;

    let records: Vec<WeatherRecord> = hourly
        .iter()
        .map(|h| WeatherRecord {
            time: h["time(UTC)"].as_str().unwrap_or("").to_string(),
            ghi: h["G(h)"].as_f64().unwrap_or(0.0),
            dni: h["Gb(n)"].as_f64().unwrap_or(0.0),
            dhi: h["Gd(h)"].as_f64().unwrap_or(0.0),
            temp_air: h["T2m"].as_f64().unwrap_or(0.0),
            wind_speed: h["WS10m"].as_f64().unwrap_or(0.0),
            pressure: h["SP"].as_f64().unwrap_or(0.0),
            relative_humidity: h["RH"].as_f64().unwrap_or(0.0),
            infrared: h["IR(h)"].as_f64(),
            wind_direction: h["WD10m"].as_f64(),
            temp_dew: None,
            albedo: None,
            precipitable_water: None,
            year: None,
            month: None,
            day: None,
            hour: None,
        })
        .collect();

    Ok(WeatherData {
        records,
        metadata: WeatherMetadata {
            latitude,
            longitude,
            elevation,
            tz_offset: None,
            name: None,
            city: None,
            state: None,
            source: None,
            months_selected,
            extra: HashMap::new(),
        },
    })
}

/// Parse PVGIS hourly radiation JSON response into `WeatherData`.
pub fn parse_pvgis_hourly_json(json_str: &str) -> Result<WeatherData, Box<dyn Error>> {
    let root: Value = serde_json::from_str(json_str)?;

    let inputs = &root["inputs"]["location"];
    let latitude = inputs["latitude"].as_f64().unwrap_or(0.0);
    let longitude = inputs["longitude"].as_f64().unwrap_or(0.0);
    let elevation = inputs["elevation"].as_f64();

    let hourly = root["outputs"]["hourly"]
        .as_array()
        .ok_or("Missing outputs.hourly in PVGIS hourly response")?;

    let records: Vec<WeatherRecord> = hourly
        .iter()
        .map(|h| WeatherRecord {
            time: h["time"].as_str().unwrap_or("").to_string(),
            ghi: h["G(h)"].as_f64().unwrap_or(0.0),
            dni: h["Gb(n)"].as_f64().unwrap_or(0.0),
            dhi: h["Gd(h)"].as_f64().unwrap_or(0.0),
            temp_air: h["T2m"].as_f64().unwrap_or(0.0),
            wind_speed: h["WS10m"].as_f64().unwrap_or(0.0),
            pressure: h["SP"].as_f64().unwrap_or(0.0),
            relative_humidity: h["RH"].as_f64().unwrap_or(0.0),
            infrared: h["IR(h)"].as_f64(),
            wind_direction: h["WD10m"].as_f64(),
            temp_dew: None,
            albedo: None,
            precipitable_water: None,
            year: None,
            month: None,
            day: None,
            hour: None,
        })
        .collect();

    Ok(WeatherData {
        records,
        metadata: WeatherMetadata {
            latitude,
            longitude,
            elevation,
            tz_offset: None,
            name: None,
            city: None,
            state: None,
            source: None,
            months_selected: None,
            extra: HashMap::new(),
        },
    })
}

/// Parse PVGIS horizon JSON response into a vector of (azimuth, elevation) pairs.
/// Azimuths are converted to pvlib convention (north=0, clockwise).
pub fn parse_pvgis_horizon_json(json_str: &str) -> Result<Vec<(f64, f64)>, Box<dyn Error>> {
    let root: Value = serde_json::from_str(json_str)?;

    let profile = root["outputs"]["horizon_profile"]
        .as_array()
        .ok_or("Missing outputs.horizon_profile in PVGIS horizon response")?;

    let mut result: Vec<(f64, f64)> = profile
        .iter()
        .map(|p| {
            let az = p["A"].as_f64().unwrap_or(0.0);
            let el = p["H_hor"].as_f64().unwrap_or(0.0);
            // PVGIS uses south=0; convert to pvlib north=0 by adding 180.
            let az_pvlib = az + 180.0;
            (az_pvlib, el)
        })
        .collect();

    // Remove the duplicate north point (360 == 0).
    result.retain(|&(az, _)| az < 360.0);

    Ok(result)
}

// ---------------------------------------------------------------------------
// TMY3 / EPW file readers
// ---------------------------------------------------------------------------

/// Read a TMY3 CSV file.
///
/// TMY3 files have two header lines: the first contains site metadata
/// (USAF, Name, State, TZ, latitude, longitude, altitude), the second
/// contains column names. Data rows follow.
///
/// Modeled after `pvlib.iotools.read_tmy3`.
pub fn read_tmy3(filepath: &str) -> Result<WeatherData, Box<dyn Error>> {
    let file = File::open(filepath)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // First line: metadata
    let meta_line = lines.next().ok_or("TMY3 file is empty")??;
    let meta_fields: Vec<&str> = meta_line.split(',').collect();
    if meta_fields.len() < 7 {
        return Err("TMY3 metadata line has fewer than 7 fields".into());
    }
    let metadata = WeatherMetadata {
        latitude: meta_fields[4].trim().parse()?,
        longitude: meta_fields[5].trim().parse()?,
        elevation: Some(meta_fields[6].trim().parse()?),
        tz_offset: Some(meta_fields[3].trim().parse()?),
        name: Some(meta_fields[1].trim().to_string()),
        city: Some(meta_fields[1].trim().to_string()),
        state: Some(meta_fields[2].trim().to_string()),
        source: Some(format!("USAF {}", meta_fields[0].trim())),
        months_selected: None,
        extra: HashMap::new(),
    };

    // Second line: column headers
    let header_line = lines.next().ok_or("TMY3 file missing column header line")??;
    let headers: Vec<String> = header_line.split(',').map(|s| s.trim().to_string()).collect();

    let idx = |name: &str| -> Result<usize, Box<dyn Error>> {
        headers.iter().position(|h| h == name)
            .ok_or_else(|| format!("TMY3 column '{}' not found", name).into())
    };
    let i_date = idx("Date (MM/DD/YYYY)")?;
    let i_time = idx("Time (HH:MM)")?;
    let i_ghi = idx("GHI (W/m^2)")?;
    let i_dni = idx("DNI (W/m^2)")?;
    let i_dhi = idx("DHI (W/m^2)")?;
    let i_temp = idx("Dry-bulb (C)")?;
    let i_dew = idx("Dew-point (C)")?;
    let i_rh = idx("RHum (%)")?;
    let i_pres = idx("Pressure (mbar)")?;
    let i_wdir = idx("Wdir (degrees)")?;
    let i_wspd = idx("Wspd (m/s)")?;
    let i_alb = idx("Alb (unitless)")?;
    let i_pwat = idx("Pwat (cm)")?;

    let mut records = Vec::new();
    for line_result in lines {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();

        let get_field = |i: usize| -> Result<&str, Box<dyn Error>> {
            fields.get(i).copied().ok_or_else(|| format!("missing TMY3 field {}", i).into())
        };

        // Parse date MM/DD/YYYY
        let date_field = get_field(i_date)?;
        let date_parts: Vec<&str> = date_field.split('/').collect();
        if date_parts.len() != 3 {
            return Err(format!("malformed TMY3 date field: {:?}", date_field).into());
        }
        let month: u32 = date_parts[0].parse()?;
        let day: u32 = date_parts[1].parse()?;
        let year: i32 = date_parts[2].parse()?;

        // Parse time HH:MM — TMY3 uses 1-24, 24:00 means midnight next day
        let time_field = get_field(i_time)?;
        let time_parts: Vec<&str> = time_field.split(':').collect();
        if time_parts.is_empty() {
            return Err(format!("malformed TMY3 time field: {:?}", time_field).into());
        }
        let raw_hour: u32 = time_parts[0].parse()?;
        let hour = raw_hour % 24;

        let parse_f64 = |i: usize| -> Result<f64, Box<dyn Error>> {
            get_field(i)?.trim().parse::<f64>().map_err(|e| e.into())
        };

        let time_str = format!("{:04}{:02}{:02}:{:02}{:02}",
            year, month, day, hour, 0);

        records.push(WeatherRecord {
            time: time_str,
            ghi: parse_f64(i_ghi)?,
            dni: parse_f64(i_dni)?,
            dhi: parse_f64(i_dhi)?,
            temp_air: parse_f64(i_temp)?,
            wind_speed: parse_f64(i_wspd)?,
            pressure: parse_f64(i_pres)?,
            relative_humidity: parse_f64(i_rh)?,
            infrared: None,
            wind_direction: Some(parse_f64(i_wdir)?),
            temp_dew: Some(parse_f64(i_dew)?),
            albedo: Some(parse_f64(i_alb)?),
            precipitable_water: Some(parse_f64(i_pwat)?),
            year: Some(year),
            month: Some(month),
            day: Some(day),
            hour: Some(hour),
        });
    }

    Ok(WeatherData { metadata, records })
}

/// Read an EPW (EnergyPlus Weather) file.
///
/// EPW files have 8 header lines (LOCATION, DESIGN CONDITIONS, etc.)
/// followed by hourly data rows. The LOCATION line provides site metadata.
///
/// Modeled after `pvlib.iotools.read_epw`.
pub fn read_epw(filepath: &str) -> Result<WeatherData, Box<dyn Error>> {
    let file = File::open(filepath)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // First line: LOCATION,city,state,country,data_type,WMO_code,lat,lon,tz,elev
    let loc_line = lines.next().ok_or("EPW file is empty")??;
    let loc_fields: Vec<&str> = loc_line.split(',').collect();
    if loc_fields.len() < 10 {
        return Err("EPW LOCATION line has fewer than 10 fields".into());
    }
    let metadata = WeatherMetadata {
        latitude: loc_fields[6].trim().parse()?,
        longitude: loc_fields[7].trim().parse()?,
        elevation: Some(loc_fields[9].trim().parse()?),
        tz_offset: Some(loc_fields[8].trim().parse()?),
        name: Some(loc_fields[1].trim().to_string()),
        city: Some(loc_fields[1].trim().to_string()),
        state: Some(loc_fields[2].trim().to_string()),
        source: Some(loc_fields[4].trim().to_string()),
        months_selected: None,
        extra: HashMap::new(),
    };

    // Skip remaining 7 header lines
    for _ in 0..7 {
        lines.next().ok_or("EPW file has fewer than 8 header lines")??;
    }

    // Data columns (0-indexed):
    //  0=year, 1=month, 2=day, 3=hour, 4=minute, 5=data_source,
    //  6=temp_air, 7=temp_dew, 8=rh, 9=pressure,
    //  10=etr, 11=etrn, 12=ghi_infrared, 13=ghi, 14=dni, 15=dhi,
    //  ...20=wind_dir, 21=wind_speed, ...28=precipitable_water,
    //  ...32=albedo
    let mut records = Vec::new();
    for line_result in lines {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 29 {
            continue;
        }

        let get_field = |i: usize| -> Result<&str, Box<dyn Error>> {
            fields.get(i).copied().ok_or_else(|| format!("missing EPW field {}", i).into())
        };
        let parse_f64 = |i: usize| -> Result<f64, Box<dyn Error>> {
            get_field(i)?.trim().parse::<f64>().map_err(|e| e.into())
        };

        let try_parse_f64 = |i: usize| -> Option<f64> {
            fields.get(i).and_then(|s| s.trim().parse::<f64>().ok())
        };

        // EPW hour is 1-24; convert to 0-23
        let raw_hour: u32 = get_field(3)?.trim().parse()?;
        let hour = if raw_hour == 0 { 0 } else { raw_hour - 1 };
        let year: i32 = get_field(0)?.trim().parse()?;
        let month: u32 = get_field(1)?.trim().parse()?;
        let day: u32 = get_field(2)?.trim().parse()?;

        let time_str = format!("{:04}{:02}{:02}:{:02}{:02}",
            year, month, day, hour, 0);

        // EPW standard columns (0-indexed):
        //  0-5: year, month, day, hour, minute, data_source
        //  6: temp_air, 7: temp_dew, 8: rh, 9: pressure
        //  10: etr, 11: etrn, 12: ghi_infrared, 13: ghi, 14: dni, 15: dhi
        //  16-19: illuminance fields, 20: wind_dir, 21: wind_speed
        //  22-27: sky cover, visibility, ceiling, weather obs/codes
        //  28: precipitable_water, 29: aod, 30: snow_depth, 31: days_since_snow
        //  32: albedo, 33: liquid_precip_depth, 34: liquid_precip_qty
        records.push(WeatherRecord {
            time: time_str,
            temp_air: parse_f64(6)?,
            wind_speed: parse_f64(21)?,
            pressure: parse_f64(9)?,
            relative_humidity: parse_f64(8)?,
            ghi: parse_f64(13)?,
            dni: parse_f64(14)?,
            dhi: parse_f64(15)?,
            infrared: Some(parse_f64(12)?),
            wind_direction: Some(parse_f64(20)?),
            temp_dew: Some(parse_f64(7)?),
            precipitable_water: try_parse_f64(28),
            albedo: try_parse_f64(32),
            year: Some(year),
            month: Some(month),
            day: Some(day),
            hour: Some(hour),
        });
    }

    Ok(WeatherData { metadata, records })
}
