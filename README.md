# pvlib-rust

A Rust port of [pvlib-python](https://github.com/pvlib/pvlib-python) — the open-source solar photovoltaic modeling library.

pvlib-rust provides the same algorithms, accuracy, and modeling capabilities as pvlib-python, with Rust's performance and safety guarantees. It covers the full PV simulation pipeline: from solar position and clear-sky irradiance through module temperature and single-diode physics to DC/AC power output.

## Performance

| Operation | pvlib-python | pvlib-rust | Speedup |
|-----------|-------------|------------|---------|
| TMY year (8760 hours) | ~2-5 seconds | ~4 ms | **500-1250x** |
| Single timestep | ~0.5 ms | ~0.5 us | **~1000x** |
| Parallelism | Manual (multiprocessing) | Automatic (rayon) | Built-in |

## Features

- **170+ public functions** across **24 modules**
- **319 tests** with end-to-end validation
- **Batch processing** — rayon-parallelized, full TMY year (8760 hours) in ~4 ms
- Full simulation pipeline via `ModelChain` (scalar) and `BatchModelChain` (time series)
- Builder pattern API for ergonomic configuration
- Multiple model options at each step (5 transposition, 6 temperature, 5 IAM, 3 inverter)
- Weather file I/O (TMY3, EPW) and PVGIS API client
- IV curve fitting and single-diode model parameter extraction (Bishop88, De Soto)
- `#[inline]` annotations on all 99 hot-path functions for optimal loop performance
- No unsafe code

## Quick Start

Add to your `Cargo.toml`:

```toml
[dependencies]
pvlib = { git = "https://github.com/p-vbordei/pvlib-rust" }
```

### Batch Simulation (Recommended for Production)

Simulate an entire year in milliseconds using rayon-parallelized batch processing:

```rust
use chrono::TimeZone;
use chrono_tz::US::Mountain;
use pvlib::batch::{BatchModelChain, WeatherSeries};
use pvlib::location::Location;
use pvlib::irradiance::DiffuseModel;

let location = Location::new(39.74, -105.18, Mountain, 1830.0, "Golden, CO");

// Builder pattern for ergonomic configuration
let mc = BatchModelChain::pvwatts(location, 30.0, 180.0, 5000.0)
    .with_gamma_pdc(-0.004)
    .with_inverter(5000.0, 0.96)
    .with_albedo(0.2)
    .with_transposition(DiffuseModel::Perez);

// Weather data — one value per hour (8760 for a full year)
let weather = WeatherSeries {
    times,        // Vec<DateTime<Tz>>
    ghi,          // Vec<f64> — global horizontal irradiance [W/m2]
    dni,          // Vec<f64> — direct normal irradiance [W/m2]
    dhi,          // Vec<f64> — diffuse horizontal irradiance [W/m2]
    temp_air,     // Vec<f64> — ambient temperature [C]
    wind_speed,   // Vec<f64> — wind speed [m/s]
    albedo: None, // Optional<Vec<f64>> — ground albedo
};

// Full year simulation — ~4ms on modern hardware
let results = mc.run(&weather).unwrap();

println!("Annual energy: {:.0} kWh", results.total_energy_wh() / 1000.0);
println!("Peak power: {:.0} W", results.peak_power());
println!("Capacity factor: {:.1}%", results.capacity_factor(5000.0) * 100.0);

// Access per-timestep results
for i in 0..results.ac_power.len() {
    // results.solar_zenith[i], results.poa_global[i],
    // results.cell_temperature[i], results.dc_power[i], results.ac_power[i]
}
```

Individual batch functions are also available for custom pipelines:

```rust
use pvlib::batch;

// All batch functions accept &[f64] and return Vec<f64>
// Automatically parallelized across CPU cores via rayon
let am = batch::airmass_relative_batch(&zenith_vec);
let (ghi, dni, dhi) = batch::ineichen_batch(&zenith_vec, &am, 3.0, 1830.0);
let (dni, dhi) = batch::erbs_batch(&ghi, &zenith, &doy, &dni_extra);
let iam = batch::iam_physical_batch(&aoi, 1.526, 4.0, 0.002);
let pac = batch::pvwatts_ac_batch(&pdc, 5000.0, 0.96, 0.9637);
```

### Single-Timestep Simulation

For real-time or embedded use cases:

```rust
use chrono::TimeZone;
use chrono_tz::US::Mountain;
use pvlib::location::Location;
use pvlib::pvsystem::{PVSystem, Array, FixedMount};
use pvlib::modelchain::{ModelChain, WeatherInput};

let location = Location::new(39.74, -105.18, Mountain, 1830.0, "Golden, CO");
let array = Array {
    mount: Box::new(FixedMount { surface_tilt: 30.0, surface_azimuth: 180.0 }),
    nameplate_dc: 5000.0,
    gamma_pdc: -0.004,
    modules_per_string: 10,
    strings: 2,
    albedo: 0.2,
};
let system = PVSystem::new(vec![array], 5000.0);

let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);

let weather = WeatherInput {
    time: Mountain.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap(),
    ghi: Some(900.0),
    dni: Some(800.0),
    dhi: Some(150.0),
    temp_air: 30.0,
    wind_speed: 3.0,
    albedo: Some(0.2),
};

let result = mc.run_model_from_weather(&weather).unwrap();
println!("DC: {:.0} W, AC: {:.0} W, Tcell: {:.1} C",
    result.dc_power, result.ac_power, result.cell_temperature);
```

### Step-by-Step Usage

Use individual module functions for maximum control:

```rust
use pvlib::{atmosphere, clearsky, irradiance, temperature, iam, inverter};

// Solar position (NREL SPA algorithm)
let location = pvlib::location::Location::new(39.74, -105.18,
    chrono_tz::US::Mountain, 1830.0, "Golden, CO");
let time = chrono_tz::US::Mountain.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap();
let solpos = location.get_solarposition(time).unwrap();

// Airmass
let am_rel = atmosphere::get_relative_airmass(solpos.zenith);
let pressure = atmosphere::alt2pres(1830.0);
let am_abs = atmosphere::get_absolute_airmass(am_rel, pressure);

// Clear sky irradiance (Ineichen model)
let cs = location.get_clearsky(time, "ineichen");

// Irradiance decomposition (GHI -> DNI + DHI)
let (dni, dhi) = irradiance::erbs(cs.ghi, solpos.zenith, 172, 1366.1);

// Angle of incidence on tilted surface
let aoi = irradiance::aoi(30.0, 180.0, solpos.zenith, solpos.azimuth);

// Incidence angle modifier (Fresnel reflection losses)
let iam_val = iam::physical(aoi, 1.526, 4.0, 0.002);

// Cell temperature (SAPM model)
let (t_cell, _) = temperature::sapm_cell_temperature(
    800.0, 30.0, 3.0, -3.56, -0.075, 3.0, 1000.0);

// AC power (PVWatts inverter)
let pac = inverter::pvwatts_ac(4000.0, 5000.0, 0.96, 0.9637);
```

## Modules

### Core Physics

| Module | Description | Key Functions |
|--------|-------------|---------------|
| `solarposition` | Solar position (SPA + analytical) | `get_solarposition`, `declination_spencer71`, `hour_angle`, `solar_zenith_analytical`, `sun_rise_set_transit_geometric` |
| `atmosphere` | Atmospheric properties | `get_relative_airmass`, `alt2pres`, `pres2alt`, `gueymard94_pw`, `kasten96_lt`, `tdew_from_rh`, `windspeed_powerlaw` |
| `clearsky` | Clear sky irradiance models | `ineichen`, `haurwitz`, `bird`, `simplified_solis`, `detect_clearsky` |
| `irradiance` | Transposition & decomposition | `aoi`, `perez`, `haydavies`, `klucher`, `reindl`, `isotropic`, `king`, `erbs`, `disc`, `dirindex`, `boland`, `dirint`, `get_total_irradiance` |
| `location` | Location with convenience methods | `get_solarposition`, `get_clearsky`, `get_airmass`, `lookup_altitude` |

### System Modeling

| Module | Description | Key Functions |
|--------|-------------|---------------|
| `temperature` | Cell/module temperature (6 models) | `sapm_cell_temperature`, `pvsyst_cell_temperature`, `faiman`, `fuentes`, `noct_sam`, `ross`, `generic_linear` |
| `iam` | Incidence angle modifiers (5 models) | `ashrae`, `physical`, `martin_ruiz`, `schlick`, `sapm`, `martin_ruiz_diffuse`, `schlick_diffuse` |
| `pvsystem` | PV system & single-diode params | `calcparams_desoto`, `calcparams_cec`, `calcparams_pvsyst`, `sapm`, `sapm_effective_irradiance` |
| `singlediode` | Single-diode equation solvers | `bishop88`, `bishop88_mpp`, `bishop88_i_from_v`, `bishop88_v_from_i`, `estimate_voc` |
| `inverter` | Inverter models (3 models + multi-MPPT) | `pvwatts_ac`, `sandia`, `adr`, `pvwatts_multi`, `sandia_multi` |
| `tracking` | Single-axis tracker with backtracking | `singleaxis`, `calc_surface_orientation`, `calc_axis_tilt`, `calc_cross_axis_tilt` |
| `modelchain` | End-to-end simulation pipeline | `run_model_from_weather`, `run_model_from_poa`, `run_model_from_effective_irradiance`, `complete_irradiance` |

### Batch Processing

| Module | Description | Key Functions |
|--------|-------------|---------------|
| `batch` | Rayon-parallelized batch operations | `BatchModelChain`, `solar_position_batch`, `ineichen_batch`, `erbs_batch`, `disc_batch`, `perez_batch`, `total_irradiance_batch`, `sapm_cell_temperature_batch`, `pvwatts_ac_batch` |

### Advanced Features

| Module | Description | Key Functions |
|--------|-------------|---------------|
| `bifacial` | Bifacial irradiance | `get_irradiance_infinite_sheds` |
| `shading` | Row-to-row shading | `masking_angle`, `ground_angle`, `masking_angle_passias`, `projected_solar_zenith_angle`, `shaded_fraction1d` |
| `soiling` | Soiling losses | `hsu`, `kimber`, `accumulation_model` |
| `snow` | Snow coverage & losses | `fully_covered_nrel`, `coverage_nrel`, `dc_loss_nrel` |
| `spectrum` | Spectral mismatch | `spectral_mismatch_modifier`, `spectral_factor_sapm`, `spectral_factor_caballero` |
| `scaling` | Geographic smoothing (WVM) | `wvm_smoothing` |
| `albedo` | Surface albedo | `inland_water_dvoracek`, `surface_albedo` |
| `pvarray` | Module efficiency models | `pvefficiency_adr`, `huld` |
| `transformer` | Transformer losses | `simple_efficiency` |

### I/O & Data

| Module | Description | Key Functions |
|--------|-------------|---------------|
| `iotools` | Weather file I/O & APIs | `read_tmy3`, `read_epw`, `get_pvgis_tmy`, `get_pvgis_hourly`, `get_pvgis_horizon`, `retrieve_sam` |
| `ivtools` | IV curve fitting | `fit_sandia_simple`, `fit_desoto`, `rectify_iv_curve` |

## ModelChain Configuration

Both `ModelChain` (scalar) and `BatchModelChain` (time series) support configurable model selection:

| Step | Available Models |
|------|-----------------|
| **Transposition** | Isotropic, Hay-Davies, Perez, Klucher, Reindl |
| **Temperature** | SAPM, PVsyst, Faiman, Fuentes, NOCT_SAM, PVWatts |
| **IAM (AOI)** | Physical (Fresnel), ASHRAE, Martin-Ruiz, SAPM, No loss |
| **DC Model** | PVWatts |
| **AC Model** | PVWatts, Sandia, ADR |
| **Losses** | PVWatts, No loss |

**Factory constructors:**
- `ModelChain::with_pvwatts()` / `BatchModelChain::pvwatts()` — PVWatts DC/AC, Physical AOI, Perez transposition
- `ModelChain::with_sapm()` — SAPM temperature, ASHRAE AOI, Hay-Davies transposition
- `ModelChain::with_config()` — fully custom model selection

## Weather Data

### Read TMY3 files

```rust
let data = pvlib::iotools::read_tmy3("weather.csv").unwrap();
println!("Location: {} ({}, {})",
    data.metadata.name.as_deref().unwrap_or("Unknown"),
    data.metadata.latitude, data.metadata.longitude);
for record in &data.records {
    println!("GHI: {} W/m2, Temp: {} C", record.ghi, record.temp_air);
}
```

### Read EPW files

```rust
let data = pvlib::iotools::read_epw("weather.epw").unwrap();
```

### Fetch from PVGIS API

```rust
// TMY data for any location worldwide
let tmy = pvlib::iotools::get_pvgis_tmy(45.0, 8.0, None, None, None).unwrap();

// Hourly radiation data
let hourly = pvlib::iotools::get_pvgis_hourly(45.0, 8.0, 2020, 2023,
    false, None, None, None).unwrap();

// Horizon profile
let horizon = pvlib::iotools::get_pvgis_horizon(45.0, 8.0).unwrap();
```

## References

This library implements algorithms from peer-reviewed solar energy literature:

- **Solar Position**: NREL Solar Position Algorithm (Reda & Andreas, 2004), Spencer (1971)
- **Clear Sky**: Ineichen & Perez (2002), Bird & Hulstrom (1981), Haurwitz (1945), Simplified Solis (Ineichen, 2008)
- **Irradiance**: Perez et al. (1990), Hay & Davies (1980), Erbs et al. (1982), DISC (Maxwell, 1987), Boland et al. (2008)
- **Temperature**: SAPM (King et al., 2004), PVsyst, Faiman (2008), Fuentes (1987), Ross (1980), NOCT SAM
- **IAM**: Martin & Ruiz (2001), De Soto et al. (2006), Schlick approximation, ASHRAE
- **Single Diode**: Bishop (1988), De Soto et al. (2006), Newton-Raphson solver
- **Inverter**: PVWatts (Dobos, 2014), Sandia (King et al., 2007), ADR (Driesse, 2023)
- **Tracking**: Anderson & Mikofski (2020), Marion & Dobos (2013)
- **Soiling**: HSU (Coello & Boyle, 2019), Kimber et al. (2006)
- **Bifacial**: Infinite sheds model (Mikofski et al., 2019)
- **Spectrum**: Caballero et al. (2018), First Solar spectral correction

## Comparison with pvlib-python

pvlib-rust covers the core simulation pipeline of pvlib-python with significant performance advantages:

| Aspect | pvlib-python | pvlib-rust |
|--------|-------------|------------|
| **Speed** | ~2-5s / TMY year | ~4ms / TMY year |
| **Parallelism** | Manual multiprocessing | Automatic via rayon |
| **Type safety** | Runtime errors | Compile-time checks |
| **Memory** | pandas DataFrame overhead | Zero-copy slices |
| **Deployment** | Requires Python runtime | Single static binary |
| **Data structures** | pandas DataFrames | Native Rust structs/Vec |
| **IO sources** | 20+ weather data adapters | TMY3, EPW, PVGIS, SAM |
| **Spectral models** | SPECTRL2 + 6 mismatch models | 3 mismatch models |

## Building & Testing

```bash
cargo build              # Build the library
cargo test               # Run all 319 tests
cargo test --release     # Run tests with optimizations (faster batch)
cargo clippy             # Lint checks
cargo doc --open         # Generate and view API documentation
cargo bench              # Run benchmarks (if configured)
```

## License

Licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

This project is a Rust port of [pvlib-python](https://github.com/pvlib/pvlib-python), developed by the pvlib community. All credit for the underlying algorithms and models goes to the original authors and the pvlib-python contributors.
