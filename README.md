# pvlib-rust

A Rust port of [pvlib-python](https://github.com/pvlib/pvlib-python) — the open-source solar photovoltaic modeling library.

pvlib-rust provides the same algorithms, accuracy, and modeling capabilities as pvlib-python, with Rust's performance and safety guarantees. It covers the full PV simulation pipeline: from solar position and clear-sky irradiance through module temperature and single-diode physics to DC/AC power output.

## Features

- **170+ public functions** across **24 modules**
- **319 tests** with end-to-end validation
- **Batch processing** — rayon-parallelized, full TMY year (8760 hours) in ~4ms
- Full simulation pipeline via `ModelChain` and `BatchModelChain`
- Multiple model options at each step (transposition, temperature, IAM, inverter)
- Weather file I/O (TMY3, EPW) and PVGIS API client
- IV curve fitting and single-diode model parameter extraction
- `#[inline]` annotations on all hot-path functions for optimal loop performance
- No unsafe code

## Quick Start

Add to your `Cargo.toml`:

```toml
[dependencies]
pvlib = { git = "https://github.com/p-vbordei/pvlib-rust" }
```

### Full System Simulation

```rust
use chrono::TimeZone;
use chrono_tz::US::Mountain;
use pvlib::location::Location;
use pvlib::pvsystem::{PVSystem, Array, FixedMount};
use pvlib::modelchain::{ModelChain, WeatherInput};

// Define location and system
let location = Location::new(39.74, -105.18, Mountain, 1830.0, "Golden, CO");
let array = Array {
    mount: Box::new(FixedMount { surface_tilt: 30.0, surface_azimuth: 180.0 }),
    nameplate_dc: 5000.0,     // 5 kW system
    gamma_pdc: -0.004,         // -0.4%/°C temperature coefficient
    modules_per_string: 10,
    strings: 2,
    albedo: 0.2,
};
let system = PVSystem::new(vec![array], 5000.0);

// Create ModelChain with PVWatts configuration
let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);

// Run simulation
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
println!("DC Power: {:.0} W", result.dc_power);
println!("AC Power: {:.0} W", result.ac_power);
println!("Cell Temp: {:.1} °C", result.cell_temperature);
```

### Batch Simulation (Production Workloads)

Simulate an entire year in milliseconds using rayon-parallelized batch processing:

```rust
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

// Load weather data (8760 hourly timesteps)
let weather = WeatherSeries { times, ghi, dni, dhi, temp_air, wind_speed, albedo: None };

// Full year simulation — ~4ms on modern hardware
let results = mc.run(&weather).unwrap();

println!("Annual energy: {:.0} kWh", results.total_energy_wh() / 1000.0);
println!("Peak power: {:.0} W", results.peak_power());
println!("Capacity factor: {:.1}%", results.capacity_factor(5000.0) * 100.0);
```

Individual batch functions are also available for custom pipelines:

```rust
use pvlib::batch;

let zenith_vec = vec![30.0; 8760];
let am = batch::airmass_relative_batch(&zenith_vec);
let (ghi, dni, dhi) = batch::ineichen_batch(&zenith_vec, &am, 3.0, 1830.0);
```

### Step-by-Step Usage

```rust
use pvlib::{solarposition, atmosphere, clearsky, irradiance, temperature, iam, inverter};

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

// Irradiance decomposition
let (dni, dhi) = irradiance::erbs(cs.ghi, solpos.zenith, 172, 1366.1);

// Angle of incidence on tilted surface
let aoi = irradiance::aoi(30.0, 180.0, solpos.zenith, solpos.azimuth);

// Incidence angle modifier
let iam_loss = iam::physical(aoi, 1.526, 4.0, 0.002);

// Cell temperature
let (t_cell, _) = temperature::sapm_cell_temperature(
    800.0, 30.0, 3.0, -3.56, -0.075, 3.0, 1000.0);
```

## Modules

### Core Physics

| Module | Description | Key Functions |
|--------|-------------|---------------|
| `solarposition` | Solar position (SPA + analytical) | `get_solarposition`, `declination_spencer71`, `hour_angle`, `solar_zenith_analytical` |
| `atmosphere` | Atmospheric properties | `get_relative_airmass`, `alt2pres`, `gueymard94_pw`, `kasten96_lt`, `windspeed_powerlaw` |
| `clearsky` | Clear sky irradiance | `ineichen`, `haurwitz`, `bird`, `simplified_solis` |
| `irradiance` | Irradiance transposition & decomposition | `aoi`, `perez`, `haydavies`, `erbs`, `disc`, `dirindex`, `get_total_irradiance` |
| `location` | Location with convenience methods | `get_solarposition`, `get_clearsky`, `get_airmass` |

### System Modeling

| Module | Description | Key Functions |
|--------|-------------|---------------|
| `temperature` | Cell/module temperature | `sapm_cell_temperature`, `pvsyst_cell_temperature`, `faiman`, `fuentes`, `noct_sam`, `ross` |
| `iam` | Incidence angle modifiers | `ashrae`, `physical`, `martin_ruiz`, `schlick`, `sapm` |
| `pvsystem` | PV system & single-diode params | `calcparams_desoto`, `calcparams_cec`, `calcparams_pvsyst`, `sapm` |
| `singlediode` | Single-diode equation solvers | `bishop88`, `bishop88_mpp`, `bishop88_i_from_v`, `i_from_v`, `v_from_i` |
| `inverter` | Inverter models | `pvwatts_ac`, `sandia`, `adr`, `pvwatts_multi`, `sandia_multi` |
| `tracking` | Single-axis tracker | `singleaxis`, `calc_surface_orientation` |
| `modelchain` | End-to-end simulation pipeline | `run_model_from_weather`, `run_model_from_poa`, `with_pvwatts`, `with_sapm` |

### Advanced Features

| Module | Description | Key Functions |
|--------|-------------|---------------|
| `bifacial` | Bifacial irradiance | `get_irradiance_infinite_sheds` |
| `shading` | Row-to-row shading | `masking_angle`, `ground_angle`, `projected_solar_zenith_angle`, `shaded_fraction1d` |
| `soiling` | Soiling losses | `hsu`, `kimber`, `accumulation_model` |
| `snow` | Snow coverage & losses | `fully_covered_nrel`, `coverage_nrel`, `dc_loss_nrel` |
| `spectrum` | Spectral mismatch | `spectral_mismatch_modifier`, `spectral_factor_sapm`, `spectral_factor_caballero` |
| `scaling` | Geographic smoothing | `wvm_smoothing` |
| `albedo` | Surface albedo | `inland_water_dvoracek`, `surface_albedo` |
| `pvarray` | Module efficiency models | `pvefficiency_adr`, `huld` |
| `transformer` | Transformer losses | `simple_efficiency` |

### I/O & Data

| Module | Description | Key Functions |
|--------|-------------|---------------|
| `iotools` | Weather data I/O | `read_tmy3`, `read_epw`, `get_pvgis_tmy`, `get_pvgis_hourly`, `retrieve_sam` |
| `ivtools` | IV curve fitting | `fit_sandia_simple`, `fit_desoto`, `rectify_iv_curve` |

## ModelChain Configuration

The `ModelChain` supports multiple model selections:

| Step | Options |
|------|---------|
| **Transposition** | Isotropic, Hay-Davies, Perez, Klucher, Reindl |
| **Temperature** | SAPM, PVsyst, Faiman, Fuentes, NOCT_SAM, PVWatts |
| **IAM (AOI)** | Physical (Fresnel), ASHRAE, Martin-Ruiz, SAPM, No loss |
| **Spectral** | No loss |
| **DC Model** | PVWatts |
| **AC Model** | PVWatts, Sandia, ADR |
| **Losses** | PVWatts, No loss |

Factory constructors:
- `ModelChain::with_pvwatts()` — PVWatts DC/AC, Physical AOI, Perez transposition
- `ModelChain::with_sapm()` — SAPM temperature, ASHRAE AOI, Hay-Davies transposition
- `ModelChain::with_config()` — fully custom model selection

## Weather Data

### Read TMY3 files

```rust
let data = pvlib::iotools::read_tmy3("weather.csv").unwrap();
println!("Location: {} ({}, {})", data.metadata.name.unwrap(),
    data.metadata.latitude, data.metadata.longitude);
for record in &data.records {
    println!("GHI: {} W/m2, Temp: {} °C", record.ghi, record.temp_air);
}
```

### Read EPW files

```rust
let data = pvlib::iotools::read_epw("weather.epw").unwrap();
```

### Fetch from PVGIS API

```rust
let tmy = pvlib::iotools::get_pvgis_tmy(45.0, 8.0, None, None, None).unwrap();
```

## References

This library implements algorithms from peer-reviewed solar energy literature:

- **Solar Position**: NREL Solar Position Algorithm (Reda & Andreas, 2004), Spencer (1971)
- **Clear Sky**: Ineichen & Perez (2002), Bird (1981), Haurwitz (1945), Simplified Solis (Ineichen, 2008)
- **Irradiance**: Perez et al. (1990), Hay & Davies (1980), Erbs et al. (1982), DISC (Maxwell, 1987)
- **Temperature**: SAPM (King et al., 2004), PVsyst, Faiman (2008), Fuentes (1987), Ross (1980)
- **IAM**: Martin & Ruiz (2001), De Soto et al. (2006), ASHRAE
- **Single Diode**: Bishop (1988), De Soto et al. (2006)
- **Inverter**: PVWatts (Dobos, 2014), Sandia (King et al., 2007), ADR (Driesse, 2023)
- **Tracking**: Marion & Dobos (2013)
- **Soiling**: HSU (Coello & Boyle, 2019), Kimber et al. (2006)
- **Bifacial**: Infinite sheds model (Mikofski et al., 2019)

## Comparison with pvlib-python

pvlib-rust covers the core simulation pipeline of pvlib-python. Some differences:

- **Data structures**: Uses Rust structs and enums instead of pandas DataFrames
- **Single timestep**: Functions operate on single values (no vectorized operations over time series yet)
- **IO tools**: Covers TMY3, EPW, and PVGIS; pvlib-python has 20+ additional data source adapters
- **No spectral irradiance model**: SPECTRL2 is not yet ported (requires extensive lookup tables)

## Building & Testing

```bash
cargo build          # Build the library
cargo test           # Run all 319 tests
cargo clippy         # Lint checks
cargo doc --open     # Generate and view documentation
```

## License

Licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

This project is a Rust port of [pvlib-python](https://github.com/pvlib/pvlib-python), developed by the pvlib community. All credit for the underlying algorithms and models goes to the original authors and the pvlib-python contributors.
