//! # pvlib-rust
//!
//! A Rust port of [pvlib-python](https://github.com/pvlib/pvlib-python), the
//! open-source solar photovoltaic modeling library. The crate covers the full
//! PV simulation pipeline: solar position, clear-sky irradiance, transposition,
//! incidence-angle modifier, cell temperature, and DC/AC power.
//!
//! ## Feature flags
//!
//! - `pvgis` (default) — enables network I/O for PVGIS and NREL SAM helpers
//!   in [`iotools`]. Pulls in `reqwest` with `rustls-tls`. Disable with
//!   `default-features = false` to build the pure-compute parts of the
//!   library on targets where blocking HTTP is not available
//!   (e.g. `wasm32-unknown-unknown`).
//!
//! ## Quick start — batch simulation
//!
//! ```no_run
//! use chrono::TimeZone;
//! use chrono_tz::US::Mountain;
//! use pvlib::batch::{BatchModelChain, WeatherSeries};
//! use pvlib::irradiance::DiffuseModel;
//! use pvlib::Location;
//!
//! let location = Location::try_new(39.74, -105.18, Mountain, 1830.0, "Golden, CO").unwrap();
//!
//! let mc = BatchModelChain::pvwatts(location, 30.0, 180.0, 5000.0)
//!     .with_gamma_pdc(-0.004)
//!     .with_inverter(5000.0, 0.96)
//!     .with_albedo(0.2)
//!     .with_transposition(DiffuseModel::Perez)
//!     .with_auto_decomposition(true)
//!     .with_system_losses(0.14);
//!
//! # let times = Vec::<chrono::DateTime<chrono_tz::Tz>>::new();
//! # let zeros = vec![0.0_f64; times.len()];
//! let weather = WeatherSeries {
//!     times,
//!     ghi: zeros.clone(), dni: zeros.clone(), dhi: zeros.clone(),
//!     temp_air: zeros.clone(), wind_speed: zeros.clone(),
//!     albedo: None,
//! };
//!
//! let results = mc.run(&weather).unwrap();
//! println!("Annual energy: {:.0} kWh", results.total_energy_wh() / 1000.0);
//! ```
//!
//! See the [README](https://github.com/p-vbordei/pvlib-rust) for the full
//! Quick Start, step-by-step usage, and the comparison with pvlib-python.
//!
//! ## Algorithm notes
//!
//! - [`irradiance::dirint`] and [`irradiance::dirint_series`] are faithful
//!   ports of `pvlib.irradiance.dirint`, including the full 6×6×7×5
//!   coefficient tensor. The scalar version uses the `delta_kt' = -1`
//!   fallback bin; call the series version for time-aware persistence.
//! - [`clearsky::detect_clearsky`] is a faithful port of the Reno–Hansen
//!   (2016) windowed 5-criterion algorithm with iterative α rescaling.
//! - [`bifacial`] exposes the real 2-D view-factor primitives
//!   (`vf_row_sky_2d`, `vf_row_sky_2d_integ`, `vf_row_ground_2d`,
//!   `vf_row_ground_2d_integ`, `unshaded_ground_fraction`,
//!   `solar_projection_tangent`) used inside
//!   `pvlib.bifacial.infinite_sheds`, plus a composed rear-irradiance
//!   helper [`bifacial::rear_irradiance_sheds`]. The front-face /
//!   shadow-fraction pipeline of the full `get_irradiance_poa` is
//!   tracked in `ROADMAP.md`.

// Private data tables for algorithms that require large coefficient sets.
mod dirint_coeffs;

pub mod location;
pub mod solarposition;
pub mod atmosphere;
pub mod clearsky;
pub mod irradiance;
pub mod temperature;
pub mod iam;
pub mod tracking;
pub mod pvsystem;
pub mod inverter;
pub mod modelchain;
pub mod bifacial;
pub mod singlediode;
pub mod iotools;
pub mod shading;
pub mod soiling;
pub mod snow;
pub mod spectrum;
pub mod scaling;
pub mod albedo;
pub mod pvarray;
pub mod transformer;
pub mod ivtools;
pub mod batch;

// ---------------------------------------------------------------------------
// Crate-root re-exports for the common entry points. Users can write
// `pvlib::Location` instead of `pvlib::location::Location`, etc.
// ---------------------------------------------------------------------------

pub use location::{Location, LocationError};
pub use batch::{BatchModelChain, WeatherSeries, SimulationSeries};
pub use modelchain::{ModelChain, ModelChainConfig, ModelChainResult, WeatherInput};
pub use pvsystem::{PVSystem, Array, Mount, FixedMount, SingleAxisTrackerMount};
