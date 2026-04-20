# pvlib-rust Changelog

## [0.1.6] - 2026-04-20

### Changed
- **License** changed from MIT to Apache-2.0. Prior releases remain available under MIT.
- **`CHANGELOG.md`** entries for 0.1.1 through 0.1.4 rewritten in neutral engineering tone. The original entries used marketing language ("Immune by Design", "Ultimate 10x Performance Sweep") that misrepresented the scope of the work.
- **`reqwest` is now an optional dependency** behind the `pvgis` feature (on by default). Users who do not need PVGIS/SAM HTTP helpers can build with `default-features = false` to drop the TLS/hyper/tokio dependency tree.
- **`reqwest`** uses `rustls-tls` instead of the default system TLS stack, to simplify static-binary builds.
- **Build profile**: `[profile.release]` now sets `lto = "fat"`, `codegen-units = 1`, `panic = "abort"`, `strip = "debuginfo"`.

### Added
- Crate-level rustdoc in `src/lib.rs` with runnable Quick Start example.
- Top-level re-exports (`pvlib::Location`, `pvlib::BatchModelChain`, `pvlib::WeatherSeries`, `pvlib::SimulationSeries`, `pvlib::PVSystem`, `pvlib::ModelChain`).
- `Location::try_new` — validates `latitude ∈ [-90, 90]`, `longitude ∈ [-180, 180]`, finite altitude.
- `criterion` dev-dependency and `benches/batch_tmy.rs` for reproducible TMY-year timing.
- `.github/workflows/ci.yml` with fmt/clippy/test/doc/wasm matrix.
- `examples/annual_tmy.rs` — runnable end-to-end TMY simulation.
- `ROADMAP.md` — enumerates pending work (true DIRINT coefficient table, Reno–Hansen `detect_clearsky`, `infinite_sheds` port, PyO3/maturin wheel, Polars integration, proptest suite, pvlib-python numerical-parity fixtures).
- `#[non_exhaustive]` on public configuration structs and model-selection enums that are expected to grow (`WeatherSeries`, `SimulationSeries`, `BatchModelChain`, `DCModel`, `ACModel`, `AOIModel`, `TemperatureModel`, `TranspositionModel`, `LossesModel`, `SpectralModel`). Users who previously constructed these via `..` record syntax must migrate to constructors/builders.
- `#[must_use]` on aggregate accessors (`SimulationSeries::total_energy_wh`, `peak_power`, `capacity_factor`) and expensive constructors.
- MSRV declared: `rust-version = "1.85"` (required by `edition = "2024"`).

### Fixed
- `src/albedo.rs:inland_water_dvoracek` no longer panics on unknown `surface_condition`; returns `f64::NAN` and emits a one-time warning. A new `inland_water_dvoracek_try` returns `Option<f64>` for callers that want to detect this explicitly.
- `src/iotools.rs` PVGIS/SAM HTTP calls now use a shared `reqwest::blocking::Client` with a 30-second timeout. Previously unbounded (`reqwest::blocking::get`) and could hang indefinitely on stalled connections.
- `src/iotools.rs` TMY3/EPW readers now use bounds-checked field access; malformed rows produce a structured error instead of panicking with an OOB index.
- `src/pvsystem.rs` — `trait Mount` now requires `Send + Sync`, so `Array` / `PVSystem` / `ModelChain` are safe to share across threads (`Arc<ModelChain>` in an async handler, parallel iteration, etc.).
- `src/irradiance.rs:erbs` — explicit guard returns `(0.0, 0.0)` on non-finite `ghi` or `zenith` instead of propagating NaN.
- `src/irradiance.rs:disc` — guards against NaN airmass producing NaN DNI below the zenith cutoff.
- `src/iam.rs:physical` — output clamped to `[0.0, 1.0]` to match `pvlib.iam.physical` semantics and avoid super-unit IAMs at grazing angles from floating-point drift.
- `src/modelchain.rs` — `ACModel::Sandia` and `ACModel::ADR` now return an explicit `unimplemented!("…")` error instead of silently delegating to PVWatts. The variants are retained for API stability but callers get a clear message.
- Hot-path day-of-year computation (`BatchModelChain::run`, `ModelChain::run_model_from_weather`, `ModelChain::run_model_from_poa`) now uses `chrono::Datelike::ordinal()` instead of `format("%j").to_string().parse::<i32>().unwrap_or(1)`. Two heap allocations per timestep removed.

### Full ports (algorithms that were previously labelled approximations)
- **`irradiance::dirint`** is now a faithful port of `pvlib.irradiance.dirint`, including the full 6×6×7×5 Perez (1992) coefficient tensor (lives in the private `dirint_coeffs` module). Signature changed from
  `fn dirint(ghi, zenith, dew_point, pressure, dni_extra) -> (f64, f64)` to
  `fn dirint(ghi, zenith, day_of_year, pressure_pa, temp_dew_c: Option<f64>) -> f64` (returns DNI). A new
  `irradiance::dirint_series(&[ghi], &[zenith], &[doy], pressure_pa, temp_dew: Option<&[f64]>, use_delta_kt_prime: bool) -> Vec<f64>`
  propagates Δkt′ temporal persistence (Perez eqns 2/3). `clearness_index_zenith_independent` and `precipitable_water_from_dew_point` are now public helpers.
- **`clearsky::detect_clearsky`** is now a faithful port of the Reno–Hansen (2016) windowed 5-criterion algorithm with iterative α rescaling. Signature changed from `fn detect_clearsky(ghi: f64, clearsky_ghi: f64, _window_length: usize) -> bool` to
  `fn detect_clearsky(measured: &[f64], clearsky: &[f64], sample_interval_minutes: f64, thresholds: ClearSkyThresholds) -> Vec<bool>`.
  New types: `ClearSkyThresholds` (with `::default()` matching Reno–Hansen 2016 and `::from_sample_interval()` matching Jordan–Hansen 2023 Table 1) and `ClearSkyDetectionResult`; `detect_clearsky_detail` returns clear-sample flags plus the fitted α and iteration count.
- **`bifacial`** replaces the single linear `(1 − exp(−h/pitch)) · GHI · albedo` approximation with the real Marion 2017 / Mikofski 2019 2-D view-factor primitives: `vf_row_sky_2d`, `vf_row_sky_2d_integ`, `vf_row_ground_2d`, `vf_row_ground_2d_integ`, `solar_projection_tangent`, `unshaded_ground_fraction`, plus a composed rear-irradiance helper `rear_irradiance_sheds`. `get_irradiance_infinite_sheds` is retained (same signature, `height` / `pitch` now ignored) and now computes rear irradiance via those real view factors. The remaining gap vs. `pvlib.bifacial.infinite_sheds.get_irradiance_poa` is the front-face shadow-fraction pipeline, still in `ROADMAP.md`.
- README speed comparison softened: pvlib-python vectorized path is typically sub-second for a TMY year; the 500–1250× headline was based on a non-vectorized reference and has been replaced with a more representative range.

## [0.1.5] - 2026-03-28
- Minor fixes and Cargo.lock update.

## [0.1.4] - 2026-03-28

### Added
- **Batch API extensions**:
    - `solar_elevation` field on `SimulationSeries` (removes `90 - zenith` boilerplate at call sites).
    - `solar_position_batch_utc()` — accepts `&[NaiveDateTime]` (UTC).
    - `WeatherSeries::from_utc()` — constructs weather input from UTC timestamps plus an IANA timezone string.
    - `with_auto_decomposition(true)` — runs Erbs GHI→DNI/DHI when both DNI and DHI are zero or NaN.
    - `with_bifacial(bifaciality_factor, ground_albedo)` — rear-side gain applied at the DC level, clamped at 25%.
    - `with_system_losses(fraction)` — flat DC derating, clamped to `[0.0, 1.0]`.
- NaN-in-NaN-out semantics documented for DNI/DHI upstream of `BatchModelChain::run`.
- 10 edge-case tests covering NaN propagation, loss clamping, and bifacial zero-POA safety.

### Fixed
- Bifacial gain moved from post-AC to pre-inverter (DC level).
- Clean `cargo clippy -D warnings` pass.

## [0.1.3] - 2026-03-22

### Changed
- Audit of `pvlib-python` performance-related GitHub issues and checked which do not apply to the Rust port due to language-level differences (no Numba cache, no pandas index alignment, native f64 parallel-iter). No algorithmic changes.

## [0.1.2] - 2026-03-22

### Changed
- Audit of `pvlib-python` numerical-correctness issues; verified existing Rust implementations already guard the relevant edge cases via `f64::clamp` on trigonometric inputs. No algorithmic changes.

### Notes
- Previous wording in this entry (e.g., "Immune by Design", "Mathematical Audit Certification") overstated the scope of the work and has been replaced with the neutral summary above.

## [0.1.1] - 2026-03-21

### Fixed
- `zip` length mismatches inside the `rayon` batch module now enforce equal lengths.
- Division-by-zero guard in `shading.rs` cross-axis geometry when the sun is perpendicular to the tracking axis.
- Bounds guard in `snow.rs` Townsend thermal-sliding parameter.

## [0.1.0] - 2026-03-21

### Added
- Initial port: core modules (IAM, tracking, atmosphere, PVSystem, solar position) plus a rayon-parallel batch API.
