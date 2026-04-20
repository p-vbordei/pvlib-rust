# pvlib-rust roadmap

Items deliberately deferred past v0.1.6. They are tracked here so the
current `src/` reflects only what has been implemented and validated, and
so downstream users can see what is genuinely on the way versus marketing.

## Numerical parity (high priority)

1. ~~True DIRINT~~ — **shipped in 0.1.6.** Full Perez (1992) 4-D
   coefficient tensor in `src/dirint_coeffs.rs`, scalar and series APIs
   in `irradiance::{dirint, dirint_series}`.
2. ~~True `detect_clearsky` (Reno–Hansen 2016)~~ — **shipped in 0.1.6.**
   Windowed 5-criterion algorithm with iterative α rescaling in
   `clearsky::{detect_clearsky, detect_clearsky_detail,
   ClearSkyThresholds}`.
3. **Full `infinite_sheds.get_irradiance_poa`** — the 2-D view-factor
   primitives (`vf_row_sky_2d*`, `vf_row_ground_2d*`,
   `unshaded_ground_fraction`, `solar_projection_tangent`) are in
   `src/bifacial.rs` as of 0.1.6. What's still missing vs.
   `pvlib.bifacial.infinite_sheds` is the front-face shadow-fraction
   pipeline (`_shaded_fraction`, `_poa_sky_diffuse_pv` with split
   shaded / unshaded integrals, multi-row ground view factor
   `vf_ground_sky_2d_integ`) that composes front POA from these primitives.
4. **Perez 1990 coefficient table** — replace the 7-digit coefficients
   currently in `irradiance::perez` with the 4-digit published values
   used by `pvlib-python` (the extra digits appear fabricated and make
   bit-identical parity impossible).
5. **pvlib-python reference fixtures** — commit a pinned CSV produced
   by a specific `pvlib-python` version (e.g. 0.11.2) for a canonical
   location + weather input, and add `tests/test_python_parity.rs` that
   diffs `ModelChain` output against it with documented tolerances.
6. **Solar position refraction correction** — wire the `spa` crate's
   `pressure` and `temperature` arguments through
   `solarposition::get_solarposition`. Currently ignored, which drifts
   the apparent zenith by up to ~0.1° at low sun.

## API / ergonomics (medium priority)

7. **Unified error type** — replace the `Box<dyn Error>` returns in
   `iotools` and the `Result<_, String>` returns in `batch` with a
   `thiserror`-based `PvlibError` enum.
8. **`#[non_exhaustive]` on `WeatherSeries` / `SimulationSeries` /
   `BatchModelChain`** — once a non-breaking constructor exists, these
   should become non-exhaustive so additional outputs or config
   options do not constitute breaking changes. Adding today would
   break every test and every user that constructs `WeatherSeries`
   via struct-literal syntax.
9. **`Box<dyn Mount>` → `enum Mount`** — the `Mount` trait object in
   `pvsystem::Array` is a public API hazard and blocks `Clone`/`Copy`.
   A closed `enum Mount { Fixed(..), SingleAxis(..) }` covers the only
   two implementations in the crate today.
10. **`BatchModelChain::run` input validation** — return a structured
    `Err` (not `assert_eq!` / `panic`) on length-mismatched weather
    vectors. Add tilt / azimuth / capacity range checks.

## Performance (medium priority)

11. **SoA output writes in `BatchModelChain::run`** — the current
    implementation collects an 11-tuple-per-row `Vec` then transposes
    with 12 sequential `iter().map().collect()` passes. Pre-sizing the
    output vectors and writing by index eliminates the intermediate
    allocation and one full pass.
12. **Cached SPA per day / per hour** — the `spa` crate recomputes
    every per-day term (nutation, obliquity, apparent sidereal time)
    on every call. Caching these across a contiguous day in the batch
    loop is the single biggest speedup available.
13. **De-duplicate `aoi()` / hoist `to_radians()`** in Perez and
    HayDavies transposition and in the batch main loop.
14. **`with_min_len` on rayon iterators** in cheap kernels
    (`airmass_relative_batch`, `aoi_batch`, `iam_ashrae_batch`) — for
    8760 items and sub-microsecond work the fork/join overhead exceeds
    the gain.

## Ecosystem (strategic)

15. **PyO3 / maturin wheel** — expose `BatchModelChain` and the
    `*_batch` functions as a `pvlib-accel` wheel on PyPI. Accept
    `numpy.ndarray` zero-copy, release the GIL around the rayon loop.
    This is the single largest adoption lever — the pvlib user base
    lives in Python.
16. **Polars / Arrow feature** — optional `polars = ["dep:polars"]`
    and `arrow = ["dep:arrow"]` features with
    `impl From<&DataFrame> for WeatherSeries` and
    `impl Into<DataFrame> for SimulationSeries`.

## Testing (medium priority)

17. **Property-based tests** — `proptest` invariants for
    `iam::physical` monotonicity, `irradiance::erbs` kt ∈ [0,1], solar
    position round-trip, inverter saturation.
18. **Edge-case suite** — polar night / polar day / equator equinox,
    DST transitions, leap-year Feb 29, year-boundary timestamps.
19. **Fuzzing** — `cargo fuzz` targets for `iotools::read_tmy3` and
    `iotools::read_epw` (classic crash-surface).
20. **Coverage gate** — `cargo llvm-cov` in CI with a floor (e.g.
    80 %); `src/ivtools.rs` (720 LOC) currently has no dedicated test
    file.
