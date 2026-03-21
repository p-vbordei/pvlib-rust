# PvLib-Rust Architecture & Development Guide

## Mission

A faithful Rust port of [pvlib-python](https://github.com/pvlib/pvlib-python) — the industry-standard open-source library for simulating photovoltaic energy systems. This library aims to provide the same algorithms, accuracy, and API ergonomics as pvlib-python, with Rust's performance and safety guarantees.

## Module Map (pvlib-python → pvlib-rust)

Each module below maps directly to a pvlib-python module. Implement them in priority order.

### Phase 1: Core Physics (Must Have)

| pvlib-python module | Rust module | Status | Description |
|---|---|---|---|
| `pvlib.location` | `location.rs` | ✅ Done | Location struct with lat/lon/tz/altitude |
| `pvlib.solarposition` | `solarposition.rs` | ✅ Started | Solar position via NREL SPA |
| `pvlib.atmosphere` | `atmosphere.rs` | ❌ TODO | Airmass, pressure, refraction |
| `pvlib.clearsky` | `clearsky.rs` | ❌ TODO | Clear-sky irradiance models |
| `pvlib.irradiance` | `irradiance.rs` | ❌ TODO | Decomposition & transposition |
| `pvlib.temperature` | `temperature.rs` | ❌ TODO | Cell/module temperature models |
| `pvlib.pvsystem` | `pvsystem.rs` | ❌ TODO | PV system definitions |
| `pvlib.inverter` | `inverter.rs` | ❌ TODO | Inverter performance models |
| `pvlib.iam` | `iam.rs` | ❌ TODO | Incidence angle modifiers |

### Phase 2: System Simulation

| pvlib-python module | Rust module | Status | Description |
|---|---|---|---|
| `pvlib.singlediode` | `singlediode.rs` | ❌ TODO | Single-diode equation solver |
| `pvlib.modelchain` | `modelchain.rs` | ❌ TODO | End-to-end simulation pipeline |
| `pvlib.tracking` | `tracking.rs` | ❌ TODO | Single/dual-axis tracker geometry |
| `pvlib.shading` | `shading.rs` | ❌ TODO | Shading loss calculations |
| `pvlib.soiling` | `soiling.rs` | ❌ TODO | Soiling loss models |
| `pvlib.snow` | `snow.rs` | ❌ TODO | Snow cover models |

### Phase 3: Advanced

| pvlib-python module | Rust module | Status | Description |
|---|---|---|---|
| `pvlib.bifacial` | `bifacial/` | ❌ TODO | Bifacial irradiance models |
| `pvlib.spectrum` | `spectrum/` | ❌ TODO | Spectral irradiance & mismatch |
| `pvlib.scaling` | `scaling.rs` | ❌ TODO | Wavelet variability scaling |
| `pvlib.iotools` | `iotools/` | ❌ TODO | Weather data file readers |

---

## Implementation Guidelines

### 1. Match pvlib-python function signatures

Every public function should mirror the pvlib-python API as closely as Rust allows:

```python
# pvlib-python
pvlib.irradiance.perez(surface_tilt, surface_azimuth, dhi, dni, dni_extra,
                        solar_zenith, solar_azimuth, airmass)
```

```rust
// pvlib-rust
pub fn perez(surface_tilt: f64, surface_azimuth: f64, dhi: f64, dni: f64,
             dni_extra: f64, solar_zenith: f64, solar_azimuth: f64,
             airmass: f64) -> PerezResult
```

### 2. Reference papers in docstrings

Every model function MUST cite the original paper. Example:

```rust
/// Estimate DNI and DHI from GHI using the Erbs (1982) model.
///
/// # References
/// Erbs, D.G., S.A. Klein, and J.A. Duffie, 1982,
/// "Estimation of the diffuse radiation fraction for hourly,
/// daily and monthly-average global radiation,"
/// Solar Energy, 28(4), pp. 293-302.
pub fn erbs(ghi: f64, zenith: f64, day_of_year: u32) -> ErbsResult { ... }
```

### 3. Validate against pvlib-python outputs

For every function, write tests that compare output against pvlib-python reference values. Generate reference values by running the Python version:

```python
import pvlib
result = pvlib.irradiance.perez(40, 180, 100, 600, 1360, 30, 180, 1.5)
print(repr(result))  # Use this as the Rust test expected value
```

### 4. Use `f64` throughout

All angles in degrees. All irradiance in W/m². All temperatures in °C. All power in W or kW (document which). All pressures in Pa.

### 5. Error handling

Use `thiserror` for library errors. Return `Result` for functions that can fail (e.g., invalid inputs). Use assertions or clamps for physically impossible values (negative irradiance, zenith > 180°).

---

## Key Algorithms to Implement (with references)

### Atmosphere (`atmosphere.rs`)

1. **Relative airmass** — Kasten & Young (1989)
   - `AM_rel = 1 / (cos(zenith_rad) + 0.50572 * (96.07995 - zenith)^(-1.6364))`
   - Clamp zenith < 90°
   - Source: pvlib/atmosphere.py `get_relative_airmass()`

2. **Absolute airmass** — altitude correction
   - `AM_abs = AM_rel * P / 101325`
   - Where P = pressure at altitude
   - Source: pvlib/atmosphere.py `get_absolute_airmass()`

3. **Pressure from altitude** — barometric formula
   - `P = 101325 * (1 - 2.25577e-5 * altitude)^5.25588`
   - Source: pvlib/atmosphere.py `alt2pres()`

### Clear Sky (`clearsky.rs`)

1. **Ineichen & Perez (2002)**
   - The standard clear-sky model. Uses Linke turbidity factor.
   - Source: pvlib/clearsky.py `ineichen()`
   - Paper: Ineichen, P. and Perez, R., 2002. "A new airmass independent formulation for the Linke turbidity coefficient"

2. **Haurwitz (1945)** — simple GHI-only model
   - `GHI = 1098 * cos(zenith) * exp(-0.057 / cos(zenith))`
   - Good for quick estimates

### Irradiance Decomposition (`irradiance.rs`)

1. **Erbs (1982)** — GHI → DNI + DHI
   - Clearness index: `kt = GHI / (ETI * cos(zenith))`
   - Piecewise polynomial for diffuse fraction
   - Source: pvlib/irradiance.py `erbs()`

2. **DIRINT (Perez 1992)** — improved decomposition
   - Uses persistence index and atmospheric moisture
   - More accurate than Erbs in tropical climates
   - Source: pvlib/irradiance.py `dirint()`

3. **Boland (2008)** — logistic regression model
   - Single-equation, continuous function
   - Source: pvlib/irradiance.py `boland()`

### POA Transposition (`irradiance.rs`)

**This is the highest-impact algorithm to get right.**

1. **Isotropic** — simplest, baseline
   - `POA_diffuse = DHI * (1 + cos(tilt)) / 2`
   - `POA_ground = GHI * albedo * (1 - cos(tilt)) / 2`
   - Source: pvlib/irradiance.py `get_total_irradiance(model='isotropic')`

2. **Perez (1990)** — industry standard, highest accuracy
   - Uses brightness coefficients (F1, F2) from 8 sky clearness bins
   - Accounts for circumsolar, horizon brightening, and isotropic diffuse
   - Source: pvlib/irradiance.py `perez()`
   - Paper: Perez et al., 1990, "Modeling daylight availability and irradiance components from direct and global irradiance"
   - **Critical: The F1/F2 coefficient table has 8 rows × 6 columns. Copy exactly from pvlib source.**

3. **Hay-Davies (1980)** — good middle ground
   - `POA_diffuse = DHI * ((1 - A) * (1 + cos(tilt)) / 2 + A * R_b)`
   - Where A = DNI/ETI (anisotropy index), R_b = cos(AOI)/cos(zenith)
   - Source: pvlib/irradiance.py `haydavies()`

4. **Klucher (1979)** — better for overcast
   - Modulation factor F = 1 - (DHI/GHI)²
   - Source: pvlib/irradiance.py `klucher()`

5. **Angle of Incidence (AOI)**
   - `cos(AOI) = cos(zenith)*cos(tilt) + sin(zenith)*sin(tilt)*cos(azimuth - surface_azimuth)`
   - Source: pvlib/irradiance.py `aoi()`

6. **Extraterrestrial irradiance**
   - Spencer (1971): `ETI = 1361 * (1 + 0.033412 * cos(2π * (day - 3) / 365))`
   - Source: pvlib/irradiance.py `get_extra_radiation()`

### Incidence Angle Modifier (`iam.rs`)

1. **Physical (Fresnel)** — glass/air interface reflection
   - Uses Snell's law and Fresnel equations
   - `n_glass` typically 1.526, `extinction_coeff * thickness` typically 0.0002
   - Source: pvlib/iam.py `physical()`

2. **ASHRAE** — simple parametric
   - `IAM = 1 - b * (1/cos(AOI) - 1)`
   - Default b = 0.05
   - Source: pvlib/iam.py `ashrae()`

3. **Martin-Ruiz (2001)** — empirical
   - Source: pvlib/iam.py `martin_ruiz()`

### Cell Temperature (`temperature.rs`)

1. **SAPM (King 2004)** — Sandia Array Performance Model
   - `T_cell = POA * exp(a + b * wind_speed) + T_ambient`
   - Pre-defined coefficients for ~20 mounting configurations
   - Source: pvlib/temperature.py `sapm_cell()`

2. **Faiman (2008)** — simplified
   - `T_cell = T_ambient + POA / (U0 + U1 * wind_speed)`
   - Source: pvlib/temperature.py `faiman()`

3. **PVSyst** — widely used in commercial tools
   - Similar to SAPM with different parameterization
   - Source: pvlib/temperature.py `pvsyst_cell()`

### PV System Power (`pvsystem.rs`)

1. **PVWatts DC** — linear model
   - `P_dc = capacity * (POA / 1000) * (1 + gamma * (T_cell - 25))`
   - Source: pvlib/pvsystem.py `pvwatts_dc()`

2. **PVWatts AC** — simple inverter
   - `P_ac = min(P_dc * eta_inv, P_ac_max)`
   - Source: pvlib/pvsystem.py `pvwatts_ac()`

3. **Sandia Inverter** — detailed inverter model
   - Uses 5 coefficients from CEC database
   - Source: pvlib/inverter.py `sandia()`

### Single Diode (`singlediode.rs`)

- Solves: `I = I_L - I_0 * (exp((V + I*R_s) / (n*N_s*V_th)) - 1) - (V + I*R_s) / R_sh`
- Requires Lambert W function (implement or use a crate)
- Source: pvlib/singlediode.py
- Papers: De Soto et al. (2006), Dobos (2012)

---

## Testing Strategy

### Unit Tests
Every function gets tests comparing against pvlib-python reference outputs.

### Integration Tests
`tests/test_modelchain.rs` — run a full simulation for a known location/date and compare total energy output against pvlib-python.

### Reference Value Generation
Create a Python script `tests/generate_references.py` that uses pvlib-python to generate JSON test fixtures:

```python
import pvlib, json
loc = pvlib.location.Location(44.4268, 26.1025, 'Europe/Bucharest', 85)
times = pd.date_range('2024-06-21', periods=24, freq='h', tz='Europe/Bucharest')
solpos = loc.get_solarposition(times)
# ... serialize to JSON for Rust tests
```

---

## Perez Coefficient Table

This is the most important constant table in the library. Copy exactly:

```
// Perez et al. (1990) Table 1 — brightness coefficients
// Columns: f11, f12, f13, f21, f22, f23
// Rows: 8 sky clearness bins (epsilon ranges)
const PEREZ_COEFFICIENTS: [[f64; 6]; 8] = [
    // epsilon: [1.000, 1.065)
    [-0.0083117, 0.5877277, -0.0620636, -0.0596012, 0.0721249, -0.0220216],
    // epsilon: [1.065, 1.230)
    [0.1299457, 0.6825954, -0.1513752, -0.0189325, 0.0659650, -0.0288748],
    // epsilon: [1.230, 1.500)
    [0.3296958, 0.4868735, -0.2210958, 0.0554140, -0.0639588, -0.0260542],
    // epsilon: [1.500, 1.950)
    [0.5682053, 0.1874990, -0.2951290, 0.1088631, -0.1519229, -0.0139754],
    // epsilon: [1.950, 2.800)
    [0.8730280, -0.3920403, -0.3616149, 0.2255647, -0.4620442,  0.0012448],
    // epsilon: [2.800, 4.500)
    [1.1326077, -1.2367284, -0.4118494, 0.2877813, -0.8230357,  0.0558225],
    // epsilon: [4.500, 6.200)
    [1.0601591, -1.5999137, -0.3589221, 0.2642124, -1.1272340,  0.1310694],
    // epsilon: [6.200, ∞)
    [0.6777470, -0.3272588, -0.2504286, 0.1561313, -1.3765031,  0.2506212],
];
```

---

## Common Pitfalls

1. **Degrees vs Radians** — pvlib-python works in degrees. So should we. Convert internally only where needed.
2. **Zenith > 90°** — Sun is below horizon. Most irradiance functions should return 0, not NaN.
3. **Division by cos(zenith)** — Always clamp `cos(zenith) >= cos(87.534°) ≈ 0.0436` to avoid infinity.
4. **Negative irradiance** — Clamp to 0. Physical impossibility.
5. **Temperature coefficient sign** — gamma is negative for silicon (~-0.004 /°C). Power decreases with temperature.
6. **Albedo range** — Should be [0, 1]. Default 0.25.
7. **Airmass at night** — Return NaN or None, not infinity.
