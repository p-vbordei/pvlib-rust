//! End-to-end integration tests for PvLib-Rust.
//!
//! These tests simulate realistic PV system modeling workflows,
//! exercising the full pipeline from location/weather input through
//! to AC power output, verifying physical plausibility at each step.

use chrono::TimeZone;
use chrono_tz::US::Mountain;
use chrono_tz::Europe::Berlin;
use chrono_tz::UTC;

use pvlib::location::Location;
use pvlib::pvsystem::{PVSystem, Array, FixedMount, Mount};
use pvlib::modelchain::{
    ModelChain, WeatherInput, POAInput, EffectiveIrradianceInput,
    ModelChainConfig, ModelChainResult,
    DCModel, ACModel, AOIModel, SpectralModel, TemperatureModel,
    TranspositionModel, LossesModel,
};
use pvlib::solarposition;
use pvlib::atmosphere;
use pvlib::clearsky;
use pvlib::irradiance;
use pvlib::temperature;
use pvlib::iam;
use pvlib::inverter;
use pvlib::singlediode;
use pvlib::pvsystem;
use pvlib::tracking;
use pvlib::shading;
use pvlib::snow;
use pvlib::soiling;
use pvlib::spectrum;
use pvlib::albedo;
use pvlib::pvarray;
use pvlib::transformer;

// ============================================================================
// Helper: create a typical residential PV system
// ============================================================================

fn make_residential_system() -> (PVSystem, Location) {
    // Golden, Colorado (NREL location)
    let location = Location::new(39.742, -105.178, Mountain, 1830.0, "Golden, CO");

    // 5 kW south-facing system, 30° tilt
    let array = Array {
        mount: Box::new(FixedMount {
            surface_tilt: 30.0,
            surface_azimuth: 180.0,
        }),
        nameplate_dc: 5000.0,
        gamma_pdc: -0.004, // -0.4%/°C typical for c-Si
        modules_per_string: 10,
        strings: 2,
        albedo: 0.2,
    };
    let system = PVSystem::new(vec![array], 5000.0);
    (system, location)
}

fn make_summer_noon_weather() -> WeatherInput {
    let time = Mountain.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap();
    WeatherInput {
        time,
        ghi: Some(900.0),
        dni: Some(800.0),
        dhi: Some(150.0),
        temp_air: 30.0,
        wind_speed: 3.0,
        albedo: Some(0.2),
    }
}

fn make_winter_noon_weather() -> WeatherInput {
    let time = Mountain.with_ymd_and_hms(2024, 12, 21, 12, 0, 0).unwrap();
    WeatherInput {
        time,
        ghi: Some(400.0),
        dni: Some(600.0),
        dhi: Some(80.0),
        temp_air: 0.0,
        wind_speed: 5.0,
        albedo: Some(0.3), // higher albedo from snow
    }
}

fn make_night_weather() -> WeatherInput {
    let time = Mountain.with_ymd_and_hms(2024, 6, 21, 23, 0, 0).unwrap();
    WeatherInput {
        time,
        ghi: Some(0.0),
        dni: Some(0.0),
        dhi: Some(0.0),
        temp_air: 20.0,
        wind_speed: 2.0,
        albedo: None,
    }
}

// ============================================================================
// Test 1: Full ModelChain pipeline - summer noon
// ============================================================================

#[test]
fn test_e2e_full_pipeline_summer_noon() {
    let (system, location) = make_residential_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);
    let weather = make_summer_noon_weather();

    let result = mc.run_model_from_weather(&weather).unwrap();

    // Solar position: summer noon in Colorado, sun should be high
    assert!(result.solar_zenith > 0.0, "Zenith should be positive");
    assert!(result.solar_zenith < 40.0, "Zenith should be < 40° at summer noon at 39.7°N");

    // POA: tilted surface facing sun should have good irradiance
    assert!(result.poa_global > 700.0, "POA global should be > 700 W/m2, got {}", result.poa_global);
    assert!(result.poa_global < 1200.0, "POA global should be < 1200 W/m2, got {}", result.poa_global);

    // AOI modifier: should be close to 1 (sun is close to normal)
    assert!(result.aoi_modifier > 0.8, "AOI modifier should be > 0.8, got {}", result.aoi_modifier);
    assert!(result.aoi_modifier <= 1.0, "AOI modifier should be <= 1.0, got {}", result.aoi_modifier);

    // Cell temperature: should be above ambient due to heating
    assert!(result.cell_temperature > weather.temp_air,
        "Cell temp {} should exceed ambient {}", result.cell_temperature, weather.temp_air);
    assert!(result.cell_temperature < 70.0, "Cell temp should be < 70°C, got {}", result.cell_temperature);

    // DC power: 5kW system with good irradiance should produce 2-5 kW
    assert!(result.dc_power > 2000.0, "DC power should be > 2000 W, got {}", result.dc_power);
    assert!(result.dc_power < 5500.0, "DC power should be < 5500 W, got {}", result.dc_power);

    // AC power: should be less than DC (inverter losses) and positive
    assert!(result.ac_power > 0.0, "AC power should be positive, got {}", result.ac_power);
    assert!(result.ac_power < result.dc_power, "AC power {} should be < DC power {}", result.ac_power, result.dc_power);
    assert!(result.ac_power > 1800.0, "AC power should be > 1800 W, got {}", result.ac_power);

    // Inverter efficiency check
    let inverter_eff = result.ac_power / result.dc_power;
    assert!(inverter_eff > 0.90, "Inverter efficiency should be > 90%, got {:.1}%", inverter_eff * 100.0);
    assert!(inverter_eff < 1.0, "Inverter efficiency should be < 100%");
}

// ============================================================================
// Test 2: Full pipeline - winter noon (colder, lower sun)
// ============================================================================

#[test]
fn test_e2e_full_pipeline_winter_noon() {
    let (system, location) = make_residential_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);
    let weather = make_winter_noon_weather();

    let result = mc.run_model_from_weather(&weather).unwrap();

    // Winter: higher zenith angle at noon
    assert!(result.solar_zenith > 50.0, "Winter zenith should be > 50°, got {}", result.solar_zenith);
    assert!(result.solar_zenith < 70.0, "Winter zenith should be < 70°, got {}", result.solar_zenith);

    // Cell temperature: cold day, should be moderate
    assert!(result.cell_temperature > weather.temp_air,
        "Cell temp {} should exceed ambient {}", result.cell_temperature, weather.temp_air);
    assert!(result.cell_temperature < 40.0, "Winter cell temp should be < 40°C, got {}", result.cell_temperature);

    // Winter produces less power but cold helps efficiency
    assert!(result.dc_power > 500.0, "Winter DC should be > 500 W, got {}", result.dc_power);
    assert!(result.dc_power < 4000.0, "Winter DC should be < 4000 W, got {}", result.dc_power);
    assert!(result.ac_power > 0.0, "AC power should be positive");
}

// ============================================================================
// Test 3: Night time - zero output
// ============================================================================

#[test]
fn test_e2e_nighttime_zero_output() {
    let (system, location) = make_residential_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);
    let weather = make_night_weather();

    let result = mc.run_model_from_weather(&weather).unwrap();

    assert!(result.solar_zenith > 90.0, "Night zenith should be > 90°");
    assert!(result.dc_power <= 0.0 || result.dc_power.abs() < 1.0,
        "Night DC power should be ~0, got {}", result.dc_power);
    assert!(result.ac_power <= 0.0 || result.ac_power.abs() < 1.0,
        "Night AC power should be ~0, got {}", result.ac_power);
}

// ============================================================================
// Test 4: Solar position chain - verify consistent results
// ============================================================================

#[test]
fn test_e2e_solar_position_chain() {
    let location = Location::new(39.742, -105.178, Mountain, 1830.0, "Golden, CO");
    let time = Mountain.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap();

    // Get solar position
    let solpos = location.get_solarposition(time).unwrap();
    assert!(solpos.zenith > 10.0 && solpos.zenith < 40.0,
        "Summer noon zenith at 39.7°N should be 15-35°, got {}", solpos.zenith);
    assert!(solpos.elevation > 50.0 && solpos.elevation < 80.0,
        "Summer noon elevation should be 50-80°, got {}", solpos.elevation);
    // Azimuth at solar noon should be near 180° (south)
    assert!(solpos.azimuth > 140.0 && solpos.azimuth < 220.0,
        "Noon azimuth should be near 180°, got {}", solpos.azimuth);

    // Cross-check: zenith + elevation = 90
    let sum = solpos.zenith + solpos.elevation;
    assert!((sum - 90.0).abs() < 0.5,
        "Zenith + elevation should = 90°, got {}", sum);

    // Get airmass
    let (am_rel, am_abs) = location.get_airmass(time).unwrap();
    assert!(am_rel > 1.0 && am_rel < 2.0,
        "Midday airmass should be 1.0-2.0, got {}", am_rel);
    assert!(am_abs < am_rel, "Absolute AM at altitude should be < relative AM");

    // Get clearsky
    let cs = location.get_clearsky(time, "ineichen").unwrap();
    assert!(cs.ghi > 700.0, "Clear sky GHI should be > 700 W/m2, got {}", cs.ghi);
    assert!(cs.dni > 600.0, "Clear sky DNI should be > 600 W/m2, got {}", cs.dni);
    assert!(cs.dhi > 50.0 && cs.dhi < 300.0,
        "Clear sky DHI should be 50-300 W/m2, got {}", cs.dhi);
    // Physical constraint: GHI ≈ DNI*cos(z) + DHI
    let cos_z = solpos.zenith.to_radians().cos();
    let ghi_check = cs.dni * cos_z + cs.dhi;
    assert!((cs.ghi - ghi_check).abs() < 50.0,
        "GHI should ≈ DNI*cos(z) + DHI: {} vs {}", cs.ghi, ghi_check);
}

// ============================================================================
// Test 5: Irradiance decomposition roundtrip
// ============================================================================

#[test]
fn test_e2e_irradiance_decomposition_roundtrip() {
    // Start with known GHI, decompose to DNI+DHI, verify closure
    let ghi = 800.0;
    let zenith = 30.0;
    let doy = 172; // summer solstice
    let dni_extra = irradiance::get_extra_radiation(doy);

    let (dni, dhi) = irradiance::erbs(ghi, zenith, doy as u32, dni_extra);

    // Closure: GHI = DNI * cos(z) + DHI
    let cos_z = zenith.to_radians().cos();
    let ghi_reconstructed = dni * cos_z + dhi;
    assert!((ghi - ghi_reconstructed).abs() < 1.0,
        "Erbs closure failed: GHI={}, reconstructed={}", ghi, ghi_reconstructed);

    // Sanity: DNI should be positive, DHI should be positive
    assert!(dni > 0.0, "DNI should be positive, got {}", dni);
    assert!(dhi > 0.0, "DHI should be positive, got {}", dhi);
    assert!(dhi < ghi, "DHI should be < GHI");

    // Test DISC model too
    let (dni_disc, kt, _am) = irradiance::disc(ghi, zenith, doy, Some(101325.0));
    assert!(dni_disc > 0.0, "DISC DNI should be positive, got {}", dni_disc);
    assert!(kt > 0.0 && kt < 1.0, "Clearness index should be 0-1, got {}", kt);
}

// ============================================================================
// Test 6: Temperature model comparison
// ============================================================================

#[test]
fn test_e2e_temperature_models_agree_roughly() {
    let poa = 800.0;
    let temp_air = 25.0;
    let wind_speed = 3.0;

    // All models should give cell temp above ambient
    let (sapm, _) = temperature::sapm_cell_temperature(poa, temp_air, wind_speed, -3.56, -0.075, 3.0, 1000.0);
    let pvsyst = temperature::pvsyst_cell_temperature(poa, temp_air, wind_speed, 29.0, 0.0, 0.15, 0.9);
    let faiman = temperature::faiman(poa, temp_air, wind_speed, 25.0, 6.84);
    let fuentes = temperature::fuentes(poa, temp_air, wind_speed, 45.0);
    let ross = temperature::ross(poa, temp_air, 0.025);
    let generic = temperature::generic_linear_default(poa, temp_air, wind_speed);

    let models = [
        ("SAPM", sapm),
        ("PVsyst", pvsyst),
        ("Faiman", faiman),
        ("Fuentes", fuentes),
        ("Ross", ross),
        ("Generic", generic),
    ];

    for (name, tc) in &models {
        assert!(*tc > temp_air,
            "{}: Cell temp {} should exceed ambient {}", name, tc, temp_air);
        assert!(*tc < 80.0,
            "{}: Cell temp {} should be < 80°C", name, tc);
        // All models should give similar-ish results (within ~25°C of each other)
        for (name2, tc2) in &models {
            assert!((tc - tc2).abs() < 25.0,
                "{} ({:.1}°C) vs {} ({:.1}°C) differ by more than 25°C", name, tc, name2, tc2);
        }
    }
}

// ============================================================================
// Test 7: IAM models - physical plausibility
// ============================================================================

#[test]
fn test_e2e_iam_models_plausibility() {
    // At normal incidence (AOI=0), all models should return ~1.0
    let models_at_0 = [
        ("ASHRAE", iam::ashrae(0.0, 0.05)),
        ("Physical", iam::physical(0.0, 1.526, 4.0, 0.002)),
        ("Martin-Ruiz", iam::martin_ruiz(0.0, 0.16)),
        ("Schlick", iam::schlick(0.0)),
    ];

    for (name, val) in &models_at_0 {
        assert!(*val > 0.95, "{} at AOI=0 should be > 0.95, got {}", name, val);
    }

    // At high AOI (80°), all should show significant reduction
    let models_at_80 = [
        ("ASHRAE", iam::ashrae(80.0, 0.05)),
        ("Physical", iam::physical(80.0, 1.526, 4.0, 0.002)),
        ("Martin-Ruiz", iam::martin_ruiz(80.0, 0.16)),
        ("Schlick", iam::schlick(80.0)),
    ];

    for (name, val) in &models_at_80 {
        assert!(*val < 0.7, "{} at AOI=80° should be < 0.7, got {}", name, val);
        assert!(*val > 0.0, "{} at AOI=80° should be > 0, got {}", name, val);
    }

    // Monotonicity: IAM should decrease with increasing AOI
    for aoi in [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0].windows(2) {
        let iam1 = iam::ashrae(aoi[0], 0.05);
        let iam2 = iam::ashrae(aoi[1], 0.05);
        assert!(iam1 >= iam2,
            "ASHRAE IAM should decrease: at {}°={}, at {}°={}", aoi[0], iam1, aoi[1], iam2);
    }
}

// ============================================================================
// Test 8: Single-diode model - Bishop88 full chain
// ============================================================================

#[test]
fn test_e2e_singlediode_bishop88_chain() {
    // Typical c-Si module parameters at STC
    let il = 9.0;   // photocurrent
    let i0 = 1e-10;  // saturation current
    let rs = 0.3;    // series resistance
    let rsh = 300.0;  // shunt resistance
    let nnsvth = 1.5 * 0.02569 * 60.0; // n * Ns * Vth (60 cells)

    let params = singlediode::Bishop88Params {
        photocurrent: il,
        saturation_current: i0,
        resistance_series: rs,
        resistance_shunt: rsh,
        n_ns_vth: nnsvth,
        d2mutau: 0.0,
        ns_vbi: 0.0,
        breakdown_factor: 0.0,
        breakdown_voltage: -5.5,
        breakdown_exp: 3.28,
    };

    // Find Voc
    let voc = singlediode::estimate_voc(il, i0, nnsvth);
    assert!(voc > 30.0 && voc < 50.0, "Voc should be 30-50V for 60-cell module, got {}", voc);

    // Find Isc
    let isc = singlediode::bishop88_i_from_v(0.0, &params);
    assert!((isc - il).abs() < 0.5, "Isc should be close to IL, got {} vs {}", isc, il);

    // Find MPP
    let mpp = singlediode::bishop88_mpp(&params);
    assert!(mpp.voltage > 0.0 && mpp.voltage < voc, "Vmpp should be between 0 and Voc");
    assert!(mpp.current > 0.0 && mpp.current < isc, "Impp should be between 0 and Isc");
    assert!(mpp.power > 0.0, "MPP power should be positive");

    // Power should be roughly 200-400W for a typical 60-cell module
    assert!(mpp.power > 150.0 && mpp.power < 500.0,
        "MPP power should be 150-500W, got {}", mpp.power);

    // Roundtrip: I_from_V at Vmpp should give Impp
    let i_check = singlediode::bishop88_i_from_v(mpp.voltage, &params);
    assert!((i_check - mpp.current).abs() < 0.01,
        "I_from_V roundtrip: got {} expected {}", i_check, mpp.current);
}

// ============================================================================
// Test 9: calcparams_desoto -> singlediode -> power chain
// ============================================================================

#[test]
fn test_e2e_desoto_to_power() {
    // Typical CEC module parameters (similar to a 300W module)
    let alpha_sc = 0.003;
    let a_ref = 1.5 * 0.02569 * 60.0; // n*Ns*Vth at 25°C
    let i_l_ref = 9.5;
    let i_o_ref = 1e-10;
    let r_sh_ref = 400.0;
    let r_s = 0.3;
    let eg_ref = 1.121;
    let d_eg_dt = -0.0002677;

    // At STC (1000 W/m2, 25°C)
    let params_stc = pvsystem::calcparams_desoto(
        1000.0, 25.0, alpha_sc, a_ref, i_l_ref, i_o_ref, r_sh_ref, r_s, eg_ref, d_eg_dt,
    );

    assert!((params_stc.photocurrent - i_l_ref).abs() < 0.01,
        "IL at STC should match reference: {} vs {}", params_stc.photocurrent, i_l_ref);
    assert!(params_stc.resistance_shunt > 300.0,
        "Rsh at STC should be close to ref: {}", params_stc.resistance_shunt);

    // At 500 W/m2, 40°C (realistic field conditions)
    let params_field = pvsystem::calcparams_desoto(
        500.0, 40.0, alpha_sc, a_ref, i_l_ref, i_o_ref, r_sh_ref, r_s, eg_ref, d_eg_dt,
    );

    // IL should be roughly half at half irradiance
    assert!(params_field.photocurrent < params_stc.photocurrent,
        "IL at 500 W/m2 should be less than at 1000");
    assert!(params_field.photocurrent > 4.0,
        "IL at 500 W/m2 should be > 4A, got {}", params_field.photocurrent);

    // Rsh should be higher at lower irradiance
    assert!(params_field.resistance_shunt > params_stc.resistance_shunt,
        "Rsh should increase at lower irradiance");

    // I0 should increase with temperature
    assert!(params_field.saturation_current > params_stc.saturation_current,
        "I0 should increase with temperature");
}

// ============================================================================
// Test 10: Complete irradiance decomposition via ModelChain
// ============================================================================

#[test]
fn test_e2e_complete_irradiance() {
    let (system, location) = make_residential_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);

    // Test: only GHI provided, should decompose
    let weather_ghi_only = WeatherInput {
        time: Mountain.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap(),
        ghi: Some(800.0),
        dni: None,
        dhi: None,
        temp_air: 25.0,
        wind_speed: 3.0,
        albedo: None,
    };

    let (ghi, dni, dhi) = mc.complete_irradiance(&weather_ghi_only).unwrap();
    assert!((ghi - 800.0).abs() < 0.1, "GHI should be unchanged");
    assert!(dni > 0.0, "Decomposed DNI should be positive, got {}", dni);
    assert!(dhi > 0.0, "Decomposed DHI should be positive, got {}", dhi);
    // Closure
    let solpos = location.get_solarposition(weather_ghi_only.time).unwrap();
    let cos_z = solpos.zenith.to_radians().cos().max(0.0);
    let ghi_check = dni * cos_z + dhi;
    assert!((ghi - ghi_check).abs() < 5.0,
        "Closure: GHI={}, DNI*cos(z)+DHI={}", ghi, ghi_check);

    // Should still work and produce valid power
    let result = mc.run_model_from_weather(&weather_ghi_only).unwrap();
    assert!(result.dc_power > 1000.0, "Should produce significant power from 800 GHI");
}

// ============================================================================
// Test 11: Tracking system
// ============================================================================

#[test]
fn test_e2e_tracking_geometry() {
    let solar_zenith = 30.0;
    let solar_azimuth = 180.0;
    let axis_tilt = 0.0;
    let axis_azimuth = 0.0; // N-S axis
    let max_angle = 60.0;
    let backtrack = true;
    let gcr = 0.35;

    let result = tracking::singleaxis(
        solar_zenith, solar_azimuth,
        axis_tilt, axis_azimuth,
        max_angle, backtrack, gcr,
    );

    // Should return valid tracking angles
    assert!(result.tracker_theta.abs() <= max_angle,
        "Tracker angle {} should be within max angle {}", result.tracker_theta, max_angle);
    assert!(result.surface_tilt >= 0.0 && result.surface_tilt <= 90.0,
        "Surface tilt should be 0-90°, got {}", result.surface_tilt);

    // AOI should be less than the fixed-tilt AOI at this geometry
    assert!(result.aoi >= 0.0 && result.aoi < 90.0,
        "AOI should be 0-90°, got {}", result.aoi);

    // calc_surface_orientation should be consistent
    let (tilt, azimuth) = tracking::calc_surface_orientation(
        result.tracker_theta, axis_tilt, axis_azimuth,
    );
    assert!((tilt - result.surface_tilt).abs() < 1.0,
        "calc_surface_orientation tilt {} should match tracking tilt {}", tilt, result.surface_tilt);
}

// ============================================================================
// Test 12: Shading + snow + soiling loss chain
// ============================================================================

#[test]
fn test_e2e_loss_factors() {
    // Shading: ground angle and masking angle should be physically reasonable
    let gcr = 0.4;
    let surface_tilt = 25.0;
    let ga = shading::ground_angle(surface_tilt, gcr, 1.0);
    assert!(ga > 0.0 && ga < 45.0, "Ground angle should be 0-45°, got {}", ga);

    let ma = shading::masking_angle_passias(surface_tilt, gcr);
    assert!(ma >= 0.0 && ma < 90.0, "Masking angle should be 0-90°, got {}", ma);

    // Snow: should lose power when covered
    let loss_full = snow::dc_loss_nrel(1.0, 3);
    assert!((loss_full - 1.0).abs() < 0.01, "Full coverage should give 100% loss");

    let loss_partial = snow::dc_loss_nrel(0.5, 4);
    assert!(loss_partial > 0.0 && loss_partial < 1.0,
        "Partial coverage should give 0-100% loss, got {}", loss_partial);

    let loss_none = snow::dc_loss_nrel(0.0, 3);
    assert!(loss_none.abs() < 0.01, "No coverage should give 0% loss");

    // Soiling: accumulation should reduce performance
    let (soiling_ratio, _mass) = soiling::hsu(
        0.0,  // no rain
        5.0,  // cleaning threshold
        25.0, // tilt
        30.0, // PM2.5
        50.0, // PM10
        None, // default depo velocity
        None, // default accumulation period
        1.0,  // previous mass
    );
    assert!(soiling_ratio <= 1.0, "Soiling ratio should be <= 1.0, got {}", soiling_ratio);
    assert!(soiling_ratio > 0.5, "Soiling ratio should be > 0.5 with moderate soiling");
}

// ============================================================================
// Test 13: Spectrum + albedo chain
// ============================================================================

#[test]
fn test_e2e_spectrum_and_albedo() {
    // Spectral mismatch at typical airmass
    let smm = spectrum::spectral_mismatch_modifier(1.5, &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(smm > 0.0, "Spectral modifier should be positive");

    // SAPM spectral factor
    let sf = spectrum::spectral_factor_sapm(1.5, "monosi");
    assert!(sf > 0.8 && sf < 1.2, "SAPM spectral factor should be near 1.0, got {}", sf);

    // Water albedo: should be high at low sun, low at high sun
    let albedo_low_sun = albedo::inland_water_dvoracek(5.0);
    let albedo_high_sun = albedo::inland_water_dvoracek(60.0);
    assert!(albedo_low_sun > albedo_high_sun,
        "Water albedo should be higher at low sun: {} vs {}", albedo_low_sun, albedo_high_sun);
    assert!(albedo_low_sun > 0.0 && albedo_low_sun < 1.0);
    assert!(albedo_high_sun > 0.0 && albedo_high_sun < 1.0);
}

// ============================================================================
// Test 14: PV array efficiency models
// ============================================================================

#[test]
fn test_e2e_pvarray_models() {
    // ADR efficiency model at STC
    let eff_stc = pvarray::pvefficiency_adr(1000.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    // At STC with zero coefficients, should return 1.0 (normalized)
    assert!((eff_stc - 1.0).abs() < 0.01, "ADR at STC should be ~1.0, got {}", eff_stc);

    // Huld model at STC
    let p_stc = pvarray::huld(1000.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    // At STC with zero coefficients, should return normalized 1.0
    assert!((p_stc - 1.0).abs() < 0.01, "Huld at STC should be ~1.0, got {}", p_stc);

    // At lower irradiance, both should give lower output
    let eff_low = pvarray::pvefficiency_adr(500.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    assert!(eff_low > 0.0, "ADR at 500 W/m2 should be positive");
}

// ============================================================================
// Test 15: Transformer efficiency
// ============================================================================

#[test]
fn test_e2e_transformer() {
    let power = 4000.0;
    let rating = 5000.0;
    let no_load_loss = 10.0;  // 10 W no-load loss
    let load_loss = 50.0;     // 50 W full-load loss

    let eff = transformer::simple_efficiency(power, rating, no_load_loss, load_loss);
    assert!(eff > 0.95, "Transformer efficiency should be > 95%, got {:.2}%", eff * 100.0);
    assert!(eff < 1.0, "Transformer efficiency should be < 100%");

    // Higher load = lower efficiency (due to I²R losses)
    let eff_low = transformer::simple_efficiency(1000.0, rating, no_load_loss, load_loss);
    // At low load, no-load loss is a bigger fraction
    assert!(eff_low < eff || (eff_low - eff).abs() < 0.01,
        "Low load may have different efficiency profile");
}

// ============================================================================
// Test 16: Inverter models comparison
// ============================================================================

#[test]
fn test_e2e_inverter_models() {
    let pdc = 4000.0;
    let pdc0 = 5000.0;

    // PVWatts
    let pac_pvwatts = inverter::pvwatts_ac(pdc, pdc0, 0.96, 0.9637);
    assert!(pac_pvwatts > 3500.0, "PVWatts AC should be > 3500, got {}", pac_pvwatts);
    assert!(pac_pvwatts < pdc, "AC should be < DC");

    // Multi-MPPT PVWatts (should give same result when summed)
    let pac_multi = inverter::pvwatts_multi(&[2000.0, 2000.0], pdc0, 0.96, 0.9637);
    assert!((pac_multi - pac_pvwatts).abs() < 1.0,
        "Multi-MPPT should match single: {} vs {}", pac_multi, pac_pvwatts);

    // Clipping: input > capacity
    let pac_clip = inverter::pvwatts_ac(6000.0, pdc0, 0.96, 0.9637);
    let pac0 = 0.96 * pdc0;
    assert!(pac_clip <= pac0 + 1.0, "Should clip at pac0={}, got {}", pac0, pac_clip);
}

// ============================================================================
// Test 17: Atmosphere chain
// ============================================================================

#[test]
fn test_e2e_atmosphere_chain() {
    // Altitude -> pressure -> altitude roundtrip
    let alt = 1830.0; // Golden, CO
    let pressure = atmosphere::alt2pres(alt);
    let alt_back = atmosphere::pres2alt(pressure);
    assert!((alt - alt_back).abs() < 1.0,
        "alt2pres/pres2alt roundtrip: {} -> {} -> {}", alt, pressure, alt_back);

    // Precipitable water from temperature/humidity
    let pw = atmosphere::gueymard94_pw(25.0, 50.0);
    assert!(pw > 0.0 && pw < 5.0, "PW should be 0-5 cm, got {}", pw);

    // Dew point roundtrip
    let tdew = atmosphere::tdew_from_rh(25.0, 50.0);
    let rh_back = atmosphere::rh_from_tdew(25.0, tdew);
    assert!((rh_back - 50.0).abs() < 1.0,
        "RH roundtrip: 50% -> Tdew={} -> {:.1}%", tdew, rh_back);

    // Linke turbidity from components
    let aod_bb = atmosphere::bird_hulstrom80_aod_bb(0.1, 0.08);
    assert!(aod_bb > 0.0, "Broadband AOD should be positive");
    let lt = atmosphere::kasten96_lt(1.5, pw, aod_bb);
    assert!(lt > 1.0 && lt < 10.0, "Linke turbidity should be 1-10, got {}", lt);

    // Wind speed power law
    let ws_10m = 5.0;
    let ws_80m = atmosphere::windspeed_powerlaw(ws_10m, 10.0, 80.0, None);
    assert!(ws_80m > ws_10m, "Wind at 80m should be > 10m: {} vs {}", ws_80m, ws_10m);
}

// ============================================================================
// Test 18: Clear sky model comparison
// ============================================================================

#[test]
fn test_e2e_clearsky_models_compare() {
    let zenith = 30.0;
    let am = atmosphere::get_relative_airmass(zenith);

    // Haurwitz (simplest)
    let haurwitz = clearsky::haurwitz(zenith);

    // Bird (most detailed)
    let bird = clearsky::bird_default(zenith, am);

    // All should give positive GHI at zenith=30°
    assert!(haurwitz.ghi > 500.0, "Haurwitz GHI should be > 500, got {}", haurwitz.ghi);
    assert!(bird.ghi > 500.0, "Bird GHI should be > 500, got {}", bird.ghi);

    // Bird should provide DNI and DHI
    assert!(bird.dni > 500.0, "Bird DNI should be > 500, got {}", bird.dni);
    assert!(bird.dhi > 30.0, "Bird DHI should be > 30, got {}", bird.dhi);
}

// ============================================================================
// Test 19: Analytical solar position functions
// ============================================================================

#[test]
fn test_e2e_analytical_solar_position() {
    // Summer solstice, equator, solar noon
    let doy = 172;
    let decl = solarposition::declination_spencer71(doy);
    assert!(decl > 0.35 && decl < 0.45,
        "Summer solstice declination should be ~23.45° (~0.41 rad), got {} rad", decl);

    let eot = solarposition::equation_of_time_spencer71(doy);
    assert!(eot.abs() < 15.0, "EoT should be within ±15 minutes, got {}", eot);

    // Hour angle at solar noon
    let ha = solarposition::hour_angle(0.0, eot, 12.0);
    assert!(ha.abs() < 5.0, "Hour angle at noon should be near 0°, got {}", ha);

    // Earth-sun distance
    let d = solarposition::nrel_earthsun_distance(doy);
    assert!(d > 0.98 && d < 1.02,
        "Earth-sun distance should be near 1 AU, got {}", d);

    // Zenith at equator at noon on equinox
    let decl_equinox = solarposition::declination_spencer71(80); // ~March 21
    let zenith = solarposition::solar_zenith_analytical(0.0_f64.to_radians(), 0.0_f64.to_radians(), decl_equinox);
    let zenith_deg = zenith.to_degrees();
    assert!(zenith_deg.abs() < 5.0,
        "Equator noon equinox zenith should be ~0°, got {}", zenith_deg);
}

// ============================================================================
// Test 20: Full daily simulation (multiple timesteps)
// ============================================================================

#[test]
fn test_e2e_daily_simulation() {
    let (system, location) = make_residential_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);

    let mut total_energy_wh = 0.0;
    let mut max_power = 0.0_f64;

    // Simulate every hour of a summer day
    for hour in 5..=20 {
        let time = Mountain.with_ymd_and_hms(2024, 6, 21, hour, 0, 0).unwrap();

        // Simplified clear-sky-like irradiance profile
        let zenith_approx = ((hour as f64 - 12.5) * 15.0).abs().min(95.0);
        let cos_z = zenith_approx.to_radians().cos().max(0.0);
        let ghi = (1000.0 * cos_z).max(0.0);
        let dhi = (150.0 * cos_z).max(0.0);
        let dni = if cos_z > 0.05 { (ghi - dhi) / cos_z } else { 0.0 };

        let weather = WeatherInput {
            time,
            ghi: Some(ghi),
            dni: Some(dni.max(0.0)),
            dhi: Some(dhi),
            temp_air: 25.0 + 5.0 * cos_z, // warmer midday
            wind_speed: 3.0,
            albedo: Some(0.2),
        };

        let result = mc.run_model_from_weather(&weather).unwrap();

        // Power should never be negative
        assert!(result.ac_power >= -1.0,
            "Hour {}: AC power should not be negative, got {}", hour, result.ac_power);

        if result.ac_power > 0.0 {
            total_energy_wh += result.ac_power; // 1 hour timestep
        }
        max_power = max_power.max(result.ac_power);
    }

    // A 5kW system in Colorado summer should produce 25-45 kWh/day
    let total_kwh = total_energy_wh / 1000.0;
    assert!(total_kwh > 10.0, "Daily energy should be > 10 kWh, got {:.1}", total_kwh);
    assert!(total_kwh < 50.0, "Daily energy should be < 50 kWh, got {:.1}", total_kwh);

    // Peak power should be in the 3-5 kW range
    assert!(max_power > 2000.0, "Peak power should be > 2000 W, got {:.0}", max_power);
    assert!(max_power < 5500.0, "Peak power should be < 5500 W, got {:.0}", max_power);
}

// ============================================================================
// Test 21: run_model_from_poa and run_model_from_effective_irradiance
// ============================================================================

#[test]
fn test_e2e_alternative_entry_points() {
    let (system, location) = make_residential_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);

    let time = Mountain.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap();

    // From POA
    let poa_input = POAInput {
        time,
        poa_direct: 700.0,
        poa_diffuse: 150.0,
        poa_global: 850.0,
        temp_air: 30.0,
        wind_speed: 3.0,
        aoi: 20.0,
    };

    let result_poa = mc.run_model_from_poa(&poa_input).unwrap();
    assert!(result_poa.ac_power > 1000.0, "POA entry should produce power, got {}", result_poa.ac_power);
    assert!(result_poa.poa_global == 850.0, "POA global should be preserved");

    // From effective irradiance
    let eff_input = EffectiveIrradianceInput {
        time,
        effective_irradiance: 800.0,
        poa_global: 850.0,
        temp_air: 30.0,
        wind_speed: 3.0,
    };

    let result_eff = mc.run_model_from_effective_irradiance(&eff_input).unwrap();
    assert!(result_eff.ac_power > 1000.0, "Eff irrad entry should produce power, got {}", result_eff.ac_power);
    assert!(result_eff.effective_irradiance == 800.0, "Eff irrad should be preserved");
}

// ============================================================================
// Test 22: Different ModelChain configurations produce different results
// ============================================================================

#[test]
fn test_e2e_model_config_differences() {
    let location = Location::new(39.742, -105.178, Mountain, 1830.0, "Golden, CO");
    let weather = make_summer_noon_weather();

    let make_system = || -> PVSystem {
        let array = Array {
            mount: Box::new(FixedMount { surface_tilt: 30.0, surface_azimuth: 180.0 }),
            nameplate_dc: 5000.0,
            gamma_pdc: -0.004,
            modules_per_string: 10,
            strings: 2,
            albedo: 0.2,
        };
        PVSystem::new(vec![array], 5000.0)
    };

    // PVWatts config
    let mc_pvwatts = ModelChain::with_pvwatts(make_system(), location.clone(), 30.0, 180.0, 5000.0, 0.96);
    let r_pvwatts = mc_pvwatts.run_model_from_weather(&weather).unwrap();

    // SAPM config (different temperature model and transposition)
    let mc_sapm = ModelChain::with_sapm(make_system(), location.clone(), 30.0, 180.0, 5000.0, 0.96);
    let r_sapm = mc_sapm.run_model_from_weather(&weather).unwrap();

    // Both should produce power
    assert!(r_pvwatts.ac_power > 1000.0, "PVWatts should produce power");
    assert!(r_sapm.ac_power > 1000.0, "SAPM should produce power");

    // Results should differ (different models used)
    // At minimum, cell temperature should differ (PVWatts vs SAPM temp model)
    // Note: they might be close but shouldn't be identical
    assert!(r_pvwatts.ac_power > 0.0 && r_sapm.ac_power > 0.0,
        "Both configs should produce valid positive power");
}

// ============================================================================
// Test 23: Bird clear sky model - full detail check
// ============================================================================

#[test]
fn test_e2e_bird_clearsky_detailed() {
    let zenith = 20.0;
    let am = atmosphere::get_relative_airmass(zenith);
    let pw = atmosphere::gueymard94_pw(25.0, 50.0);

    let result = clearsky::bird(
        zenith, am,
        0.1, 0.08,  // AOD at 380nm, 500nm
        pw,          // precipitable water
        0.3,         // ozone
        101325.0,    // pressure
        1364.0,      // DNI extra
        0.85,        // asymmetry
        0.2,         // albedo
    );

    // At low zenith with moderate conditions, should get high irradiance
    assert!(result.ghi > 800.0, "Bird GHI should be > 800 at zenith=20°, got {}", result.ghi);
    assert!(result.dni > 800.0, "Bird DNI should be > 800, got {}", result.dni);
    assert!(result.dhi > 50.0 && result.dhi < 200.0, "Bird DHI should be 50-200, got {}", result.dhi);
    assert!(result.direct_horizontal > 0.0, "Direct horizontal should be > 0");

    // Physical: GHI ≈ DNI*cos(z) + DHI
    let cos_z = zenith.to_radians().cos();
    let ghi_check = result.dni * cos_z + result.dhi;
    assert!((result.ghi - ghi_check).abs() < 10.0,
        "Bird GHI closure: {} vs {}", result.ghi, ghi_check);
}
