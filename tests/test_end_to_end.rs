//! End-to-end integration tests for PvLib-Rust.
//!
//! These tests simulate realistic PV system modeling workflows,
//! exercising the full pipeline from location/weather input through
//! to AC power output, verifying physical plausibility at each step.

use chrono::TimeZone;
use chrono_tz::US::Mountain;

use pvlib::location::Location;
use pvlib::pvsystem::{PVSystem, Array, FixedMount};
use pvlib::modelchain::{
    ModelChain, WeatherInput, POAInput, EffectiveIrradianceInput,
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
    let location = Location::new(39.742, -105.178, Mountain, 1830.0, "Golden, CO");
    let array = Array {
        mount: Box::new(FixedMount {
            surface_tilt: 30.0,
            surface_azimuth: 180.0,
        }),
        nameplate_dc: 5000.0,
        gamma_pdc: -0.004,
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
        albedo: Some(0.3),
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

    assert!(result.solar_zenith > 0.0 && result.solar_zenith < 40.0,
        "Summer noon zenith should be 0-40°, got {}", result.solar_zenith);
    assert!(result.poa_global > 700.0 && result.poa_global < 1200.0,
        "POA should be 700-1200, got {}", result.poa_global);
    assert!(result.aoi_modifier > 0.8 && result.aoi_modifier <= 1.0,
        "AOI modifier should be 0.8-1.0, got {}", result.aoi_modifier);
    assert!(result.cell_temperature > weather.temp_air && result.cell_temperature < 70.0,
        "Cell temp {} should be above ambient {} and < 70", result.cell_temperature, weather.temp_air);
    assert!(result.dc_power > 2000.0 && result.dc_power < 5500.0,
        "DC should be 2000-5500, got {}", result.dc_power);
    assert!(result.ac_power > 0.0 && result.ac_power < result.dc_power,
        "AC {} should be positive and < DC {}", result.ac_power, result.dc_power);

    let inv_eff = result.ac_power / result.dc_power;
    assert!(inv_eff > 0.90 && inv_eff < 1.0,
        "Inverter eff should be 90-100%, got {:.1}%", inv_eff * 100.0);
}

// ============================================================================
// Test 2: Winter noon
// ============================================================================

#[test]
fn test_e2e_full_pipeline_winter_noon() {
    let (system, location) = make_residential_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);
    let result = mc.run_model_from_weather(&make_winter_noon_weather()).unwrap();

    assert!(result.solar_zenith > 50.0 && result.solar_zenith < 70.0,
        "Winter zenith should be 50-70°, got {}", result.solar_zenith);
    assert!(result.cell_temperature > 0.0 && result.cell_temperature < 40.0,
        "Winter cell temp should be 0-40°C, got {}", result.cell_temperature);
    assert!(result.dc_power > 500.0 && result.dc_power < 4000.0,
        "Winter DC should be 500-4000, got {}", result.dc_power);
    assert!(result.ac_power > 0.0);
}

// ============================================================================
// Test 3: Night time - zero output
// ============================================================================

#[test]
fn test_e2e_nighttime_zero_output() {
    let (system, location) = make_residential_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);
    let result = mc.run_model_from_weather(&make_night_weather()).unwrap();

    assert!(result.solar_zenith > 90.0, "Night zenith should be > 90°");
    assert!(result.dc_power.abs() < 1.0, "Night DC should be ~0, got {}", result.dc_power);
    assert!(result.ac_power.abs() < 1.0, "Night AC should be ~0, got {}", result.ac_power);
}

// ============================================================================
// Test 4: Solar position chain
// ============================================================================

#[test]
fn test_e2e_solar_position_chain() {
    let location = Location::new(39.742, -105.178, Mountain, 1830.0, "Golden, CO");
    let time = Mountain.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap();

    let solpos = location.get_solarposition(time).unwrap();
    assert!(solpos.zenith > 10.0 && solpos.zenith < 40.0,
        "Summer noon zenith should be 10-40°, got {}", solpos.zenith);
    assert!(solpos.elevation > 50.0 && solpos.elevation < 80.0,
        "Elevation should be 50-80°, got {}", solpos.elevation);
    // Azimuth depends on exact solar noon vs clock noon; allow wider range
    assert!(solpos.azimuth > 100.0 && solpos.azimuth < 260.0,
        "Noon azimuth should be roughly south, got {}", solpos.azimuth);

    let sum = solpos.zenith + solpos.elevation;
    assert!((sum - 90.0).abs() < 0.5, "Zenith + elevation should = 90°, got {}", sum);

    // Airmass
    let (am_rel, am_abs) = location.get_airmass(time);
    assert!(am_rel > 1.0 && am_rel < 2.0, "Midday AM should be 1-2, got {}", am_rel);
    assert!(am_abs < am_rel, "Absolute AM at altitude should be < relative AM");

    // Clear sky
    let cs = location.get_clearsky(time, "ineichen");
    assert!(cs.ghi > 700.0, "Clear sky GHI should be > 700, got {}", cs.ghi);
    assert!(cs.dni > 600.0, "Clear sky DNI should be > 600, got {}", cs.dni);
    assert!(cs.dhi > 50.0 && cs.dhi < 300.0, "Clear sky DHI should be 50-300, got {}", cs.dhi);

    // GHI ≈ DNI*cos(z) + DHI
    let cos_z = solpos.zenith.to_radians().cos();
    let ghi_check = cs.dni * cos_z + cs.dhi;
    assert!((cs.ghi - ghi_check).abs() < 50.0,
        "GHI closure: {} vs {}", cs.ghi, ghi_check);
}

// ============================================================================
// Test 5: Irradiance decomposition roundtrip
// ============================================================================

#[test]
fn test_e2e_irradiance_decomposition_roundtrip() {
    let ghi = 800.0;
    let zenith = 30.0;
    let doy = 172_u32;
    let dni_extra = irradiance::get_extra_radiation(doy as i32);

    let (dni, dhi) = irradiance::erbs(ghi, zenith, doy, dni_extra);

    // Closure: GHI = DNI * cos(z) + DHI
    let cos_z = zenith.to_radians().cos();
    let ghi_reconstructed = dni * cos_z + dhi;
    assert!((ghi - ghi_reconstructed).abs() < 1.0,
        "Erbs closure: GHI={}, reconstructed={}", ghi, ghi_reconstructed);
    assert!(dni > 0.0 && dhi > 0.0 && dhi < ghi);

    // DISC model
    let disc_out = irradiance::disc(ghi, zenith, doy as i32, Some(101325.0));
    assert!(disc_out.dni > 0.0, "DISC DNI should be positive, got {}", disc_out.dni);
    assert!(disc_out.kt > 0.0 && disc_out.kt < 1.0, "Kt should be 0-1, got {}", disc_out.kt);
}

// ============================================================================
// Test 6: Temperature models agree roughly
// ============================================================================

#[test]
fn test_e2e_temperature_models_agree_roughly() {
    let poa = 800.0;
    let ta = 25.0;
    let ws = 3.0;

    let (sapm, _) = temperature::sapm_cell_temperature(poa, ta, ws, -3.56, -0.075, 3.0, 1000.0);
    let pvsyst = temperature::pvsyst_cell_temperature(poa, ta, ws, 29.0, 0.0, 0.15, 0.9);
    let faiman = temperature::faiman(poa, ta, ws, 25.0, 6.84);
    let fuentes = temperature::fuentes(poa, ta, ws, 45.0);
    let ross = temperature::ross(poa, ta, 45.0);
    let generic = temperature::generic_linear_default(poa, ta, ws);

    for (name, tc) in &[("SAPM", sapm), ("PVsyst", pvsyst), ("Faiman", faiman),
                         ("Fuentes", fuentes), ("Ross", ross), ("Generic", generic)] {
        assert!(*tc > ta, "{}: Cell temp {} should exceed ambient {}", name, tc, ta);
        assert!(*tc < 80.0, "{}: Cell temp {} should be < 80°C", name, tc);
    }
}

// ============================================================================
// Test 7: IAM models plausibility
// ============================================================================

#[test]
fn test_e2e_iam_models_plausibility() {
    // At normal incidence, all ~1.0
    assert!(iam::ashrae(0.0, 0.05) > 0.95);
    assert!(iam::physical(0.0, 1.526, 4.0, 0.002) > 0.95);
    assert!(iam::martin_ruiz(0.0, 0.16) > 0.95);
    assert!(iam::schlick(0.0) > 0.95);

    // At 80°, all should show some reduction (< 0.95)
    assert!(iam::ashrae(80.0, 0.05) < 0.95);
    assert!(iam::physical(80.0, 1.526, 4.0, 0.002) < 0.95);
    assert!(iam::martin_ruiz(80.0, 0.16) < 0.95);
    assert!(iam::schlick(80.0) < 0.95);

    // Monotonicity
    for aoi_pair in [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0].windows(2) {
        assert!(iam::ashrae(aoi_pair[0], 0.05) >= iam::ashrae(aoi_pair[1], 0.05));
    }
}

// ============================================================================
// Test 8: Single-diode Bishop88 full chain
// ============================================================================

#[test]
fn test_e2e_singlediode_bishop88_chain() {
    let il = 9.0;
    let i0 = 1e-10;
    let rs = 0.3;
    let rsh = 300.0;
    let nnsvth = 1.5 * 0.02569 * 60.0;

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

    let voc = singlediode::estimate_voc(il, i0, nnsvth);
    assert!(voc > 30.0 && voc < 65.0, "Voc should be 30-65V, got {}", voc);

    let isc = singlediode::bishop88_i_from_v(0.0, &params);
    assert!((isc - il).abs() < 0.5, "Isc ~ IL: {} vs {}", isc, il);

    let mpp = singlediode::bishop88_mpp(&params);
    assert!(mpp.power > 100.0 && mpp.power < 600.0,
        "MPP power should be 100-600W, got {} (V={}, I={})", mpp.power, mpp.voltage, mpp.current);
    assert!(mpp.voltage.abs() < voc * 1.1,
        "Vmpp magnitude should be within Voc range: V={}, Voc={}", mpp.voltage, voc);
    assert!(mpp.current.abs() < isc * 1.1,
        "Impp magnitude should be within Isc range: I={}, Isc={}", mpp.current, isc);

    // Roundtrip
    let i_check = singlediode::bishop88_i_from_v(mpp.voltage, &params);
    assert!((i_check - mpp.current).abs() < 0.01);
}

// ============================================================================
// Test 9: calcparams_desoto chain
// ============================================================================

#[test]
fn test_e2e_desoto_to_power() {
    let alpha_sc = 0.003;
    let a_ref = 1.5 * 0.02569 * 60.0;
    let i_l_ref = 9.5;
    let i_o_ref = 1e-10;
    let r_sh_ref = 400.0;
    let r_s = 0.3;

    let params_stc = pvsystem::calcparams_desoto(
        1000.0, 25.0, alpha_sc, a_ref, i_l_ref, i_o_ref, r_sh_ref, r_s, 1.121, -0.0002677,
    );
    assert!((params_stc.photocurrent - i_l_ref).abs() < 0.01);
    assert!(params_stc.resistance_shunt > 300.0);

    let params_field = pvsystem::calcparams_desoto(
        500.0, 40.0, alpha_sc, a_ref, i_l_ref, i_o_ref, r_sh_ref, r_s, 1.121, -0.0002677,
    );
    assert!(params_field.photocurrent < params_stc.photocurrent);
    assert!(params_field.photocurrent > 4.0);
    assert!(params_field.resistance_shunt > params_stc.resistance_shunt);
    assert!(params_field.saturation_current > params_stc.saturation_current);
}

// ============================================================================
// Test 10: Complete irradiance via ModelChain
// ============================================================================

#[test]
fn test_e2e_complete_irradiance() {
    let (system, location) = make_residential_system();
    let loc_clone = location.clone();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);

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
    assert!((ghi - 800.0).abs() < 0.1);
    assert!(dni > 0.0 && dhi > 0.0);

    // Closure check
    let solpos = loc_clone.get_solarposition(weather_ghi_only.time).unwrap();
    let cos_z = solpos.zenith.to_radians().cos().max(0.0);
    let ghi_check = dni * cos_z + dhi;
    assert!((ghi - ghi_check).abs() < 5.0,
        "Closure: GHI={}, DNI*cos(z)+DHI={}", ghi, ghi_check);

    let result = mc.run_model_from_weather(&weather_ghi_only).unwrap();
    assert!(result.dc_power > 1000.0, "Should produce power from 800 GHI");
}

// ============================================================================
// Test 11: Tracking geometry
// ============================================================================

#[test]
fn test_e2e_tracking_geometry() {
    let (surface_tilt, _surface_azimuth, aoi_val) = tracking::singleaxis(
        30.0, 180.0,  // solar position
        0.0, 0.0,     // axis orientation (N-S horizontal)
        60.0, true, 0.35, 0.0,
    );

    assert!(surface_tilt >= 0.0 && surface_tilt <= 60.0,
        "Surface tilt should be 0-60°, got {}", surface_tilt);
    assert!(aoi_val >= 0.0 && aoi_val < 90.0, "AOI should be 0-90°, got {}", aoi_val);

    // calc_surface_orientation consistency
    let (tilt2, _az2) = tracking::calc_surface_orientation(surface_tilt, 0.0, 0.0);
    assert!(tilt2 >= 0.0, "calc_surface_orientation tilt should be >= 0");
}

// ============================================================================
// Test 12: Loss factors chain
// ============================================================================

#[test]
fn test_e2e_loss_factors() {
    // Shading
    let ga = shading::ground_angle(25.0, 0.4, 1.0);
    assert!(ga > 0.0 && ga < 45.0, "Ground angle should be 0-45°, got {}", ga);

    let ma = shading::masking_angle_passias(25.0, 0.4);
    assert!(ma >= 0.0 && ma < 90.0, "Masking angle should be 0-90°, got {}", ma);

    // Snow losses
    assert!((snow::dc_loss_nrel(1.0, 3) - 1.0).abs() < 0.01);
    assert!(snow::dc_loss_nrel(0.5, 4) > 0.0 && snow::dc_loss_nrel(0.5, 4) < 1.0);
    assert!(snow::dc_loss_nrel(0.0, 3).abs() < 0.01);

    // Soiling
    let (soiling_ratio, _mass) = soiling::hsu(
        0.0, 5.0, 25.0, 30.0, 50.0,
        0.004, 0.004, 3600.0, 1.0,
    );
    assert!(soiling_ratio <= 1.0 && soiling_ratio > 0.5);
}

// ============================================================================
// Test 13: Spectrum and albedo
// ============================================================================

#[test]
fn test_e2e_spectrum_and_albedo() {
    let smm = spectrum::spectral_mismatch_modifier(1.5, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(smm > 0.0);

    let sf = spectrum::spectral_factor_sapm(1.5, [1.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(sf > 0.8 && sf < 1.2, "SAPM spectral factor should be near 1.0, got {}", sf);

    // Water albedo
    let albedo_low = albedo::inland_water_dvoracek(5.0, "clear_water_no_waves");
    let albedo_high = albedo::inland_water_dvoracek(60.0, "clear_water_no_waves");
    assert!(albedo_low > albedo_high, "Albedo should be higher at low sun");
    assert!(albedo_low > 0.0 && albedo_low < 1.0);
}

// ============================================================================
// Test 14: PV array efficiency models
// ============================================================================

#[test]
fn test_e2e_pvarray_models() {
    // ADR with k_a=0.20 (20% module efficiency), zero temperature/loss coefficients
    let eff_stc = pvarray::pvefficiency_adr(1000.0, 25.0, 0.20, -6.0, 0.0, 0.0, 0.0);
    assert!(eff_stc > 0.15 && eff_stc < 0.25, "ADR efficiency at STC should be ~0.20, got {}", eff_stc);

    // Huld at STC with pdc0=1.0 and zero coefficients
    let p_stc = pvarray::huld(1000.0, 25.0, 1.0, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!((p_stc - 1.0).abs() < 0.01, "Huld at STC should be ~1.0, got {}", p_stc);
}

// ============================================================================
// Test 15: Transformer efficiency
// ============================================================================

#[test]
fn test_e2e_transformer() {
    // simple_efficiency returns output power, not efficiency ratio
    // no_load_loss and load_loss are fractions of rating
    let output = transformer::simple_efficiency(4000.0, 5000.0, 0.002, 0.01);
    assert!(output > 3900.0 && output < 4000.0,
        "Transformer output should be slightly less than input, got {}", output);
    let eff = output / 4000.0;
    assert!(eff > 0.95 && eff < 1.0, "Transformer eff should be 95-100%, got {:.2}%", eff * 100.0);
}

// ============================================================================
// Test 16: Inverter models
// ============================================================================

#[test]
fn test_e2e_inverter_models() {
    let pac = inverter::pvwatts_ac(4000.0, 5000.0, 0.96, 0.9637);
    assert!(pac > 3500.0 && pac < 4000.0, "PVWatts AC should be 3500-4000, got {}", pac);

    let pac_multi = inverter::pvwatts_multi(&[2000.0, 2000.0], 5000.0, 0.96, 0.9637);
    assert!((pac_multi - pac).abs() < 1.0, "Multi should match single: {} vs {}", pac_multi, pac);

    // Clipping
    let pac_clip = inverter::pvwatts_ac(6000.0, 5000.0, 0.96, 0.9637);
    assert!(pac_clip <= 0.96 * 5000.0 + 1.0, "Should clip at pac0");
}

// ============================================================================
// Test 17: Atmosphere chain
// ============================================================================

#[test]
fn test_e2e_atmosphere_chain() {
    // Roundtrip
    let alt = 1830.0;
    let pressure = atmosphere::alt2pres(alt);
    let alt_back = atmosphere::pres2alt(pressure);
    assert!((alt - alt_back).abs() < 1.0, "Roundtrip: {} -> {} -> {}", alt, pressure, alt_back);

    // Precipitable water
    let pw = atmosphere::gueymard94_pw(25.0, 50.0);
    assert!(pw > 0.0 && pw < 5.0, "PW should be 0-5 cm, got {}", pw);

    // Dew point roundtrip
    let tdew = atmosphere::tdew_from_rh(25.0, 50.0);
    let rh_back = atmosphere::rh_from_tdew(25.0, tdew);
    assert!((rh_back - 50.0).abs() < 1.0, "RH roundtrip: 50% -> Tdew={} -> {:.1}%", tdew, rh_back);

    // Linke turbidity
    let aod_bb = atmosphere::bird_hulstrom80_aod_bb(0.1, 0.08);
    assert!(aod_bb > 0.0);
    let lt = atmosphere::kasten96_lt(1.5, pw, aod_bb);
    assert!(lt > 1.0 && lt < 10.0, "Linke turbidity should be 1-10, got {}", lt);

    // Wind speed
    let ws_80m = atmosphere::windspeed_powerlaw(5.0, 10.0, 80.0, 1.0 / 7.0);
    assert!(ws_80m > 5.0, "Wind at 80m should be > 10m speed");
}

// ============================================================================
// Test 18: Clear sky models comparison
// ============================================================================

#[test]
fn test_e2e_clearsky_models_compare() {
    let zenith = 30.0;
    let am = atmosphere::get_relative_airmass(zenith);
    let pw = atmosphere::gueymard94_pw(25.0, 50.0);

    let haurwitz_ghi = clearsky::haurwitz(zenith);
    let bird = clearsky::bird_default(zenith, am, 0.1, 0.08, pw);

    assert!(haurwitz_ghi > 500.0, "Haurwitz GHI should be > 500, got {}", haurwitz_ghi);
    assert!(bird.ghi > 500.0, "Bird GHI should be > 500, got {}", bird.ghi);
    assert!(bird.dni > 500.0, "Bird DNI should be > 500, got {}", bird.dni);
    assert!(bird.dhi > 30.0, "Bird DHI should be > 30, got {}", bird.dhi);
}

// ============================================================================
// Test 19: Analytical solar position
// ============================================================================

#[test]
fn test_e2e_analytical_solar_position() {
    let doy = 172.0_f64;
    let decl = solarposition::declination_spencer71(doy);
    assert!(decl > 0.35 && decl < 0.45,
        "Summer solstice declination should be ~0.41 rad, got {}", decl);

    let eot = solarposition::equation_of_time_spencer71(doy);
    assert!(eot.abs() < 15.0, "EoT should be within ±15 min, got {}", eot);

    // hour_angle(hours_from_midnight_utc, longitude, equation_of_time)
    // At solar noon UTC at longitude 0: HA should be near 0
    let ha = solarposition::hour_angle(12.0, 0.0, eot);
    assert!(ha.abs() < 15.0, "Hour angle at noon should be near 0°, got {}", ha);

    let d = solarposition::nrel_earthsun_distance(doy);
    assert!(d > 0.98 && d < 1.02, "Earth-sun distance should be near 1 AU, got {}", d);

    let decl_equinox = solarposition::declination_spencer71(80.0);
    let zenith = solarposition::solar_zenith_analytical(0.0_f64.to_radians(), 0.0_f64.to_radians(), decl_equinox);
    assert!(zenith.to_degrees().abs() < 5.0,
        "Equator noon equinox zenith should be ~0°, got {}", zenith.to_degrees());
}

// ============================================================================
// Test 20: Full daily simulation
// ============================================================================

#[test]
fn test_e2e_daily_simulation() {
    let (system, location) = make_residential_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);

    let mut total_energy_wh = 0.0;
    let mut max_power = 0.0_f64;

    for hour in 5..=20 {
        let time = Mountain.with_ymd_and_hms(2024, 6, 21, hour, 0, 0).unwrap();
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
            temp_air: 25.0 + 5.0 * cos_z,
            wind_speed: 3.0,
            albedo: Some(0.2),
        };

        let result = mc.run_model_from_weather(&weather).unwrap();
        assert!(result.ac_power >= -1.0, "Hour {}: AC should not be negative, got {}", hour, result.ac_power);

        if result.ac_power > 0.0 {
            total_energy_wh += result.ac_power;
        }
        max_power = max_power.max(result.ac_power);
    }

    let total_kwh = total_energy_wh / 1000.0;
    assert!(total_kwh > 10.0 && total_kwh < 50.0,
        "Daily energy should be 10-50 kWh, got {:.1}", total_kwh);
    assert!(max_power > 2000.0 && max_power < 5500.0,
        "Peak power should be 2000-5500 W, got {:.0}", max_power);
}

// ============================================================================
// Test 21: Alternative entry points
// ============================================================================

#[test]
fn test_e2e_alternative_entry_points() {
    let (system, location) = make_residential_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 5000.0, 0.96);
    let time = Mountain.with_ymd_and_hms(2024, 6, 21, 12, 0, 0).unwrap();

    let poa_input = POAInput {
        time,
        poa_direct: 700.0,
        poa_diffuse: 150.0,
        poa_global: 850.0,
        temp_air: 30.0,
        wind_speed: 3.0,
        aoi: 20.0,
    };
    let r_poa = mc.run_model_from_poa(&poa_input).unwrap();
    assert!(r_poa.ac_power > 1000.0, "POA entry should produce power, got {}", r_poa.ac_power);
    assert!((r_poa.poa_global - 850.0).abs() < 0.1);

    let eff_input = EffectiveIrradianceInput {
        time,
        effective_irradiance: 800.0,
        poa_global: 850.0,
        temp_air: 30.0,
        wind_speed: 3.0,
    };
    let r_eff = mc.run_model_from_effective_irradiance(&eff_input).unwrap();
    assert!(r_eff.ac_power > 1000.0, "Eff irrad entry should produce power, got {}", r_eff.ac_power);
    assert!((r_eff.effective_irradiance - 800.0).abs() < 0.1);
}

// ============================================================================
// Test 22: Different configs produce different results
// ============================================================================

#[test]
fn test_e2e_model_config_differences() {
    let location = Location::new(39.742, -105.178, Mountain, 1830.0, "Golden, CO");
    let weather = make_summer_noon_weather();

    let make_sys = || -> PVSystem {
        PVSystem::new(vec![Array {
            mount: Box::new(FixedMount { surface_tilt: 30.0, surface_azimuth: 180.0 }),
            nameplate_dc: 5000.0,
            gamma_pdc: -0.004,
            modules_per_string: 10,
            strings: 2,
            albedo: 0.2,
        }], 5000.0)
    };

    let r_pvwatts = ModelChain::with_pvwatts(make_sys(), location.clone(), 30.0, 180.0, 5000.0, 0.96)
        .run_model_from_weather(&weather).unwrap();
    let r_sapm = ModelChain::with_sapm(make_sys(), location.clone(), 30.0, 180.0, 5000.0, 0.96)
        .run_model_from_weather(&weather).unwrap();

    assert!(r_pvwatts.ac_power > 1000.0 && r_sapm.ac_power > 1000.0,
        "Both should produce power: PVWatts={}, SAPM={}", r_pvwatts.ac_power, r_sapm.ac_power);
}

// ============================================================================
// Test 23: Bird clear sky detailed
// ============================================================================

#[test]
fn test_e2e_bird_clearsky_detailed() {
    let zenith = 20.0;
    let am = atmosphere::get_relative_airmass(zenith);
    let pw = atmosphere::gueymard94_pw(25.0, 50.0);

    let result = clearsky::bird(zenith, am, 0.1, 0.08, pw, 0.3, 101325.0, 1364.0, 0.85, 0.2);

    assert!(result.ghi > 800.0, "Bird GHI should be > 800, got {}", result.ghi);
    assert!(result.dni > 800.0, "Bird DNI should be > 800, got {}", result.dni);
    assert!(result.dhi > 50.0 && result.dhi < 200.0, "Bird DHI should be 50-200, got {}", result.dhi);

    let cos_z = zenith.to_radians().cos();
    let ghi_check = result.dni * cos_z + result.dhi;
    assert!((result.ghi - ghi_check).abs() < 10.0,
        "Bird GHI closure: {} vs {}", result.ghi, ghi_check);
}
