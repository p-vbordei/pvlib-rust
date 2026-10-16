// Regression tests for the pvwatts_ac pdc0 bug: ModelChain::calc_ac_power and
// BatchModelChain::run used to pass the DC nameplate as `pdc0` to
// `inverter::pvwatts_ac`, so AC clipping happened at the DC capacity instead
// of the inverter capacity. The `inverter_pac0` / `inverter_capacity` fields
// were dead code. Wrong for any DC/AC ratio != 1.
//
// Reference plant: Chisineu-Cris, Romania — 54.1 MW DC / 49.57 MW AC
// (DC/AC ratio 1.0914). With eta_inv_nom = 0.96 the correct AC ceiling is
// 0.96 * 49.57e6 = 47.5872 MW; the old (buggy) ceiling was
// 0.96 * 54.1e6 = 51.936 MW.

use chrono::TimeZone;
use chrono_tz::Europe::Bucharest;
use pvlib::batch;
use pvlib::inverter;
use pvlib::location::Location;
use pvlib::modelchain::{EffectiveIrradianceInput, ModelChain};
use pvlib::pvsystem::{Array, FixedMount, PVSystem};

const DC_NAMEPLATE: f64 = 54.1e6; // W
const INVERTER_PDC0: f64 = 49.57e6; // W, plant AC rating used as pvwatts_ac pdc0
const ETA_INV_NOM: f64 = 0.96;
/// AC ceiling per the pvwatts_ac formula: pac0 = eta_inv_nom * pdc0.
const PAC0_CORRECT: f64 = ETA_INV_NOM * INVERTER_PDC0; // 47.5872 MW
const PAC0_BUGGY: f64 = ETA_INV_NOM * DC_NAMEPLATE; // 51.936 MW

fn chisineu() -> Location {
    Location::new(46.52, 21.52, Bucharest, 100.0, "Chisineu-Cris")
}

fn chisineu_system() -> PVSystem {
    let mount = Box::new(FixedMount {
        surface_tilt: 25.0,
        surface_azimuth: 180.0,
    });
    let array = Array {
        mount,
        nameplate_dc: DC_NAMEPLATE,
        gamma_pdc: -0.004,
        modules_per_string: 26,
        strings: 1,
        albedo: 0.2,
    };
    PVSystem::new(vec![array], DC_NAMEPLATE)
}

fn chisineu_modelchain(inverter_pac0: f64) -> ModelChain {
    ModelChain::new(
        chisineu_system(),
        chisineu(),
        25.0,
        180.0,
        inverter_pac0,
        ETA_INV_NOM,
    )
}

/// Strong-irradiance input that drives DC power well above the inverter
/// rating, so the AC output must sit exactly on the clipping ceiling.
fn strong_sun_input() -> EffectiveIrradianceInput {
    EffectiveIrradianceInput {
        time: Bucharest.with_ymd_and_hms(2025, 6, 21, 13, 0, 0).unwrap(),
        effective_irradiance: 1100.0,
        poa_global: 1100.0,
        temp_air: 5.0,
        wind_speed: 2.0,
    }
}

#[test]
fn test_modelchain_clips_at_inverter_ac_capacity() {
    let mc = chisineu_modelchain(INVERTER_PDC0);
    let result = mc
        .run_model_from_effective_irradiance(&strong_sun_input())
        .unwrap();

    // DC power is far above the inverter rating, so the inverter must clip.
    assert!(
        result.dc_power > INVERTER_PDC0,
        "test premise: dc_power {} must exceed inverter rating {}",
        result.dc_power,
        INVERTER_PDC0
    );
    // Clipped exactly at pac0 = eta_inv_nom * inverter_pac0 ...
    assert!(
        (result.ac_power - PAC0_CORRECT).abs() < 1.0,
        "ac_power {} should clip at eta * inverter capacity = {}",
        result.ac_power,
        PAC0_CORRECT
    );
    // ... NOT at the DC-nameplate-derived ceiling (the old bug).
    assert!(
        result.ac_power < PAC0_BUGGY - 1.0e6,
        "ac_power {} must be well below the buggy DC-nameplate ceiling {}",
        result.ac_power,
        PAC0_BUGGY
    );
}

#[test]
fn test_modelchain_falls_back_to_dc_nameplate_when_inverter_unset() {
    // inverter_pac0 <= 0.0 means "not set": pre-fix behavior is preserved,
    // pdc0 falls back to the DC nameplate total.
    let mc = chisineu_modelchain(0.0);
    let result = mc
        .run_model_from_effective_irradiance(&strong_sun_input())
        .unwrap();

    let expected = inverter::pvwatts_ac(result.dc_power, DC_NAMEPLATE, ETA_INV_NOM, 0.9637);
    assert!(
        (result.ac_power - expected).abs() < 1e-6,
        "fallback ac_power {} should equal pvwatts_ac with DC nameplate pdc0 = {}",
        result.ac_power,
        expected
    );
    // At this DC power the fallback clips at the (old) nameplate ceiling.
    assert!((result.ac_power - PAC0_BUGGY).abs() < 1.0);
}

/// One clear-sky summer-noon timestep at Chisineu with irradiance high enough
/// to drive DC power above the inverter rating.
fn strong_sun_weather() -> batch::WeatherSeries {
    batch::WeatherSeries {
        times: vec![Bucharest.with_ymd_and_hms(2025, 6, 21, 13, 0, 0).unwrap()],
        ghi: vec![1000.0],
        dni: vec![950.0],
        dhi: vec![120.0],
        temp_air: vec![0.0],
        wind_speed: vec![2.0],
        albedo: None,
    }
}

#[test]
fn test_batch_modelchain_clips_at_inverter_ac_capacity() {
    let chain = batch::BatchModelChain::pvwatts(chisineu(), 25.0, 180.0, DC_NAMEPLATE)
        .with_gamma_pdc(-0.004)
        .with_inverter(INVERTER_PDC0, ETA_INV_NOM);

    let result = chain.run(&strong_sun_weather()).unwrap();

    assert!(
        result.dc_power[0] > INVERTER_PDC0,
        "test premise: dc_power {} must exceed inverter rating {}",
        result.dc_power[0],
        INVERTER_PDC0
    );
    assert!(
        (result.ac_power[0] - PAC0_CORRECT).abs() < 1.0,
        "batch ac_power {} should clip at eta * inverter capacity = {}",
        result.ac_power[0],
        PAC0_CORRECT
    );
    assert!(
        result.ac_power[0] < PAC0_BUGGY - 1.0e6,
        "batch ac_power {} must be well below the buggy DC-nameplate ceiling {}",
        result.ac_power[0],
        PAC0_BUGGY
    );
}

#[test]
fn test_batch_default_inverter_preserves_dc_nameplate_behavior() {
    // Without with_inverter(), inverter_capacity defaults to the DC nameplate,
    // so the AC output matches pvwatts_ac with pdc0 = system_capacity_dc
    // (pre-fix behavior preserved for users who never set the inverter).
    let chain = batch::BatchModelChain::pvwatts(chisineu(), 25.0, 180.0, DC_NAMEPLATE)
        .with_gamma_pdc(-0.004);

    let result = chain.run(&strong_sun_weather()).unwrap();
    let expected = inverter::pvwatts_ac(result.dc_power[0], DC_NAMEPLATE, 0.96, 0.9637);
    assert!(
        (result.ac_power[0] - expected).abs() < 1e-6,
        "default batch ac_power {} should equal pvwatts_ac with DC nameplate pdc0 = {}",
        result.ac_power[0],
        expected
    );
}

#[test]
fn test_modelchain_and_batch_agree_on_clipped_output() {
    // Both pipelines, given DC power above the inverter rating, must clip at
    // the identical ceiling pac0 = eta_inv_nom * inverter_pdc0.
    let mc = chisineu_modelchain(INVERTER_PDC0);
    let mc_result = mc
        .run_model_from_effective_irradiance(&strong_sun_input())
        .unwrap();

    let chain = batch::BatchModelChain::pvwatts(chisineu(), 25.0, 180.0, DC_NAMEPLATE)
        .with_gamma_pdc(-0.004)
        .with_inverter(INVERTER_PDC0, ETA_INV_NOM);
    let batch_result = chain.run(&strong_sun_weather()).unwrap();

    assert!(mc_result.dc_power > INVERTER_PDC0);
    assert!(batch_result.dc_power[0] > INVERTER_PDC0);
    assert!(
        (mc_result.ac_power - batch_result.ac_power[0]).abs() < 1e-6,
        "ModelChain ({}) and BatchModelChain ({}) must clip at the same ceiling",
        mc_result.ac_power,
        batch_result.ac_power[0]
    );
}
