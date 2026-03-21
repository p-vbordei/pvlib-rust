use chrono::TimeZone;
use chrono_tz::US::Eastern;
use pvlib::location::Location;
use pvlib::pvsystem::{PVSystem, Array, FixedMount};
use pvlib::modelchain::{
    ModelChain, WeatherInput, POAInput, EffectiveIrradianceInput,
};

/// Helper: create a standard test system (5 kW, Tucson-like location).
fn make_test_system() -> (PVSystem, Location) {
    let mount = Box::new(FixedMount { surface_tilt: 30.0, surface_azimuth: 180.0 });
    let array = Array {
        mount,
        nameplate_dc: 5000.0,
        gamma_pdc: -0.004,
        modules_per_string: 10,
        strings: 1,
        albedo: 0.25,
    };
    let system = PVSystem::new(vec![array], 5000.0);
    let location = Location::new(32.2, -110.9, Eastern, 700.0, "Tucson");
    (system, location)
}

/// Helper: create a noon timestamp on a summer day.
fn summer_noon() -> chrono::DateTime<chrono_tz::Tz> {
    Eastern.with_ymd_and_hms(2020, 6, 15, 12, 0, 0).unwrap()
}

// ---------------------------------------------------------------------------
// with_pvwatts factory + run_model_from_weather
// ---------------------------------------------------------------------------

#[test]
fn test_pvwatts_factory_run_model_from_weather() {
    let (system, location) = make_test_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 4800.0, 0.96);

    let weather = WeatherInput {
        time: summer_noon(),
        ghi: Some(800.0),
        dni: Some(600.0),
        dhi: Some(200.0),
        temp_air: 30.0,
        wind_speed: 2.0,
        albedo: Some(0.25),
    };

    let result = mc.run_model_from_weather(&weather).expect("run_model_from_weather failed");

    assert!(result.solar_zenith >= 0.0 && result.solar_zenith < 90.0,
        "Zenith should be daytime, got {}", result.solar_zenith);
    assert!(result.airmass > 0.0, "Airmass should be positive");
    assert!(result.poa_global > 0.0, "POA global should be positive");
    assert!(result.aoi_modifier > 0.0 && result.aoi_modifier <= 1.0,
        "AOI modifier should be in (0, 1], got {}", result.aoi_modifier);
    assert!(result.effective_irradiance > 0.0, "Effective irradiance should be positive");
    assert!(result.cell_temperature > 30.0, "Cell temp should exceed ambient");
    assert!(result.dc_power > 0.0, "DC power should be positive");
    assert!(result.ac_power > 0.0, "AC power should be positive");
    assert!(result.ac_power <= result.dc_power,
        "AC power should not exceed DC power");
}

#[test]
fn test_pvwatts_nighttime() {
    let (system, location) = make_test_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 4800.0, 0.96);

    let weather = WeatherInput {
        time: Eastern.with_ymd_and_hms(2020, 6, 15, 0, 0, 0).unwrap(),
        ghi: Some(0.0),
        dni: Some(0.0),
        dhi: Some(0.0),
        temp_air: 20.0,
        wind_speed: 1.0,
        albedo: None,
    };

    let result = mc.run_model_from_weather(&weather).expect("nighttime run failed");
    assert_eq!(result.dc_power, 0.0, "DC power at night should be 0");
    assert_eq!(result.ac_power, 0.0, "AC power at night should be 0");
}

// ---------------------------------------------------------------------------
// with_sapm factory + run_model_from_weather
// ---------------------------------------------------------------------------

#[test]
fn test_sapm_factory_run_model_from_weather() {
    let (system, location) = make_test_system();
    let mc = ModelChain::with_sapm(system, location, 30.0, 180.0, 4800.0, 0.96);

    let weather = WeatherInput {
        time: summer_noon(),
        ghi: Some(900.0),
        dni: Some(700.0),
        dhi: Some(150.0),
        temp_air: 25.0,
        wind_speed: 3.0,
        albedo: Some(0.2),
    };

    let result = mc.run_model_from_weather(&weather).expect("SAPM run failed");

    assert!(result.poa_global > 0.0);
    assert!(result.dc_power > 0.0, "SAPM DC power should be positive");
    assert!(result.ac_power > 0.0, "SAPM AC power should be positive");
    // SAPM temperature model should give different cell temp than PVWatts
    // (both should be above ambient though)
    assert!(result.cell_temperature > 25.0, "Cell temp above ambient");
}

// ---------------------------------------------------------------------------
// run_model_from_poa
// ---------------------------------------------------------------------------

#[test]
fn test_run_model_from_poa() {
    let (system, location) = make_test_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 4800.0, 0.96);

    let poa_input = POAInput {
        time: summer_noon(),
        poa_direct: 700.0,
        poa_diffuse: 200.0,
        poa_global: 900.0,
        temp_air: 28.0,
        wind_speed: 2.0,
        aoi: 20.0,
    };

    let result = mc.run_model_from_poa(&poa_input).expect("run_model_from_poa failed");

    assert!(result.poa_global > 0.0);
    assert!(result.effective_irradiance > 0.0);
    assert!(result.dc_power > 0.0);
    assert!(result.ac_power > 0.0);
    // AOI modifier for Physical model at 20 degrees should be close to 1
    assert!(result.aoi_modifier > 0.95, "IAM at 20 deg AOI should be high");
}

// ---------------------------------------------------------------------------
// run_model_from_effective_irradiance
// ---------------------------------------------------------------------------

#[test]
fn test_run_model_from_effective_irradiance() {
    let (system, location) = make_test_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 4800.0, 0.96);

    let eff_input = EffectiveIrradianceInput {
        time: summer_noon(),
        effective_irradiance: 850.0,
        poa_global: 900.0,
        temp_air: 30.0,
        wind_speed: 1.5,
    };

    let result = mc.run_model_from_effective_irradiance(&eff_input)
        .expect("run_model_from_effective_irradiance failed");

    assert!((result.effective_irradiance - 850.0).abs() < 1e-6,
        "Effective irradiance should match input");
    assert!(result.dc_power > 0.0);
    assert!(result.ac_power > 0.0);
    assert!(result.cell_temperature > 30.0);
}

// ---------------------------------------------------------------------------
// complete_irradiance
// ---------------------------------------------------------------------------

#[test]
fn test_complete_irradiance_all_provided() {
    let (system, location) = make_test_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 4800.0, 0.96);

    let weather = WeatherInput {
        time: summer_noon(),
        ghi: Some(800.0),
        dni: Some(600.0),
        dhi: Some(200.0),
        temp_air: 25.0,
        wind_speed: 2.0,
        albedo: None,
    };

    let (ghi, dni, dhi) = mc.complete_irradiance(&weather).unwrap();
    assert!((ghi - 800.0).abs() < 1e-6);
    assert!((dni - 600.0).abs() < 1e-6);
    assert!((dhi - 200.0).abs() < 1e-6);
}

#[test]
fn test_complete_irradiance_from_ghi_only() {
    let (system, location) = make_test_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 4800.0, 0.96);

    let weather = WeatherInput {
        time: summer_noon(),
        ghi: Some(800.0),
        dni: None,
        dhi: None,
        temp_air: 25.0,
        wind_speed: 2.0,
        albedo: None,
    };

    let (ghi, dni, dhi) = mc.complete_irradiance(&weather).unwrap();
    assert!((ghi - 800.0).abs() < 1e-6, "GHI should be preserved");
    assert!(dni >= 0.0, "DNI from Erbs should be non-negative");
    assert!(dhi >= 0.0, "DHI from Erbs should be non-negative");
    // Closure check: GHI ~ DNI * cos(zenith) + DHI
    // (approximate since we don't know exact zenith here)
    assert!(dni + dhi > 0.0, "DNI + DHI should be positive for positive GHI");
}

#[test]
fn test_complete_irradiance_from_dni_dhi() {
    let (system, location) = make_test_system();
    let mc = ModelChain::with_pvwatts(system, location, 30.0, 180.0, 4800.0, 0.96);

    let weather = WeatherInput {
        time: summer_noon(),
        ghi: None,
        dni: Some(600.0),
        dhi: Some(200.0),
        temp_air: 25.0,
        wind_speed: 2.0,
        albedo: None,
    };

    let (ghi, dni, dhi) = mc.complete_irradiance(&weather).unwrap();
    assert!((dni - 600.0).abs() < 1e-6);
    assert!((dhi - 200.0).abs() < 1e-6);
    assert!(ghi > 200.0, "GHI should be at least DHI");
}

// ---------------------------------------------------------------------------
// Legacy run_model still works
// ---------------------------------------------------------------------------

#[test]
fn test_legacy_run_model() {
    let (system, location) = make_test_system();
    let mc = ModelChain::new(system, location, 30.0, 180.0, 4800.0, 0.96);

    let result = mc.run_model(summer_noon(), 800.0, 600.0, 200.0, 25.0, 2.0)
        .expect("legacy run_model failed");

    assert!(result.poa_global > 0.0);
    assert!(result.dc_power > 0.0);
    assert!(result.ac_power > 0.0);
}

// ---------------------------------------------------------------------------
// Model config enum coverage
// ---------------------------------------------------------------------------

#[test]
fn test_different_transposition_models() {
    use pvlib::modelchain::{ModelChainConfig, DCModel, ACModel, AOIModel, SpectralModel, TemperatureModel, TranspositionModel, LossesModel};

    let models = [
        TranspositionModel::Isotropic,
        TranspositionModel::HayDavies,
        TranspositionModel::Klucher,
        TranspositionModel::Reindl,
        TranspositionModel::Perez,
    ];

    for transposition in &models {
        let (system, location) = make_test_system();
        let config = ModelChainConfig {
            dc_model: DCModel::PVWatts,
            ac_model: ACModel::PVWatts,
            aoi_model: AOIModel::NoLoss,
            spectral_model: SpectralModel::NoLoss,
            temperature_model: TemperatureModel::PVWatts,
            transposition_model: *transposition,
            losses_model: LossesModel::NoLoss,
        };
        let mc = ModelChain::with_config(system, location, 30.0, 180.0, 4800.0, 0.96, config);

        let weather = WeatherInput {
            time: summer_noon(),
            ghi: Some(800.0),
            dni: Some(600.0),
            dhi: Some(200.0),
            temp_air: 25.0,
            wind_speed: 2.0,
            albedo: Some(0.25),
        };

        let result = mc.run_model_from_weather(&weather)
            .expect(&format!("Failed with transposition {:?}", transposition));
        assert!(result.dc_power > 0.0,
            "DC power should be positive with {:?}", transposition);
    }
}

#[test]
fn test_different_temperature_models() {
    use pvlib::modelchain::{ModelChainConfig, DCModel, ACModel, AOIModel, SpectralModel, TemperatureModel, TranspositionModel, LossesModel};

    let models = [
        TemperatureModel::PVWatts,
        TemperatureModel::SAPM,
        TemperatureModel::PVSyst,
        TemperatureModel::Faiman,
        TemperatureModel::Fuentes,
        TemperatureModel::NOCT_SAM,
    ];

    for temp_model in &models {
        let (system, location) = make_test_system();
        let config = ModelChainConfig {
            dc_model: DCModel::PVWatts,
            ac_model: ACModel::PVWatts,
            aoi_model: AOIModel::ASHRAE,
            spectral_model: SpectralModel::NoLoss,
            temperature_model: *temp_model,
            transposition_model: TranspositionModel::Isotropic,
            losses_model: LossesModel::NoLoss,
        };
        let mc = ModelChain::with_config(system, location, 30.0, 180.0, 4800.0, 0.96, config);

        let weather = WeatherInput {
            time: summer_noon(),
            ghi: Some(800.0),
            dni: Some(600.0),
            dhi: Some(200.0),
            temp_air: 25.0,
            wind_speed: 2.0,
            albedo: Some(0.25),
        };

        let result = mc.run_model_from_weather(&weather)
            .expect(&format!("Failed with temp model {:?}", temp_model));
        assert!(result.cell_temperature > 25.0,
            "Cell temp should exceed ambient with {:?}, got {}", temp_model, result.cell_temperature);
    }
}
