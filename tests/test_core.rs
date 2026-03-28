use chrono_tz::US::Pacific;
use pvlib::location::Location;
use pvlib::atmosphere::{get_relative_airmass, alt2pres, get_absolute_airmass};
use pvlib::clearsky::haurwitz;
use pvlib::irradiance::{aoi, isotropic, erbs, get_extra_radiation,
    aoi_projection, beam_component, get_ground_diffuse, poa_components,
    get_total_irradiance, DiffuseModel, disc, erbs_driesse, king, dirindex,
    perez, perez_driesse};
use pvlib::temperature::sapm_cell_temperature;
use pvlib::singlediode::{i_from_v, v_from_i};
use pvlib::iam::ashrae;
use pvlib::pvsystem::{calcparams_desoto, calcparams_cec, calcparams_pvsyst,
    sapm as sapm_model, sapm_effective_irradiance, SAPMParams};
use pvlib::inverter::pvwatts_ac;

#[test]
fn test_location_creation() {
    let loc = Location::new(32.2, -110.9, Pacific, 700.0, "Tucson");
    assert_eq!(loc.name, "Tucson");
    assert_eq!(loc.latitude, 32.2);
}

#[test]
fn test_airmass() {
    let am_rel = get_relative_airmass(0.0); // Zenith = 0
    assert!((am_rel - 1.0).abs() < 0.05, "AM relative at zenith should be ~ 1.0");

    let press = alt2pres(0.0);
    assert!((press - 101325.0).abs() < 1.0, "Sea level pressure should be 101325 Pa");

    let am_abs = get_absolute_airmass(am_rel, press);
    assert!((am_abs - am_rel).abs() < 0.01, "AM absolute should equal AM relative at sea level");
}

#[test]
fn test_clearsky_haurwitz() {
    let ghi = haurwitz(0.0); // Zenith 0
    assert!(ghi > 1000.0 && ghi < 1100.0, "GHI at zenith should be high");

    let ghi_down = haurwitz(90.0);
    assert_eq!(ghi_down, 0.0, "GHI when down should be 0");
}

#[test]
fn test_irradiance_aoi() {
    // horizontal surface, sun at zenith = 0 -> aoi = 0
    let aoi_val = aoi(0.0, 180.0, 0.0, 180.0);
    assert!((aoi_val - 0.0).abs() < 1e-6);

    // surface tilted 30 deg south, sun at zenith 30 deg south (azimuth 180) => aoi = 0
    let aoi_val2 = aoi(30.0, 180.0, 30.0, 180.0);
    assert!((aoi_val2 - 0.0).abs() < 1e-4);
}

#[test]
fn test_temperature_sapm() {
    let (t_cell, t_mod) = sapm_cell_temperature(1000.0, 20.0, 1.0, -3.56, -0.075, 3.0, 1000.0);
    assert!(t_mod > 20.0, "Module temp should be higher than ambient");
    assert!(t_cell > t_mod, "Cell temp should be higher than module temp");
}

#[test]
fn test_iam_ashrae() {
    let iam_0 = ashrae(0.0, 0.05);
    assert!((iam_0 - 1.0).abs() < 1e-6);

    let iam_60 = ashrae(60.0, 0.05);
    // 1 - 0.05 * (1/cos(60) - 1) = 1 - 0.05 * (2 - 1) = 0.95
    assert!((iam_60 - 0.95).abs() < 1e-6);
}

#[test]
fn test_pvsystem_and_inverter() {
    use pvlib::pvsystem::{PVSystem, Array, FixedMount};
    let mount = Box::new(FixedMount { surface_tilt: 20.0, surface_azimuth: 180.0 });
    let array = Array { 
        mount, nameplate_dc: 5000.0, gamma_pdc: -0.004, 
        modules_per_string: 10, strings: 1, albedo: 0.2 
    };
    let system = PVSystem::new(vec![array], 5000.0);
    let dc_power = system.get_dc_power_total(1000.0, 25.0);
    assert!((dc_power - 5000.0).abs() < 1e-4, "DC power at STC should match nameplate");

    let ac_power = pvwatts_ac(dc_power, 5000.0, 0.96, 0.9637);
    assert!((ac_power - 4800.0).abs() < 1.0, "AC power should be capped at inverter rating");
}

#[test]
fn test_irradiance_erbs() {
    let (dni, dhi) = erbs(1000.0, 30.0, 1, 1361.0);
    assert!(dni > 0.0 && dhi > 0.0, "DNI and DHI should be positive for valid GHI");
    assert!((dni * 30.0_f64.to_radians().cos() + dhi - 1000.0).abs() < 10.0, "Closure GHI = DNI*cosZ + DHI");
}

#[test]
fn test_singlediode() {
    // Basic test: Short circuit current (V=0) should be approx light current I_L
    let i_sc = i_from_v(0.0, 5.0, 1e-9, 0.1, 1000.0, 1.5);
    assert!((i_sc - 5.0).abs() < 0.1, "I_sc should be very close to light current");
    
    // Open circuit voltage (I=0)
    let v_oc = v_from_i(0.0, 5.0, 1e-9, 0.1, 1000.0, 1.5);
    assert!(v_oc > 20.0, "V_oc should be a reasonable positive voltage");
}

#[test]
fn test_aoi_projection() {
    // Horizontal surface, sun at zenith -> cos(0) = 1
    let proj = aoi_projection(0.0, 180.0, 0.0, 180.0);
    assert!((proj - 1.0).abs() < 1e-6);

    // Surface tilted 30S, sun at zenith 30 from south -> cos(0) = 1
    let proj2 = aoi_projection(30.0, 180.0, 30.0, 180.0);
    assert!((proj2 - 1.0).abs() < 1e-4);

    // Steep tilt, sun behind -> negative projection
    // Surface tilted 80 toward south (az=180), sun at zenith 80 from north (az=0)
    let proj3 = aoi_projection(80.0, 180.0, 80.0, 0.0);
    assert!(proj3 < 0.0, "Sun behind steeply tilted surface should give negative projection");
}

#[test]
fn test_beam_component() {
    let beam = beam_component(0.0, 180.0, 30.0, 180.0, 1000.0);
    // cos(30 deg) * 1000 ~ 866
    assert!((beam - 866.025).abs() < 1.0);

    // Sun behind steeply tilted surface -> beam should be 0
    let beam2 = beam_component(80.0, 180.0, 80.0, 0.0, 1000.0);
    assert_eq!(beam2, 0.0);
}

#[test]
fn test_get_ground_diffuse() {
    // Horizontal surface -> (1 - cos(0))/2 = 0
    let gd = get_ground_diffuse(0.0, 1000.0, 0.25);
    assert!((gd - 0.0).abs() < 1e-6);

    // Vertical surface -> (1 - cos(90))/2 = 0.5
    let gd2 = get_ground_diffuse(90.0, 1000.0, 0.25);
    assert!((gd2 - 125.0).abs() < 1.0);

    // 30 degree tilt
    let gd3 = get_ground_diffuse(30.0, 1000.0, 0.25);
    let expected = 1000.0 * 0.25 * (1.0 - 30.0_f64.to_radians().cos()) * 0.5;
    assert!((gd3 - expected).abs() < 1e-6);
}

#[test]
fn test_poa_components() {
    let result = poa_components(30.0, 1000.0, 100.0, 50.0);
    let expected_direct = 1000.0 * 30.0_f64.to_radians().cos();
    assert!((result.poa_direct - expected_direct).abs() < 1e-6);
    assert!((result.poa_diffuse - 150.0).abs() < 1e-6);
    assert!((result.poa_global - (expected_direct + 150.0)).abs() < 1e-6);
    assert_eq!(result.poa_sky_diffuse, 100.0);
    assert_eq!(result.poa_ground_diffuse, 50.0);

    // AOI > 90 -> direct = 0
    let result2 = poa_components(100.0, 1000.0, 100.0, 50.0);
    assert_eq!(result2.poa_direct, 0.0);
}

#[test]
fn test_get_total_irradiance_isotropic() {
    let result = get_total_irradiance(
        30.0, 180.0, 30.0, 180.0,
        800.0, 1000.0, 200.0,
        0.25, DiffuseModel::Isotropic,
        None, None,
    );
    assert!(result.poa_global > 0.0);
    assert!(result.poa_direct > 0.0);
    assert!(result.poa_diffuse > 0.0);
    // Direct should be close to DNI since AOI ~ 0
    assert!((result.poa_direct - 800.0).abs() < 1.0);
}

#[test]
fn test_disc() {
    let result = disc(1000.0, 30.0, 180, Some(101325.0));
    assert!(result.dni > 0.0, "DNI should be positive");
    assert!(result.kt > 0.0 && result.kt <= 1.0, "Kt should be in (0, 1]");
    assert!(result.airmass > 0.0, "Airmass should be positive");

    // Night: zenith > 87
    let night = disc(0.0, 90.0, 180, Some(101325.0));
    assert_eq!(night.dni, 0.0);
}

#[test]
fn test_erbs_driesse() {
    let result = erbs_driesse(1000.0, 30.0, 180);
    assert!(result.dni > 0.0, "DNI should be positive");
    assert!(result.dhi > 0.0, "DHI should be positive");
    assert!(result.kt > 0.0 && result.kt <= 1.0, "Kt should be in (0, 1]");
    // Check closure: GHI ~ DNI * cos(z) + DHI
    let closure = result.dni * 30.0_f64.to_radians().cos() + result.dhi;
    assert!((closure - 1000.0).abs() < 10.0, "Closure should hold: GHI ~ DNI*cosZ + DHI");

    // Night
    let night = erbs_driesse(0.0, 90.0, 180);
    assert_eq!(night.dni, 0.0);
}

#[test]
fn test_calcparams_desoto() {
    // Typical c-Si module parameters at STC
    let params = calcparams_desoto(
        1000.0, 25.0,       // effective_irradiance, temp_cell (STC)
        0.003,              // alpha_sc
        1.6,                // a_ref
        6.0,                // I_L_ref
        1e-10,              // I_o_ref
        300.0,              // R_sh_ref
        0.5,                // R_s
        1.121,              // EgRef (c-Si)
        -0.0002677,         // dEgdT
    );
    // At STC, photocurrent ~ I_L_ref
    assert!((params.photocurrent - 6.0).abs() < 1e-6, "IL at STC should equal I_L_ref");
    // At STC, Rsh = R_sh_ref
    assert!((params.resistance_shunt - 300.0).abs() < 1e-6, "Rsh at STC should equal R_sh_ref");
    // At STC, Rs = R_s
    assert!((params.resistance_series - 0.5).abs() < 1e-10);
    // At STC, nNsVth = a_ref (since Tcell = Tref)
    assert!((params.n_ns_vth - 1.6).abs() < 1e-6, "nNsVth at STC should equal a_ref");
    // I0 at STC should equal I_o_ref
    assert!((params.saturation_current - 1e-10).abs() < 1e-15, "I0 at STC should equal I_o_ref");
}

#[test]
fn test_calcparams_desoto_temperature() {
    // Test at higher temperature: 50C
    let params = calcparams_desoto(
        1000.0, 50.0,
        0.003, 1.6, 6.0, 1e-10, 300.0, 0.5, 1.121, -0.0002677,
    );
    // nNsVth scales with temperature
    let expected_nns = 1.6 * (50.0 + 273.15) / (25.0 + 273.15);
    assert!((params.n_ns_vth - expected_nns).abs() < 1e-6);
    // IL increases slightly with temperature
    let expected_il = 6.0 + 0.003 * 25.0;
    assert!((params.photocurrent - expected_il).abs() < 1e-6);
    // I0 increases with temperature
    assert!(params.saturation_current > 1e-10, "I0 should increase with temperature");
}

#[test]
fn test_calcparams_desoto_low_irradiance() {
    let params = calcparams_desoto(
        200.0, 25.0,
        0.003, 1.6, 6.0, 1e-10, 300.0, 0.5, 1.121, -0.0002677,
    );
    // IL scales with irradiance
    assert!((params.photocurrent - 6.0 * 200.0 / 1000.0).abs() < 1e-6);
    // Rsh = R_sh_ref * 1000/200 = 1500
    assert!((params.resistance_shunt - 1500.0).abs() < 1e-6);
}

#[test]
fn test_calcparams_cec() {
    // CEC adjusts alpha_sc by (1 - Adjust/100)
    let adjust = 10.0; // 10% reduction
    let params_cec = calcparams_cec(
        1000.0, 50.0,
        0.003, 1.6, 6.0, 1e-10, 300.0, 0.5, adjust, 1.121, -0.0002677,
    );
    let params_desoto = calcparams_desoto(
        1000.0, 50.0,
        0.003 * 0.9, 1.6, 6.0, 1e-10, 300.0, 0.5, 1.121, -0.0002677,
    );
    assert!((params_cec.photocurrent - params_desoto.photocurrent).abs() < 1e-10);
    assert!((params_cec.saturation_current - params_desoto.saturation_current).abs() < 1e-20);
}

#[test]
fn test_calcparams_pvsyst() {
    let params = calcparams_pvsyst(
        1000.0, 25.0,
        0.003,          // alpha_sc
        1.1,            // gamma_ref
        0.0005,         // mu_gamma
        6.0,            // I_L_ref
        1e-10,          // I_o_ref
        300.0,          // R_sh_ref
        3000.0,         // R_sh_0
        0.5,            // R_s
        60,             // cells_in_series
        1.121,          // EgRef
    );
    // At STC, photocurrent ~ I_L_ref
    assert!((params.photocurrent - 6.0).abs() < 1e-6);
    assert!(params.resistance_shunt > 0.0);
    assert!(params.n_ns_vth > 0.0);
    assert!(params.saturation_current > 0.0);
}

#[test]
fn test_sapm() {
    let module = SAPMParams {
        isco: 5.0, impo: 4.5, voco: 40.0, vmpo: 32.0,
        aisc: 0.0005, aimp: -0.0003,
        bvoco: -0.12, mbvoc: 0.0, bvmpo: -0.14, mbvmp: 0.0,
        n: 1.2, cells_in_series: 60,
        c0: 1.0, c1: 0.0, c2: 1.0, c3: 0.0,
        a0: 0.92, a1: 0.058, a2: -0.0088, a3: 0.00065, a4: -0.000018,
        b0: 1.0, b1: -0.002438, b2: 3.103e-4, b3: -1.246e-5,
        b4: 2.112e-7, b5: -1.359e-9,
        fd: 1.0,
    };

    let result = sapm_model(1000.0, 25.0, &module);
    // At STC with C0=1, C1=0: i_mp ~ impo * 1 * Ee * (1+0) = 4.5
    assert!((result.i_sc - 5.0).abs() < 0.01, "Isc at STC should be ~Isco");
    assert!((result.i_mp - 4.5).abs() < 0.01, "Imp at STC should be ~Impo");
    assert!(result.v_oc > 0.0);
    assert!(result.v_mp > 0.0);
    assert!(result.p_mp > 0.0);
    // Voc at STC: Voco + Ns*delta*ln(1) + Bvoco*0 = Voco
    assert!((result.v_oc - 40.0).abs() < 0.1, "Voc at STC should be ~Voco");
}

#[test]
fn test_sapm_effective_irradiance() {
    let module = SAPMParams {
        isco: 5.0, impo: 4.5, voco: 40.0, vmpo: 32.0,
        aisc: 0.0005, aimp: -0.0003,
        bvoco: -0.12, mbvoc: 0.0, bvmpo: -0.14, mbvmp: 0.0,
        n: 1.2, cells_in_series: 60,
        c0: 1.0, c1: 0.0, c2: 1.0, c3: 0.0,
        a0: 0.92, a1: 0.058, a2: -0.0088, a3: 0.00065, a4: -0.000018,
        b0: 1.0, b1: -0.002438, b2: 3.103e-4, b3: -1.246e-5,
        b4: 2.112e-7, b5: -1.359e-9,
        fd: 1.0,
    };

    let ee = sapm_effective_irradiance(800.0, 200.0, 1.5, 20.0, &module);
    assert!(ee > 0.0, "Effective irradiance should be positive");
    // Should be less than poa_direct + poa_diffuse due to losses
    assert!(ee < 1000.0 + 50.0, "Effective irradiance should be reasonable");
}

#[test]
fn test_king_diffuse() {
    // Horizontal surface: king should equal isotropic (second term vanishes)
    let king_horiz = king(0.0, 200.0, 1000.0, 30.0);
    let iso_horiz = isotropic(0.0, 200.0);
    assert!((king_horiz - iso_horiz).abs() < 1e-6,
        "King on horizontal should match isotropic");

    // Tilted surface: result should be positive and differ from isotropic
    let king_tilt = king(30.0, 200.0, 1000.0, 30.0);
    assert!(king_tilt > 0.0, "King diffuse should be positive");

    // At high zenith, the horizon brightening term (0.012*z - 0.04) grows
    let king_high_z = king(30.0, 200.0, 1000.0, 60.0);
    assert!(king_high_z > 0.0);

    // Vertical surface gets more ground/horizon contribution
    let king_vert = king(90.0, 200.0, 1000.0, 30.0);
    assert!(king_vert > 0.0);
}

#[test]
fn test_king_nonnegative() {
    // King should never return negative
    let result = king(45.0, 10.0, 50.0, 85.0);
    assert!(result >= 0.0, "King should be non-negative");
}

#[test]
fn test_dirindex() {
    // Under clear sky conditions (ghi == ghi_clearsky), dirindex should
    // return approximately dni_clearsky
    let ghi_clear = 800.0;
    let dni_clear = 700.0;
    let zenith = 30.0;

    let dni = dirindex(ghi_clear, ghi_clear, dni_clear, zenith, 180, Some(101325.0));
    // When ghi == ghi_clearsky, dirint(ghi) == dirint(ghi_clear), so ratio = 1
    assert!((dni - dni_clear).abs() < 1.0,
        "Under clear sky, dirindex should return dni_clearsky, got {}", dni);

    // Under cloudy conditions (ghi < ghi_clear), DNI should be reduced
    let dni_cloudy = dirindex(400.0, ghi_clear, dni_clear, zenith, 180, Some(101325.0));
    assert!(dni_cloudy < dni_clear, "Cloudy conditions should reduce DNI");
    assert!(dni_cloudy >= 0.0, "DNI should be non-negative");
}

#[test]
fn test_dirindex_night() {
    let dni = dirindex(0.0, 0.0, 0.0, 95.0, 180, Some(101325.0));
    assert_eq!(dni, 0.0, "Nighttime dirindex should be 0");
}

#[test]
fn test_perez_driesse_basic() {
    // Perez-Driesse should produce positive diffuse for a tilted surface
    let dni_extra = get_extra_radiation(180);
    let result = perez_driesse(
        30.0, 180.0, 200.0, 800.0, dni_extra, 30.0, 180.0, None,
    );
    assert!(result > 0.0, "Perez-Driesse diffuse should be positive, got {}", result);
}

#[test]
fn test_perez_driesse_horizontal() {
    // On a horizontal surface, diffuse should approximate DHI (isotropic-like)
    let dni_extra = get_extra_radiation(180);
    let result = perez_driesse(
        0.0, 180.0, 200.0, 800.0, dni_extra, 30.0, 180.0, None,
    );
    // For horizontal surface: cos(tilt)=1, sin(tilt)=0
    // term1 = (1-F1)*1, term2 = F1*cos(0)/cos(z), term3 = 0
    // result should be close to DHI
    assert!(result > 0.0);
    assert!((result - 200.0).abs() < 50.0,
        "Perez-Driesse on horizontal should be near DHI, got {}", result);
}

#[test]
fn test_perez_driesse_vs_perez_agreement() {
    // Perez-Driesse should produce results close to the original Perez model
    let dni_extra = get_extra_radiation(180);
    let surface_tilt = 30.0;
    let surface_azimuth = 180.0;
    let solar_zenith = 40.0;
    let solar_azimuth = 180.0;
    let dhi = 150.0;
    let dni = 700.0;
    let am = pvlib::atmosphere::get_relative_airmass(solar_zenith);
    let aoi_val = aoi(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth);

    let perez_result = perez(
        surface_tilt, surface_azimuth, dhi, dni, dni_extra,
        solar_zenith, solar_azimuth, am, aoi_val,
    );
    let pd_result = perez_driesse(
        surface_tilt, surface_azimuth, dhi, dni, dni_extra,
        solar_zenith, solar_azimuth, Some(am),
    );

    // They should agree within ~20% for typical conditions
    let diff_pct = if perez_result > 0.0 {
        ((pd_result - perez_result) / perez_result).abs() * 100.0
    } else {
        0.0
    };
    assert!(diff_pct < 25.0,
        "Perez-Driesse and Perez should be similar: PD={}, P={}, diff={}%",
        pd_result, perez_result, diff_pct);
}

#[test]
fn test_perez_driesse_nonnegative() {
    let dni_extra = get_extra_radiation(180);
    // Edge case: very low irradiance
    let result = perez_driesse(45.0, 180.0, 5.0, 10.0, dni_extra, 85.0, 180.0, None);
    assert!(result >= 0.0, "Perez-Driesse should be non-negative");

    // Night case
    let night = perez_driesse(45.0, 180.0, 0.0, 0.0, dni_extra, 95.0, 180.0, None);
    assert!(night >= 0.0, "Perez-Driesse at night should be non-negative");
}
