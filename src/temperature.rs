/// Sandia Array Performance Model (SAPM) for cell and module temperature.
/// 
/// # Arguments
/// * `poa_global` - Total incident irradiance on the module in W/m^2.
/// * `temp_air` - Ambient air temperature in Celsius.
/// * `wind_speed` - Wind speed in m/s at standard 10m height.
/// * `a` - SAPM parameter `a` (empirical).
/// * `b` - SAPM parameter `b` (empirical).
/// * `delta_t` - SAPM parameter `delta_t` (temperature difference between cell and module back surface).
/// * `irrad_ref` - Reference irradiance, typically 1000.0 W/m^2.
/// 
/// # Returns
/// A tuple of `(temp_cell, temp_module)` in Celsius.
#[inline]
pub fn sapm_cell_temperature(poa_global: f64, temp_air: f64, wind_speed: f64, a: f64, b: f64, delta_t: f64, irrad_ref: f64) -> (f64, f64) {
    if poa_global <= 0.0 {
        return (temp_air, temp_air);
    }
    
    let temp_module = temp_air + poa_global * (a + b * wind_speed).exp();
    let temp_cell = temp_module + (poa_global / irrad_ref) * delta_t;
    
    (temp_cell, temp_module)
}

/// PVsyst cell temperature model.
/// 
/// # Arguments
/// * `poa_global` - Total incident irradiance on the module in W/m^2.
/// * `temp_air` - Ambient air temperature in Celsius.
/// * `wind_speed` - Wind speed in m/s.
/// * `u_c` - Constant heat transfer component (W/(m^2 K)).
/// * `u_v` - Convective heat transfer component (W/(m^3 s K)).
/// * `module_efficiency` - Module efficiency as a decimal (e.g. 0.15 for 15%).
/// * `alpha_absorption` - Absorption coefficient (e.g., 0.9).
/// 
/// # Returns
/// Cell temperature in Celsius.
#[inline]
pub fn pvsyst_cell_temperature(poa_global: f64, temp_air: f64, wind_speed: f64, u_c: f64, u_v: f64, module_efficiency: f64, alpha_absorption: f64) -> f64 {
    if poa_global <= 0.0 {
        return temp_air;
    }
    let h_total = (u_c + u_v * wind_speed).max(0.01);
    temp_air + (alpha_absorption * poa_global * (1.0 - module_efficiency)) / h_total
}

/// Faiman (2008) cell temperature model.
/// `T_cell = T_ambient + POA / (U0 + U1 * wind_speed)`
/// 
/// # References
/// Faiman, D., 2008, "Assessing the outdoor operating temperature of photovoltaic modules," 
/// Progress in Photovoltaics 16(4), pp. 307-315.
/// 
/// # Arguments
/// * `poa_global` - Total incident irradiance on the module in W/m^2.
/// * `temp_air` - Ambient air temperature in Celsius.
/// * `wind_speed` - Wind speed in m/s.
/// * `u0` - Constant heat transfer coefficient (typically 25.0).
/// * `u1` - Convective heat transfer coefficient (typically 6.84).
#[inline]
pub fn faiman(poa_global: f64, temp_air: f64, wind_speed: f64, u0: f64, u1: f64) -> f64 {
    if poa_global <= 0.0 {
        return temp_air;
    }
    let h_total = (u0 + u1 * wind_speed).max(0.01);
    temp_air + poa_global / h_total
}

/// Fuentes (1987) module temperature model.
/// 
/// Accounts for thermal capacitance and inertia of PV modules.
/// 
/// # References
/// Fuentes, M.K., 1987. "A simplified thermal model for flat-plate photovoltaic arrays".
#[inline]
pub fn fuentes(poa_global: f64, temp_air: f64, wind_speed: f64, inoct: f64) -> f64 {
    // Simplified equilibrium approximation (actual is differential eq over time)
    if poa_global <= 0.0 { return temp_air; }
    let h = 1.2 + 0.8 * wind_speed; // rough convective cooling
    temp_air + poa_global * (inoct - 20.0) / (800.0 * h)
}

/// Ross module temperature model.
///
/// # References
/// Ross, R. G., 1980. "Flat-Plate Photovoltaic Array Design Optimization".
#[inline]
pub fn ross(poa_global: f64, temp_air: f64, inoct: f64) -> f64 {
    if poa_global <= 0.0 {
        return temp_air;
    }
    let k = (inoct - 20.0) / 800.0;
    temp_air + k * poa_global
}

/// SAPM module back-surface temperature model.
///
/// T_module = poa_global * exp(a + b * wind_speed) + temp_air
///
/// # Arguments
/// * `poa_global` - Total incident irradiance [W/m^2].
/// * `temp_air` - Ambient dry bulb temperature [C].
/// * `wind_speed` - Wind speed at 10m height [m/s].
/// * `a` - SAPM parameter a (default -3.56 for glass/polymer open rack).
/// * `b` - SAPM parameter b (default -0.0750 for glass/polymer open rack).
///
/// # Returns
/// Module back-surface temperature in Celsius.
///
/// # References
/// King, D. et al, 2004, "Sandia Photovoltaic Array Performance Model", SAND Report 3535.
#[inline]
pub fn sapm_module(poa_global: f64, temp_air: f64, wind_speed: f64, a: f64, b: f64) -> f64 {
    poa_global * (a + b * wind_speed).exp() + temp_air
}

/// SAPM module temperature with default parameters (glass/polymer, open rack).
#[inline]
pub fn sapm_module_default(poa_global: f64, temp_air: f64, wind_speed: f64) -> f64 {
    sapm_module(poa_global, temp_air, wind_speed, -3.56, -0.0750)
}

/// Calculate cell temperature from module back-surface temperature using SAPM.
///
/// T_cell = module_temperature + poa_global / irrad_ref * delta_t
///
/// # Arguments
/// * `module_temperature` - Module back-surface temperature [C].
/// * `poa_global` - Total incident irradiance [W/m^2].
/// * `delta_t` - Temperature difference between cell and module back [C].
/// * `irrad_ref` - Reference irradiance, default 1000 W/m^2.
///
/// # Returns
/// Cell temperature in Celsius.
///
/// # References
/// King, D. et al, 2004, "Sandia Photovoltaic Array Performance Model", SAND Report 3535.
#[inline]
pub fn sapm_cell_from_module(module_temperature: f64, poa_global: f64, delta_t: f64, irrad_ref: f64) -> f64 {
    module_temperature + poa_global / irrad_ref * delta_t
}

/// SAPM cell from module with default parameters (delta_t=3, irrad_ref=1000).
#[inline]
pub fn sapm_cell_from_module_default(module_temperature: f64, poa_global: f64) -> f64 {
    sapm_cell_from_module(module_temperature, poa_global, 3.0, 1000.0)
}

/// Adjustment to NOCT for mounting standoff distance.
///
/// Piecewise function matching SAM implementation.
fn adj_for_mounting_standoff(x: f64) -> f64 {
    if x <= 0.0 {
        0.0
    } else if x < 0.5 {
        18.0
    } else if x < 1.5 {
        11.0
    } else if x < 2.5 {
        6.0
    } else if x <= 3.5 {
        2.0
    } else {
        0.0
    }
}

/// NOCT cell temperature model from the System Advisor Model (SAM).
///
/// Note: This function assumes effective_irradiance == poa_global, i.e., no
/// spectral or IAM losses are applied to the irradiance used for temperature
/// calculation. This matches the common use case where effective irradiance
/// is not separately tracked for the thermal model. In pvlib-python, an
/// optional `effective_irradiance` parameter allows decoupling these; this
/// implementation does not yet support that.
///
/// # Arguments
/// * `poa_global` - Total incident irradiance [W/m^2].
/// * `temp_air` - Ambient dry bulb temperature [C].
/// * `wind_speed` - Wind speed [m/s].
/// * `noct` - Nominal operating cell temperature [C].
/// * `module_efficiency` - Module efficiency at reference conditions.
/// * `transmittance_absorptance` - Combined tau*alpha coefficient (default 0.9).
/// * `array_height` - Height above ground: 1 or 2 stories (default 1).
/// * `mount_standoff` - Distance between module and mounting surface in inches (default 4.0).
///
/// # Returns
/// Cell temperature in Celsius.
///
/// # References
/// Gilman, P. et al, 2018, "SAM Photovoltaic Model Technical Reference Update", NREL/TP-6A20-67399.
#[inline]
pub fn noct_sam(
    poa_global: f64,
    temp_air: f64,
    wind_speed: f64,
    noct: f64,
    module_efficiency: f64,
    transmittance_absorptance: f64,
    array_height: u32,
    mount_standoff: f64,
) -> f64 {
    let wind_adj = match array_height {
        1 => 0.51 * wind_speed,
        2 => 0.61 * wind_speed,
        _ => 0.51 * wind_speed, // default to height 1
    };

    let noct_adj = noct + adj_for_mounting_standoff(mount_standoff);
    let tau_alpha = transmittance_absorptance;

    let cell_temp_init = poa_global / 800.0 * (noct_adj - 20.0);
    let heat_loss = 1.0 - module_efficiency / tau_alpha;
    let wind_loss = 9.5 / (5.7 + 3.8 * wind_adj);

    temp_air + cell_temp_init * heat_loss * wind_loss
}

/// NOCT SAM with default parameters (transmittance_absorptance=0.9, array_height=1, mount_standoff=4.0).
#[inline]
pub fn noct_sam_default(
    poa_global: f64,
    temp_air: f64,
    wind_speed: f64,
    noct: f64,
    module_efficiency: f64,
) -> f64 {
    noct_sam(poa_global, temp_air, wind_speed, noct, module_efficiency, 0.9, 1, 4.0)
}

/// Generic linear cell temperature model.
///
/// T_cell = temp_air + poa_global * (absorptance - module_efficiency) / (u_c + u_v * wind_speed)
///
/// # Arguments
/// * `poa_global` - Total incident irradiance [W/m^2].
/// * `temp_air` - Ambient dry bulb temperature [C].
/// * `wind_speed` - Wind speed at 10m height [m/s].
/// * `u_c` - Combined heat transfer coefficient at zero wind [(W/m^2)/C].
/// * `u_v` - Wind influence on heat transfer [(W/m^2)/C/(m/s)].
/// * `module_efficiency` - Module electrical efficiency.
/// * `absorptance` - Light absorptance of the module.
///
/// # Returns
/// Cell temperature in Celsius.
///
/// # References
/// Driesse, A. et al, 2022, "PV Module Operating Temperature Model Equivalence
/// and Parameter Translation", IEEE PVSC.
#[inline]
pub fn generic_linear(
    poa_global: f64,
    temp_air: f64,
    wind_speed: f64,
    u_c: f64,
    u_v: f64,
    module_efficiency: f64,
    absorptance: f64,
) -> f64 {
    let heat_input = poa_global * (absorptance - module_efficiency);
    let total_loss_factor = u_c + u_v * wind_speed;
    temp_air + heat_input / total_loss_factor
}

/// Generic linear model with default parameters (u_c=29.0, u_v=0.0, efficiency=0.15, absorptance=0.9).
#[inline]
pub fn generic_linear_default(poa_global: f64, temp_air: f64, wind_speed: f64) -> f64 {
    generic_linear(poa_global, temp_air, wind_speed, 29.0, 0.0, 0.15, 0.9)
}

