pub trait Mount {
    fn get_surface_tilt(&self) -> f64;
    fn get_surface_azimuth(&self) -> f64;
}

#[derive(Debug, Clone)]
pub struct FixedMount {
    pub surface_tilt: f64,
    pub surface_azimuth: f64,
}

impl Mount for FixedMount {
    fn get_surface_tilt(&self) -> f64 { self.surface_tilt }
    fn get_surface_azimuth(&self) -> f64 { self.surface_azimuth }
}

#[derive(Debug, Clone)]
pub struct SingleAxisTrackerMount {
    pub axis_tilt: f64,
    pub axis_azimuth: f64,
    pub max_angle: f64,
    pub backtrack: bool,
    pub gcr: f64,
}

impl Mount for SingleAxisTrackerMount {
    fn get_surface_tilt(&self) -> f64 { self.axis_tilt } // Dynamic in reality, simplified here
    fn get_surface_azimuth(&self) -> f64 { self.axis_azimuth }
}

/// Represents a subset of a PV system strings with identical orientation
pub struct Array {
    pub mount: Box<dyn Mount>,
    pub nameplate_dc: f64,
    pub gamma_pdc: f64,
    pub modules_per_string: u32,
    pub strings: u32,
    pub albedo: f64,
}

pub struct PVSystem {
    pub arrays: Vec<Array>,
    pub inverter_capacity: f64,
}

impl PVSystem {
    pub fn new(arrays: Vec<Array>, inverter_capacity: f64) -> Self {
        Self { arrays, inverter_capacity }
    }

    #[inline]
    pub fn get_nameplate_dc_total(&self) -> f64 {
        self.arrays.iter().map(|a| a.nameplate_dc).sum()
    }

    #[inline]
    pub fn get_dc_power_total(&self, poa_global: f64, temp_cell: f64) -> f64 {
        let mut total_power = 0.0;
        for array in &self.arrays {
             let pdc = array.nameplate_dc * (poa_global / 1000.0) * (1.0 + array.gamma_pdc * (temp_cell - 25.0));
             total_power += pdc.max(0.0);
        }
        total_power
    }
}

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant in eV/K
const KB_EV: f64 = 8.617333262e-05;
/// Boltzmann constant in J/K
const KB_J: f64 = 1.380649e-23;
/// Elementary charge in coulombs
const Q_C: f64 = 1.602176634e-19;

// ---------------------------------------------------------------------------
// Single-diode model parameters output
// ---------------------------------------------------------------------------

/// Five parameters for the single-diode equation.
#[derive(Debug, Clone, Copy)]
pub struct SDMParams {
    /// Light-generated (photo) current [A]
    pub photocurrent: f64,
    /// Diode saturation current [A]
    pub saturation_current: f64,
    /// Series resistance [ohm]
    pub resistance_series: f64,
    /// Shunt resistance [ohm]
    pub resistance_shunt: f64,
    /// Product n * Ns * Vth at operating conditions [V]
    pub n_ns_vth: f64,
}

/// Calculate five single-diode model parameters using the De Soto et al. model.
///
/// # Parameters
/// - `effective_irradiance`: Irradiance converted to photocurrent [W/m2]
/// - `temp_cell`: Average cell temperature [C]
/// - `alpha_sc`: Short-circuit current temperature coefficient [A/C]
/// - `a_ref`: n * Ns * Vth at reference conditions [V]
/// - `i_l_ref`: Photo current at reference conditions [A]
/// - `i_o_ref`: Diode saturation current at reference conditions [A]
/// - `r_sh_ref`: Shunt resistance at reference conditions [ohm]
/// - `r_s`: Series resistance at reference conditions [ohm]
/// - `eg_ref`: Bandgap energy at reference temperature [eV] (default 1.121 for c-Si)
/// - `d_eg_dt`: Temperature dependence of bandgap [1/K] (default -0.0002677)
///
/// # References
/// W. De Soto et al., 2006, "Improvement and validation of a model for
/// photovoltaic array performance", Solar Energy, vol 80, pp. 78-88.
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn calcparams_desoto(
    effective_irradiance: f64,
    temp_cell: f64,
    alpha_sc: f64,
    a_ref: f64,
    i_l_ref: f64,
    i_o_ref: f64,
    r_sh_ref: f64,
    r_s: f64,
    eg_ref: f64,
    d_eg_dt: f64,
) -> SDMParams {
    let irrad_ref = 1000.0;
    let temp_ref = 25.0;
    let tref_k = temp_ref + 273.15;
    let tcell_k = temp_cell + 273.15;

    let eg = eg_ref * (1.0 + d_eg_dt * (tcell_k - tref_k));

    let n_ns_vth = a_ref * (tcell_k / tref_k);

    let photocurrent = effective_irradiance / irrad_ref
        * (i_l_ref + alpha_sc * (tcell_k - tref_k));

    let saturation_current = i_o_ref
        * (tcell_k / tref_k).powi(3)
        * (eg_ref / (KB_EV * tref_k) - eg / (KB_EV * tcell_k)).exp();

    let resistance_shunt = if effective_irradiance > 0.0 {
        r_sh_ref * (irrad_ref / effective_irradiance)
    } else {
        f64::INFINITY
    };

    SDMParams {
        photocurrent,
        saturation_current,
        resistance_series: r_s,
        resistance_shunt,
        n_ns_vth,
    }
}

/// Calculate five single-diode model parameters using the CEC model.
///
/// The CEC model differs from De Soto by applying an adjustment factor
/// to the short-circuit current temperature coefficient.
///
/// # Parameters
/// Same as [`calcparams_desoto`] plus:
/// - `adjust`: Adjustment to alpha_sc temperature coefficient [percent]
///
/// # References
/// A. Dobos, 2012, "An Improved Coefficient Calculator for the California
/// Energy Commission 6 Parameter Photovoltaic Module Model", Journal of
/// Solar Energy Engineering, vol 134.
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn calcparams_cec(
    effective_irradiance: f64,
    temp_cell: f64,
    alpha_sc: f64,
    a_ref: f64,
    i_l_ref: f64,
    i_o_ref: f64,
    r_sh_ref: f64,
    r_s: f64,
    adjust: f64,
    eg_ref: f64,
    d_eg_dt: f64,
) -> SDMParams {
    calcparams_desoto(
        effective_irradiance,
        temp_cell,
        alpha_sc * (1.0 - adjust / 100.0),
        a_ref,
        i_l_ref,
        i_o_ref,
        r_sh_ref,
        r_s,
        eg_ref,
        d_eg_dt,
    )
}

/// Calculate five single-diode model parameters using the PVsyst v6 model.
///
/// # Parameters
/// - `effective_irradiance`: Irradiance converted to photocurrent [W/m2]
/// - `temp_cell`: Average cell temperature [C]
/// - `alpha_sc`: Short-circuit current temperature coefficient [A/C]
/// - `gamma_ref`: Diode ideality factor at reference [unitless]
/// - `mu_gamma`: Temperature coefficient for diode ideality factor [1/K]
/// - `i_l_ref`: Photo current at reference conditions [A]
/// - `i_o_ref`: Saturation current at reference conditions [A]
/// - `r_sh_ref`: Shunt resistance at reference conditions [ohm]
/// - `r_sh_0`: Shunt resistance at zero irradiance [ohm]
/// - `r_s`: Series resistance [ohm]
/// - `cells_in_series`: Number of cells in series
/// - `eg_ref`: Bandgap energy at reference temperature [eV]
///
/// # References
/// K. Sauer, T. Roessler, C. W. Hansen, 2015, "Modeling the Irradiance and
/// Temperature Dependence of Photovoltaic Modules in PVsyst",
/// IEEE Journal of Photovoltaics v5(1).
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn calcparams_pvsyst(
    effective_irradiance: f64,
    temp_cell: f64,
    alpha_sc: f64,
    gamma_ref: f64,
    mu_gamma: f64,
    i_l_ref: f64,
    i_o_ref: f64,
    r_sh_ref: f64,
    r_sh_0: f64,
    r_s: f64,
    cells_in_series: u32,
    eg_ref: f64,
) -> SDMParams {
    let irrad_ref = 1000.0;
    let temp_ref = 25.0;
    let r_sh_exp: f64 = 5.5;
    let tref_k = temp_ref + 273.15;
    let tcell_k = temp_cell + 273.15;

    // gamma adjusted for temperature
    let gamma = gamma_ref + mu_gamma * (temp_cell - temp_ref);

    // nNsVth = gamma * k/q * Ns * Tcell_K
    let n_ns_vth = gamma * KB_J / Q_C * (cells_in_series as f64) * tcell_k;

    // Photocurrent
    let photocurrent = effective_irradiance / irrad_ref
        * (i_l_ref + alpha_sc * (tcell_k - tref_k));

    // Saturation current (PVsyst uses q*Eg/(k*gamma) formulation)
    let saturation_current = i_o_ref
        * (tcell_k / tref_k).powi(3)
        * ((Q_C * eg_ref) / (KB_J * gamma) * (1.0 / tref_k - 1.0 / tcell_k)).exp();

    // Shunt resistance: PVsyst exponential model
    let rsh_tmp = (r_sh_ref - r_sh_0 * (-r_sh_exp).exp()) / (1.0 - (-r_sh_exp).exp());
    let rsh_base = rsh_tmp.max(0.0);
    let resistance_shunt = rsh_base
        + (r_sh_0 - rsh_base) * (-r_sh_exp * effective_irradiance / irrad_ref).exp();

    SDMParams {
        photocurrent,
        saturation_current,
        resistance_series: r_s,
        resistance_shunt,
        n_ns_vth,
    }
}

// ---------------------------------------------------------------------------
// SAPM (Sandia Array Performance Model)
// ---------------------------------------------------------------------------

/// SAPM module parameters.
#[derive(Debug, Clone, Copy)]
pub struct SAPMParams {
    /// Short-circuit current at reference [A]
    pub isco: f64,
    /// Max-power current at reference [A]
    pub impo: f64,
    /// Open-circuit voltage at reference [V]
    pub voco: f64,
    /// Max-power voltage at reference [V]
    pub vmpo: f64,
    /// Isc temperature coefficient [1/C]
    pub aisc: f64,
    /// Imp temperature coefficient [1/C]
    pub aimp: f64,
    /// Voc temperature coefficient [V/C]
    pub bvoco: f64,
    /// Irradiance dependence for BetaVoc [V/C]
    pub mbvoc: f64,
    /// Vmp temperature coefficient [V/C]
    pub bvmpo: f64,
    /// Irradiance dependence for BetaVmp [V/C]
    pub mbvmp: f64,
    /// Diode factor [unitless]
    pub n: f64,
    /// Number of cells in series
    pub cells_in_series: u32,
    /// Empirical coefficients C0..C3
    pub c0: f64,
    pub c1: f64,
    pub c2: f64,
    pub c3: f64,
    /// Airmass coefficients A0..A4
    pub a0: f64,
    pub a1: f64,
    pub a2: f64,
    pub a3: f64,
    pub a4: f64,
    /// AOI coefficients B0..B5
    pub b0: f64,
    pub b1: f64,
    pub b2: f64,
    pub b3: f64,
    pub b4: f64,
    pub b5: f64,
    /// Fraction of diffuse irradiance used
    pub fd: f64,
}

/// SAPM output: key points on the I-V curve.
#[derive(Debug, Clone, Copy)]
pub struct SAPMOutput {
    /// Short-circuit current [A]
    pub i_sc: f64,
    /// Current at max-power point [A]
    pub i_mp: f64,
    /// Open-circuit voltage [V]
    pub v_oc: f64,
    /// Voltage at max-power point [V]
    pub v_mp: f64,
    /// Power at max-power point [W]
    pub p_mp: f64,
}

/// Sandia PV Array Performance Model (SAPM).
///
/// Generates 5 points on a PV module's I-V curve according to
/// SAND2004-3535.
///
/// # Parameters
/// - `effective_irradiance`: Irradiance reaching module cells [W/m2]
/// - `temp_cell`: Cell temperature [C]
/// - `module`: SAPM module parameters
///
/// # References
/// King, D. et al, 2004, "Sandia Photovoltaic Array Performance Model",
/// SAND Report 3535, Sandia National Laboratories.
#[inline]
pub fn sapm(effective_irradiance: f64, temp_cell: f64, module: &SAPMParams) -> SAPMOutput {
    let irradiance_ref = 1000.0;
    let temperature_ref = 25.0;

    let ee = effective_irradiance / irradiance_ref;
    let dt = temp_cell - temperature_ref;

    let delta = module.n * KB_J * (temp_cell + 273.15) / Q_C;
    let cells = module.cells_in_series as f64;

    let bvmpo = module.bvmpo + module.mbvmp * (1.0 - ee);
    let bvoco = module.bvoco + module.mbvoc * (1.0 - ee);

    let log_ee = if ee > 0.0 { ee.ln() } else { f64::NEG_INFINITY };

    let i_sc = module.isco * ee * (1.0 + module.aisc * dt);

    let i_mp = module.impo * (module.c0 * ee + module.c1 * ee.powi(2))
        * (1.0 + module.aimp * dt);

    let v_oc = (module.voco + cells * delta * log_ee + bvoco * dt).max(0.0);

    let v_mp = (module.vmpo
        + module.c2 * cells * delta * log_ee
        + module.c3 * cells * (delta * log_ee).powi(2)
        + bvmpo * dt)
        .max(0.0);

    let p_mp = i_mp * v_mp;

    SAPMOutput { i_sc, i_mp, v_oc, v_mp, p_mp }
}

/// SAPM spectral factor: fourth-degree polynomial in airmass.
///
/// Calculates the spectral mismatch factor f1 for the SAPM model.
fn sapm_spectral_factor(airmass_absolute: f64, module: &SAPMParams) -> f64 {
    let am = airmass_absolute;
    let f1 = module.a0
        + module.a1 * am
        + module.a2 * am.powi(2)
        + module.a3 * am.powi(3)
        + module.a4 * am.powi(4);

    if f1.is_nan() { 0.0 } else { f1.max(0.0) }
}

/// Calculate SAPM effective irradiance.
///
/// Accounts for spectral and angle-of-incidence losses using the SAPM model:
/// `Ee = f1(AM) * (Eb * f2(AOI) + fd * Ed)`
///
/// # Parameters
/// - `poa_direct`: Direct irradiance on the module [W/m2]
/// - `poa_diffuse`: Diffuse irradiance on the module [W/m2]
/// - `airmass_absolute`: Absolute airmass [unitless]
/// - `aoi_val`: Angle of incidence [degrees]
/// - `module`: SAPM module parameters
///
/// # References
/// King, D. et al, 2004, "Sandia Photovoltaic Array Performance Model",
/// SAND2004-3535, Sandia National Laboratories.
#[inline]
pub fn sapm_effective_irradiance(
    poa_direct: f64,
    poa_diffuse: f64,
    airmass_absolute: f64,
    aoi_val: f64,
    module: &SAPMParams,
) -> f64 {
    let f1 = sapm_spectral_factor(airmass_absolute, module);
    let f2 = crate::iam::sapm(aoi_val, module.b0, module.b1, module.b2, module.b3, module.b4, module.b5);

    f1 * (poa_direct * f2 + module.fd * poa_diffuse)
}
