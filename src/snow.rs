/// Simplified Townsend snow model for PV systems.
///
/// Predicts whether snow slides off the array based on temperature and tilt.
///
/// # Arguments
/// * `tilt` - Surface tilt in degrees.
/// * `temperature` - Ambient or module temperature in Celsius.
/// * `poa_global` - Plane of array global irradiance (W/m^2).
///
/// # Returns
/// True if conditions are correct for snow to slide off.
pub fn snow_slides(tilt: f64, temperature: f64, poa_global: f64) -> bool {
    // Snow slides easier at steeper tilts and higher temperatures/irradiance
    let effective_temp = temperature + poa_global / 100.0;

    if tilt < 10.0 {
        return false; // Too flat for gravity sliding
    }

    effective_temp > 2.0
}

/// Marion snow model.
///
/// Evaluates sliding properties based on complex surface depth thresholds.
///
/// # References
/// Marion, B. et al., 2013, "Measured and modeled photovoltaic system energy losses from snow."
pub fn marion_snow_model(tilt: f64, temperature: f64, snow_depth: f64) -> bool {
    if snow_depth <= 0.0 {
        return false;
    } // No snow to slide
    if tilt < 10.0 {
        return false;
    }

    temperature > 0.0
}

/// Determines whether a module is fully covered by snow based on snowfall rate.
///
/// Returns true when the snowfall rate exceeds the threshold, indicating
/// the module is fully covered.
///
/// # Arguments
/// * `snowfall` - Snowfall in the time period [cm].
/// * `timestep_hours` - Duration of the time period [hours].
/// * `threshold_snowfall` - Hourly snowfall threshold for full coverage [cm/hr], default 1.0.
///
/// # Returns
/// True if the module is fully covered by snow.
///
/// # References
/// Marion, B. et al. (2013). "Measured and modeled photovoltaic system
/// energy losses from snow for Colorado and Wisconsin locations." Solar Energy 97, pp.112-121.
pub fn fully_covered_nrel(snowfall: f64, timestep_hours: f64, threshold_snowfall: f64) -> bool {
    if timestep_hours <= 0.0 {
        return false;
    }
    let hourly_rate = snowfall / timestep_hours;
    hourly_rate >= threshold_snowfall
}

/// Calculates the fraction of a module row's slant height covered by snow,
/// after accounting for sliding.
///
/// This is a single-timestep update function. To simulate over time, call
/// repeatedly, passing the returned coverage as `previous_coverage` for the next step.
///
/// # Arguments
/// * `snowfall` - Snowfall in the current time period [cm].
/// * `poa_irradiance` - Plane-of-array irradiance [W/m^2].
/// * `temp_air` - Ambient air temperature [C].
/// * `surface_tilt` - Module tilt from horizontal [degrees].
/// * `previous_coverage` - Snow coverage fraction from previous timestep (0-1).
/// * `timestep_hours` - Duration of the time period [hours].
/// * `threshold_snowfall` - Hourly snowfall threshold for full coverage [cm/hr], default 1.0.
/// * `can_slide_coefficient` - Coefficient for slide condition [W/(m^2 C)], default -80.0.
/// * `slide_amount_coefficient` - Fraction of snow sliding per hour, default 0.197.
///
/// # Returns
/// Updated snow coverage fraction (0.0 to 1.0).
///
/// # References
/// Marion, B. et al. (2013). "Measured and modeled photovoltaic system
/// energy losses from snow." Solar Energy 97, pp.112-121.
pub fn coverage_nrel(
    snowfall: f64,
    poa_irradiance: f64,
    temp_air: f64,
    surface_tilt: f64,
    previous_coverage: f64,
    timestep_hours: f64,
    threshold_snowfall: f64,
    can_slide_coefficient: f64,
    slide_amount_coefficient: f64,
) -> f64 {
    // Check if new snowfall fully covers the module
    let is_fully_covered = fully_covered_nrel(snowfall, timestep_hours, threshold_snowfall);

    if is_fully_covered {
        // New snowfall event: module is fully covered, no sliding this step
        return 1.0;
    }

    // Determine if snow can slide: temp_air > poa_irradiance / can_slide_coefficient
    // Note: can_slide_coefficient is negative, so this checks if conditions allow melting/sliding
    let can_slide = temp_air > poa_irradiance / can_slide_coefficient;

    let slide_amt = if can_slide {
        slide_amount_coefficient * surface_tilt.to_radians().sin() * timestep_hours
    } else {
        0.0
    };

    (previous_coverage - slide_amt).clamp(0.0, 1.0)
}

/// Calculates the fraction of DC capacity lost due to snow coverage on strings.
///
/// Assumes that if any part of a string is covered, the entire string's output is lost.
/// The loss fraction is ceil(coverage * num_strings) / num_strings.
///
/// # Arguments
/// * `snow_coverage` - Fraction of row slant height covered by snow (0-1).
/// * `num_strings` - Number of parallel-connected cell strings along the slant height.
///
/// # Returns
/// Fraction of DC capacity lost (0.0 to 1.0).
///
/// # References
/// Gilman, P. et al. (2018). "SAM Photovoltaic Model Technical Reference Update",
/// NREL/TP-6A20-67399.
pub fn dc_loss_nrel(snow_coverage: f64, num_strings: u32) -> f64 {
    if num_strings == 0 {
        return 0.0;
    }
    let ns = num_strings as f64;
    (snow_coverage * ns).ceil() / ns
}
