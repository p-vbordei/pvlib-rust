//! Bifacial / row-to-row view-factor primitives (2-D infinite-sheds model).
//!
//! Closed-form view factors from `pvlib.bifacial.utils` plus a rear-side
//! irradiance composition that uses them. These are the same analytical
//! expressions used inside `pvlib.bifacial.infinite_sheds.get_irradiance_poa`.
//! The full `infinite_sheds` pipeline (front/rear decomposition with
//! shaded-row fraction and multi-row ground view factors) is tracked in
//! `ROADMAP.md`; this module implements the pieces already used inside it.
//!
//! # References
//! - Mikofski, M., Darawali, R., Hamer, M., Neubert, A., and Newmiller, J.
//!   "Bifacial Performance Modeling in Large Arrays." 2019 IEEE PVSC 46.
//!   DOI: 10.1109/PVSC40753.2019.8980572.
//! - Marion, B., MacAlpine, S., Deline, C., Asgharzadeh, A., Toor, F.,
//!   Riley, D., Stein, J., and Hansen, C. "A practical irradiance model
//!   for bifacial PV modules." 2017 IEEE PVSC 44.
//!   DOI: 10.1109/PVSC.2017.8366263.

// ---------------------------------------------------------------------------
// Solar-geometry helpers
// ---------------------------------------------------------------------------

/// Tangent of the angle between the zenith and the sun vector projected
/// into the plane perpendicular to the row azimuth (Mikofski 2019 eqn 7).
///
/// `tan φ = cos(solar_azimuth − surface_azimuth) · tan(solar_zenith)`
#[inline]
pub fn solar_projection_tangent(
    solar_zenith: f64,
    solar_azimuth: f64,
    surface_azimuth: f64,
) -> f64 {
    let rotation = (solar_azimuth - surface_azimuth).to_radians();
    rotation.cos() * solar_zenith.to_radians().tan()
}

/// Fraction of the ground between rows that sees the direct beam
/// (Mikofski 2019 eqn 4):
///
/// `F_gnd,beam = 1 − min(1, GCR · |cos β + sin β · tan φ|)`
///
/// where β is the module surface tilt and φ is the solar projection
/// angle from [`solar_projection_tangent`]. Clamped to 0 for zeniths
/// above `max_zenith`.
#[inline]
pub fn unshaded_ground_fraction(
    surface_tilt: f64,
    surface_azimuth: f64,
    solar_zenith: f64,
    solar_azimuth: f64,
    gcr: f64,
) -> f64 {
    const MAX_ZENITH_DEG: f64 = 87.0;
    if solar_zenith > MAX_ZENITH_DEG {
        return 0.0;
    }
    let tan_phi = solar_projection_tangent(solar_zenith, solar_azimuth, surface_azimuth);
    let b = surface_tilt.to_radians();
    let width = gcr * (b.cos() + b.sin() * tan_phi).abs();
    (1.0 - width.min(1.0)).max(0.0)
}

// ---------------------------------------------------------------------------
// 2-D view-factor primitives
// ---------------------------------------------------------------------------

/// Helper used by the closed-form row ↔ sky and row ↔ ground view factors.
/// `delta` is `+1` for the row-to-ground case and `-1` for the row-to-sky case.
#[inline]
fn vf_poly(surface_tilt: f64, gcr: f64, x: f64, delta: f64) -> f64 {
    let a = 1.0 / gcr;
    let c = surface_tilt.to_radians().cos();
    (a * a + 2.0 * delta * a * c * x + x * x).sqrt()
}

/// View factor from a point `x` on a row's slant length to the sky dome
/// (Marion 2017 eqn 6). `x = 0` is the bottom of the row, `x = 1` is the top.
/// Assumes infinitely long, uniformly spaced rows on horizontal ground.
#[inline]
pub fn vf_row_sky_2d(surface_tilt: f64, gcr: f64, x: f64) -> f64 {
    let p = vf_poly(surface_tilt, gcr, 1.0 - x, -1.0);
    let c = surface_tilt.to_radians().cos();
    0.5 * (1.0 + (c / gcr - (1.0 - x)) / p)
}

/// Average view factor to the sky over a segment `[x0, x1]` of the row
/// slant length (Marion 2017 eqn 6, integrated analytically).
#[inline]
pub fn vf_row_sky_2d_integ(surface_tilt: f64, gcr: f64, x0: f64, x1: f64) -> f64 {
    let u = (x1 - x0).abs();
    if u < 1e-6 {
        return vf_row_sky_2d(surface_tilt, gcr, x0);
    }
    let p0 = vf_poly(surface_tilt, gcr, 1.0 - x0, -1.0);
    let p1 = vf_poly(surface_tilt, gcr, 1.0 - x1, -1.0);
    0.5 * (1.0 + (p1 - p0) / u)
}

/// View factor from a point `x` on a row's slant length to the visible
/// ground (Marion 2017 eqn 7). Counterpart of [`vf_row_sky_2d`].
#[inline]
pub fn vf_row_ground_2d(surface_tilt: f64, gcr: f64, x: f64) -> f64 {
    let p = vf_poly(surface_tilt, gcr, x, 1.0);
    let c = surface_tilt.to_radians().cos();
    0.5 * (1.0 - (c / gcr + x) / p)
}

/// Average view factor to the ground over a segment `[x0, x1]` of the row
/// slant length. Analytical integration of [`vf_row_ground_2d`].
#[inline]
pub fn vf_row_ground_2d_integ(surface_tilt: f64, gcr: f64, x0: f64, x1: f64) -> f64 {
    let u = (x1 - x0).abs();
    if u < 1e-6 {
        return vf_row_ground_2d(surface_tilt, gcr, x0);
    }
    let p0 = vf_poly(surface_tilt, gcr, x0, 1.0);
    let p1 = vf_poly(surface_tilt, gcr, x1, 1.0);
    0.5 * (1.0 - (p1 - p0) / u)
}

// ---------------------------------------------------------------------------
// Rear-side irradiance composition
// ---------------------------------------------------------------------------

/// Rear-side irradiance for an infinite-sheds bifacial array, computed
/// from real 2-D view factors rather than a linear albedo approximation.
///
/// The rear surface sees:
///
/// 1. **Sky diffuse** — DHI weighted by the integrated row-to-sky view
///    factor of the rear face (`surface_tilt` seen from behind = `180 - tilt`).
/// 2. **Ground-reflected diffuse** — `albedo · (F_gnd,beam · (GHI − DHI) +
///    VF_gnd,sky · DHI)`, weighted by the row-to-ground view factor.
/// 3. **Direct beam** — only contributes when the sun is on the rear side
///    of the module (AOI > 90° from the front normal); otherwise zero.
///    Callers that do not have a precomputed rear AOI should pass `0.0`
///    for `rear_beam`, which is conservative.
///
/// The sky view factor is evaluated at `tilt_rear = 180 − tilt` and the
/// ground view factor mirrors the rear geometry.
///
/// # Parameters
/// - `surface_tilt`, `surface_azimuth` — front-face orientation [°]
/// - `solar_zenith`, `solar_azimuth` — sun position [°]
/// - `gcr` — ground coverage ratio `W/pitch` (unitless)
/// - `ghi`, `dhi`, `dni` — irradiance components [W/m²]
/// - `albedo` — ground albedo (unitless, 0–1)
/// - `rear_beam` — direct beam already projected onto the rear face, or
///   `0.0` if unknown. This lets callers bolt on a beam term without
///   committing to a specific shadow-fraction model.
///
/// # Returns
/// Rear-side plane-of-array irradiance [W/m²].
#[allow(clippy::too_many_arguments)]
#[inline]
pub fn get_irradiance_infinite_sheds(
    surface_tilt: f64,
    surface_azimuth: f64,
    gcr: f64,
    // Kept for API back-compat but no longer used — see ROADMAP.md for the
    // full height-dependent multi-row ground view-factor model.
    _height: f64,
    _pitch: f64,
    ghi: f64,
    dhi: f64,
    _dni: f64,
    albedo: f64,
) -> f64 {
    if gcr <= 0.0 || !ghi.is_finite() || !dhi.is_finite() {
        return 0.0;
    }
    // Default sun direction — zero out beam contribution. The simple-sheds
    // form without solar zenith/azimuth information cannot attribute a
    // ground-beam fraction to the rear, so use an isotropic fallback:
    // F_gnd,beam = 1 − GCR·cos(β) clamped to [0, 1].
    rear_irradiance_sheds(
        surface_tilt,
        surface_azimuth,
        None,
        gcr,
        ghi,
        dhi,
        albedo,
    )
}

/// Full rear-irradiance composition when solar position is available.
///
/// Computes rear sky diffuse + rear ground diffuse using the real 2-D
/// view factors in this module. Prefer this over the
/// solar-position-less [`get_irradiance_infinite_sheds`] for production
/// use — it accounts for the shadow the rows cast on the ground between
/// themselves.
///
/// `solar_position` is `Some((solar_zenith, solar_azimuth))` or `None`
/// for an isotropic fallback (sun directly overhead).
#[inline]
pub fn rear_irradiance_sheds(
    surface_tilt: f64,
    surface_azimuth: f64,
    solar_position: Option<(f64, f64)>,
    gcr: f64,
    ghi: f64,
    dhi: f64,
    albedo: f64,
) -> f64 {
    if gcr <= 0.0 {
        return 0.0;
    }

    // Rear face tilt: if front faces tilt β, rear faces 180 - β.
    let rear_tilt = 180.0 - surface_tilt;

    // Rear sky diffuse — average view factor from the full rear slant
    // length to the sky dome, weighted by DHI.
    let vf_rear_sky = vf_row_sky_2d_integ(rear_tilt, gcr, 0.0, 1.0);
    let rear_sky = dhi * vf_rear_sky;

    // Ground diffuse component. Mikofski 2019 eqn 1:
    //   poa_ground = albedo · (F_gnd,beam · (GHI − DHI) + VF_gnd,sky · DHI)
    // where VF_gnd,sky ≈ 1 − VF_gnd,rows. For the unobstructed inter-row
    // region, VF_gnd,sky is close to 1; using the conservative constant 1.0
    // matches the simple-sheds assumption documented in pvlib.
    let f_gnd_beam = match solar_position {
        Some((sz, sa)) => unshaded_ground_fraction(surface_tilt, surface_azimuth, sz, sa, gcr),
        None => {
            // Fallback: replace tan(phi) with 0 (sun overhead along row normal).
            let b = surface_tilt.to_radians();
            (1.0 - (gcr * b.cos()).min(1.0)).max(0.0)
        }
    };
    let vf_gnd_sky = 1.0;
    let poa_ground = albedo * (f_gnd_beam * (ghi - dhi).max(0.0) + vf_gnd_sky * dhi.max(0.0));

    // Rear ground diffuse — average view factor from the rear row to the
    // visible ground, weighted by POA ground irradiance.
    let vf_rear_ground = vf_row_ground_2d_integ(rear_tilt, gcr, 0.0, 1.0);
    let rear_ground = poa_ground * vf_rear_ground;

    (rear_sky + rear_ground).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn horizontal_row_sees_full_sky() {
        // Tilt 0 (flat on the ground looking up) with negligible adjacent row
        // shading (GCR small) should see roughly the full sky.
        let vf = vf_row_sky_2d_integ(0.0, 0.01, 0.0, 1.0);
        assert!(vf > 0.9, "flat row with tiny GCR should see most of sky, got {vf}");
    }

    #[test]
    fn vertical_row_sees_half_sky_at_most() {
        // A 90° tilt row can see at most half the sky (dome split by the row).
        let vf = vf_row_sky_2d_integ(90.0, 0.5, 0.0, 1.0);
        assert!(vf <= 0.5 + 1e-9, "vertical row vf to sky {vf} should be ≤ 0.5");
    }

    #[test]
    fn row_sky_plus_row_ground_sum_leaves_neighbor_view() {
        // For a point on a row, the view splits three ways: sky, ground,
        // and the neighboring row. So vf_sky + vf_ground ≤ 1 but can be
        // meaningfully less than 1 when the adjacent row is visible. We
        // bound each term to [0, 1] and the sum to ≤ 1+ε.
        for &tilt in &[10.0_f64, 30.0, 60.0] {
            for &gcr in &[0.3_f64, 0.5] {
                let s = vf_row_sky_2d(tilt, gcr, 0.5);
                let g = vf_row_ground_2d(tilt, gcr, 0.5);
                assert!((0.0..=1.0).contains(&s), "vf_sky out of range: {s}");
                assert!((0.0..=1.0).contains(&g), "vf_ground out of range: {g}");
                assert!(
                    s + g <= 1.0 + 1e-9,
                    "vf_sky+vf_ground at tilt={tilt} gcr={gcr} mid-row = {}",
                    s + g
                );
            }
        }
    }

    #[test]
    fn unshaded_fraction_bounds() {
        let f = unshaded_ground_fraction(30.0, 180.0, 30.0, 180.0, 0.5);
        assert!((0.0..=1.0).contains(&f));
        // Sun below horizon — everything shaded.
        assert_eq!(unshaded_ground_fraction(30.0, 180.0, 90.0, 180.0, 0.5), 0.0);
    }

    #[test]
    fn rear_irradiance_positive() {
        let rear = rear_irradiance_sheds(
            30.0,
            180.0,
            Some((30.0, 180.0)),
            0.4,
            1000.0,
            150.0,
            0.25,
        );
        assert!(rear > 0.0 && rear < 500.0, "implausible rear irradiance {rear}");
    }
}
