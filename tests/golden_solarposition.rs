//! Gate-G1 golden parity: `solarposition::get_solarposition` /
//! `get_solarposition_apparent` vs pinned pvlib-python `spa_python` fixtures.
//!
//! See `tests/golden/mod.rs` for the shared helper, fixture schemas, and the
//! azimuth-convention reconciliation. Skips cleanly if `PVLIB_GOLDEN_DIR` is
//! unset.

#[path = "golden/mod.rs"]
mod golden;

use pvlib::location::Location;
use pvlib::solarposition::{get_solarposition, get_solarposition_apparent};

fn locations_by_slug(dir: &std::path::Path) -> std::collections::HashMap<String, Location> {
    golden::load_locations(dir)
        .into_iter()
        .map(|(slug, m)| {
            let loc = Location::new(m.lat, m.lon, chrono_tz::UTC, m.alt, &m.slug);
            (slug, loc)
        })
        .collect()
}

/// True (non-refraction-corrected) position: Rust `get_solarposition` vs
/// pvlib's `zenith` / `elevation` / `azimuth` columns (geometric, not the
/// `apparent_*` / `ephem_*` columns).
#[test]
fn true_position_matches_pvlib() {
    let Some(dir) = golden::fixture_dir() else {
        eprintln!("{}", golden::SKIP_MESSAGE);
        return;
    };
    let locations = locations_by_slug(&dir);

    let mut day = golden::TrackerSet::new();
    let mut night = golden::TrackerSet::new();
    let mut overall = golden::TrackerSet::new();

    let mut rdr = golden::open_csv_gz(&dir, "solarposition.csv.gz");
    for rec in rdr.deserialize::<golden::SolarPositionRow>() {
        let row = rec.expect("solarposition.csv.gz row");
        let loc = locations
            .get(&row.slug)
            .unwrap_or_else(|| panic!("unknown slug {} (not in locations.json)", row.slug));
        let t = golden::parse_ts_z_tz(&row.timestamp_utc);
        let sp = get_solarposition(loc, t).expect("get_solarposition");

        let zen_err = (sp.zenith - row.zenith).abs();
        let elev_err = (sp.elevation - row.elevation).abs();
        let az_err = golden::azimuth_diff(sp.azimuth, row.azimuth);

        let is_day = row.apparent_elevation > 0.0;
        let tracker = if is_day { &mut day } else { &mut night };
        let slug = row.slug.clone();
        let ts = row.timestamp_utc.clone();
        tracker.record(&format!("{slug}|zenith"), zen_err, || {
            format!(
                "slug={slug} ts={ts} col=zenith got={:.6} want={:.6}",
                sp.zenith, row.zenith
            )
        });
        tracker.record(&format!("{slug}|elevation"), elev_err, || {
            format!(
                "slug={slug} ts={ts} col=elevation got={:.6} want={:.6}",
                sp.elevation, row.elevation
            )
        });
        tracker.record(&format!("{slug}|azimuth"), az_err, || {
            format!(
                "slug={slug} ts={ts} col=azimuth got={:.6} want={:.6}",
                sp.azimuth, row.azimuth
            )
        });

        if is_day {
            let combined = zen_err.max(elev_err).max(az_err);
            overall.record("all", combined, || {
                format!(
                    "slug={slug} ts={ts} zen_err={zen_err:.6} elev_err={elev_err:.6} az_err={az_err:.6}"
                )
            });
        }
    }

    day.report("golden_solarposition: TRUE position, daytime (apparent_elevation > 0)");
    night.report(
        "golden_solarposition: TRUE position, nighttime (informational only, not asserted)",
    );
    overall.report(
        "golden_solarposition: TRUE position, per-row combined max-component error, daytime",
    );

    let max_day = day.global_max();
    assert!(
        max_day <= golden::POSITION_TOL_DEG,
        "true solar position exceeded {} deg tolerance during daytime: max={:.6} \
         (see per-location/per-column report above for the worst offenders)",
        golden::POSITION_TOL_DEG,
        max_day
    );
}

/// Apparent (refraction-corrected) position: Rust `get_solarposition_apparent`
/// (defaults: pressure 101325 Pa, temperature 12 C — matching
/// `spa_python`'s defaults used to generate the fixture) vs pvlib's
/// `apparent_zenith` / `apparent_elevation` columns.
///
/// This is the load-bearing check for Gate G1. As of the `solar-positioning`
/// crate swap (`docs/decisions/spa-solar-position.md` in solar-engine),
/// pvlib-rust computes this via the full Reda & Andreas NREL SPA
/// (±0.0003 deg), replacing the earlier `spa` crate dependency, which —
/// despite the name — implemented PSA (Blanco-Muriel et al. 2001), documented
/// as "accurate to within 0.5 minutes of arc for the year 1999 and
/// following". This test runs the current backend against all of 2025 (or
/// the 12-day subset) to confirm the swap holds the 0.01 deg Gate G1 budget
/// with wide margin, not just barely.
#[test]
fn apparent_position_matches_pvlib() {
    let Some(dir) = golden::fixture_dir() else {
        eprintln!("{}", golden::SKIP_MESSAGE);
        return;
    };
    let locations = locations_by_slug(&dir);

    let mut day = golden::TrackerSet::new();
    let mut night = golden::TrackerSet::new();
    let mut overall = golden::TrackerSet::new();

    let mut rdr = golden::open_csv_gz(&dir, "solarposition.csv.gz");
    for rec in rdr.deserialize::<golden::SolarPositionRow>() {
        let row = rec.expect("solarposition.csv.gz row");
        let loc = locations
            .get(&row.slug)
            .unwrap_or_else(|| panic!("unknown slug {} (not in locations.json)", row.slug));
        let t = golden::parse_ts_z_tz(&row.timestamp_utc);
        // None, None -> pvlib spa_python defaults (101325 Pa, 12 C), matching
        // how generate.py called solarposition.spa_python(times, lat, lon, altitude=alt).
        let sp =
            get_solarposition_apparent(loc, t, None, None).expect("get_solarposition_apparent");

        let zen_err = (sp.apparent_zenith - row.apparent_zenith).abs();
        let elev_err = (sp.apparent_elevation - row.apparent_elevation).abs();

        let is_day = row.apparent_elevation > 0.0;
        let tracker = if is_day { &mut day } else { &mut night };
        let slug = row.slug.clone();
        let ts = row.timestamp_utc.clone();
        tracker.record(&format!("{slug}|apparent_zenith"), zen_err, || {
            format!(
                "slug={slug} ts={ts} col=apparent_zenith got={:.6} want={:.6}",
                sp.apparent_zenith, row.apparent_zenith
            )
        });
        tracker.record(&format!("{slug}|apparent_elevation"), elev_err, || {
            format!(
                "slug={slug} ts={ts} col=apparent_elevation got={:.6} want={:.6}",
                sp.apparent_elevation, row.apparent_elevation
            )
        });

        if is_day {
            let combined = zen_err.max(elev_err);
            overall.record("all", combined, || {
                format!("slug={slug} ts={ts} zen_err={zen_err:.6} elev_err={elev_err:.6}")
            });
        }
    }

    day.report("golden_solarposition: APPARENT position, daytime (apparent_elevation > 0)");
    night.report(
        "golden_solarposition: APPARENT position, nighttime (informational only, not asserted)",
    );
    overall.report(
        "golden_solarposition: APPARENT position, per-row combined max-component error, daytime",
    );

    let max_day = day.global_max();
    assert!(
        max_day <= golden::POSITION_TOL_DEG,
        "apparent solar position exceeded {} deg tolerance during daytime: max={:.6} \
         (see per-location/per-column report above for the worst offenders -- \
         this is the PSA-vs-2025-drift question from docs/PLAN.md Etapa 1)",
        golden::POSITION_TOL_DEG,
        max_day
    );
}

/// Azimuth convention reconciliation, verified numerically (not just
/// asserted in a doc comment): at each location's daily solar transit
/// (minimum zenith of the day), the sun should be roughly due south
/// (azimuth ~ 180 deg) for these northern-hemisphere fleet locations, in
/// BOTH pvlib's fixture and pvlib-rust's output. If pvlib-rust's SPA backend
/// used a different zero point (e.g. "0 = South" instead of "0 = North"),
/// this would show transit azimuths near 0 instead of 180 and this test
/// would fail loudly.
#[test]
fn azimuth_zero_point_matches_pvlib_north_convention() {
    let Some(dir) = golden::fixture_dir() else {
        eprintln!("{}", golden::SKIP_MESSAGE);
        return;
    };
    let locations = locations_by_slug(&dir);

    // Find the fixture row with minimum zenith (closest to transit) per slug.
    let mut best_per_slug: std::collections::HashMap<String, golden::SolarPositionRow> =
        std::collections::HashMap::new();
    let mut rdr = golden::open_csv_gz(&dir, "solarposition.csv.gz");
    for rec in rdr.deserialize::<golden::SolarPositionRow>() {
        let row = rec.expect("solarposition.csv.gz row");
        best_per_slug
            .entry(row.slug.clone())
            .and_modify(|cur| {
                if row.zenith < cur.zenith {
                    *cur = row.clone();
                }
            })
            .or_insert(row);
    }

    assert!(
        !best_per_slug.is_empty(),
        "no rows read from solarposition.csv.gz"
    );
    let mut slugs: Vec<&String> = best_per_slug.keys().collect();
    slugs.sort();
    for slug in slugs {
        let row = &best_per_slug[slug];
        let loc = locations.get(slug).unwrap();
        let t = golden::parse_ts_z_tz(&row.timestamp_utc);
        let sp = get_solarposition(loc, t).unwrap();
        eprintln!(
            "  transit check {slug}: fixture azimuth={:.2} rust azimuth={:.2} (zenith={:.2})",
            row.azimuth, sp.azimuth, row.zenith
        );
        assert!(
            (row.azimuth - 180.0).abs() < 60.0,
            "{slug}: fixture azimuth at daily transit should be near 180 (south) for a \
             northern-hemisphere location under the 'degrees clockwise from North' \
             convention, got {}",
            row.azimuth
        );
        assert!(
            (sp.azimuth - 180.0).abs() < 60.0,
            "{slug}: pvlib-rust azimuth at daily transit should be near 180 (south), got {} \
             -- indicates pvlib-rust's SPA backend's azimuth zero point does NOT match pvlib's",
            sp.azimuth
        );
    }
}
