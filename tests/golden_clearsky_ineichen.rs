//! Gate-G1 golden parity: the Ineichen clearsky chain, replicated exactly
//! per MANIFEST.json ("clearsky_ineichen.csv.gz" -> "chain"):
//!
//!   am_rel = atmosphere::get_relative_airmass(apparent_zenith)   // Kasten-Young 1989
//!   pressure = atmosphere::alt2pres(altitude)
//!   am_abs = atmosphere::get_absolute_airmass(am_rel, pressure)
//!   dni_extra = irradiance::get_extra_radiation(day_of_year)     // Spencer 1971
//!   clearsky::ineichen(apparent_zenith, am_abs, linke_turbidity, altitude, dni_extra)
//!
//! `apparent_zenith` is this fixture's own input, joined from
//! `solarposition.csv.gz` by `(slug, timestamp_utc)` since
//! `clearsky_ineichen.csv.gz` doesn't carry its own solar-angle column (see
//! generate.py: both come from the same `sp` DataFrame).

#[path = "golden/mod.rs"]
mod golden;

use pvlib::{atmosphere, clearsky, irradiance};

#[test]
fn ineichen_chain_matches_pvlib() {
    let Some(dir) = golden::fixture_dir() else {
        eprintln!("{}", golden::SKIP_MESSAGE);
        return;
    };
    let locations = golden::load_locations(&dir);
    let angles = golden::load_solar_angles(&dir);

    let mut stats = golden::TrackerSet::new();
    let mut night_violations: Vec<String> = Vec::new();
    let mut tol_violations: Vec<String> = Vec::new();
    let mut n_tol_violations = 0usize;

    let mut rdr = golden::open_csv_gz(&dir, "clearsky_ineichen.csv.gz");
    for rec in rdr.deserialize::<golden::ClearskyRow>() {
        let row = rec.expect("clearsky_ineichen.csv.gz row");
        let loc = locations
            .get(&row.slug)
            .unwrap_or_else(|| panic!("unknown slug {} (not in locations.json)", row.slug));
        let sa = golden::lookup_angles(&angles, &row.slug, &row.timestamp_utc);

        let am_rel = atmosphere::get_relative_airmass(sa.apparent_zenith);
        let pressure = atmosphere::alt2pres(loc.alt);
        let am_abs = atmosphere::get_absolute_airmass(am_rel, pressure);
        let ts = golden::parse_ts_z(&row.timestamp_utc);
        let dni_extra = irradiance::get_extra_radiation(golden::day_of_year(ts) as i32);

        let cs = clearsky::ineichen(
            sa.apparent_zenith,
            am_abs,
            row.linke_turbidity,
            loc.alt,
            dni_extra,
        );

        // Fixture rows are exactly 0.0 at night (asserted at fixture-generation
        // time too); pvlib-rust's ineichen() short-circuits to exact 0.0 for
        // zenith >= 90 deg, so this should hold trivially -- verify it does.
        let is_night = row.ghi == 0.0 && row.dni == 0.0 && row.dhi == 0.0;
        if is_night && (cs.ghi != 0.0 || cs.dni != 0.0 || cs.dhi != 0.0) {
            night_violations.push(format!(
                "slug={} ts={} TL={} got=({:.6},{:.6},{:.6}) want=(0,0,0)",
                row.slug, row.timestamp_utc, row.linke_turbidity, cs.ghi, cs.dni, cs.dhi
            ));
        }

        let slug = row.slug.clone();
        let ts_s = row.timestamp_utc.clone();
        let tl = row.linke_turbidity;
        for (col, got, want) in [
            ("ghi", cs.ghi, row.ghi),
            ("dni", cs.dni, row.dni),
            ("dhi", cs.dhi, row.dhi),
        ] {
            let err = (got - want).abs();
            stats.record(&format!("{slug}|{col}"), err, || {
                format!("slug={slug} ts={ts_s} TL={tl} col={col} got={got:.6} want={want:.6}")
            });
            let tol = golden::irradiance_tol(want);
            if err > tol {
                n_tol_violations += 1;
                if tol_violations.len() < 20 {
                    tol_violations.push(format!(
                        "slug={slug} ts={ts_s} TL={tl} col={col} got={got:.6} want={want:.6} \
                         err={err:.6} tol={tol:.6}"
                    ));
                }
            }
        }
    }

    stats.report(
        "golden_clearsky_ineichen: ghi/dni/dhi vs pvlib (irradiance tolerance: max(1 W/m^2, 0.5%))",
    );
    if !night_violations.is_empty() {
        eprintln!(
            "  NIGHT VIOLATIONS ({} rows nonzero when fixture says exactly 0):",
            night_violations.len()
        );
        for v in night_violations.iter().take(10) {
            eprintln!("    {v}");
        }
    }
    if n_tol_violations > 0 {
        eprintln!(
            "  TOLERANCE VIOLATIONS: {n_tol_violations} (row, component) pairs exceeded max(1 W/m^2, 0.5% of want); first {}:",
            tol_violations.len()
        );
        for v in &tol_violations {
            eprintln!("    {v}");
        }
    }

    assert!(
        night_violations.is_empty(),
        "{} rows: ineichen produced a nonzero component at night where pvlib is exactly 0",
        night_violations.len()
    );
    assert_eq!(
        n_tol_violations, 0,
        "{n_tol_violations} (row, component) pairs exceeded the irradiance tolerance -- see report above"
    );
}
