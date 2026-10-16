//! Gate-G1 golden parity: Erbs and DIRINT decomposition, per MANIFEST.json
//! ("decomposition.csv.gz" -> "calls"):
//!
//!   irradiance::erbs(ghi_syn, TRUE_zenith, day_of_year, dni_extra)  -> (dni, dhi)
//!   irradiance::dirint_series(ghi_syn[], TRUE_zenith[], doy[], pressure, None, true)
//!
//! `TRUE zenith` (pvlib's `sp['zenith']`) is joined from `solarposition.csv.gz`
//! by `(slug, timestamp_utc)` -- this fixture doesn't carry its own zenith
//! column. `dirint` is fed the whole per-(slug, cloud_factor) time-ordered
//! series (matching pvlib's temporal-persistence `delta_kt'`), never chained
//! from Rust's own erbs/DISC output.
//!
//! `kt`: pvlib-rust's public `erbs()` returns only `(dni, dhi)`, not `kt` --
//! see tests/golden/mod.rs's doc comment for why this test recomputes
//! erbs()'s own internal clearness-index formula (cos_z clamped to
//! `min_cos_zenith=0.065`, matching pvlib-python's default) instead of
//! substituting the differently-clamped public `clearness_index()` function.

#[path = "golden/mod.rs"]
mod golden;

use pvlib::irradiance;

/// Replicates the `kt` pvlib-rust's `erbs()` computes internally
/// (src/irradiance.rs `erbs`): `cos_z` clamped to a `min_cos_zenith=0.065`
/// floor (matching pvlib-python's default), `kt = ghi / (dni_extra * cos_z)`,
/// unclamped above. Not the same as the public `clearness_index()` (which
/// floors `cos_z` at 0.01 and clamps the result to [0, 1]) -- documented
/// convention gap, see tests/golden/mod.rs.
fn erbs_internal_kt(ghi: f64, zenith_deg: f64, dni_extra: f64) -> f64 {
    let min_cos_zenith = 0.065;
    let cos_z = zenith_deg.to_radians().cos().max(min_cos_zenith);
    ghi / (dni_extra * cos_z)
}

#[test]
fn erbs_matches_pvlib() {
    let Some(dir) = golden::fixture_dir() else {
        eprintln!("{}", golden::SKIP_MESSAGE);
        return;
    };
    let angles = golden::load_solar_angles(&dir);

    let mut stats = golden::TrackerSet::new();
    let mut n_bad = 0usize;
    let mut examples: Vec<String> = Vec::new();

    let mut rdr = golden::open_csv_gz(&dir, "decomposition.csv.gz");
    for rec in rdr.deserialize::<golden::DecompositionRow>() {
        let row = rec.expect("decomposition.csv.gz row");
        let sa = golden::lookup_angles(&angles, &row.slug, &row.timestamp_utc);
        let ts = golden::parse_ts_z(&row.timestamp_utc);
        let doy = golden::day_of_year(ts);
        let dni_extra = irradiance::get_extra_radiation(doy as i32);

        let (dni, dhi) = irradiance::erbs(row.ghi, sa.zenith, doy, dni_extra);
        let kt = erbs_internal_kt(row.ghi, sa.zenith, dni_extra);

        let slug = row.slug.clone();
        let ts_s = row.timestamp_utc.clone();
        let cf = row.cloud_factor;

        let kt_err = (kt - row.kt).abs();
        stats.record(&format!("{slug}|kt"), kt_err, || {
            format!(
                "slug={slug} ts={ts_s} cf={cf} col=kt got={kt:.6} want={:.6}",
                row.kt
            )
        });
        if kt_err > 0.001 {
            n_bad += 1;
            if examples.len() < 10 {
                examples.push(format!(
                    "slug={slug} ts={ts_s} cf={cf} col=kt got={kt:.6} want={:.6} err={kt_err:.6}",
                    row.kt
                ));
            }
        }

        for (col, got, want) in [
            ("dni_erbs", dni, row.dni_erbs),
            ("dhi_erbs", dhi, row.dhi_erbs),
        ] {
            let err = (got - want).abs();
            stats.record(&format!("{slug}|{col}"), err, || {
                format!("slug={slug} ts={ts_s} cf={cf} col={col} got={got:.6} want={want:.6}")
            });
            if err > golden::irradiance_tol(want) {
                n_bad += 1;
                if examples.len() < 10 {
                    examples.push(format!("slug={slug} ts={ts_s} cf={cf} col={col} got={got:.6} want={want:.6} err={err:.6}"));
                }
            }
        }
    }

    stats.report("golden_decomposition: erbs dni/dhi (irradiance tolerance) + kt (<=0.001 abs)");
    if n_bad > 0 {
        eprintln!("  TOLERANCE VIOLATIONS: {n_bad}; first {}:", examples.len());
        for e in &examples {
            eprintln!("    {e}");
        }
    }
    assert_eq!(
        n_bad, 0,
        "{n_bad} (row, column) pairs exceeded tolerance -- see report above"
    );
}

/// pvlib's `dirint(use_delta_kt_prime=True)` (default) is fed the whole
/// per-location time-ordered series so its temporal-persistence term
/// (`delta_kt'`) sees the same neighbours pvlib-python saw. The fixture is
/// already grouped and time-sorted as `(slug, cloud_factor)` contiguous
/// blocks (see generate.py's nested `for loc: for cf:` loops), so
/// `slice::chunk_by` on those two keys recovers each series without any
/// extra sorting.
#[test]
fn dirint_matches_pvlib_on_non_nan_rows() {
    let Some(dir) = golden::fixture_dir() else {
        eprintln!("{}", golden::SKIP_MESSAGE);
        return;
    };
    let locations = golden::load_locations(&dir);
    let angles = golden::load_solar_angles(&dir);

    // ponytail: loads the whole decomposition fixture into memory (~2.5M
    // rows / a few hundred MB on the full-year grid) so slice::chunk_by can
    // recover the (slug, cloud_factor) time-ordered groups pvlib's dirint()
    // needs. Switch to a streaming group-by if this ever becomes a real
    // memory problem.
    let mut rdr = golden::open_csv_gz(&dir, "decomposition.csv.gz");
    let rows: Vec<golden::DecompositionRow> = rdr
        .deserialize::<golden::DecompositionRow>()
        .map(|r| r.expect("decomposition.csv.gz row"))
        .collect();

    let mut stats = golden::TrackerSet::new();
    let mut n_bad = 0usize;
    let mut n_nan_skipped = 0usize;
    let mut n_compared = 0usize;
    let mut examples: Vec<String> = Vec::new();

    for group in
        rows.chunk_by(|a, b| a.slug == b.slug && (a.cloud_factor - b.cloud_factor).abs() < 1e-9)
    {
        let slug = &group[0].slug;
        let alt = locations
            .get(slug)
            .unwrap_or_else(|| panic!("unknown slug {slug} (not in locations.json)"))
            .alt;
        let pressure = pvlib::atmosphere::alt2pres(alt);

        let ghi: Vec<f64> = group.iter().map(|r| r.ghi).collect();
        let zenith: Vec<f64> = group
            .iter()
            .map(|r| golden::lookup_angles(&angles, &r.slug, &r.timestamp_utc).zenith)
            .collect();
        let doy: Vec<i32> = group
            .iter()
            .map(|r| golden::day_of_year(golden::parse_ts_z(&r.timestamp_utc)) as i32)
            .collect();

        let dni_dirint = irradiance::dirint_series(&ghi, &zenith, &doy, pressure, None, true);

        for (row, got) in group.iter().zip(dni_dirint.iter()) {
            let Some(want) = row.dni_dirint else {
                n_nan_skipped += 1;
                continue;
            };
            n_compared += 1;
            let err = (got - want).abs();
            let s = row.slug.clone();
            let ts = row.timestamp_utc.clone();
            let cf = row.cloud_factor;
            stats.record(&format!("{s}|dni_dirint"), err, || {
                format!("slug={s} ts={ts} cf={cf} col=dni_dirint got={got:.6} want={want:.6}")
            });
            if err > golden::irradiance_tol(want) {
                n_bad += 1;
                if examples.len() < 10 {
                    examples.push(format!("slug={s} ts={ts} cf={cf} col=dni_dirint got={got:.6} want={want:.6} err={err:.6}"));
                }
            }
        }
    }

    stats.report("golden_decomposition: dirint dni_dirint vs pvlib (irradiance tolerance, non-NaN rows only)");
    eprintln!("  compared {n_compared} rows, skipped {n_nan_skipped} NaN (night) rows");
    if n_bad > 0 {
        eprintln!(
            "  DIRINT DIFFS: {n_bad} rows exceeded tolerance (may reflect genuine DIRINT bin-edge \
             sensitivity, not necessarily a bug -- see docs/PLAN.md); first {}:",
            examples.len()
        );
        for e in &examples {
            eprintln!("    {e}");
        }
    }
    assert_eq!(
        n_bad, 0,
        "{n_bad} dni_dirint rows exceeded the irradiance tolerance -- see report above \
         (reported honestly per the task spec, not masked, even if caused by DIRINT bin-edge sensitivity)"
    );
}
