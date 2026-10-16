//! Gate-G1 golden parity: isotropic transposition chain, per
//! MANIFEST.json ("transposition_isotropic.csv.gz" -> "calls"):
//!
//!   aoi              = irradiance::aoi(tilt, saz, solar_zenith, solar_azimuth)
//!   poa_beam         = irradiance::beam_component(tilt, saz, solar_zenith, solar_azimuth, dni)
//!   poa_sky_diffuse  = irradiance::isotropic(tilt, dhi)
//!   poa_ground_diffuse = irradiance::get_ground_diffuse(tilt, ghi, albedo)
//!   poa_global       = poa_beam + poa_sky_diffuse + poa_ground_diffuse
//!
//! `solar_zenith` = SPA apparent_zenith, `solar_azimuth` = SPA azimuth (the
//! pvlib ModelChain convention) -- joined from `solarposition.csv.gz` by
//! `(slug, timestamp_utc)` since this fixture doesn't carry its own angle
//! columns. Every other input (surface_tilt, surface_azimuth, albedo, ghi,
//! dni, dhi) is this fixture's own column.

#[path = "golden/mod.rs"]
mod golden;

use pvlib::irradiance;

#[test]
fn isotropic_transposition_matches_pvlib() {
    let Some(dir) = golden::fixture_dir() else {
        eprintln!("{}", golden::SKIP_MESSAGE);
        return;
    };
    let angles = golden::load_solar_angles(&dir);

    let mut stats = golden::TrackerSet::new();
    let mut n_bad = 0usize;
    let mut examples: Vec<String> = Vec::new();

    let mut rdr = golden::open_csv_gz(&dir, "transposition_isotropic.csv.gz");
    for rec in rdr.deserialize::<golden::TranspositionRow>() {
        let row = rec.expect("transposition_isotropic.csv.gz row");
        let sa = golden::lookup_angles(&angles, &row.slug, &row.timestamp_utc);
        let (solar_zenith, solar_azimuth) = (sa.apparent_zenith, sa.azimuth);

        let aoi = irradiance::aoi(
            row.surface_tilt,
            row.surface_azimuth,
            solar_zenith,
            solar_azimuth,
        );
        let poa_beam = irradiance::beam_component(
            row.surface_tilt,
            row.surface_azimuth,
            solar_zenith,
            solar_azimuth,
            row.dni,
        );
        let poa_sky_diffuse = irradiance::isotropic(row.surface_tilt, row.dhi);
        let poa_ground_diffuse =
            irradiance::get_ground_diffuse(row.surface_tilt, row.ghi, row.albedo);
        let poa_global = poa_beam + poa_sky_diffuse + poa_ground_diffuse;

        let slug = row.slug.clone();
        let ts = row.timestamp_utc.clone();
        // Same numeric tolerance rule (max(1, 0.5%)) applied to every
        // quantity including `aoi`, per the task spec -- aoi is in degrees,
        // not W/m^2, but the instruction is explicit that the irradiance
        // rule reuses this same numeric formula everywhere in this test.
        for (col, got, want) in [
            ("aoi", aoi, row.aoi),
            ("poa_beam", poa_beam, row.poa_beam),
            ("poa_sky_diffuse", poa_sky_diffuse, row.poa_sky_diffuse),
            (
                "poa_ground_diffuse",
                poa_ground_diffuse,
                row.poa_ground_diffuse,
            ),
            ("poa_global", poa_global, row.poa_global),
        ] {
            let err = (got - want).abs();
            stats.record(&format!("{slug}|{col}"), err, || {
                format!(
                    "slug={slug} ts={ts} tilt={} saz={} alb={} col={col} got={got:.6} want={want:.6}",
                    row.surface_tilt, row.surface_azimuth, row.albedo
                )
            });
            if err > golden::irradiance_tol(want) {
                n_bad += 1;
                if examples.len() < 10 {
                    examples.push(format!(
                        "slug={slug} ts={ts} col={col} got={got:.6} want={want:.6} err={err:.6}"
                    ));
                }
            }
        }
    }

    stats.report(
        "golden_transposition: aoi/poa_beam/poa_sky_diffuse/poa_ground_diffuse/poa_global vs pvlib",
    );
    if n_bad > 0 {
        eprintln!("  TOLERANCE VIOLATIONS: {n_bad}; first {}:", examples.len());
        for e in &examples {
            eprintln!("    {e}");
        }
    }
    assert_eq!(
        n_bad, 0,
        "{n_bad} (row, column) pairs exceeded the irradiance tolerance -- see report above"
    );
}
