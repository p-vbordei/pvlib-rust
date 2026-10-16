//! Gate-G1 golden parity: PVWatts DC + AC, per MANIFEST.json
//! ("pvwatts.csv.gz" -> "calls").
//!
//! pvlib-rust has no standalone `pvwatts_dc` free function -- the classic
//! PVWatts DC formula lives in `pvsystem::PVSystem::get_dc_power_total`
//! (`nameplate_dc * (poa_global/1000) * (1 + gamma_pdc*(temp_cell-25))`,
//! summed over a system's arrays, clamped to >= 0 per array). A one-array
//! `PVSystem` with `nameplate_dc = pdc0_dc` reproduces
//! `pvsystem.pvwatts_dc(g, tc, pdc0, gamma_pdc)` exactly for this fixture's
//! grid (the >= 0 clamp never triggers: the grid's most extreme
//! gamma_pdc/temp_cell combination is `1 + (-0.005)*(70-25) = 0.775 > 0`).
//!
//! AC uses the free function `inverter::pvwatts_ac`, fed the fixture's own
//! `pdc` column (pvlib's DC output) and `pdc0_ac_r<ratio>` column as inputs
//! -- this isolates the AC-clipping primitive itself (independent of the
//! ModelChain/BatchModelChain pdc0 wiring bug fixed in commit b86b7ec) and
//! is exactly what "nails the AC-clipping plateau" per MANIFEST.json.
//! Ratio tags are discovered from the CSV header rather than hardcoded, so
//! this doesn't silently stop testing ratios if the fixture is regenerated
//! with different `%.4f` rounding.

#[path = "golden/mod.rs"]
mod golden;

use pvlib::inverter;
use pvlib::pvsystem::{Array, FixedMount, PVSystem};

const ETA_INV_NOM: f64 = 0.9637;
const ETA_INV_REF: f64 = 0.9637;

#[test]
fn pvwatts_dc_and_ac_match_pvlib() {
    let Some(dir) = golden::fixture_dir() else {
        eprintln!("{}", golden::SKIP_MESSAGE);
        return;
    };

    let mut rdr = golden::open_csv_gz(&dir, "pvwatts.csv.gz");
    let headers = rdr.headers().expect("pvwatts.csv.gz headers").clone();
    let idx = |name: &str| {
        headers
            .iter()
            .position(|h| h == name)
            .unwrap_or_else(|| panic!("pvwatts.csv.gz missing column {name}"))
    };
    let (i_g, i_tc, i_pdc0, i_gamma, i_pdc) = (
        idx("g_poa"),
        idx("temp_cell"),
        idx("pdc0_dc"),
        idx("gamma_pdc"),
        idx("pdc"),
    );
    // Discover per-ratio (pdc0_ac_r<tag>, pac_r<tag>) column pairs from the
    // header rather than hardcoding the ratio tags.
    let ratio_cols: Vec<(String, usize, usize)> = headers
        .iter()
        .filter_map(|h| h.strip_prefix("pdc0_ac_r"))
        .map(|tag| {
            let tag = tag.to_string();
            (
                tag.clone(),
                idx(&format!("pdc0_ac_r{tag}")),
                idx(&format!("pac_r{tag}")),
            )
        })
        .collect();
    assert!(
        !ratio_cols.is_empty(),
        "no pdc0_ac_r* / pac_r* column pairs found in pvwatts.csv.gz header"
    );
    eprintln!(
        "  discovered dc/ac ratio tags: {:?}",
        ratio_cols.iter().map(|(t, ..)| t).collect::<Vec<_>>()
    );

    let mut dc_stats = golden::TrackerSet::new();
    let mut ac_stats = golden::TrackerSet::new();
    let mut n_bad = 0usize;
    let mut examples: Vec<String> = Vec::new();

    for rec in rdr.records() {
        let rec = rec.expect("pvwatts.csv.gz record");
        let g_poa: f64 = rec[i_g].parse().unwrap();
        let temp_cell: f64 = rec[i_tc].parse().unwrap();
        let pdc0_dc: f64 = rec[i_pdc0].parse().unwrap();
        let gamma_pdc: f64 = rec[i_gamma].parse().unwrap();
        let pdc_want: f64 = rec[i_pdc].parse().unwrap();

        let array = Array {
            mount: Box::new(FixedMount {
                surface_tilt: 0.0,
                surface_azimuth: 0.0,
            }),
            nameplate_dc: pdc0_dc,
            gamma_pdc,
            modules_per_string: 1,
            strings: 1,
            albedo: 0.0,
        };
        let system = PVSystem::new(vec![array], 0.0);
        let pdc_got = system.get_dc_power_total(g_poa, temp_cell);

        let dc_err = (pdc_got - pdc_want).abs();
        let key = format!("gamma={gamma_pdc}|pdc");
        dc_stats.record(&key, dc_err, || {
            format!(
                "g_poa={g_poa} temp_cell={temp_cell} gamma_pdc={gamma_pdc} col=pdc got={pdc_got:.9} want={pdc_want:.9}"
            )
        });
        if dc_err > golden::power_tol(pdc_want) {
            n_bad += 1;
            if examples.len() < 10 {
                examples.push(format!(
                    "col=pdc g_poa={g_poa} temp_cell={temp_cell} gamma_pdc={gamma_pdc} \
                     got={pdc_got:.9} want={pdc_want:.9} err={dc_err:.9}"
                ));
            }
        }

        // AC: feed the fixture's OWN pdc (pvlib's DC output), not pdc_got --
        // isolates the AC primitive per the task's "never chain Rust outputs
        // into the next stage" rule.
        for (tag, i_pdc0ac, i_pac) in &ratio_cols {
            let pdc0_ac: f64 = rec[*i_pdc0ac].parse().unwrap();
            let pac_want: f64 = rec[*i_pac].parse().unwrap();
            let pac_got = inverter::pvwatts_ac(pdc_want, pdc0_ac, ETA_INV_NOM, ETA_INV_REF);
            let ac_err = (pac_got - pac_want).abs();
            let key = format!("ratio={tag}|pac");
            ac_stats.record(&key, ac_err, || {
                format!(
                    "g_poa={g_poa} temp_cell={temp_cell} gamma_pdc={gamma_pdc} ratio={tag} \
                     pdc={pdc_want:.9} pdc0_ac={pdc0_ac:.9} col=pac got={pac_got:.9} want={pac_want:.9}"
                )
            });
            if ac_err > golden::power_tol(pac_want) {
                n_bad += 1;
                if examples.len() < 10 {
                    examples.push(format!(
                        "col=pac ratio={tag} g_poa={g_poa} temp_cell={temp_cell} \
                         got={pac_got:.9} want={pac_want:.9} err={ac_err:.9}"
                    ));
                }
            }
        }
    }

    dc_stats.report("golden_pvwatts: pvwatts_dc (via PVSystem::get_dc_power_total) vs pvlib (0.1% power tolerance)");
    ac_stats.report("golden_pvwatts: pvwatts_ac (via inverter::pvwatts_ac, fed fixture's own pdc) vs pvlib (0.1% power tolerance)");
    if n_bad > 0 {
        eprintln!("  TOLERANCE VIOLATIONS: {n_bad}; first {}:", examples.len());
        for e in &examples {
            eprintln!("    {e}");
        }
    }
    assert_eq!(
        n_bad, 0,
        "{n_bad} (row, column) pairs exceeded the 0.1% power tolerance -- see report above"
    );
}
