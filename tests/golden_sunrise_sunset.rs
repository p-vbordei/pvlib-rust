//! Gate-G1 golden parity: `solarposition::sun_rise_set_transit_geometric`
//! (Spencer71 declination + equation of time) vs pinned pvlib-python
//! `sun_rise_set_transit_spa` (full NREL SPA) fixtures.
//!
//! Tolerance per docs/PLAN.md Etapa 1 / MANIFEST.json "day_night_boundary":
//! sunrise/sunset/transit may differ by up to one 15-min slot (900 s)
//! between the simplified geometric method and the full SPA method (and,
//! separately, between SPA and ephem's refraction model) -- documented, not
//! masked. ephem comparison is reported only, never asserted.

#[path = "golden/mod.rs"]
mod golden;

use chrono::{Datelike, NaiveDate, TimeZone};
use pvlib::solarposition::{
    declination_spencer71, equation_of_time_spencer71, sun_rise_set_transit_geometric,
};

/// Rust's geometric sunrise/sunset/transit, expressed as a UTC `DateTime`,
/// for a given UTC calendar date and location. `None` on polar day/night
/// (not expected for this fleet's latitude band).
fn rust_rise_set_transit(
    date: NaiveDate,
    lat: f64,
    lon: f64,
) -> Option<(
    chrono::DateTime<chrono::Utc>,
    chrono::DateTime<chrono::Utc>,
    chrono::DateTime<chrono::Utc>,
)> {
    let doy = date.ordinal() as f64;
    let declination = declination_spencer71(doy);
    let eot = equation_of_time_spencer71(doy);
    let rst = sun_rise_set_transit_geometric(lat, lon, declination, eot, 0.0)?;
    let midnight = chrono::Utc.from_utc_datetime(&date.and_hms_opt(0, 0, 0).unwrap());
    let to_dt = |hours: f64| {
        midnight + chrono::Duration::milliseconds((hours * 3_600_000.0).round() as i64)
    };
    Some((to_dt(rst.sunrise), to_dt(rst.sunset), to_dt(rst.transit)))
}

#[test]
fn geometric_sunrise_sunset_matches_spa_within_one_slot() {
    let Some(dir) = golden::fixture_dir() else {
        eprintln!("{}", golden::SKIP_MESSAGE);
        return;
    };
    let locations = golden::load_locations(&dir);

    let mut spa = golden::TrackerSet::new(); // asserted vs SPA
    let mut ephem = golden::TrackerSet::new(); // reported only, vs ephem

    let mut rdr = golden::open_csv_gz(&dir, "sunrise_sunset.csv.gz");
    let mut n_polar = 0usize;
    for rec in rdr.deserialize::<golden::SunriseSunsetRow>() {
        let row = rec.expect("sunrise_sunset.csv.gz row");
        let loc = locations
            .get(&row.slug)
            .unwrap_or_else(|| panic!("unknown slug {} (not in locations.json)", row.slug));
        let date = NaiveDate::parse_from_str(&row.date, "%Y-%m-%d")
            .unwrap_or_else(|e| panic!("bad date {:?}: {e}", row.date));

        let Some((rise, set, transit)) = rust_rise_set_transit(date, loc.lat, loc.lon) else {
            n_polar += 1;
            continue;
        };

        let want_rise = golden::parse_rfc3339_utc(&row.sunrise_spa);
        let want_set = golden::parse_rfc3339_utc(&row.sunset_spa);
        let want_transit = golden::parse_rfc3339_utc(&row.transit_spa);

        let secs = |a: chrono::DateTime<chrono::Utc>, b: chrono::DateTime<chrono::Utc>| {
            ((a - b).num_milliseconds() as f64 / 1000.0).abs()
        };

        let slug = row.slug.clone();
        let date_s = row.date.clone();
        spa.record(&format!("{slug}|sunrise"), secs(rise, want_rise), || {
            format!("slug={slug} date={date_s} col=sunrise got={rise} want={want_rise}")
        });
        spa.record(&format!("{slug}|sunset"), secs(set, want_set), || {
            format!("slug={slug} date={date_s} col=sunset got={set} want={want_set}")
        });
        spa.record(
            &format!("{slug}|transit"),
            secs(transit, want_transit),
            || format!("slug={slug} date={date_s} col=transit got={transit} want={want_transit}"),
        );

        // Informational only: ephem uses its own default-refraction horizon,
        // never asserted against.
        if !row.sunrise_ephem.is_empty() {
            let want_rise_ephem = golden::parse_rfc3339_utc(&row.sunrise_ephem);
            ephem.record(&format!("{slug}|sunrise_vs_ephem"), secs(rise, want_rise_ephem), || {
                format!("slug={slug} date={date_s} col=sunrise_vs_ephem got={rise} want={want_rise_ephem}")
            });
        }
        if !row.sunset_ephem.is_empty() {
            let want_set_ephem = golden::parse_rfc3339_utc(&row.sunset_ephem);
            ephem.record(&format!("{slug}|sunset_vs_ephem"), secs(set, want_set_ephem), || {
                format!("slug={slug} date={date_s} col=sunset_vs_ephem got={set} want={want_set_ephem}")
            });
        }
    }

    spa.report("golden_sunrise_sunset: geometric vs SPA (asserted, +/- 900s slot rule)");
    ephem.report("golden_sunrise_sunset: geometric vs ephem (informational only, not asserted)");
    if n_polar > 0 {
        eprintln!(
            "  ({n_polar} rows skipped: sun_rise_set_transit_geometric returned None -- polar day/night)"
        );
    }

    let max_err = spa.global_max();
    assert!(
        max_err <= golden::SUNRISE_SUNSET_TOL_SECONDS,
        "geometric sunrise/sunset/transit exceeded the +/-{} s (one 15-min slot) tolerance \
         vs SPA: max={:.1} s (see report above)",
        golden::SUNRISE_SUNSET_TOL_SECONDS,
        max_err
    );
}
