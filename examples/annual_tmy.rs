//! End-to-end batch simulation of a full TMY year.
//!
//! Synthesises a cloudless-day GHI curve for every hour of 2024, runs it
//! through `BatchModelChain`, and prints aggregate energy statistics.
//! Useful as a sanity check and a quick benchmark (`cargo run --release
//! --example annual_tmy`).
//!
//! This example is intentionally self-contained — it does not hit the
//! network and does not require an EPW/TMY3 file on disk.

use chrono::{Duration, TimeZone};
use chrono_tz::US::Mountain;
use pvlib::batch::{BatchModelChain, WeatherSeries};
use pvlib::irradiance::DiffuseModel;
use pvlib::Location;

fn main() {
    // Golden, Colorado (NREL HQ).
    let location =
        Location::try_new(39.74, -105.18, Mountain, 1830.0, "Golden, CO").expect("valid location");

    let mc = BatchModelChain::pvwatts(location.clone(), 30.0, 180.0, 5_000.0)
        .with_gamma_pdc(-0.004)
        .with_inverter(5_000.0, 0.96)
        .with_albedo(0.2)
        .with_transposition(DiffuseModel::Perez)
        .with_auto_decomposition(true)
        .with_system_losses(0.14);

    // One hourly timestamp for every hour of 2024 (8784 — leap year).
    let start = Mountain.with_ymd_and_hms(2024, 1, 1, 0, 0, 0).unwrap();
    let n_hours = 366 * 24;
    let times: Vec<_> = (0..n_hours)
        .map(|h| start + Duration::hours(h))
        .collect();

    // Synthesise a very rough GHI profile: a half-sinusoid peaking at solar
    // noon, scaled by a crude seasonal factor. This is NOT a realistic TMY
    // curve — it is a smoke-test signal for the pipeline.
    let ghi: Vec<f64> = times
        .iter()
        .map(|t| {
            use chrono::{Datelike, Timelike};
            let doy = t.ordinal() as f64;
            let hod = t.hour() as f64 + t.minute() as f64 / 60.0;
            let seasonal = (std::f64::consts::TAU * (doy - 80.0) / 365.25).sin() * 0.3 + 0.9;
            let diurnal = ((hod - 12.0) * std::f64::consts::PI / 12.0).cos().max(0.0);
            (950.0 * seasonal * diurnal).max(0.0)
        })
        .collect();

    let zeros = vec![0.0; n_hours as usize];
    let weather = WeatherSeries {
        times,
        ghi,
        dni: zeros.clone(),
        dhi: zeros.clone(),
        temp_air: vec![15.0; n_hours as usize],
        wind_speed: vec![2.0; n_hours as usize],
        albedo: None,
    };

    let t0 = std::time::Instant::now();
    let results = mc.run(&weather).expect("batch run");
    let elapsed = t0.elapsed();

    println!("pvlib-rust — annual TMY simulation");
    println!("  timesteps:       {}", results.ac_power.len());
    println!("  wall clock:      {:.2?}", elapsed);
    println!("  annual energy:   {:.1} kWh", results.total_energy_wh() / 1000.0);
    println!("  peak AC power:   {:.0} W", results.peak_power());
    println!(
        "  capacity factor: {:.1}%",
        results.capacity_factor(5000.0) * 100.0
    );
    println!("  nan timesteps:   {}", results.nan_count());
}
