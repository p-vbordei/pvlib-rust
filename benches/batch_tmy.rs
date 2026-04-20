//! Criterion benchmark: full `BatchModelChain::run` on a TMY-year-sized
//! weather series (8760 hourly rows).
//!
//! Run with `cargo bench --bench batch_tmy`. Reports wall-clock per iteration
//! and confidence intervals on the stable rustc toolchain.

use chrono::{Duration, TimeZone};
use chrono_tz::US::Mountain;
use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use pvlib::batch::{BatchModelChain, WeatherSeries};
use pvlib::irradiance::DiffuseModel;
use pvlib::Location;

fn build_weather(n_hours: usize) -> WeatherSeries {
    let start = Mountain.with_ymd_and_hms(2024, 1, 1, 0, 0, 0).unwrap();
    let times: Vec<_> = (0..n_hours as i64)
        .map(|h| start + Duration::hours(h))
        .collect();
    let ghi: Vec<f64> = times
        .iter()
        .map(|t| {
            use chrono::{Datelike, Timelike};
            let doy = t.ordinal() as f64;
            let hod = t.hour() as f64;
            let seasonal = (std::f64::consts::TAU * (doy - 80.0) / 365.25).sin() * 0.3 + 0.9;
            let diurnal = ((hod - 12.0) * std::f64::consts::PI / 12.0).cos().max(0.0);
            (950.0 * seasonal * diurnal).max(0.0)
        })
        .collect();
    WeatherSeries {
        times,
        ghi,
        dni: vec![0.0; n_hours],
        dhi: vec![0.0; n_hours],
        temp_air: vec![15.0; n_hours],
        wind_speed: vec![2.0; n_hours],
        albedo: None,
    }
}

fn bench_batch_tmy(c: &mut Criterion) {
    let location =
        Location::try_new(39.74, -105.18, Mountain, 1830.0, "Golden, CO").unwrap();
    let mc = BatchModelChain::pvwatts(location, 30.0, 180.0, 5_000.0)
        .with_gamma_pdc(-0.004)
        .with_inverter(5_000.0, 0.96)
        .with_albedo(0.2)
        .with_transposition(DiffuseModel::Perez)
        .with_auto_decomposition(true)
        .with_system_losses(0.14);

    let weather = build_weather(8760);

    let mut group = c.benchmark_group("batch_tmy");
    group.throughput(Throughput::Elements(8760));
    group.bench_function("8760_hours_pvwatts_perez", |b| {
        b.iter(|| {
            let r = mc.run(&weather).expect("batch run");
            criterion::black_box(r.total_energy_wh())
        });
    });
    group.finish();
}

criterion_group!(benches, bench_batch_tmy);
criterion_main!(benches);
