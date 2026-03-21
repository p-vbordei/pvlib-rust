use pvlib::iam::*;

// --- physical ---

#[test]
fn test_physical_normal_incidence() {
    // At normal incidence (aoi=0), IAM should be very close to 1.0
    let iam = physical(0.0, 1.526, 4.0, 0.002);
    assert!((iam - 1.0).abs() < 0.01, "Physical IAM at normal incidence should be ~1.0, got {}", iam);
}

#[test]
fn test_physical_high_aoi() {
    // At high AOI, IAM should decrease
    let iam_0 = physical(0.0, 1.526, 4.0, 0.002);
    let iam_60 = physical(60.0, 1.526, 4.0, 0.002);
    let iam_80 = physical(80.0, 1.526, 4.0, 0.002);
    assert!(iam_60 < iam_0, "IAM should decrease with AOI");
    assert!(iam_80 < iam_60, "IAM should decrease further at 80 deg");
    assert!(iam_80 > 0.0, "IAM should still be positive at 80 deg");
}

#[test]
fn test_physical_at_90() {
    let iam = physical(90.0, 1.526, 4.0, 0.002);
    assert_eq!(iam, 0.0, "Physical IAM at 90 deg should be 0");
}

#[test]
fn test_physical_nan() {
    let iam = physical(f64::NAN, 1.526, 4.0, 0.002);
    assert!(iam.is_nan(), "Physical IAM for NaN input should be NaN");
}

// --- schlick ---

#[test]
fn test_schlick_normal_incidence() {
    let iam = schlick(0.0);
    assert!((iam - 1.0).abs() < 1e-10, "Schlick at AOI=0 should be 1.0");
}

#[test]
fn test_schlick_at_90() {
    let iam = schlick(90.0);
    assert_eq!(iam, 0.0, "Schlick at AOI=90 should be 0");
}

#[test]
fn test_schlick_at_60() {
    // iam = 1 - (1 - cos(60))^5 = 1 - (1 - 0.5)^5 = 1 - 0.03125 = 0.96875
    let iam = schlick(60.0);
    assert!((iam - 0.96875).abs() < 1e-6, "Schlick at 60 deg should be 0.96875, got {}", iam);
}

#[test]
fn test_schlick_decreases_with_aoi() {
    let iam_20 = schlick(20.0);
    let iam_50 = schlick(50.0);
    let iam_80 = schlick(80.0);
    assert!(iam_20 > iam_50, "Schlick should decrease with AOI");
    assert!(iam_50 > iam_80, "Schlick should decrease further");
}

// --- schlick_diffuse ---

#[test]
fn test_schlick_diffuse_horizontal() {
    // For horizontal surface (tilt=0), sky IAM should be high
    let (iam_sky, _iam_ground) = schlick_diffuse(0.001); // avoid exact zero
    assert!(iam_sky > 0.9, "Sky IAM for horizontal should be high, got {}", iam_sky);
}

#[test]
fn test_schlick_diffuse_tilted() {
    // Compare with Python reference values for tilt=20
    let (iam_sky, iam_ground) = schlick_diffuse(20.0);
    assert!((iam_sky - 0.9625).abs() < 0.001,
        "Sky IAM at 20 deg tilt should be ~0.9625, got {}", iam_sky);
    assert!((iam_ground - 0.627).abs() < 0.01,
        "Ground IAM at 20 deg tilt should be ~0.627, got {}", iam_ground);
}

#[test]
fn test_schlick_diffuse_vertical() {
    // At 90 deg tilt, sky and ground should be equal
    let (iam_sky, iam_ground) = schlick_diffuse(90.0);
    assert!((iam_sky - iam_ground).abs() < 0.01,
        "At vertical, sky and ground IAM should be similar: sky={}, ground={}", iam_sky, iam_ground);
}

// --- martin_ruiz_diffuse ---

#[test]
fn test_martin_ruiz_diffuse_tilted() {
    let (iam_sky, iam_ground) = martin_ruiz_diffuse(30.0, 0.16);
    assert!(iam_sky > 0.0 && iam_sky <= 1.0, "Sky IAM should be in (0,1], got {}", iam_sky);
    assert!(iam_ground > 0.0 && iam_ground <= 1.0, "Ground IAM should be in (0,1], got {}", iam_ground);
}

#[test]
fn test_martin_ruiz_diffuse_vertical_symmetry() {
    // At 90 deg, sky and ground should be equal
    let (iam_sky, iam_ground) = martin_ruiz_diffuse(90.0, 0.16);
    assert!((iam_sky - iam_ground).abs() < 0.01,
        "At vertical, sky and ground should be similar: sky={}, ground={}", iam_sky, iam_ground);
}

#[test]
fn test_martin_ruiz_diffuse_complementary() {
    // iam_sky at tilt=30 should equal iam_ground at tilt=150
    let (iam_sky_30, _) = martin_ruiz_diffuse(30.0, 0.16);
    let (_, iam_ground_150) = martin_ruiz_diffuse(150.0, 0.16);
    assert!((iam_sky_30 - iam_ground_150).abs() < 0.01,
        "Complementary property: sky@30={} should equal ground@150={}", iam_sky_30, iam_ground_150);
}

// --- existing functions ---

#[test]
fn test_ashrae_basic() {
    let iam = ashrae(0.0, 0.05);
    assert!((iam - 1.0).abs() < 1e-6);
    let iam_90 = ashrae(90.0, 0.05);
    assert_eq!(iam_90, 0.0);
}

#[test]
fn test_martin_ruiz_basic() {
    let iam = martin_ruiz(0.0, 0.16);
    assert!((iam - 1.0).abs() < 0.01, "Martin-Ruiz at AOI=0 should be ~1.0, got {}", iam);
    let iam_90 = martin_ruiz(90.0, 0.16);
    assert_eq!(iam_90, 0.0);
}
