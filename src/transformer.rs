/// Calculate transformer output power accounting for efficiency losses.
///
/// Uses a quadratic model where load losses scale with the square of output power.
///
/// # Arguments
/// * `input_power` - Real AC power input to the transformer in W.
/// * `transformer_rating` - Nominal output power of the transformer in VA.
/// * `no_load_loss` - Constant losses as fraction of rating (0-1).
/// * `load_loss` - Load-dependent losses as fraction of rating (0-1).
///
/// # Returns
/// Real AC power output in W.
#[inline]
pub fn simple_efficiency(
    input_power: f64,
    transformer_rating: f64,
    no_load_loss: f64,
    load_loss: f64,
) -> f64 {
    let input_normalized = input_power / transformer_rating;

    let a = load_loss;
    let b = 1.0;
    let c = no_load_loss - input_normalized;

    // Alternative quadratic formula to avoid divide-by-zero when a == 0
    let disc = (b * b - 4.0 * a * c).sqrt();
    let output_normalized = 2.0 * c / (-b - disc);

    output_normalized * transformer_rating
}
