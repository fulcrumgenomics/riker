/// Divide two `u64` values; return 0.0 if the denominator is zero.
#[must_use]
pub fn safe_div(num: u64, denom: u64) -> f64 {
    if denom == 0 { 0.0 } else { num as f64 / denom as f64 }
}

/// Divide two `f64` values; return 0.0 if the denominator is zero or NaN.
#[must_use]
pub fn safe_div_f(num: f64, denom: f64) -> f64 {
    if denom == 0.0 || !denom.is_finite() { 0.0 } else { num / denom }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    // ── safe_div ────────────────────────────────────────────────────────────

    #[test]
    fn test_safe_div_normal() {
        assert_eq!(safe_div(10, 5), 2.0);
    }

    #[test]
    fn test_safe_div_zero_denom() {
        assert_eq!(safe_div(10, 0), 0.0);
    }

    #[test]
    fn test_safe_div_zero_both() {
        assert_eq!(safe_div(0, 0), 0.0);
    }

    // ── safe_div_f ──────────────────────────────────────────────────────────

    #[test]
    fn test_safe_div_f_normal() {
        assert_eq!(safe_div_f(10.0, 5.0), 2.0);
    }

    #[test]
    fn test_safe_div_f_zero_denom() {
        assert_eq!(safe_div_f(10.0, 0.0), 0.0);
    }

    #[test]
    fn test_safe_div_f_nan_denom() {
        assert_eq!(safe_div_f(10.0, f64::NAN), 0.0);
    }

    #[test]
    fn test_safe_div_f_inf_denom() {
        assert_eq!(safe_div_f(10.0, f64::INFINITY), 0.0);
    }
}
