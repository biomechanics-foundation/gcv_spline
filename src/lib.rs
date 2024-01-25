use crate::woltring::{gcvspl::fit_gcv_spline, splder::evaluate_spline};

pub mod spline;

pub use spline::GcvSpline;

pub mod woltring;

pub mod prelude {
    pub use crate::spline::GcvSpline;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sin_eval() {
        let knots: Vec<f64> = (0..=100).map(|e| e as f64).collect();
        let data = knots.iter().map(|e| (e * 0.01).sin()).collect();
        let weights = vec![1.0; knots.len()];

        let coefs = fit_gcv_spline(&knots, &data, &weights, 3, 0.0).unwrap();
        let value = evaluate_spline(0, 3, 50.5, &knots, &coefs, 0);
        assert!((value - 0.505_f64.sin()).abs() < 1e-15)
    }

    #[test]
    fn test_sin_derivative() {
        let knots: Vec<f64> = (0..=100).map(|e| e as f64).collect();
        let data = knots.iter().map(|e| (e * 0.01).sin()).collect();
        let weights = vec![1.0; knots.len()];

        let coefs = fit_gcv_spline(&knots, &data, &weights, 3, 0.0).unwrap();
        let value = evaluate_spline(1, 3, 50.5, &knots, &coefs, 0);
        assert!((value - 0.01 * 0.505_f64.cos()).abs() < 1e-15)
    }

    #[test]
    fn test_sin_second_derivative() {
        let knots: Vec<f64> = (0..=100).map(|e| e as f64).collect();
        let data = knots.iter().map(|e| (e * 0.01).sin()).collect();
        let weights = vec![1.0; knots.len()];

        let coefs = fit_gcv_spline(&knots, &data, &weights, 3, 0.0).unwrap();
        let value = evaluate_spline(2, 3, 50.5, &knots, &coefs, 0);
        assert!((value + 0.01 * 0.01 * 0.505_f64.sin()).abs() < 1e-15)
    }
}
