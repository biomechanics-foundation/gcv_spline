use crate::woltring::gcvspl::fit_gcv_spline;
use crate::woltring::splder::evaluate_spline;

pub mod gcv_spline;

mod woltring;


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sin_eval() {
        let knots: Vec<f64> = (0..=100).map(|e| e as f64).collect();
        let data = knots.clone().iter().map(|e| (e * 0.01).sin()).collect();
        let weights = vec![1.0; knots.len()];

        let coefs = fit_gcv_spline(&knots, &data, &weights, 3, 0.0).unwrap();
        let value = evaluate_spline(0, 3, 50.5, &knots, &coefs, 0);
        assert!((value - 0.505_f64.sin()).abs() < 1e-12)
    }

    #[test]
    fn test_sin_derivative() {
        let knots: Vec<f64> = (0..=100).map(|e| e as f64).collect();
        let data = knots.clone().iter().map(|e| (e * 0.01).sin()).collect();
        let weights = vec![1.0; knots.len()];

        let coefs = fit_gcv_spline(&knots, &data, &weights, 3, 0.0).unwrap();
        let value = evaluate_spline(1, 3, 50.5, &knots, &coefs, 0);
        assert!((value - 0.01 * 0.505_f64.cos()).abs() < 1e-12)
    }

    #[test]
    fn test_sin_second_derivative() {
        let knots: Vec<f64> = (0..=100).map(|e| e as f64).collect();
        let data = knots.clone().iter().map(|e| (e * 0.01).sin()).collect();
        let weights = vec![1.0; knots.len()];

        let coefs = fit_gcv_spline(&knots, &data, &weights, 3, 0.0).unwrap();
        let value = evaluate_spline(2, 3, 50.5, &knots, &coefs, 0);
        assert!((value + 0.01 * 0.01 * 0.505_f64.sin()).abs() < 1e-12)
    }
}
