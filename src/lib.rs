//! # gcv_spline - Generalized Cross-Validated Splines for Interpolation and Derivation in Pure Rust
//! This crate implements the GCV spline, a versatile and easy-to-use spline structure for interpolating data at
//! unknown points and taking accurate derivatives of smooth data.
//!
//! GCV splines were first developed by Herman J. Woltring. The modules inside the private woltring module are based on
//! his [FORTRAN package](https://isbweb.org/software/sigproc/gcvspl/gcvspl.f) and a
//! [C translation](https://isbweb.org/software/sigproc/gcvspl/Twisk/gcvspl.c) by D. Twisk. Comments from these
//! versions are included in this implementation.
//!
//! # Examples
//! ```
//! use gcv_spline::GcvSpline;
//!
//! // Points describe the function y = x**2
//! let time = vec![0., 1., 3., 4., 5., 6.];
//! let values = vec![0., 1., 9., 16., 25., 36.];
//!
//! let spline = GcvSpline::from_data(&time, &values).unwrap();
//! // Interpolate at a missing point.
//! // This should be close to the expected value of 4.
//! assert!((spline.single_point(2.) - 4.).abs() < 1e-12);
//! // Take derivatives at the interpolated point.
//! // First derivative should be close to 4.
//! assert!((spline.point_derivative(2., 1) - 4.).abs() < 1e-12);
//! // Second derivative should be close to 2.
//! assert!((spline.point_derivative(2., 2) - 2.).abs() < 1e-12);
//! ```

pub mod spline;
pub mod woltring;

pub use crate::spline::GcvSpline;
pub use crate::woltring::support::FittingError;

#[cfg(test)]
mod tests {
    use crate::woltring::{gcvspl::fit_gcv_spline, splder::evaluate_spline};
    use crate::GcvSpline;

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

    #[test]
    fn test_default() {
        let spline = GcvSpline::new();
        assert_eq!(spline.knots(), vec![0., 1.]);
    }
}
