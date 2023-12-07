use crate::gcvspl::fit_gcv_spline;
use crate::splc::fit_spline_coefficients_with_stats;

pub mod gcv_spline;

mod basis;
mod support;
mod prep;
mod bandet;
mod bansol;
mod trinv;
mod splc;
mod gcvspl;
mod search;
mod splder;


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vec() {
        let result = gcv_spline::test_vec(-3.3, 5f64);
        assert_eq!(result, vec![-3.3, 5.0])
    }
}
