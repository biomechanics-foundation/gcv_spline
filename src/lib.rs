pub mod gcv_spline;

mod basis;
mod support;
mod prep;
mod bandet;
mod bansol;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vec() {
        let result = gcv_spline::test_vec(-3.3, 5f64);
        assert_eq!(result, vec![-3.3, 5.0])
    }
}
