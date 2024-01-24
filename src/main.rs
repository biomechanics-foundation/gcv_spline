use crate::gcvspl::fit_gcv_spline;
use crate::splder::evaluate_spline;

pub mod gcv_spline;

mod basis;
mod support;
mod prep;
mod bandet;
mod bansol;
mod trinv;
mod splc;
mod gcvspl;
mod splder;
mod search;


pub fn main() {
    let knots: Vec<f64> = (0..=100).map(|e| e as f64).collect();
    let data = knots.clone().iter().map(|e| (e * 0.01).sin()).collect();
    let weights = vec![1.0; knots.len()];

    let coefs = fit_gcv_spline(&knots, &data, &weights, 3, 0.0).unwrap();
    //println!("Number of coefs: {}", coefs.len());
    println!("Coefs: {:?}", coefs);

    let value = evaluate_spline(0, 2, 50.5, &knots, &coefs, 0);
    println!("Evaluation: {}", value);
    //println!("{}", 0.01 * 0.505_f64.cos())
}