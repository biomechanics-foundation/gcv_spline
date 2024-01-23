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

    let (coefs, variance, stats) = fit_gcv_spline(&knots, &data, &weights, 3, 1.0).unwrap();
    //println!("Number of coefs: {}", coefs.len());
    println!("Coefs: {:?}", coefs);
    println!("Variance: {}", variance);
    println!("Stats: {:?}", stats);

    let value = evaluate_spline(0, 3, 50.5, &knots, &coefs, 0);
    println!("Evaluation: {}", value)
}