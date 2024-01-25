use gcv_spline::{GcvSpline, FittingError};

fn main() {
    let vec: Vec<f64> = vec![5., 6., 2., 4., 3., 6., 8., 9.];
    let time: Vec<f64> = vec![0., 1., 2., 3., 4., 5., 6., 7.];

    let spline = GcvSpline::from_data(&time, &vec).unwrap();
    let evaluated: Vec<f32> = spline.points(&time).iter().map(|e| *e as f32).collect();
    dbg!(evaluated[1] as f64);
}