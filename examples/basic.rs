use gcv_spline::GcvSpline;

fn main() {
    let vec: Vec<f64> = vec![5., 6., 2., 4., 3., 6., 8., 9.];
    let time: Vec<f64> = vec![0., 1., 2., 3., 4., 5., 6., 7.];

    let spline = GcvSpline::from_data(&time, &vec).unwrap();

    let new_time = vec![0., 2.5, 3., 7.];
    let evaluated = spline.points(&new_time);

    dbg!(evaluated);
}