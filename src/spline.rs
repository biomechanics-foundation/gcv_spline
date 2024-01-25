use crate::woltring::gcvspl::fit_gcv_spline;
use crate::woltring::splder::evaluate_spline;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GcvSpline {
    knots: Vec<f64>,
    coefficients: Vec<f64>,
    half_order: usize,
}

impl GcvSpline {
    pub fn new() -> Self {
        GcvSpline {
            knots: vec![0.; 2],
            coefficients: vec![0.; 2],
            half_order: 1
        }
    }

    pub fn from_data(time: &Vec<f64>, data: &Vec<f64>) -> Self {
        Self::from_full_parameters(time, data, &vec![1.; time.len()], 3, 0.)
    }

    pub fn from_data_and_half_order(time: &Vec<f64>, data: &Vec<f64>, half_order: usize) -> Self {
        Self::from_full_parameters(time, data, &vec![1.; time.len()], half_order, 0.)
    }

    pub fn from_data_half_order_and_smoothing(time: &Vec<f64>, data: &Vec<f64>, half_order: usize, error_variance: f64)
        -> Self {
        Self::from_full_parameters(time, data, &vec![1.; time.len()], half_order, error_variance)
    }

    pub fn from_full_parameters(time: &Vec<f64>, data: &Vec<f64>, weights: &Vec<f64>, half_order: usize,
                                error_variance: f64) -> Self {
        GcvSpline {
            knots: time.clone(),
            coefficients: fit_gcv_spline(time, data, weights, half_order, error_variance)?,
            half_order
        }
    }

    pub fn single_point(&self, point: f64) -> f64{
        self.point_derivative(point, 0)
    }

    pub fn points(&self, points: Vec<f64>) -> Vec<f64> {
        let result = vec![0.; points.len()];
        for i in 0 .. points.len() {
            result[i] = self.single_point(points[i]);
        }
        result
    }

    pub fn point_derivative(&self, point: f64, derivative_order: usize) -> f64 {
        let end = *self.knots.last().expect("Time cannot be empty");
        let start = *self.knots.first().expect("Time cannot be empty");
        let knot_guess = ((point - start) / (end - start) * *self.knots.len() as f64) as usize;

        evaluate_spline(derivative_order, self.half_order, point, &self.knots, &self.coefficients, knot_guess)
    }

    pub fn derivative(&self, points: Vec<f64>, derivative_order: usize) -> Vec<f64> {
        let result = vec![0.; points.len()];
        for i in 0 .. points.len() {
            result[i] = self.point_derivative(points[i], derivative_order);
        }
        result
    }

    pub fn first_derivative(&self, points: Vec<f64>) -> Vec<f64> {
        self.derivative(points, 1)
    }

    pub fn second_derivative(&self, points: Vec<f64>) -> Vec<f64> {
        self.derivative(points, 2)
    }

    pub fn third_derivative(&self, points: Vec<f64>) -> Vec<f64> {
        self.derivative(points, 3)
    }

    pub fn time(&self) -> Vec<f64> {
        self.knots.clone()
    }

    pub fn knots(&self) -> Vec<f64> {
        self.time()
    }
}

impl Default for GcvSpline {
    fn default() -> Self {
        GcvSpline::new()
    }
}
