use crate::woltring::gcvspl::fit_gcv_spline;
use crate::woltring::splder::evaluate_spline;
use crate::woltring::support::FittingError;

/// Represents a GCV spline fitted to provided data.
#[derive(Clone, Debug, PartialEq)]
pub struct GcvSpline {
    knots: Vec<f64>,
    coefficients: Vec<f64>,
    half_order: usize,
}

impl GcvSpline {

    /// Fits a GcvSpline from user-provided time (knots) and data vectors. This method uses default values for
    /// half-order, error variance, and weights, which are generally applicable. This method should be used in most
    /// cases.
    pub fn from_data(time: &Vec<f64>, data: &Vec<f64>) -> Result<Self, FittingError> {
        Self::from_full_parameters(time, data, &vec![1.; time.len()], 3, 0.)
    }

    /// Fits a GcvSpline from user-provided time (knots) and data vectors and a half-order. A half-order *m* will
    /// produce a GcvSpline with a degree of 2 * *m* - 1.
    pub fn from_data_and_half_order(time: &Vec<f64>, data: &Vec<f64>, half_order: usize) -> Result<Self, FittingError> {
        Self::from_full_parameters(time, data, &vec![1.; time.len()], half_order, 0.)
    }

    /// Fits a GcvSpline from user-provided time (knots) and data vectors, a half-order, and an error variance.
    /// Increasing the provided error variance will introduce smoothing into the spline. With higher error variance
    /// values, input data points will not be matched as closely. Generally, it is recommended to low-pass filter noisy
    /// input data before fitting with an error variance of 0 rather than using this type of smoothing.
    pub fn from_data_half_order_and_smoothing(time: &Vec<f64>, data: &Vec<f64>, half_order: usize, error_variance: f64)
        -> Result<Self, FittingError> {
        Self::from_full_parameters(time, data, &vec![1.; time.len()], half_order, error_variance)
    }

    /// Fits a GcvSpline from user-provided time (knots) and data vectors, a half-order, an error variance, and a
    /// vector of weights. The weights define how important individual fitting points are.
    pub fn from_full_parameters(time: &Vec<f64>, data: &Vec<f64>, weights: &Vec<f64>, half_order: usize,
                                error_variance: f64) -> Result<Self, FittingError> {
        let coefficients = fit_gcv_spline(time, data, weights, half_order, error_variance)?;
        Ok(GcvSpline {
            knots: time.clone(),
            coefficients,
            half_order
        })
    }

    /// Creates a GcvSpline with default values. This does not describe any user-provided data.
    pub fn new() -> Self {
        GcvSpline {
            knots: vec![0., 1.],
            coefficients: vec![0.; 2],
            half_order: 1
        }
    }

    /// Evaluates a GCV spline at a single point.
    pub fn single_point(&self, point: f64) -> f64{
        self.point_derivative(point, 0)
    }

    /// Evaluates a GCV spline at a set of points.
    pub fn points(&self, points: &Vec<f64>) -> Vec<f64> {
        let mut result = vec![0.; points.len()];
        for i in 0 .. points.len() {
            result[i] = self.single_point(points[i]);
        }
        result
    }

    /// Evaluates a derivative of a given order at a single point.
    pub fn point_derivative(&self, point: f64, derivative_order: usize) -> f64 {
        let end = *self.knots.last().expect("Time cannot be empty");
        let start = *self.knots.first().expect("Time cannot be empty");
        let knot_guess = ((point - start) / (end - start) * self.knots.len() as f64) as usize;

        evaluate_spline(derivative_order, self.half_order, point, &self.knots, &self.coefficients, knot_guess)
    }

    /// Evaluates a derivative of a given order at a set of points.
    pub fn derivative(&self, points: &Vec<f64>, derivative_order: usize) -> Vec<f64> {
        let mut result = vec![0.; points.len()];
        for i in 0 .. points.len() {
            result[i] = self.point_derivative(points[i], derivative_order);
        }
        result
    }

    /// Evaluates the first derivative at a set of points.
    pub fn first_derivative(&self, points: &Vec<f64>) -> Vec<f64> {
        self.derivative(points, 1)
    }

    /// Evaluates the second derivative at a set of points.
    pub fn second_derivative(&self, points: &Vec<f64>) -> Vec<f64> {
        self.derivative(points, 2)
    }

    /// Evaluates the third derivative at a set of points.
    pub fn third_derivative(&self, points: &Vec<f64>) -> Vec<f64> {
        self.derivative(points, 3)
    }

    /// Returns a copy of the time vector used to fit the GCV spline.
    pub fn time(&self) -> Vec<f64> {
        self.knots.clone()
    }

    /// Returns a copy of the knots vector used to fit the GCV spline.
    pub fn knots(&self) -> Vec<f64> {
        self.time()
    }
}

impl Default for GcvSpline {
    fn default() -> Self {
        GcvSpline::new()
    }
}
