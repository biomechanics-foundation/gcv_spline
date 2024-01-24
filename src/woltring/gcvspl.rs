use crate::woltring::basis::create_basis;
use crate::woltring::prep::create_weighted_matrix;
use crate::woltring::splc::fit_spline_coefficients_with_stats;
use crate::woltring::support::{check_increasing, check_order, check_vector_length, FittingError};

pub fn fit_gcv_spline(knots: &Vec<f64>, data: &Vec<f64>, weight_factors: &Vec<f64>,
                      half_order: usize, error_variance: f64)
        -> Result<Vec<f64>, FittingError> {
    let num_knots = knots.len();
    check_order(half_order, num_knots)?;
    check_increasing(knots)?;
    check_vector_length(data, num_knots)?;
    check_vector_length(weight_factors, num_knots)?;

    let smoothing_ratio = 2.0;
    let tau = 1.618033983;
    let epsilon = 1E-15;
    let tolerance = 1E-6;
    let mut coefficients = vec![0.0; num_knots];
    let mut stats = vec![0.0; 6];

    // Compute design matrices and norms
    let (spline_tableau, basis_l1_norm): (Vec<f64>, f64) = create_basis(half_order, &knots)?;
    let (weighted_matrix, mut weighted_matrix_norm): (Vec<f64>, f64)
        = create_weighted_matrix(half_order, &knots, &weight_factors)?;
    weighted_matrix_norm /= basis_l1_norm;

    // Store temporary GCV function values
    let (mut gcv_f1, mut gcv_f2, mut gcv_f3, mut gcv_f4): (f64, f64, f64, f64);
    let (mut smoothing_1, mut smoothing_2, mut smoothing_3, mut smoothing_4): (f64, f64, f64, f64);
    let mut traced_matrix = vec![0.0; weighted_matrix.len()];
    // Zero variance case
    if error_variance == 0.0 {
        smoothing_1 = 0.0; //let (g, c, s, tm)
        let _gcv_f1 = fit_spline_coefficients_with_stats(
            half_order, &data, &weight_factors, error_variance, smoothing_1, epsilon, &spline_tableau,
            &weighted_matrix, weighted_matrix_norm, &mut coefficients, &mut stats,
            &mut traced_matrix)?;
    } else {
        let mut solved = false;
        smoothing_1 = 1.0 / weighted_matrix_norm;
        smoothing_2 = smoothing_1 * smoothing_ratio;
        gcv_f2 = fit_spline_coefficients_with_stats(
            half_order, &data, &weight_factors, error_variance, smoothing_2, epsilon, &spline_tableau,
            &weighted_matrix, weighted_matrix_norm, &mut coefficients, &mut stats,
            &mut traced_matrix)?;
        gcv_f1 = fit_spline_coefficients_with_stats(
            half_order, &data, &weight_factors, error_variance, smoothing_1, epsilon, &spline_tableau,
            &weighted_matrix, weighted_matrix_norm, &mut coefficients, &mut stats,
            &mut traced_matrix)?;
        while gcv_f1 <= gcv_f2 && !solved {
            if stats[3] <= 0.0 {
                solved = true;
            } else {
                smoothing_2 = smoothing_1;
                gcv_f2 = gcv_f1;
                smoothing_1 /= smoothing_ratio;
                gcv_f1 = fit_spline_coefficients_with_stats(
                    half_order, &data, &weight_factors, error_variance, smoothing_1, epsilon, &spline_tableau,
                    &weighted_matrix, weighted_matrix_norm, &mut coefficients, &mut stats,
                    &mut traced_matrix)?;
            }
        }
        if !solved {
            smoothing_3 = smoothing_2 * smoothing_ratio;
            gcv_f3 = fit_spline_coefficients_with_stats(
                half_order, &data, &weight_factors, error_variance, smoothing_3, epsilon, &spline_tableau,
                &weighted_matrix, weighted_matrix_norm, &mut coefficients, &mut stats,
                &mut traced_matrix)?;
            while gcv_f3 <= gcv_f2 && !solved {
                if stats[3] >= 1.0 {
                    solved = true;
                } else {
                    // smoothing_2 = smoothing_3;
                    gcv_f2 = gcv_f3;
                    smoothing_3 *= smoothing_ratio;
                    gcv_f3 = fit_spline_coefficients_with_stats(
                        half_order, &data, &weight_factors, error_variance, smoothing_3, epsilon, &spline_tableau,
                        &weighted_matrix, weighted_matrix_norm, &mut coefficients, &mut stats,
                        &mut traced_matrix)?;
                }
            }
            if !solved {
                smoothing_2 = smoothing_3;
                // gcv_f2 = gcv_f3;
                let mut alpha = (smoothing_2 - smoothing_1) / tau;
                smoothing_4 = smoothing_1 + alpha;
                smoothing_3 = smoothing_2 - alpha;
                gcv_f3 = fit_spline_coefficients_with_stats(
                    half_order, &data, &weight_factors, error_variance, smoothing_3, epsilon, &spline_tableau,
                    &weighted_matrix, weighted_matrix_norm, &mut coefficients, &mut stats,
                    &mut traced_matrix)?;
                gcv_f4 = fit_spline_coefficients_with_stats(
                    half_order, &data, &weight_factors, error_variance, smoothing_4, epsilon, &spline_tableau,
                    &weighted_matrix, weighted_matrix_norm, &mut coefficients, &mut stats,
                    &mut traced_matrix)?;
                while !solved {
                    if gcv_f3 <= gcv_f4 {
                        smoothing_2 = smoothing_4;
                        // gcv_f2 = gcv_f4;
                        let error = (smoothing_2 - smoothing_1) / (smoothing_1 + smoothing_2);
                        if error * error + 1.0 == 1.0 || error <= tolerance {
                            solved = true;
                        } else {
                            smoothing_4 = smoothing_3;
                            gcv_f4 = gcv_f3;
                            alpha = alpha / tau;
                            smoothing_3 = smoothing_2 - alpha;
                            gcv_f3 = fit_spline_coefficients_with_stats(
                                half_order, &data, &weight_factors, error_variance, smoothing_3, epsilon, &spline_tableau,
                                &weighted_matrix, weighted_matrix_norm, &mut coefficients, &mut stats,
                                &mut traced_matrix)?;
                        }
                    } else {
                        smoothing_1 = smoothing_3;
                        // gcv_f1 = gcv_f3;
                        let error = (smoothing_2 - smoothing_1) / (smoothing_1 + smoothing_2);
                        if error * error + 1.0 <= 1.0 || error <= tolerance {
                            solved = true;
                        } else {
                            smoothing_3 = smoothing_4;
                            gcv_f3 = gcv_f4;
                            alpha /= tau;
                            smoothing_4 = smoothing_1 + alpha;
                            gcv_f4 = fit_spline_coefficients_with_stats(
                                half_order, &data, &weight_factors, error_variance, smoothing_4, epsilon, &spline_tableau,
                                &weighted_matrix, weighted_matrix_norm, &mut coefficients, &mut stats,
                                &mut traced_matrix)?;
                        }
                    }
                }
                smoothing_1 = 0.5 * (smoothing_1 + smoothing_2);
                let _gcv_f1 = fit_spline_coefficients_with_stats(
                    half_order, &data, &weight_factors, error_variance, smoothing_1, epsilon, &spline_tableau,
                    &weighted_matrix, weighted_matrix_norm, &mut coefficients, &mut stats,
                    &mut traced_matrix)?;
            }
        }
    }

    Ok(coefficients)
}