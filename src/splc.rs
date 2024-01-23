use crate::bandet::{consume_and_decompose};
use crate::bansol::solve_decomposed_system;
use crate::support::{check_order, check_vector_length, FittingError};
use crate::trinv::trace_inverse;

pub fn fit_spline_coefficients_with_stats(half_order: usize, data: &Vec<f64>,
                                          weight_factors: &Vec<f64>, variance: f64,
                                          real_smoothing: f64, tolerance: f64,
                                          spline_tableau: &Vec<f64>, weighted_tableau: &Vec<f64>,
                                          weighted_norm: f64)
        -> Result<(f64, Vec<f64>, Vec<f64>, Vec<f64>), FittingError> {
    let num_knots: usize = spline_tableau.len() / (2 * half_order - 1);
    check_order(half_order, num_knots)?;
    check_vector_length(data, num_knots)?;
    check_vector_length(weight_factors, num_knots)?;

    let mut smoothing = real_smoothing;
    let mut stats = vec![0.0; 6];
    let mut inverted_weighted_matrix = vec![0.0; weighted_tableau.len()];
    let splc: f64;

    stats[3] = real_smoothing / (1.0 + real_smoothing);

    // Pseudo least squares polynomial if smoothing parameter is too large
    if real_smoothing * weighted_norm * tolerance > 1.0 {
        stats[3] = 1.0;
        smoothing = 1.0 / (tolerance * weighted_norm);
    }

    // Pseudo interpolation if smoothing parameter is too small
    if real_smoothing * weighted_norm < tolerance {
        smoothing = tolerance / weighted_norm;
        stats[3] = 0.0;
    }

    // Calculate inverted weighted matrix
    for knot_index in 1 ..= num_knots {
        let lower_bound = -1 * std::cmp::min(half_order, knot_index - 1) as i32;
        let upper_bound = std::cmp::min(half_order, num_knots - knot_index) as i32;

        for inner in lower_bound ..= upper_bound {
            let index = ((knot_index as i32 - 1) * (half_order as i32 * 2 + 1) + inner + half_order as i32) as usize;
            let index_b = ((knot_index as i32 - 1) * (half_order as i32 * 2 - 1) + inner + half_order as i32 - 1) as usize;
            if inner.abs() as usize == half_order {
                inverted_weighted_matrix[index] = smoothing * weighted_tableau[index];
            } else {
                inverted_weighted_matrix[index] = spline_tableau[index_b] + smoothing * weighted_tableau[index];
            }
        }
    }

    // Solve matrix system inverted_weighted_matrix * coefficients = data,
    // evaluate TRACE[spline_tableau * inverted_weighted_matrix**-1]
    let decomposed_weighted_matrix = consume_and_decompose(inverted_weighted_matrix, half_order)?;
    let coefficients = solve_decomposed_system(&decomposed_weighted_matrix, data, half_order)?;
    let (traced_matrix, trace): (Vec<f64>, f64) = trace_inverse(weighted_tableau, decomposed_weighted_matrix, half_order)?;
    stats[2] = trace;
    let normalized_trace = trace / num_knots as f64;

    // Compute mean squared weighted residual
    let mut residual = 0.0;
    for knot_index in 1 ..= num_knots {
        let mut point = -data[knot_index - 1];

        let lower_bound = -1 * std::cmp::min(half_order - 1, knot_index - 1) as i32;
        let upper_bound = std::cmp::min(half_order - 1, num_knots - knot_index) as i32;

        for inner in lower_bound ..= upper_bound {
            let index = ((knot_index as i32 - 1) * (half_order as i32 * 2 - 1) + inner + half_order as i32 - 1) as usize;
            point += spline_tableau[index]
                * coefficients[(knot_index as i32 + inner - 1) as usize];
        }
        residual += point * point * weight_factors[knot_index - 1];
    }
    residual /= num_knots as f64;

    let estimated_variance = residual / normalized_trace; // Estimated variance
    stats[5] = estimated_variance;
    stats[0] = estimated_variance / normalized_trace; // GCV function value
    stats[1] = residual; // mean squared residual

    if variance < 0.0 {
        // Unknown variance: GCV
        stats[4] = estimated_variance - residual;
        splc = stats[0];
    } else {
        // Known variance: estimated mean squared error
        stats[4] = residual - variance * (2.0 * normalized_trace - 1.0);
        splc = stats[4];
    }

    Ok((splc, coefficients, stats, traced_matrix))
}