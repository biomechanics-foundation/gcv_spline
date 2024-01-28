use num_traits::Float;
use crate::woltring::bandet::{consume_and_decompose};
use crate::woltring::bansol::solve_decomposed_system;
use crate::woltring::support::{check_order, check_vector_length, FittingError};
use crate::woltring::trinv::trace_inverse;

pub(crate) fn fit_spline_coefficients_with_stats<T: Float>(half_order: usize, data: &Vec<T>,
                                          weight_factors: &Vec<T>, variance: T,
                                          real_smoothing: T, tolerance: T,
                                          spline_tableau: &Vec<T>, weighted_tableau: &Vec<T>,
                                          weighted_norm: T, coefs_current: &mut Vec<T>,
                                          stats_current: &mut Vec<T>,
                                          traced_current: &mut Vec<T> )
                                          -> Result<T, FittingError> {
    let num_knots: usize = spline_tableau.len() / (2 * half_order - 1);
    check_order(half_order, num_knots)?;
    check_vector_length(data, num_knots)?;
    check_vector_length(weight_factors, num_knots)?;

    let mut smoothing = real_smoothing;
    let mut stats = vec![T::from(0.).expect("Cannot convert to type from f64"); 6];
    let mut inverted_weighted_matrix = vec![T::from(0.).expect("Cannot convert to type from f64");
                                            weighted_tableau.len()];
    let splc: T;

    stats[3] = real_smoothing / (T::from(1.).expect("Cannot convert to type from f64") + real_smoothing);

    // Pseudo least squares polynomial if smoothing parameter is too large
    if real_smoothing * weighted_norm * tolerance > T::from(1.).expect("Cannot convert to type from f64") {
        stats[3] = T::from(1.).expect("Cannot convert to type from f64");
        smoothing = T::from(1.).expect("Cannot convert to type from f64") / (tolerance * weighted_norm);
    }

    // Pseudo interpolation if smoothing parameter is too small
    if real_smoothing * weighted_norm < tolerance {
        smoothing = tolerance / weighted_norm;
        stats[3] = T::from(0.).expect("Cannot convert to type from f64");
    }

    // Calculate inverted weighted matrix
    for knot_index in 1 ..= num_knots {
        let lower_bound = -1 * half_order.min(knot_index - 1) as i32;
        let upper_bound = half_order.min(num_knots - knot_index) as i32;

        for inner in lower_bound ..= upper_bound {
            let index = ((knot_index as i32 - 1) * (half_order as i32 * 2 + 1) + inner +
                half_order as i32) as usize;
            let index_b = ((knot_index as i32 - 1) * (half_order as i32 * 2 - 1) + inner +
                half_order as i32 - 1) as usize;
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
    let (traced_matrix, trace): (Vec<T>, T) = trace_inverse(weighted_tableau, decomposed_weighted_matrix, half_order)?;
    stats[2] = trace;
    let normalized_trace = trace / T::from(num_knots).expect("Cannot convert to type from usize");

    // Compute mean squared weighted residual
    let mut residual = T::from(0.).expect("Cannot convert to type from f64");
    for knot_index in 1 ..= num_knots {
        let mut point = -data[knot_index - 1];

        let lower_bound = -1 * (half_order - 1).min(knot_index - 1) as i32;
        let upper_bound = (half_order - 1).min(num_knots - knot_index) as i32;

        for inner in lower_bound ..= upper_bound {
            let index = ((knot_index as i32 - 1) *
                (half_order as i32 * 2 - 1) + inner + half_order as i32 - 1) as usize;
            point = point + (spline_tableau[index]
                * coefficients[(knot_index as i32 + inner - 1) as usize]);
        }
        residual = residual + (point * point * weight_factors[knot_index - 1]);
    }
    residual = residual / T::from(num_knots).expect("Cannot convert to type from usize");

    let estimated_variance = residual / normalized_trace; // Estimated variance
    stats[5] = estimated_variance;
    stats[0] = estimated_variance / normalized_trace; // GCV function value
    stats[1] = residual; // mean squared residual

    if variance < T::from(0.).expect("Cannot convert to type from f64") {
        // Unknown variance: GCV
        stats[4] = estimated_variance - residual;
        splc = stats[0];
    } else {
        // Known variance: estimated mean squared error
        stats[4] = residual - variance * (T::from(2.).expect("Cannot convert to type from f64")
            * normalized_trace - T::from(1.).expect("Cannot convert to type from f64"));
        splc = stats[4];
    }

    *coefs_current = coefficients;
    *stats_current = stats;
    *traced_current = traced_matrix;
    Ok(splc)
}