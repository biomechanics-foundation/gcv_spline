use num_traits::Float;
use crate::woltring::support::{check_increasing, check_order, check_vector_length, FittingError};

pub(crate) fn create_weighted_matrix<T: Float>(half_order: usize, knots: &Vec<T>, weights_diagonal: &Vec<T>)
    -> Result<(Vec<T>, T), FittingError> {
    let num_knots = knots.len();
    check_increasing(knots)?;
    check_order(half_order, num_knots)?;
    check_vector_length(weights_diagonal, num_knots)?;

    let mut weighted_matrix = vec![T::from(0.).expect("Cannot convert to type from f64");
                                   ((2 * half_order) + 1) * num_knots];

    // Calculate factor
    let mut factor_1 = T::from(-1.).expect("Cannot convert to type from f64");
    if half_order != 1 {
        for idx in 2 ..= half_order {
            factor_1 = factor_1 * -T::from(idx).expect("Cannot convert to type from usize");
        }
        for idx in half_order + 1 ..= 2 * half_order - 1 {
            factor_1 = factor_1 * T::from(idx).expect("Cannot convert to type from usize");
        }
    }

    // Evaluate unweighted matrix by column
    let mut index_1 = 1;
    let mut index_2 = half_order;
    let mut order_index = half_order + 1;
    for knot_index in 1 ..= num_knots {
        let mut increment = 2 * half_order + 1;
        let factor;
        if knot_index > num_knots - half_order {
            factor_1 = -factor_1;
            factor = factor_1;
        } else {
            if knot_index < half_order + 1 {
                increment = 1;
                factor = factor_1;
            } else {
                factor = factor_1 * (knots[knot_index + half_order - 1] - knots[knot_index - half_order - 1]);
            }
        }
        if knot_index > half_order + 1 {
            index_1 += 1;
        }
        if index_2 < num_knots {
            index_2 += 1;
        }
        let mut matrix_index = order_index;

        let mut matrix_factor = factor;
        let mut knot_value = knots[index_1 - 1];
        for idx in index_1 + 1 ..=  index_2 {
            matrix_factor = matrix_factor / (knot_value - knots[idx - 1]);
        }
        weighted_matrix[matrix_index - 1] = matrix_factor;
        matrix_index += half_order * 2;
        if index_1 + 1 <= index_2 - 1 {
            for outer in index_1 + 1 ..= index_2 - 1 {
                matrix_factor = factor;
                knot_value = knots[outer - 1];
                for inner in index_1 ..= outer - 1 {
                    matrix_factor = matrix_factor / (knot_value - knots[inner - 1]);
                }
                for inner in outer + 1 ..= index_2 {
                    matrix_factor = matrix_factor / (knot_value - knots[inner - 1]);
                }
                weighted_matrix[matrix_index - 1] = matrix_factor;
                matrix_index += half_order * 2;
            }
        }
        matrix_factor = factor;
        knot_value = knots[index_2 - 1];
        for idx in index_1 ..= index_2 - 1 {
            matrix_factor = matrix_factor / (knot_value - knots[idx - 1]);
        }
        weighted_matrix[matrix_index - 1] = matrix_factor;
        order_index += increment;
    }

    // Compute weighted matrix and L1 norm
    let mut matrix_index = 0;
    let mut matrix_norm = T::from(0.).expect("Cannot convert to type from f64");
    for outer in 1 ..= num_knots {
        let current_weight = weights_diagonal[outer - 1];
        for _inner in 1 ..= half_order * 2 + 1 {
            matrix_index += 1;
            weighted_matrix[matrix_index - 1] = weighted_matrix[matrix_index - 1] / current_weight;
            matrix_norm = matrix_norm + weighted_matrix[matrix_index - 1].abs();
        }
    }
    matrix_norm = matrix_norm / T::from(num_knots).expect("Cannot convert to type from usize");

    Ok((weighted_matrix, matrix_norm))
}