use num_traits::Float;
use crate::woltring::support::{check_matrix_size, check_order, FittingError};

pub(crate) fn trace_inverse<T: Float>(basis_tableau: &Vec<T>, mut decomp_matrix: Vec<T>, half_order: usize)
    -> Result<(Vec<T>, T), FittingError> {
    let num_knots = decomp_matrix.len() / (2 * half_order + 1);
    check_order(half_order, num_knots)?;
    check_matrix_size(basis_tableau, &decomp_matrix)?;

    let mut inversion;

    // Invert decomp_matrix and store in place
    decomp_matrix[(num_knots - 1) * (half_order * 2 + 1) + half_order] =
        T::from(1.).expect("Cannot convert to type from f64") / decomp_matrix[(num_knots - 1) * (half_order * 2 + 1) + half_order]; // n-th pivot
    for knot_index in (1 ..= num_knots - 1).rev() {
        let order_index = half_order.min(num_knots - knot_index);

        inversion = T::from(1.).expect("Cannot convert to type from f64") / decomp_matrix[(knot_index - 1) * (half_order * 2 + 1) + half_order]; // i-th pivot

        // Save i-th column of L, row of U and normalize U
        for idx in 1 ..= order_index {
            decomp_matrix[(num_knots - 1) * (half_order * 2 + 1) + idx + half_order] =
                decomp_matrix[(knot_index - 1) * (half_order * 2 + 1) + idx + half_order] * inversion;
            decomp_matrix[half_order - idx] = decomp_matrix[(idx + knot_index - 1) * (half_order * 2 + 1) - idx + half_order];
        }
        inversion = inversion + inversion;

        // Invert around i-th pivot
        for outer in (1 ..= order_index).rev() {
            let mut invert_upper = T::from(0.).expect("Cannot convert to type from f64");
            let mut invert_lower = T::from(0.).expect("Cannot convert to type from f64");
            for inner in 1 ..= order_index {
                invert_upper = invert_upper - (decomp_matrix[(num_knots - 1) * (half_order * 2 + 1) + inner + half_order] *
                    decomp_matrix[(knot_index + inner - 1) * (half_order * 2 + 1) + outer - inner + half_order]);
                invert_lower = invert_lower - (decomp_matrix[half_order - inner] *
                    decomp_matrix[(knot_index + outer - 1) * (half_order * 2 + 1) + inner - outer + half_order]);
            }
            decomp_matrix[(knot_index - 1) * (half_order * 2 + 1) + outer + half_order] = invert_upper;
            decomp_matrix[(outer + knot_index - 1) * (half_order * 2 + 1) - outer + half_order] = invert_lower;
            inversion = inversion - (decomp_matrix[(num_knots - 1) * (half_order * 2 + 1) + outer + half_order]
                * invert_lower + decomp_matrix[half_order - outer] * invert_upper);
        }
        decomp_matrix[(knot_index - 1) * (half_order * 2 + 1) + half_order] = inversion / T::from(2.).expect("Cannot convert to type from f64");
    }

    // Trace and zero portions of inverted matrix
    let mut trace = T::from(0.).expect("Cannot convert to type from f64");
    for knot_index in 1 ..= num_knots {
        let lower_bound = half_order.min(knot_index - 1);
        let upper_bound = half_order.min(num_knots - knot_index);
        for idx in lower_bound ..= upper_bound {
            trace = trace + (basis_tableau[(knot_index - 1) * (half_order * 2 + 1) + idx + half_order] *
                decomp_matrix[(idx + knot_index - 1) * (half_order * 2 + 1) - idx + half_order]);
        }
    }
    for order_index in 1 ..= half_order {
        decomp_matrix[(num_knots - 1) * (half_order * 2 + 1) + order_index + half_order] = T::from(0.).expect("Cannot convert to type from f64");
        decomp_matrix[half_order - order_index] = T::from(0.).expect("Cannot convert to type from f64");
    }

    Ok((decomp_matrix, trace))
}