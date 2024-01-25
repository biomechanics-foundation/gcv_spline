use crate::woltring::support::{check_order, FittingError};

pub(crate) fn consume_and_decompose(mut matrix: Vec<f64>, half_order: usize) -> Result<Vec<f64>, FittingError> {
    let num_knots = matrix.len() / (2 * half_order + 1);
    check_order(half_order, num_knots)?;

    for knot_index in 1 ..= num_knots {
        let mut decomp_inner = matrix[(knot_index - 1) * (half_order * 2 + 1) + half_order];
        let order_index = std::cmp::min(half_order, knot_index - 1);

        if order_index >= 1 {
            for idx in 1 ..= order_index {
                decomp_inner -= matrix[(knot_index - 1) * (half_order * 2 + 1) - idx + half_order] *
                    matrix[(knot_index - idx - 1) * (half_order * 2 + 1) + idx + half_order];
            }
            matrix[(knot_index - 1) * (half_order * 2 + 1) + half_order] = decomp_inner;
        }

        let outer_limit = std::cmp::min(half_order, num_knots - knot_index);

        if outer_limit >= 1 {
            for outer in 1 ..= outer_limit {
                let mut decomp_lower = matrix[(knot_index + outer - 1) * (half_order * 2 + 1) - outer + half_order];

                let inner_limit = std::cmp::min(half_order - outer, knot_index - 1);

                if inner_limit >= 1 {
                    let mut decomp_upper = matrix[(knot_index - 1) * (half_order * 2 + 1) + outer + half_order];
                    for inner in 1 ..= inner_limit {
                        decomp_upper -= matrix[(knot_index - 1) * (half_order * 2 + 1) - inner + half_order] *
                            matrix[(knot_index - inner - 1) * (half_order * 2 + 1) + outer + inner + half_order];
                        decomp_lower -= matrix[(outer + knot_index - 1) * (half_order * 2 + 1) - outer - inner + half_order] *
                            matrix[(knot_index - inner - 1) * (half_order * 2 + 1) + inner + half_order];
                    }
                    matrix[(knot_index - 1) * (half_order * 2 + 1) + outer + half_order] = decomp_upper;
                }
                matrix[(knot_index + outer - 1) * (half_order * 2 + 1) - outer + half_order] = decomp_lower / decomp_inner;
            }
        }
    }

    Ok(matrix)
}