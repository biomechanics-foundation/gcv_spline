use crate::woltring::support::{check_order, FittingError};

pub(crate) fn solve_decomposed_system(matrix: &Vec<f64>, rhs: &Vec<f64>, half_order: usize) -> Result<Vec<f64>, FittingError> {
    let num_knots = matrix.len() / (2 * half_order + 1);
    check_order(half_order, num_knots)?;
    let mut solution = vec![0.0; num_knots];

    if half_order == 0 {
        for idx in 1 ..= num_knots {
            solution[idx - 1] = rhs[idx - 1] / matrix[(idx - 1) * (half_order * 2 + 1) + half_order];
        }
    } else {
        solution[0] = rhs[0];
        for outer in 2 ..= num_knots {
            let inner_limit = std::cmp::min(half_order, outer - 1);
            let mut solution_point = rhs[outer - 1];
            for inner in 1 ..= inner_limit {
                solution_point -= matrix[(outer - 1) * (half_order * 2 + 1) - inner + half_order] *
                    solution[outer - inner - 1];
            }
            solution[outer - 1] = solution_point;
        }

        solution[num_knots - 1] /= matrix[(num_knots - 1) * (half_order * 2 + 1) + half_order];
        for outer in (1 .. num_knots).rev() {
            let inner_limit = std::cmp::min(half_order, num_knots - outer);
            let mut solution_point = solution[outer - 1];
            for inner in 1 ..= inner_limit {
                solution_point -= matrix[(outer - 1) * (half_order * 2 + 1) + inner + half_order] *
                    solution[outer + inner - 1];
            }
            solution[outer - 1] = solution_point / matrix[(outer - 1) * (half_order * 2 + 1) + half_order];
        }
    }

    Ok(solution)
}