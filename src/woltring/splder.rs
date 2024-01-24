use crate::woltring::search::find_knot_interval;

pub fn evaluate_spline(derivative_order: usize, half_order: usize, point: f64, knots: &Vec<f64>,
                       coefficients: &Vec<f64>, knot_guess: usize) -> f64 {
    let num_knots = knots.len();

    // Derivatives of order >= 2 * half_order are always zero
    let order = half_order as i32 * 2 - derivative_order as i32;
    if order < 1 {
        return 0.0;
    }

    // Search for interval value
    let knot_interval = find_knot_interval(&knots, point, knot_guess);

    // Initialize parameters and first row of B-spline coefficients tableau
    let mut tableau = vec![0.0; 2 * half_order];
    let mut lower_index = knot_interval as i32 + 1;
    let upper_index = knot_interval as i32 + half_order as i32 * 2;
    let mut inner_index = num_knots - 2 * half_order;

    for index in lower_index ..= upper_index {
        if index >= half_order as i32 + 1 && index <= num_knots as i32 + half_order as i32 {
            tableau[(index - knot_interval as i32 - 1) as usize] = coefficients[(index - half_order as i32 - 1) as usize];
        } else {
            tableau[(index - knot_interval as i32 - 1) as usize] = 0.0;
        }
    }

    // The following loop computes differences of the B-spline coefficients. If the value of the
    // spline is required, differencing is not necessary.
    if derivative_order > 0 {
        lower_index -= half_order as i32 * 2;
        let index_bound = half_order as i32 * 2 - knot_interval as i32;
        for der_index in 1 ..= derivative_order {
            lower_index += 1;
            inner_index += 1;
            let idx_1 = 1.max(lower_index) as usize;
            let idx_2 = knot_interval.min(inner_index);
            let mut idx = idx_2 + 1;
            if idx_1 <= idx_2 {
                for _ in idx_1 ..= idx_2 {
                    idx -= 1;
                    let work_idx = (index_bound + idx as i32) as usize;
                    tableau[work_idx - 1] = (tableau[work_idx - 1] - tableau[work_idx - 2])
                        / (knots[idx + half_order * 2 - der_index - 1] - knots[idx - 1]);
                }
            }
            if lower_index < 1 {
                idx = (index_bound + 1) as usize;
                if der_index as i32 + 1 <= index_bound {
                    for _ in der_index as i32 + 1 ..= index_bound {
                        idx -= 1;
                        tableau[idx - 1] = -1.0 * tableau[idx - 2];
                    }
                }
            }
        }
        for idx in 1 ..= order {
            tableau[idx as usize - 1] = tableau[idx as usize + derivative_order - 1];
        }
    }

    let mut solution;

    // Compute lower half of the evaluation tableau
    if order - 1 >= 1 { // Tableau is ready if derivative order == 2 * half_order - 1
        for idx in 1 ..= order as usize - 1 {
            let order_idx = (num_knots as i32 - order) as usize + idx;
            let mut working_idx = order as usize;
            let mut knot_idx = knot_interval;

            // Right hand splines
            if knot_interval >= order_idx + 1 {
                for _ in order_idx + 1 ..= knot_interval {
                    tableau[working_idx - 1] = tableau[working_idx - 2] + (point - knots[knot_idx - 1]) * tableau[working_idx - 1];
                    knot_idx -= 1;
                    working_idx -= 1;
                }
            }

            // Middle B-splines
            //lk1i = knot_interval - order + 1 + idx
            let idx_1 = 1.max(knot_interval as i32 - order + 1 + idx as i32);
            let idx_2 = knot_interval.min(order_idx) as i32;
            if idx_1 <= idx_2 {
                for _ in idx_1 ..= idx_2 {
                    let working_knot = knots[(knot_idx as i32 + order - idx as i32 - 1) as usize];
                    solution = tableau[working_idx - 1];
                    tableau[working_idx - 1] = solution + (working_knot - point)
                        * (tableau[working_idx - 2] - solution) / (working_knot - knots[knot_idx - 1]);
                    working_idx -= 1;
                    knot_idx -= 1;
                }
            }

            // Left hand B-splines
            if knot_interval as i32 - order + 1 + idx as i32 <= 0 {
                knot_idx = (order - idx as i32) as usize;
                for _ in 1 ..= 1 - (knot_interval as i32 - order + 1 + idx as i32) {
                    tableau[working_idx - 1] = tableau[working_idx - 1]
                        + (knots[knot_idx - 1] - point) * tableau[working_idx - 2];
                    knot_idx -= 1;
                    working_idx -= 1;
                }
            }
        }
    }

    // Compute the return value
    solution = tableau[order as usize - 1];

    // Multiply with factorial of derivative order > 0
    if derivative_order > 0 {
        for idx in order ..= 2 * half_order as i32 - 1 {
            solution *= idx as f64;
        }
    }

    solution
}