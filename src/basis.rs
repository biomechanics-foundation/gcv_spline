use crate::support::{check_increasing, check_order, FittingError};

pub fn create_basis(half_order: usize, knots: &Vec<f64>) -> Result<(Vec<f64>, f64), FittingError> {
    let num_knots = knots.len();
    check_increasing(knots)?;
    check_order(half_order, num_knots)?;

    let spline_order = half_order * 2 - 1;

    // Linear case (half order and spline order = 1)
    if half_order == 1 {
        return Ok((vec![1.0; num_knots], 1.0));
    }

    let mut spline_tableau = vec![0.0; num_knots * spline_order];
    // General case
    for knot_index in 1 ..= num_knots {
        let mut working_vec = vec![0.0; 2 * half_order];

        // First row
        if knot_index != 1 && knot_index != num_knots {
            working_vec[2 * half_order - 2] = 1.0 / (knots[knot_index] - knots[knot_index - 2]);
        } else {
            working_vec[2 * half_order - 2] = 1.0;
        }

        // Further rows
        let current_knot = knots[knot_index - 1];
        for tableau_index in 3 ..= half_order * 2 {
            let mut working_index = half_order + 1 - tableau_index;
            let mut v_working = working_vec[working_index + half_order - 1];

            if knot_index < tableau_index {
                // Left spline
                for spline_index in knot_index + 1 ..= tableau_index {
                    let u_working = v_working;
                    v_working = working_vec[working_index + half_order];
                    working_vec[working_index + half_order - 1] =
                        u_working + (knots[spline_index - 1] - current_knot) * v_working;
                    working_index += 1;
                }
            }

            let lower_bound = std::cmp::max(knot_index - tableau_index + 1, 1);
            let upper_bound = std::cmp::min(knot_index - 1, num_knots - tableau_index);
            if lower_bound <= upper_bound {
                // Ordinary splines
                if tableau_index < half_order * 2 {
                    for spline_index in lower_bound..=upper_bound {
                        let y_working = knots[tableau_index + spline_index - 1];
                        let u_working = v_working;
                        v_working = working_vec[working_index + half_order];
                        working_vec[working_index + half_order - 1] =
                            u_working + (v_working - u_working) * (y_working - current_knot)
                                / (y_working - knots[spline_index - 1]);
                        working_index += 1;
                    }
                } else {
                    for spline_index in lower_bound..=upper_bound {
                        let u_working = v_working;
                        v_working = working_vec[working_index + half_order];
                        working_vec[working_index + half_order - 1] =
                            (current_knot - knots[spline_index - 1]) * u_working +
                                (knots[tableau_index + spline_index - 1] - current_knot) * v_working;
                        working_index += 1;
                    }
                }
            }

            if num_knots - tableau_index + 1 < knot_index {
                // Right splines
                for spline_index in num_knots - tableau_index + 1 ..= knot_index - 1 {
                    let u_working = v_working;
                    v_working = working_vec[working_index + half_order];
                    working_vec[working_index + half_order - 1] =
                        (current_knot - knots[spline_index - 1]) * u_working + v_working;
                    working_index += 1;
                }
            }
        }
        for spline_index in 0 .. spline_order {
            spline_tableau[(knot_index - 1) * spline_order + spline_index] = working_vec[spline_index];
        }
    }

    let basis_l1_norm: f64 = spline_tableau.iter().map(|element| element.abs()).sum::<f64>() / num_knots as f64;

    Ok((spline_tableau, basis_l1_norm))
}