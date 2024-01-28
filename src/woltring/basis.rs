use num_traits::Float;
use crate::woltring::support::{check_increasing, check_order, FittingError};

pub(crate) fn create_basis<T: Float>(half_order: usize, knots: &Vec<T>)
    -> Result<(Vec<T>, T), FittingError> {
    let num_knots = knots.len();
    check_increasing(knots)?;
    check_order(half_order, num_knots)?;

    let spline_order = half_order * 2 - 1;

    // Linear case (half order and spline order = 1)
    if half_order == 1 {
        return Ok((vec![T::from(1.).expect("Cannot convert to type from f64"); num_knots],
                   T::from(1.).expect("Cannot convert to type from f64")));
    }

    let mut spline_tableau =
        vec![T::from(0.).expect("Cannot convert to type from f64"); num_knots * spline_order];
    // General case
    for knot_index in 1 ..= num_knots {
        let mut working_vec =
            vec![T::from(0.).expect("Cannot convert to type from f64"); 2 * half_order];

        // First row
        if knot_index != 1 && knot_index != num_knots {
            working_vec[2 * half_order - 2] = T::from(1.).expect("Cannot convert to type from f64")
                / (knots[knot_index] - knots[knot_index - 2]);
        } else {
            working_vec[2 * half_order - 2] = T::from(1.).expect("Cannot convert to type from f64");
        }

        // Further rows
        let current_knot = knots[knot_index - 1];
        for tableau_index in 3 ..= half_order * 2 {
            let mut working_index = half_order as i32 + 1 - tableau_index as i32;
            let mut v_working = working_vec[(working_index + half_order as i32 - 1) as usize];

            if knot_index < tableau_index {
                // Left spline
                for spline_index in knot_index + 1 ..= tableau_index {
                    let u_working = v_working;
                    v_working = working_vec[(working_index + half_order as i32) as usize];
                    working_vec[(working_index + half_order as i32 - 1) as usize] =
                        u_working + (knots[spline_index - 1] - current_knot) * v_working;
                    working_index += 1;
                }
            }

            let lower_bound: i32;
            if knot_index as i32 - tableau_index as i32 + 1 < 0 {
                lower_bound = 1;
            } else {
                lower_bound = (knot_index as i32 - tableau_index as i32 + 1).max(1);
            }
            let upper_bound = (knot_index - 1).min(num_knots - tableau_index);
            if lower_bound <= upper_bound as i32 {
                // Ordinary splines
                if tableau_index < half_order * 2 {
                    for spline_index in lower_bound..=upper_bound as i32 {
                        let y_working = knots[(tableau_index as i32 + spline_index - 1) as usize];
                        let u_working = v_working;
                        v_working = working_vec[(working_index + half_order as i32) as usize];
                        working_vec[(working_index + half_order as i32 - 1) as usize] =
                            u_working + (v_working - u_working) * (y_working - current_knot)
                                / (y_working - knots[(spline_index - 1) as usize]);
                        working_index += 1;
                    }
                } else {
                    for spline_index in lower_bound..=upper_bound as i32 {
                        let u_working = v_working;
                        v_working = working_vec[(working_index + half_order as i32) as usize];
                        working_vec[(working_index + half_order as i32 - 1) as usize] =
                            (current_knot - knots[(spline_index - 1) as usize]) * u_working +
                                (knots[(tableau_index as i32 + spline_index - 1) as usize] - current_knot) * v_working;
                        working_index += 1;
                    }
                }
            }

            if num_knots - tableau_index + 1 < knot_index {
                // Right splines
                for spline_index in num_knots - tableau_index + 1 ..= knot_index - 1 {
                    let u_working = v_working;
                    v_working = working_vec[(working_index + half_order as i32) as usize];
                    working_vec[(working_index + half_order as i32 - 1) as usize] =
                        (current_knot - knots[spline_index - 1]) * u_working + v_working;
                    working_index += 1;
                }
            }
        }
        for spline_index in 0 .. spline_order {
            spline_tableau[(knot_index - 1) * spline_order + spline_index] = working_vec[spline_index];
        }
    }

    let abs_vec = spline_tableau.iter().map(|element| element.abs()).collect::<Vec<T>>();
    let mut abs_iter = abs_vec.iter();

    let mut basis_l1_norm = T::from(0.).expect("Cannot convert to type from f64");

    while let Some(item) = abs_iter.next() {
        basis_l1_norm = basis_l1_norm + item.clone();
    }
    basis_l1_norm = basis_l1_norm / T::from(num_knots).expect("Cannot convert to usize from type");

    Ok((spline_tableau, basis_l1_norm))
}