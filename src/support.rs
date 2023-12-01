#[derive(Debug)]
pub enum FittingError {
    NotEnoughKnotsForOrder(String),
    VectorLengthMismatch(String),
    KnotsNotStrictlyIncreasing(String),
    InsufficientKnots(String),
}

pub fn check_order(half_order: usize, num_knots: usize) -> Result<(), FittingError> {
    if num_knots < 2 * half_order {
        return Err(FittingError::NotEnoughKnotsForOrder(
            String::from(format!(
                "At least {} knots needed, {} provided", 2 * half_order, num_knots
            )
            )
        ));
    }
    Ok(())
}

pub fn check_vector_length(weights_diagonal: &Vec<f64>, num_knots: usize) -> Result<(), FittingError> {
    if weights_diagonal.len() != num_knots {
        return Err(FittingError::VectorLengthMismatch(
            String::from(format!(
                "Vector length {} must match number of knots, {}",
                weights_diagonal.len(), num_knots
            ))
        ));
    }
    Ok(())
}

pub fn check_increasing(knots: &Vec<f64>) -> Result<(), FittingError> {
    if knots.len() < 2 {
        return Err(FittingError::InsufficientKnots(String::from("At least 2 knots are needed")));
    }

    let mut knots_iter = knots.iter();
    let mut previous = knots_iter.next().unwrap();
    while let Some(next) = knots_iter.next() {
        if previous >= next {
            return Err(FittingError::KnotsNotStrictlyIncreasing(
                String::from("Knots must be strictly increasing")
            ));
        }
        previous = next;
    }
    Ok(())
}
