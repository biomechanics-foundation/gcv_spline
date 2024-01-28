use std::{error::Error, fmt};
use num_traits::Float;

/// Reports input errors that prevent fitting a GCV spline.
#[derive(Debug)]
pub enum FittingError {
    NotEnoughKnotsForOrder(String),
    VectorLengthMismatch(String),
    KnotsNotStrictlyIncreasing(String),
    InsufficientKnots(String),
    MatrixMismatch(String),
}

impl Error for FittingError {}
impl fmt::Display for FittingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FittingError: {:?}", self)
    }
}

pub(crate) fn check_order(half_order: usize, num_knots: usize) -> Result<(), FittingError> {
    if num_knots < 2 * half_order {
        return Err(FittingError::NotEnoughKnotsForOrder(
            String::from(format!(
                "At least {} knots (time points) needed, {} provided", 2 * half_order, num_knots
            )
            )
        ));
    }
    Ok(())
}

pub(crate) fn check_vector_length<T: Float>(vector: &Vec<T>, num_knots: usize) -> Result<(), FittingError> {
    if vector.len() != num_knots {
        return Err(FittingError::VectorLengthMismatch(
            String::from(format!(
                "Vector length {} must match number of knots (time points), {}",
                vector.len(), num_knots
            ))
        ));
    }
    Ok(())
}

pub(crate) fn check_increasing<T: Float>(knots: &Vec<T>) -> Result<(), FittingError> {
    if knots.len() < 2 {
        return Err(FittingError::InsufficientKnots(String::from("At least 2 knots (time points) are needed")));
    }

    let mut knots_iter = knots.iter();
    let mut previous = knots_iter.next().unwrap();
    while let Some(next) = knots_iter.next() {
        if *previous >= *next {
            return Err(FittingError::KnotsNotStrictlyIncreasing(
                String::from("Knots must be strictly increasing")
            ));
        }
        previous = next;
    }
    Ok(())
}

pub(crate) fn check_matrix_size<T: Float>(matrix_1: &Vec<T>, matrix_2: &Vec<T>) -> Result<(), FittingError> {
    if matrix_1.len() != matrix_2.len() {
        return Err(FittingError::MatrixMismatch(
            String::from(format!("Matrix size mismatch: {} and {}", matrix_1.len(), matrix_2.len()))
        ));
    }
    Ok(())
}
