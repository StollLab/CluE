use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use clue_oxide::clue_errors::CluEError;

pub struct PyCluEError(pub CluEError);

impl From<PyCluEError> for PyErr {
  fn from(error: PyCluEError) -> Self {
    PyValueError::new_err(format!("{}",error.0))
  }
}

impl From<CluEError> for PyCluEError {
  fn from(other: CluEError) -> Self {
    Self(other)
  }
}

