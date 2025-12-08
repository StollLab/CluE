use pyo3::prelude::*;

use numpy::PyArray;

use ndarray::Ix1;

use num_complex::Complex;

use std::collections::HashMap;

use crate::py_clue_errors::PyCluEError;
use crate::PyConfig;
use crate::dict_to_toml_string;




#[pyfunction]
pub fn run(config: HashMap::<String,PyObject>)
    -> Result<
         (Py<PyArray<f64,Ix1>>,Py<PyArray<Complex::<f64>,Ix1>>),
         PyCluEError
       >
{

  let toml_string = dict_to_toml_string(config)?;
  let pyconfig = PyConfig::from_input(toml_string)?;
  let (time_axis, signal)  = pyconfig.run()?;

  Ok( (time_axis, signal) )
}
