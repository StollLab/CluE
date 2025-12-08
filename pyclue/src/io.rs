use pyo3::prelude::*;

use numpy::PyArray;

use clue_oxide::{
  cluster::read_clusters::read_cluster_file, 
//  io::load_auxiliary_signals,
};
use clue_oxide::signal::load_csv_to_vec_signals;

use crate::py_clue_errors::PyCluEError;
use crate::py_cluster::PyClusterSet;
use crate::py_config::PyConfig;
use crate::py_structure::PyStructure;

use ndarray::{Array1,Ix1};
use num_complex::Complex;


//------------------------------------------------------------------------------
#[pyfunction]
pub fn load_time_axis(filename: String)
    -> Result< Py<PyArray<f64,Ix1>>  ,PyCluEError>
{
  let s = Array1::from_vec(clue_oxide::io::read_time_axis(&filename)?);
  Ok(Python::with_gil(|py|{
      PyArray::from_owned_array(py, s).unbind()
  }))
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn load_signal(filename: String)
    -> Result< Py<PyArray<Complex<f64>,Ix1>>  ,PyCluEError>
{
  let signals = load_csv_to_vec_signals(&filename)?;
  let n = signals.len();
  let s = Array1::from_vec(signals[n-1].data.clone());

  Ok(Python::with_gil(|py|{
      PyArray::from_owned_array(py, s).unbind()
  }))

}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn load_auxiliary_signals(
    path_to_files: String, 
    cluster_file: String,
    pystructure: &PyStructure,
    pyconfig: &PyConfig) 
    -> Result<PyClusterSet,PyCluEError>
{
  let mut cluster_set = read_cluster_file(
      &cluster_file, 
      &pystructure.structure)?;

  clue_oxide::io::load_auxiliary_signals(
      &mut cluster_set, &pyconfig.config, &path_to_files)?;

   Ok(PyClusterSet{cluster_set})
}
//------------------------------------------------------------------------------

