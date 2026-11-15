use pyo3::prelude::*;

use numpy::PyArray;

use clue_oxide::{
  cluster::read_clusters::read_cluster_file, 
//  io::load_auxiliary_signals,
};
use clue_oxide::clue_errors::CluEError;
use clue_oxide::signal::{Signal,load_csv_to_vec_signals};

use crate::py_clue_errors::PyCluEError;
use crate::py_cluster::PyClusterSet;
use crate::py_config::PyConfig;
use crate::py_structure::PyStructure;

use ndarray::{Array1,Ix1};
use num_complex::Complex;


//------------------------------------------------------------------------------
/// This function returns the time columns `tau` and `tau2` of a signal file
/// as two arrays, in the unit of time of the run that wrote the file.
/// It returns an error if the file has no time columns.
#[pyfunction]
pub fn load_tau_axes(filename: String)
    -> Result< (Py<PyArray<f64,Ix1>>, Py<PyArray<f64,Ix1>>) ,PyCluEError>
{
  let (tau_axes, _signals) = load_csv_to_vec_signals(&filename)?;
  let Some(tau_axes) = tau_axes else{
    return Err(CluEError::CannotReadFile(filename,
          "the file has no time columns \"tau\" and \"tau2\"".to_string())
        .into());
  };
  let tau = Array1::from_vec(tau_axes.tau);
  let tau2 = Array1::from_vec(tau_axes.tau2);
  Ok(Python::with_gil(|py|{
      (PyArray::from_owned_array(py, tau).unbind(),
       PyArray::from_owned_array(py, tau2).unbind())
  }))
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn load_signal(filename: String)
    -> Result< Py<PyArray<Complex<f64>,Ix1>>  ,PyCluEError>
{
  // The last signal column; the time columns are not signals.
  let s = Array1::from_vec(Signal::read_from_csv(&filename)?.data);

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

