use pyo3::prelude::*;

use clue_oxide::{
  cluster::read_clusters::read_cluster_file, 
  io::load_auxiliary_signals,
};

use crate::py_clue_errors::PyCluEError;
use crate::py_cluster::PyClusterSet;
use crate::py_config::PyConfig;
use crate::py_signal::PySignals;
use crate::py_structure::PyStructure;

use num_complex::Complex;


//------------------------------------------------------------------------------
#[pyfunction]
pub fn read_time_axis(filename: String)
    -> Result< Vec::<f64>  ,PyCluEError>
{
  Ok(clue_oxide::io::read_time_axis(&filename)?)
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn read_signal(filename: String)
    -> Result< Vec< Vec::<Complex<f64>> > ,PyCluEError>
{
  let signals = PySignals::read_from_csv(filename)?;

  Ok(signals.get_data())
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn read_auxiliary_signals(
    path_to_files: String, 
    cluster_file: String,
    pystructure: &PyStructure,
    pyconfig: &PyConfig) 
    -> Result<PyClusterSet,PyCluEError>
{
  let mut cluster_set = read_cluster_file(
      &cluster_file, 
      &pystructure.structure)?;

   load_auxiliary_signals(&mut cluster_set, &pyconfig.config, &path_to_files)?;

   Ok(PyClusterSet{cluster_set})
}
//------------------------------------------------------------------------------

