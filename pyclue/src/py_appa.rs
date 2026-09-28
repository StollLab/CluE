use pyo3::prelude::*;

use clue_oxide::cluster_methods::appa;

#[pyfunction]
pub fn appa_hahn_modulation_depth(delta_hf: f64,b:f64) -> f64
{
  appa::appa_hahn_modulation_depth(delta_hf,b)
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn appa_hahn_frequency(delta_hf: f64,b:f64) -> f64
{
  appa::appa_hahn_frequency(delta_hf,b)
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn appa_hahn_fourth_order_coefficient(delta_hf: f64,b:f64) -> f64
{
  appa::appa_hahn_fourth_order_coefficient(delta_hf,b)
}
