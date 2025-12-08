use pyo3::prelude::*;

use clue_oxide::signal::calculate_analytic_restricted_2cluster_signals as r2cce;

#[pyfunction]
pub fn hahn_three_spin_modulation_depth(delta_hf: f64,b:f64) -> f64
{
  r2cce::hahn_three_spin_modulation_depth(delta_hf, b)
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn hahn_three_spin_modulation_frequency(delta_hf: f64,b:f64) -> f64
{
  r2cce::hahn_three_spin_modulation_frequency(delta_hf,b)
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn hahn_three_spin_fourth_order_coefficient(delta_hf: f64,b:f64) -> f64
{
  r2cce::hahn_three_spin_fourth_order_coefficient(delta_hf,b)
}
