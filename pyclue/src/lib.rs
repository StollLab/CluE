pub mod analytic_restricted_2cluster;
pub mod clue_toml;
pub mod io;
pub mod lebedev;
pub mod py_clue_errors;
pub mod py_cluster;
pub mod py_config;
pub mod py_exchange_groups;
pub mod py_particle;
pub mod py_physical_constants;
pub mod py_signal;
pub mod py_spin_operators;
pub mod py_structure;
pub mod py_tensors;
pub mod run_clue;

use analytic_restricted_2cluster::*;
use clue_toml::*;
use io::*;
use lebedev::*;
use py_cluster::PyCluster;
use py_config::*;
use py_exchange_groups::*;
use py_signal::*;
use py_particle::*;
use py_physical_constants::*;
use py_spin_operators::*;
use py_structure::*;
use py_tensors::*;
use run_clue::*;

#[allow(unused)]
fn main() {
use pyo3::prelude::*;
use pyo3::exceptions::PyTypeError;
use std::iter::zip;

use std::collections::HashMap;

use num_complex::Complex;

use clue_oxide as clue;



/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule(name = "clue_oxide")]
fn clue_odide(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPhysicalConstants>()?;
    m.add_class::<PyCluster>()?;
    m.add_class::<PyConfig>()?;
    m.add_class::<PyElement>()?;
    m.add_class::<PyExchangeGroupManager>()?;
    m.add_class::<PyIntegrationGrid>()?;
    m.add_class::<PyIsotope>()?;
    m.add_class::<PyParticle>()?;
    m.add_class::<PySignal>()?;
    m.add_class::<PySignals>()?;
    m.add_class::<PyStructure>()?;
    m.add_class::<PyHamiltonianTensors>()?;
    m.add_function( wrap_pyfunction!(run,m)? )?;
    m.add_function( wrap_pyfunction!(build_lebedev,m)? )?;
    m.add_function( wrap_pyfunction!(dict_to_toml_string,m)? )?;
    m.add_function( wrap_pyfunction!(load_time_axis,m)? )?;
    m.add_function( wrap_pyfunction!(load_signal,m)? )?;
    m.add_function( wrap_pyfunction!(load_auxiliary_signals,m)? )?;
    m.add_function( wrap_pyfunction!(hahn_three_spin_modulation_depth,m)? )?;
    m.add_function( wrap_pyfunction!(hahn_three_spin_modulation_frequency,m)? )?;
    m.add_function( 
        wrap_pyfunction!(hahn_three_spin_fourth_order_coefficient,m)? )?;
    m.add_function( wrap_pyfunction!(spin_x,m)? )?;
    m.add_function( wrap_pyfunction!(spin_y,m)? )?;
    m.add_function( wrap_pyfunction!(spin_z,m)? )?;
    m.add_function( wrap_pyfunction!(spin_plus,m)? )?;
    m.add_function( wrap_pyfunction!(spin_minus,m)? )?;
    m.add_function( wrap_pyfunction!(spin_squared,m)? )?;
    Ok(())
}
}

