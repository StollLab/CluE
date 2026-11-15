use pyo3::prelude::*;

use numpy::PyArray;

use ndarray::{Array1,Ix1};

use num_complex::Complex;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

use std::collections::HashMap;

use clue_oxide::{
  config::{
    Config, 
    command_line_input::CommandLineInput,
    },
  info,
};
use clue_oxide::io::FromTOMLString;
use crate::py_clue_errors::PyCluEError;
use crate::dict_to_toml_string;

#[pyclass(name = "Config")]
#[derive(Debug)]
pub struct PyConfig{
  pub config: Config,
  pub rng: ChaCha20Rng,
}

#[pymethods]
impl PyConfig{
 
  //----------------------------------------------------------------------------
  pub fn db_print(&self){
    println!("{:?}",self);
  }
  //----------------------------------------------------------------------------
  /// This function runs the simulation and returns `(tau, tau2, signal)`:
  /// the time columns of the run, in units.time, and the highest-order 
  /// signal, with one entry per time point each.
  pub fn run(&self)
      -> Result<
          (Py<PyArray<f64,Ix1>>,Py<PyArray<f64,Ix1>>,
           Py<PyArray<Complex::<f64>,Ix1>>),
          PyCluEError
         >
  {
    // Close the output log (output_log.txt) before an error is returned, so
    // that the next run starts a new log.
    let result = clue_oxide::run(self.config.clone());
    clue_oxide::output_log::finish();
    let (tau_axes,signal) = result?;

    let tau = Array1::from_vec(tau_axes.tau);
    let tau2 = Array1::from_vec(tau_axes.tau2);

    let signal = Array1::from_vec(signal);

    let (tau, tau2) = Python::with_gil(|py|{
        (PyArray::from_owned_array(py, tau).unbind(),
         PyArray::from_owned_array(py, tau2).unbind())
    });
    let signal = Python::with_gil(|py|{
        PyArray::from_owned_array(py, signal).unbind()
    });
    Ok( (tau, tau2, signal) )
  }
  //----------------------------------------------------------------------------
  #[staticmethod]
  pub fn from_dict(config: HashMap::<String,PyObject>) -> Result<Self,PyCluEError>
  {
    let toml_string = dict_to_toml_string(config)?;
    Self::from_input(toml_string)
  }
  //----------------------------------------------------------------------------
  #[staticmethod]
  pub fn from_input(input: String) -> Result<Self,PyCluEError>
  {
  
    let config = if std::path::Path::new(&input).exists(){
      Config::from_toml_file(&input)
    }else{
      Config::from_toml_string(&input)
    }?;

    let rng = match config.rng_seed{
      Some(seed) => ChaCha20Rng::seed_from_u64(seed),
      None => ChaCha20Rng::from_entropy(),
    };

    Ok(Self{
      config,
      rng,
    })
  }
  //----------------------------------------------------------------------------
  #[staticmethod]
  fn from_command_line_input(args: String) -> Result<Self,PyCluEError>{

    let mut input = Vec::<String>::new();
    input.push("clue_oxide".to_string());

    for arg in args.split(" "){
      input.push(arg.to_string());
    }

    let input = CommandLineInput::new(input)?;

    // Decide what information to display.
    if input.show_help{
      info::help::print_help();
    }

    if input.show_license{
      info::license::print_license();
    }

    if input.show_warranty{
      info::warranty::print_warranty();
    }

    if input.show_version{
      info::version::print_version();
    }

    info::title::print_title();

    let mut config = Config::read_input(input)?;

    config.set_defaults()?;

    config.set_tau_axis()?;

    let rng = match config.rng_seed{
      Some(seed) => ChaCha20Rng::seed_from_u64(seed),
      None => ChaCha20Rng::from_entropy(),
    };

    Ok(Self{
      config,
      rng,
    })
  }
  //----------------------------------------------------------------------------
}

