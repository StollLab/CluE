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
  pub fn run(&self)
      -> Result<
          (Py<PyArray<f64,Ix1>>,Py<PyArray<Complex::<f64>,Ix1>>),
          PyCluEError
         >
  {
    let (time_axis,signal) = clue_oxide::run(self.config.clone())?;

    let time_axis = Array1::from_vec(time_axis);

    let signal = Array1::from_vec(signal);

    let time_axis = Python::with_gil(|py|{
        PyArray::from_owned_array(py, time_axis).unbind()
    });
    let signal = Python::with_gil(|py|{
        PyArray::from_owned_array(py, signal).unbind()
    });
    Ok( (time_axis, signal) )
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

    if input.show_warrenty{
      info::warrenty::print_warrenty();
    }

    if input.show_version{
      info::version::print_version();
    }

    if input.show_title{
      info::title::print_title();
    }

    let mut config = Config::read_input(input)?;

    config.set_defaults()?;

    config.construct_time_axis()?;

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

