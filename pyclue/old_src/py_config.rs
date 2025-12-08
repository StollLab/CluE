use pyo3::prelude::*;

use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

use clue_oxide::{
  config::{Config, command_line_input::CommandLineInput},
  info,
};

use crate::py_clue_errors::PyCluEError;

#[pyclass(name = "Config")]
#[derive(Debug)]
pub struct PyConfig{
  pub config: Config,
  pub rng: ChaCha20Rng,
}

#[pymethods]
impl PyConfig{
 
  fn db_print(&self){
    println!("{:?}",self);
  }
  //----------------------------------------------------------------------------
  #[staticmethod]
  fn from_input(args: String) -> Result<Self,PyCluEError>{

    let mut input = Vec::<String>::new();
    input.push("clue_oxixe".to_string());

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
  pub fn get_input_structure_file(&self) -> Option<String>{
    self.config.input_structure_file.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_input_structure_file(&mut self, value:  Option<String>){
    self.config.input_structure_file = value;
  }
  //----------------------------------------------------------------------------
  pub fn get_save_dir(&self) -> Option<String>{
    self.config.save_name.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_save_dir(&mut self, value:  Option<String>){
    self.config.save_name = value;
  }
  //----------------------------------------------------------------------------
  pub fn get_system_name(&self) -> Option<String>{
    self.config.system_name.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_system_name(&mut self, value:  Option<String>){
    self.config.system_name = value;
  }
  //----------------------------------------------------------------------------
}

