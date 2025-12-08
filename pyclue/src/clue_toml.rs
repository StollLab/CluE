use pyo3::prelude::*;

use std::collections::HashMap;

use clue_oxide::{
  clue_errors::CluEError,
};
use crate::py_clue_errors::PyCluEError;
//use crate::PyConfig;

use pyo3::ffi::c_str;


#[pyfunction]
pub fn dict_to_toml_string(dict: HashMap::<String,PyObject>) 
    -> Result<String,PyCluEError>
{
  Ok(Python::with_gil(|py| {
    let Ok(clue_string) = PyModule::from_code(
        py,
        c_str!(r##"
import toml

class CluE_String:
  data = ''

  def write(self,s):
    self.data = f'{self.data}{s}'

  def __repr__(self):
    return self.data

def dict_to_toml(x):    
  toml_str = CluE_String()

  toml.dump(x,toml_str)

  return str(toml_str)
"##
          ),
        c_str!("clue_string.py"),
        c_str!("clue_string"),
    ) else{
      return Err(CluEError::Error("invalid python code".to_string()));    
    };

    let attr = match clue_string.getattr("dict_to_toml"){
      Ok(a) => a,
      Err(err) => return Err(CluEError::Error(err.to_string())),
    };

    let call = match attr.call1((dict,)){
      Ok(a) => a,
      Err(err) => return Err(CluEError::Error(err.to_string())),
    };

    let toml_string: String = match call.extract(){
      Ok(a) => a,
      Err(err) => return Err(CluEError::Error(err.to_string())),
    };

    Ok(toml_string)
  })?)
}
