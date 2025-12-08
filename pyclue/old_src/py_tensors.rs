use pyo3::prelude::*;

use clue_oxide::quantum::tensors::HamiltonianTensors;

use crate::py_clue_errors::PyCluEError;
use crate::py_config::PyConfig;
use crate::py_structure::PyStructure;

#[pyclass(name = "HamiltonianTensors")]
#[derive(Debug)]
pub struct PyHamiltonianTensors{
  pub tensors: HamiltonianTensors,
}

#[pymethods]
impl PyHamiltonianTensors{
  #[staticmethod]
  fn generate(pystructure: &PyStructure, pyconfig: &PyConfig) 
    -> Result<Self, PyCluEError>
  {
   
    let tensors = HamiltonianTensors::generate(
        &pystructure.structure, &pyconfig.config)?;

    Ok(Self{tensors})
  }
  //----------------------------------------------------------------------------
  fn len(&self) -> usize{
    self.tensors.len()
  }
  //----------------------------------------------------------------------------
  fn is_empty(&self) -> bool{
    self.tensors.is_empty()
  }
  //----------------------------------------------------------------------------
  fn save(&self, filename: String, pystructure: &PyStructure)
    -> Result<(),PyCluEError>
  {
    Ok(self.tensors.save(&filename, &pystructure.structure)?)
  }
  //----------------------------------------------------------------------------
  fn get_one_spin_tensor(&self, n: usize) -> Vec::<f64>{
    match self.tensors.spin1_tensors.get(n){
      Some(ten) =>  vec![ten.x(),ten.y(),ten.z()],
      None => Vec::<f64>::new(),
    }
  }
  //----------------------------------------------------------------------------
  fn get_two_spin_tensor(&self, m: usize, n: usize) -> Vec::<Vec::<f64>>{

    match self.tensors.spin2_tensors.get(m,n){
      Some(ten) =>  vec![
                      vec![ten.xx(),ten.xy(),ten.xz()],
                      vec![ten.yx(),ten.yy(),ten.yz()],
                      vec![ten.zx(),ten.zy(),ten.zz()],
                    ],
      None => Vec::<Vec::<f64>>::new(),
    }

  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  fn db_print(&self){
    println!("{:?}",self);
  }
}
