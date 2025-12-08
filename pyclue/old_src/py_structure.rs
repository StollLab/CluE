use pyo3::prelude::*;

use clue_oxide::structure::Structure;
use clue_oxide::clue_errors::CluEError;

use crate::py_config::PyConfig;
use crate::py_clue_errors::PyCluEError;
use crate::py_particle::PyParticle;

#[pyclass(name = "Structure")]
#[derive(Debug)]
pub struct PyStructure{
  pub structure: Structure,
}

#[pymethods]
impl PyStructure{
  /*
  #[new]
  fn new() -> Self{
    Self{
      structure: Structure::new(),
    }
  }
  */
  //----------------------------------------------------------------------------
  #[staticmethod]
  fn build_structure(pyconfig: &mut PyConfig) -> Result<Self,PyCluEError>
  {
    let structure = Structure::build_structure(
        &mut pyconfig.rng, &pyconfig.config)?;
   
    Ok(Self{
      structure,
    })
  }
  //----------------------------------------------------------------------------
  fn get_reference_index_of_nth_active(&self, n: usize) 
      -> Result<usize,PyCluEError>
  {
    Ok(self.structure.get_reference_index_of_nth_active(n)?)
  }
  //----------------------------------------------------------------------------
  pub fn get_average_detected_spin_position(&self) -> Result<[f64;3],PyCluEError>{
    match &self.structure.detected_particle {
      Some(det_particle) => {
        let r_ave = det_particle.get_average_spin_position()?;
        Ok(r_ave.elements.clone())
      },
      None => Err(PyCluEError::from(CluEError::NoDetectedSpinNotSet))
    }
  }
  //----------------------------------------------------------------------------
  pub fn get_bath_particles(&self) -> Vec::<PyParticle>{
    self.structure.bath_particles.iter()
      .map(|p| PyParticle{particle: p.clone()})
      .collect::<Vec::<PyParticle>>()
  }
  //----------------------------------------------------------------------------
  pub fn number(&self) -> usize{
    self.structure.number()
  }
  //----------------------------------------------------------------------------
  fn db_print(&self) {
    println!("{:?}",self);
  }
}

