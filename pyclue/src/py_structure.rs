use pyo3::prelude::*;

use clue_oxide::structure::Structure;
use clue_oxide::clue_errors::CluEError;
use clue_oxide::space_3d::Vector3D;
use clue_oxide::cluster::adjacency::AdjacencyList;
use clue_oxide::structure::particle::Particle;

use crate::py_adjacency::PyAdjacencyList;
use crate::py_config::PyConfig;
use crate::py_clue_errors::PyCluEError;
use crate::py_particle::PyParticle;
use crate::PyExchangeGroupManager;

#[pyclass(name = "Structure")]
#[derive(Debug)]
pub struct PyStructure{
  pub structure: Structure,
}

#[pymethods]
impl PyStructure{
  #[new]
//  #[staticmethod]
  fn new(
      bath_particles: Vec::<PyParticle>,
      connections: PyAdjacencyList,
      cell_edges: Vec::<Vec::<f64>>,
      ) -> Self
  {
   
    let bath_particles: Vec::<Particle> = bath_particles.iter()
        .map(|p| p.particle.clone()).collect();

    let connections: AdjacencyList = connections.list;

    let mut cell_offsets = Vec::<Vector3D>::new();

    for edge in cell_edges.iter(){
      let a = Vector3D::from([edge[0],edge[1],edge[2]]);
      cell_offsets.push(a);
    } 

    Self{
      structure: Structure::new(
        bath_particles,
        connections,
        cell_offsets,        
      ),
    }
  }
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
  pub fn get_exchange_groups(&self) -> Option<PyExchangeGroupManager>{
    if let Some(exchange_group_manager) = self.structure.exchange_groups.clone()
    {
      return Some(PyExchangeGroupManager{exchange_group_manager});
    }
    None
  }
  //----------------------------------------------------------------------------
  pub fn number(&self) -> usize{
    self.structure.number()
  }
  //----------------------------------------------------------------------------
  pub fn write_pdb(&self,file_name: String) -> Result<(),PyCluEError>{
    Ok(self.structure.write_pdb(&file_name)?)
  }
  //----------------------------------------------------------------------------
  pub fn write_gro(&self,file_name: String) -> Result<(),PyCluEError>{
    Ok(self.structure.write_gro(&file_name)?)
  }
  //----------------------------------------------------------------------------
  fn db_print(&self) {
    println!("{:?}",self);
  }
}

