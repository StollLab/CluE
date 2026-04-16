use pyo3::prelude::*;

use clue_oxide::cluster::adjacency::AdjacencyList;
use crate::py_clue_errors::PyCluEError;

#[pyclass(name = "AdjacencyList")]
#[derive(Debug,Clone)]
pub struct PyAdjacencyList{
  pub list: AdjacencyList,
}

#[pymethods]
impl PyAdjacencyList{
  //----------------------------------------------------------------------------
  #[staticmethod]
  pub fn from_array(array: Vec::<Vec::<usize>>) -> Result<Self,PyCluEError> {
    let n_vertices = array.len();
    let mut list = AdjacencyList::with_capacity(n_vertices);

    for (m,neighbors) in array.iter().enumerate(){
      for n in neighbors.iter(){
        list.connect(m,*n);
      }
    }  
    Ok(Self{list})
  }
  //----------------------------------------------------------------------------
  pub fn len(&self) -> usize{
    self.list.len()
  }
  //----------------------------------------------------------------------------
  pub fn is_empty(&self) -> bool{ self.list.is_empty() }
  //----------------------------------------------------------------------------
  #[staticmethod]
  pub fn with_capacity(n: usize) -> Self{
    let list = AdjacencyList::with_capacity(n);
    Self{list}
  }
  //----------------------------------------------------------------------------
  pub fn are_connected(&self, m: usize, n: usize) -> bool{
    self.list.are_connected(m,n)
  }
  //----------------------------------------------------------------------------
  pub fn connect(&mut self, m: usize, n: usize){
    self.list.connect(m,n);
  }
  //----------------------------------------------------------------------------
  pub fn is_active(&self, n: usize) -> bool {
    self.list.is_active(n)
  }
  //----------------------------------------------------------------------------
  pub fn activate(&mut self, n: usize){
    self.list.activate(n);
  }
  //----------------------------------------------------------------------------

}

