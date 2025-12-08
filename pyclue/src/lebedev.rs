use pyo3::prelude::*;
use numpy::PyArray;

use ndarray::{Array, Ix2,Ix1};

use crate::py_clue_errors::PyCluEError;
use clue_oxide::integration_grid::IntegrationGrid;

#[pyclass(name = "IntegrationGrid")]
pub struct PyIntegrationGrid{
  grid: IntegrationGrid,
}
#[pymethods]
impl PyIntegrationGrid{
  #[staticmethod]
  pub fn lebedev(n: usize) -> Result<Self,PyCluEError>{
  
    let grid = IntegrationGrid::lebedev(n)?;

    Ok(Self{grid})
  }
  //----------------------------------------------------------------------------
  pub fn remove_3d_hemisphere(&self) -> Self{
    let grid = self.grid.remove_3d_hemisphere();
    Self{grid}
  }
  //----------------------------------------------------------------------------
  pub fn unpack(&self) 
    -> Result<(Py<PyArray<f64,Ix2>>,Py<PyArray<f64,Ix1>>),PyCluEError>
  {
    let dim = self.grid.dim;
    let n = self.grid.len();

    let points = match Array::from_shape_vec((n, dim), self.grid.points.clone()){
      Ok(p) => p,
     Err(err) => panic!("{}",err),
   };
   let weights = match Array::from_shape_vec((n), self.grid.weights.clone()){
     Ok(p) => p,
     Err(err) => panic!("{}",err),
   };

    Ok((
     Python::with_gil(|py|{
       PyArray::from_owned_array(py, points).unbind()
     }), 
     Python::with_gil(|py|{
       PyArray::from_owned_array(py, weights).unbind()
     }) 
    ))
  }
  //----------------------------------------------------------------------------
}


#[pyfunction]
pub fn build_lebedev(n: usize)
-> Result<(Py<PyArray<f64,Ix2>>,Py<PyArray<f64,Ix1>>),PyCluEError>
{
  let grid = IntegrationGrid::lebedev(n)?;

  let points = match Array::from_shape_vec((n, 3), grid.points){
    Ok(p) => p,
    Err(err) => panic!("{}",err),
  };
  let weights = match Array::from_shape_vec((n), grid.weights){
    Ok(p) => p,
    Err(err) => panic!("{}",err),
  };

  Ok((
  Python::with_gil(|py|{
    PyArray::from_owned_array(py, points).unbind()
  }), 
  Python::with_gil(|py|{
    PyArray::from_owned_array(py, weights).unbind()
  }) 
  ))
}
