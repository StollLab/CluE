use pyo3::prelude::*;

use numpy::PyArray;

use clue_oxide::structure::exchange_groups::{
  GetCentroid,
  ExchangeGroupManager,
};

use ndarray::{Array2,Ix2};

#[pyclass(name = "ExchangeGroupManager")]
#[derive(Debug,Clone)]
pub struct PyExchangeGroupManager{
  pub exchange_group_manager: ExchangeGroupManager,
}

#[pymethods]
impl PyExchangeGroupManager{
  pub fn get_centroids(&self) -> Py<PyArray<f64,Ix2>>{
  
    let n_grps = self.exchange_group_manager.exchange_groups.len();
    let mut centroids = Array2::zeros((3, n_grps));

    for (ii,ex_grp) in self.exchange_group_manager.exchange_groups
      .iter().enumerate()
    {
      let r = ex_grp.centroid();
      centroids[[0,ii]] = r.x();
      centroids[[1,ii]] = r.y();
      centroids[[2,ii]] = r.z();
    }
    Python::with_gil(|py|{
        PyArray::from_owned_array_bound(py, centroids).unbind()
    })
  }
}
