use pyo3::prelude::*;

use clue_oxide::cluster::Cluster;
use clue_oxide::cluster::find_clusters::ClusterSet;

use crate::py_clue_errors::PyCluEError;
use crate::PyStructure;

use num_complex::Complex;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass(name = "Cluster")]
pub struct PyCluster{
  pub cluster: Cluster,
}

#[pymethods]
impl PyCluster{
  #[new]
  #[pyo3(signature = (vertices) )]
  pub fn from(vertices: Vec::<usize>) -> Self{
    PyCluster{
      cluster: Cluster::from(vertices),
    }
  }
  //----------------------------------------------------------------------------
  pub fn get_vertices(&self) -> Vec::<usize>{
    self.cluster.vertices().clone()
  }
  //----------------------------------------------------------------------------
  pub fn get_signal(&self) -> Result<Vec::<Complex::<f64>>,PyCluEError>{
    
    match &self.cluster.signal{
      Ok(Some(sig)) => Ok(sig.data.clone()),
      Ok(None) => Ok(Vec::<Complex::<f64>>::new()),
      Err(err) => Err(PyCluEError(err.clone()))
    }
  } 
  //----------------------------------------------------------------------------
  pub fn to_header(&self, pystructure: &PyStructure) ->
      Result<String,PyCluEError>
  {
    Ok(self.cluster.to_header(&pystructure.structure)?)
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass(name = "ClusterSet")]
pub struct PyClusterSet{
  pub cluster_set: ClusterSet,
}

#[pymethods]
impl PyClusterSet{
  pub fn len(&self) -> usize{
    self.cluster_set.len()
  } 
  //----------------------------------------------------------------------------
  pub fn is_empty(&self) -> bool{
    self.cluster_set.is_empty()
  } 
  //----------------------------------------------------------------------------
  pub fn get_cluster(&self, vertices: Vec::<usize>) -> Option<PyCluster>
  {
    let size_idx = vertices.len() - 1;

    let index: usize = match self.cluster_set.cluster_indices[size_idx]
        .get(&vertices)
    {
      Some(idx) => *idx,
      None => return None,
    };
   
    let cluster = self.cluster_set.clusters[size_idx][index].clone();

    Some(PyCluster{cluster})

  }
  //----------------------------------------------------------------------------
  pub fn extract_clusters(&self) -> Vec::< Vec::<PyCluster> >
  { 
    let mut pyclusters = Vec::< Vec::<PyCluster> >::with_capacity(self.len());

    for size_idx in 0..self.len(){

      let n_clusters = self.cluster_set.clusters[size_idx].len();
      pyclusters.push( Vec::<PyCluster>::with_capacity(n_clusters) );

      for cluster in self.cluster_set.clusters[size_idx].iter(){
        pyclusters[size_idx].push( PyCluster{cluster: cluster.clone()} );  
      }
    }

    pyclusters
  }
  //----------------------------------------------------------------------------
  pub fn extract_clusters_of_size(&self, size: usize) -> Vec::<PyCluster> 
  { 

    let size_idx = size - 1;

    let n_clusters = self.cluster_set.clusters[size_idx].len();
    let mut pyclusters =  Vec::<PyCluster>::with_capacity(n_clusters) ;

    for cluster in self.cluster_set.clusters[size_idx].iter(){
      pyclusters.push( PyCluster{cluster: cluster.clone()} );  
    }

    pyclusters
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

