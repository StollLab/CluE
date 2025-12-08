use crate::clue_errors::CluEError;
use crate::io::FromTOMLString;

use serde::{Serialize,Deserialize};


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,Default,PartialEq,Serialize,Deserialize)]
pub struct ClusterTOML{
  pub number_clusters: Option<Vec::<usize>>,
  pub clusters: Option<Vec::<Vec::<Vec::<usize>>>>
}

impl FromTOMLString for ClusterTOML{
  fn from_toml_string(toml_str: &str) -> Result<Self,CluEError>{
    let decoded: Result<ClusterTOML,_> = toml::from_str(toml_str);
    match decoded{
      Ok(mut cluster_toml) => {
        cluster_toml.set_number_clusters();
        cluster_toml.validate()?;
        Ok(cluster_toml)
      },
      Err(err) => Err(CluEError::ClusterTOMLCannotRead( format!("{}",err) )),
    }
  }
}


impl ClusterTOML{
  fn set_number_clusters(&mut self){
    let Some(clusters) = &self.clusters else{
      return;
    };

    let mut number_clusters = Vec::<usize>::with_capacity(clusters.len());
    for clusters_of_size in clusters.iter(){
      number_clusters.push(clusters_of_size.len());
    }

    self.number_clusters = Some(number_clusters);
  }
  //----------------------------------------------------------------------------
  fn validate(&self) -> Result<(),CluEError>{
      
    let Some(clusters) = &self.clusters else{
      return Err(CluEError::ClusterTOMLNoClusters);
    };
    let Some(number_clusters) = &self.number_clusters else{
      return Err(CluEError::ClusterTOMLNoNumberClusters);
    };

    if number_clusters.len() != clusters.len(){
      return Err(CluEError::ClusterTOMLIncorrectNumberOfClusters);
    }

    for (ii,clusters_of_size) in clusters.iter().enumerate(){
      if number_clusters[ii] != clusters_of_size.len(){
        return Err(CluEError::ClusterTOMLIncorrectNumberOfClusters);
      }
    }

    Ok(())
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//==============================================================================
#[cfg(test)]
mod tests{
  use super::*;
 
  #[allow(non_snake_case)]
  #[test]
  fn test_ClusterTOML_from_string(){
    let toml_str = r##"
      number_clusters = [3,3,1]

      clusters = [
        [[1],[2],[3]],
        [[1,2],[1,3],[2,3]],
        [[1,2,3]]
      ]
    "##;
    let cluster_toml = ClusterTOML::from_toml_string(toml_str).unwrap();

    assert_eq!(cluster_toml.number_clusters,Some(vec![3,3,1]));

    let clusters = cluster_toml.clusters.unwrap();
    assert_eq!(clusters.len(),3);

    assert_eq!(clusters[0],vec![ vec![1],vec![2],vec![3]]);
    assert_eq!(clusters[1],vec![ vec![1,2],vec![1,3],vec![2,3]]);
  }
}

