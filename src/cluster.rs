use crate::clue_errors::*;
use crate::signal::Signal;
use crate::structure::Structure;

pub mod adjacency;
pub mod build_adjacency_list;
pub mod cluster_set;
pub mod cluster_toml;
pub mod connected_subgraphs;
pub mod find_clusters;
pub mod find_isolated_subgraphs;
pub mod read_clusters;
pub mod get_subclusters;
pub mod methyl_clusters;
pub mod partition;
pub mod unit_of_clustering;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `Cluster` contains both the vertices of the cluster and the simulated 
/// signal. 
#[derive(Debug,Clone,PartialEq)]
pub struct Cluster{
  pub vertices: Vec::<usize>,
  pub signal: Result<Option<Signal>,CluEError>,
}

//------------------------------------------------------------------------------
impl Cluster{
  /// This function build a `Cluster` from a list of vertices.
  pub fn from(vertices: Vec::<usize>) -> Self{
    Cluster{
      vertices,
      signal: Ok(None),
    }
  }
}
//------------------------------------------------------------------------------
impl std::fmt::Display for Cluster {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    // "[1,2,3,4]"
    let cluster = &self.vertices;

    if cluster.is_empty(){
      return write!(f,"[]");
    }

    let mut string = format!("[{}", cluster[0]);

    for vertex in cluster.iter().skip(1){
      string = format!("{},{}",string,vertex);
    } 
    write!(f,"{}]",string)
  }
}
//------------------------------------------------------------------------------
impl Cluster{
  // This function returns the cluster size. 
  fn len(&self) -> usize{
    self.vertices.len()
  }
  //----------------------------------------------------------------------------
  /*
  // This function translates `Cluster` to a `Result<String>.
  // The difference between this function and `to_string()` is that this 
  // function need a `&Structure` to convert the vertices from internal 
  // indices to be consistant across output files.
  fn to_string_result(&self,structure: &Structure) -> Result<String,CluEError>{
    // "[1,2,3,4]"
    let cluster = &self.vertices;

    if cluster.is_empty(){
      return Ok("[]".to_string());
    }

    let mut string = format!("[{}", 
        structure.get_reference_index_of_nth_active(cluster[0])? 
        );

    for &vertex in cluster.iter().skip(1){
      string = format!("{},{}",string,
          structure.get_reference_index_of_nth_active(vertex)?
          );
    } 
    Ok(format!("{}]",string))
  }
  */
  //----------------------------------------------------------------------------
  /// This function returns a reference to the cluster vertices.
  pub fn vertices(&self) -> &Vec::<usize>{
    &self.vertices
  } 
  //----------------------------------------------------------------------------
  pub fn reference_vertices(&self,structure_opt: Option<&Structure>)
    -> Result<Vec::<usize>,CluEError>
  {

    let mut vertices = self.vertices.clone();
    match structure_opt{
      Some(structure) => {
        for v in vertices.iter_mut(){
          *v = structure.get_reference_index_of_nth_active(*v)?;
        }
      },
      None => (),
    }

    Ok(vertices)
  }
  //----------------------------------------------------------------------------
  /// This function generates header for the cluster in csv files.
  pub fn to_header(&self, structure: &Structure) -> Result<String,CluEError> {
    // "clu_1_2_3_4"
    let cluster = &self.vertices;

    if cluster.is_empty(){
      return Ok("clu".to_string());
    }

    let mut string = format!("clu_{}", 
        structure.get_reference_index_of_nth_active(cluster[0])? 
        );

    for &vertex in cluster.iter().skip(1){
      string = format!("{}_{}",string,
          structure.get_reference_index_of_nth_active(vertex)?
          );
    } 

    Ok(string)
  }
  //----------------------------------------------------------------------------
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



