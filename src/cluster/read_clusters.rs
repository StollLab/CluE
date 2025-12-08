use crate::clue_errors::CluEError;
use crate::cluster::Cluster;
use crate::cluster::cluster_set::ClusterSet;
use crate::structure::Structure;

use std::io::BufRead;
use substring::Substring;

//------------------------------------------------------------------------------
/// This function reads a cluster file and build the corresponding `ClusterSet`.
pub fn read_cluster_file(filename: &str, structure: &Structure) 
  -> Result<ClusterSet,CluEError>
{
  // DEPRECATED PATH
  if filename.ends_with(".txt"){
    return read_cluster_txt(filename, structure);
  }

  let toml_string: String = match std::fs::read_to_string(filename){
    Ok(lines) => lines,
    Err(_) => return Err(CluEError::CannotOpenFile(filename.to_string())),
  };

  let cluster_set = ClusterSet::from_toml_str(&toml_string, Some(structure))?;

  Ok(cluster_set)
}
//------------------------------------------------------------------------------
// DEPRECATED
fn read_cluster_txt(filename: &str, structure: &Structure) 
  -> Result<ClusterSet,CluEError>
{
  // Get number of lines.
  let n_lines = {   
    let Ok(file) = std::fs::File::open(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()));
    };

    std::io::BufReader::new(file).lines().count()
  }; 

  // Get lines.
  let lines_iter = {
    let Ok(file) = std::fs::File::open(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()));
    };

    std::io::BufReader::new(file).lines()
  };

  let mut lines = Vec::<String>::with_capacity(n_lines);

  // Allocate space for the cluster data.
  let mut clusters = initialize_clusters(&lines,filename)?;

  for line_result in lines_iter {
    // Extract the line.
    let Ok(line) = line_result else{
      return Err(CluEError::CannotOpenFile(filename.to_string()));
    };
    if line.is_empty(){
      continue;
    }
    lines.push(line);
  }

  set_clusters_from_lines(&mut clusters, lines,structure)?;

  Ok(ClusterSet::from(clusters))
}
//------------------------------------------------------------------------------
fn set_clusters_from_lines(
    clusters: &mut [Vec<Cluster>],
    lines: Vec::<String>, 
    structure: &Structure,
    ) -> Result<(),CluEError>
{
  let mut cluster_size = 0;

  // Loop through lines.
  for line in lines.iter() {

    // Check line type.
    if line.contains('#'){ 
      if !line.contains("cluster_size") || !line.contains('['){
        continue;
      }
      cluster_size = read_cluster_size_line(line)?;
      continue;
    }
    
    let cluster_reference_indices = parse_cluster_line(line)?;

    let mut cluster_indices 
      = Vec::<usize>::with_capacity(cluster_reference_indices.len());

      for &ref_idx in cluster_reference_indices.iter(){
        let n = structure.get_nth_active_from_reference_index(ref_idx)?;
        cluster_indices.push(n);
      }

    let cluster = Cluster::from(cluster_indices);

    let size_idx = if cluster_size == 0{
      cluster.len() - 1
    }else{
      cluster_size - 1
    };

    clusters[size_idx].push(cluster)

  }

  Ok(())
}
//------------------------------------------------------------------------------
fn read_cluster_size_line(line: &str) -> Result<usize,CluEError>
{
  let Some(idx_cluster_size) = line.find("cluster_size") else{
    return Err(CluEError::ClusterLineFormatError(line.to_string()));
  };
  let Some(idx_equals) = line.find('=') else{
    return Err(CluEError::ClusterLineFormatError(line.to_string()));
  };
  let Some(idx_close_bracket) = line.find(']') else{
    return Err(CluEError::ClusterLineFormatError(line.to_string()));
  };

  if idx_cluster_size >= idx_equals{
    return Err(CluEError::ClusterLineFormatError(line.to_string()));
  }
  if idx_equals >= idx_close_bracket{
    return Err(CluEError::ClusterLineFormatError(line.to_string()));
  }

  let cluster_size_str = line.substring(idx_equals+1, idx_close_bracket)
      .split_whitespace().collect::<Vec<_>>().join(" ");

  let Ok(cluster_size) = cluster_size_str.parse::<usize>() else{
      return Err(CluEError::ClusterLineFormatError(line.to_string()));
  };

  Ok(cluster_size)
}
//------------------------------------------------------------------------------
// This function reads lines that specify single clusters.
fn parse_cluster_line(line: &str) 
  -> Result<Vec::<usize>,CluEError>
{

  let Some(open_bracket) = line.find('[') else{
    return Err(CluEError::ClusterLineFormatError(line.to_string()));
  };
  let Some(close_bracket) = line.find(']') else{
    return Err(CluEError::ClusterLineFormatError(line.to_string()));
  };

  if open_bracket + 1 >= close_bracket {
    return Err(CluEError::ClusterLineFormatError(line.to_string()));
  }

  let cluster_str = line.substring(open_bracket+1,close_bracket);
  let n: usize = cluster_str.split(',').map(|_| 1).sum();

  let mut cluster = Vec::<usize>::with_capacity(n);
  for s in cluster_str.split(','){
    let Ok(index) = s.parse::<usize>() else{
      return Err(CluEError::ClusterLineFormatError(line.to_string()));
    };
    cluster.push(index);
  }

  Ok(cluster)
}
//------------------------------------------------------------------------------
// This function uses the meta data in the cluster files to initialize 
// the clusters' `Vec::<Vec::<Cluster>` to the correct size.
fn initialize_clusters(lines: &[String],filename: &str) 
  -> Result<Vec::<Vec::<Cluster>>,CluEError>
{


  for line in lines.iter() {

    if line.contains("#[clusters"){
      let number_clusters = parse_number_clusters(line, filename)?;
      let mut clusters = Vec::<Vec::<Cluster>>::with_capacity(
          number_clusters.len());

      for &n in number_clusters.iter(){
        clusters.push(Vec::<Cluster>::with_capacity(n));
      }

      return Ok(clusters);
    }

  }

  Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()))
}
//------------------------------------------------------------------------------
// This function reads the meta data in the cluster files to learn how many 
// clusters of each size there are.
fn parse_number_clusters(line: &str, filename: &str) 
  -> Result<Vec::<usize>,CluEError>
{
  
  let new_line = remove_char(line,' ');
  
  let Some(idx) = new_line.find("number_clusters") else{
    return Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()));
  };

  let new_line = new_line.substring(idx,new_line.len());

  let Some(open_bracket) = new_line.find('[') else{
    return Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()));
  };
  let Some(close_bracket) = new_line.find(']') else{
    return Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()));
  };


  if open_bracket + 1 >= close_bracket {
    return Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()));
  }
  let new_line = new_line.substring(open_bracket + 1,close_bracket) ;

  let max_cluster_size: usize = new_line.split(',').map(|_| 1).sum();
  let mut n_clusters = Vec::<usize>::with_capacity(max_cluster_size);

  for s in new_line.split(','){
    let Ok(n) = s.parse::<usize>() else {
      return Err(CluEError::ClusterFileContainsNoHeader(filename.to_string()));
    };
    n_clusters.push(n);
  }

  Ok(n_clusters)
}
//------------------------------------------------------------------------------
// This function takes a `&str` and builds a `String` identical except 
// all instances of a specified `char` removed.
fn remove_char(s: &str, c: char) -> String {
  let mut out_str = s.trim().to_owned();
  out_str.retain(|ch| {ch != c });
  out_str
}
//------------------------------------------------------------------------------


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  use crate::structure::particle::Particle;
  use crate::cluster::adjacency::AdjacencyList;
  use crate::space_3d::Vector3D;
  use crate::elements::Element;
  //----------------------------------------------------------------------------
  #[test]
  fn test_set_clusters_from_lines(){

    let lines = vec![
      "#[clusters, number_clusters = [11,23,59] ]".to_string(),
      "#[cluster_size = 1]".to_string(),
      "[3,4,5]".to_string(),
      "[6]".to_string(),
      "[7]".to_string(),
      "#[cluster_size = 2]".to_string(),
      "[3,4,5,6]".to_string(),
      "[3,4,5,7]".to_string(),
      "#[cluster_size = 3]".to_string(),
      "[3,4,5,6,7]".to_string(),
    ];

    let mut clusters = vec![
      Vec::<Cluster>::with_capacity(3),
      Vec::<Cluster>::with_capacity(2),
      Vec::<Cluster>::with_capacity(1),
    ];

    let bath_particles = vec![
      Particle::new(Element::Carbon, 0.0,0.0,1.0),
      Particle::new(Element::Carbon, 0.0,0.0,10.0),
      Particle::new(Element::Hydrogen,  1.0,0.0,0.0),
      Particle::new(Element::Hydrogen, -0.5,0.87,0.0),
      Particle::new(Element::Hydrogen, -0.5,-0.87,0.0),
      Particle::new(Element::Hydrogen, 0.0,0.0,5.0),
      Particle::new(Element::Hydrogen, 0.0,0.0,-5.0),
    ];
    let mut structure = Structure::new(
        bath_particles,
        AdjacencyList::with_capacity(2),
        Vec::<Vector3D>::new(),
        );

    structure.map_nth_active_to_reference_indices();

    set_clusters_from_lines(&mut clusters, lines, &structure).unwrap();

    assert_eq!(*clusters[0][0].vertices(),vec![3,4,5]);
    assert_eq!(*clusters[0][1].vertices(),vec![6]);
    assert_eq!(*clusters[0][2].vertices(),vec![7]);

    assert_eq!(*clusters[1][0].vertices(),vec![3,4,5,6]);
    assert_eq!(*clusters[1][1].vertices(),vec![3,4,5,7]);

    assert_eq!(*clusters[2][0].vertices(),vec![3,4,5,6,7]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_initialize_clusters(){
    let lines = vec![
      "#[clusters, number_clusters = [11,23,59] ]".to_string(),
      "#[cluster_size = 1]".to_string(),
      "[3,4,5]".to_string(),
      "[6]".to_string(),
      "[7]".to_string(),
      "#[cluster_size = 2]".to_string(),
      "[3,4,5,6]".to_string(),
      "[3,4,5,7]".to_string(),
      "#[cluster_size = 3]".to_string(),
      "[3,4,5,6,7]".to_string(),
    ];

    let clusters = initialize_clusters(&lines,"filename").unwrap();
    assert_eq!(clusters.len(),3);
    assert!(clusters[0].capacity() >= 3);
    assert!(clusters[1].capacity() >= 2);
    assert!(clusters[2].capacity() >= 1);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_cluster_line(){
    let line = "[7,13,23,42]";
    let value = parse_cluster_line(&line).unwrap();
    assert_eq!(value,vec![7,13,23,42]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_number_clusters(){
    let line = "#[clusters, number_clusters = [4276,32395] ]".to_string();
    let value = parse_number_clusters(&line, &"<file>").unwrap();
    assert_eq!(value, vec![4276,32395]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_remove_char(){
    let s = "#[clusters, number_clusters = [4276,32395] ]".to_string();
    let value = remove_char(&s,' ');
    let expected = "#[clusters,number_clusters=[4276,32395]]".to_string();
    assert_eq!(value,expected);
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
