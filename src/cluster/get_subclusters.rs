
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This function takes a list of vertices and returns all possible proper 
/// subsets of the the vertices.
pub fn build_subclusters(cluster: &[usize]) -> Vec::<Vec::<usize>>{

  let cluster_size = cluster.len();
  
  if cluster_size == 0 {
    // There are no proper subclusters.
    return Vec::<Vec::<usize>>::new();
  }
  
  let n_subclusters = 2_usize.pow( cluster_size as u32);
  
  let mut subclusters = Vec::<Vec::<usize>>::with_capacity(n_subclusters-2);
  for _ii in 0..n_subclusters-2{
    subclusters.push(Vec::<usize>::with_capacity(cluster.len()));
  }

  for (clu_idx, vertex) in cluster.iter().enumerate(){
    let subcluster_selection 
      = include_vertex_in_which_subcluster(
          clu_idx, cluster_size);

    for idx in subcluster_selection.iter(){
      subclusters[*idx].push(*vertex);
    }
    
  }

  subclusters
}
//------------------------------------------------------------------------------
// This function returns a `Vec::<usize>` of indices, where each index
// identifies a subcluster that contains the specified vertex.
// The vertex is specified by its index (`idx`) within the cluster of size
// `cluster_size`.
//
// For example, the set {a,b,c} has 6 (from 2^3 - 2) non-empty proper subsets:
// {b,c}, {a,c}, {c}, {a,b}, {b}, {a}.
// Given this ordering of these sets in a list, the sets can be specified
// by listing all the sets that contain a as [1,3,5], contain b as [0,3,4],
// and contain c as [0,1,2].  
// The elements of the set are arbitrary; only the index (`idx`) of the element
// within the set is important.
// 011, 101, 001, 110, 010, 100
fn include_vertex_in_which_subcluster(idx: usize,cluster_size: usize) 
  -> Vec::<usize>
{
  let n: usize = 2_usize.pow(cluster_size as u32);

  let increment: usize = 2_usize.pow(idx as u32);

  let mut do_push = false;
  let mut indices = Vec::<usize>::with_capacity(n);
  for ii in 0..n{
    if ii % increment == 0 { do_push = !do_push;} 
    if ii == 0 {continue;}
    if do_push {
      indices.push(ii-1);
    }
  }

  indices
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_include_vertex_in_which_subcluster(){
    let expected = vec![
      vec![1,3,5],
      vec![0,3,4],
      vec![0,1,2],
    ];
    let cluster_size = 3;
    for ii in 0..cluster_size{
      let indices 
        = include_vertex_in_which_subcluster(ii,cluster_size);
      assert_eq!(indices,expected[ii]);
    }
  }
  //----------------------------------------------------------------------------

  #[test]
  fn test_build_subclusters(){
    let subclusters = build_subclusters(&[1,2,3]);
    assert_eq!(subclusters,vec![
        vec![2, 3], vec![1, 3], vec![3], vec![1, 2], vec![2], vec![1]] );

    let subclusters = build_subclusters(&[1]);
    assert!(subclusters.is_empty());

    let subclusters = build_subclusters(&[]);
    assert!(subclusters.is_empty());
  }
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
