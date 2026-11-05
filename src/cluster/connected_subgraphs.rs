use crate::cluster::adjacency::AdjacencyList;
use crate::cluster::find_isolated_subgraphs::find_nearly_isolated_components;

/// This fuctions takes an AdjacencyList and returns a list of sets,
/// as well as a list of map indicating where vertex is.
/// Each set contains the indices of a connected subgraph (a connected 
/// component).  The sets are ordered by their smallest vertex, and the 
/// vertices of each set are sorted.  A vertex without an entry in 
/// `adjacency_list` forms its own set.
pub fn separate_into_connected_subgraphs(adjacency_list: &AdjacencyList)
  -> ( Vec::<Vec::<usize>>, Vec::<usize> )
{
  // Every vertex is connected to every vertex of its component, even when
  // a later vertex is what joins two earlier ones (e.g. 0-2-1).
  let subgraphs = find_nearly_isolated_components(adjacency_list, None);

  let mut subgraph_ids = vec![0; adjacency_list.len()];
  for (id, subgraph) in subgraphs.iter().enumerate(){
    for &v in subgraph.iter(){
      subgraph_ids[v] = id;
    }
  }

  (subgraphs,subgraph_ids)

}

#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_separate_into_connected_subgraphs(){

    // 0-1-5
    // 2-3-4
    // 6
    // 7-8
    let mut graph = AdjacencyList::with_capacity(9);
    for v in 0..9{
      graph.connect(v,v);
    }
    graph.connect(0,1);
    graph.connect(0,5);
    graph.connect(2,3);
    graph.connect(3,4);
    graph.connect(7,8);

    let expected_subgraphs = vec![
      vec![0,1,5],
      vec![2,3,4],
      vec![6],
      vec![7,8],
    ];

    let expected_ids = vec![0,0,1,1,1,0,2,3,3];

    let (subgraphs,ids) = separate_into_connected_subgraphs(&graph);

    assert_eq!(ids,expected_ids);
    assert_eq!(subgraphs,expected_subgraphs);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_separate_into_connected_subgraphs_joins_earlier_groups(){
    // 0-2-1: vertex 2 joins vertices 0 and 1, which have no earlier bond.
    let mut graph = AdjacencyList::with_capacity(3);
    graph.connect(0,2);
    graph.connect(1,2);
    let (subgraphs,ids) = separate_into_connected_subgraphs(&graph);
    assert_eq!(subgraphs, vec![vec![0,1,2]]);
    assert_eq!(ids, vec![0,0,0]);

    // Ring closed by the last vertex, plus a separate pair: 
    // 0-3, 1-4, 2-4, 3-5, 4-5 (one ring/tree component), 6-7.
    let mut graph = AdjacencyList::with_capacity(8);
    for (a,b) in [(0,3),(1,4),(2,4),(3,5),(4,5),(6,7)]{
      graph.connect(a,b);
    }
    let (subgraphs,ids) = separate_into_connected_subgraphs(&graph);
    assert_eq!(subgraphs, vec![vec![0,1,2,3,4,5], vec![6,7]]);
    assert_eq!(ids, vec![0,0,0,0,0,0,1,1]);

    let graph = AdjacencyList::with_capacity(0);
    let (subgraphs,ids) = separate_into_connected_subgraphs(&graph);
    assert!(subgraphs.is_empty());
    assert!(ids.is_empty());
  }
}
