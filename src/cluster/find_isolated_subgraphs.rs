use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;

//------------------------------------------------------------------------------
/// This function splits a graph, defined by an `AdjacencyList`, into 
/// isotated subsraphs.  Each subgraph is itself connected, but not connected
/// to any other subgraphs.
/// The set of subgraphs partition the graph.
pub fn find_isolated_subgraphs(graph: &AdjacencyList)
  -> Result<Vec::<AdjacencyList>,CluEError>
{
  find_nearly_isolated_subgraphs(graph, None)
}
//------------------------------------------------------------------------------
/// This function splits a graph, defined by an `AdjacencyList`, into nearly 
/// isotated subsraphs.  Each subgraph is itself connected, but not connected
/// to any other subgraphs, except possibly through the root node.
/// If G is the graph and r is the root node, then set of subgraphs
/// partitions G\{r}.  If `root_node == None`, then then set of subgraphs
/// partitions G.
pub fn find_nearly_isolated_subgraphs(
    graph: &AdjacencyList, root_node: Option<usize>)
  -> Result<Vec::<AdjacencyList>,CluEError>
{
  // Each returned AdjacencyList has graph.len() entries, so this function
  // needs O(number of subgraphs * graph.len()) memory.  Use 
  // `find_nearly_isolated_components()` when only the vertices are needed.
  let components = find_nearly_isolated_components(graph, root_node);

  let mut subgraphs = Vec::<AdjacencyList>::with_capacity(components.len());

  for component in components.iter(){
    let mut subgraph = AdjacencyList::with_capacity(graph.len());

    for &v in component.iter(){
      subgraph.activate(v);

      let Some(neighbors) = graph.get_neighbors(v) else{
        continue;
      };
      // Edges to the root node are not part of any subgraph.
      for &w in neighbors.iter(){
        if w > v && Some(w) != root_node{
          subgraph.connect(v,w);
        }
      }
    }

    subgraphs.push(subgraph);
  }

  Ok(subgraphs)
}
//------------------------------------------------------------------------------
/// This function returns the vertex sets of the nearly isolated subgraphs of
/// `graph` (see `find_nearly_isolated_subgraphs()`): the connected components
/// of G\{r}, where r is `root_node` (if any).  The components are ordered by
/// their smallest vertex, and the vertices of each component are sorted.
/// Every vertex other than the root is in exactly one component, including 
/// vertices without an entry in `graph`, which form single-vertex 
/// components.
/// The search takes O(V + E) time and O(V) memory.
pub fn find_nearly_isolated_components(
    graph: &AdjacencyList, root_node: Option<usize>) -> Vec::<Vec::<usize>>
{
  let n_vertices = graph.len();

  let mut is_assigned = vec![false; n_vertices];
  if let Some(root) = root_node{
    if root < n_vertices{
      is_assigned[root] = true;
    }
  }

  let mut components = Vec::<Vec::<usize>>::new();
  let mut to_visit = Vec::<usize>::new();

  for start in 0..n_vertices{
    if is_assigned[start]{
      continue;
    }

    // Depth-first search from `start`, never passing through the root.
    let mut component = vec![start];
    is_assigned[start] = true;
    to_visit.push(start);

    while let Some(v) = to_visit.pop(){
      let Some(neighbors) = graph.get_neighbors(v) else{
        continue;
      };
      for &w in neighbors.iter(){
        if is_assigned[w]{
          continue;
        }
        is_assigned[w] = true;
        component.push(w);
        to_visit.push(w);
      }
    }

    component.sort();
    components.push(component);
  }

  components
}
//------------------------------------------------------------------------------


#[cfg(test)]
mod tests{
  use super::*;

  //----------------------------------------------------------------------------
  //  1---2   0   4
  //       \ / \ /
  //        3   5---6
  //
  //  7 (isolated)
  fn get_graph_for_tests() -> AdjacencyList{
    let mut graph = AdjacencyList::with_capacity(8);
    graph.connect(1,2);
    graph.connect(2,3);
    graph.connect(0,3);
    graph.connect(0,5);
    graph.connect(4,5);
    graph.connect(5,6);
    graph.connect(7,7);
    graph
  }
  //----------------------------------------------------------------------------
  // This function returns the sorted vertex lists of the subgraphs.
  fn get_sorted_vertex_sets(subgraphs: &[AdjacencyList]) -> Vec::<Vec::<usize>>{
    let mut sets = subgraphs.iter().map(|g| g.get_active_vertices())
      .collect::<Vec::<Vec::<usize>>>();
    sets.sort();
    sets
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_find_nearly_isolated_subgraphs(){
    let graph = get_graph_for_tests();

    let subgraphs = find_nearly_isolated_subgraphs(&graph, Some(0)).unwrap();
    assert_eq!(get_sorted_vertex_sets(&subgraphs),
        vec![vec![1,2,3], vec![4,5,6], vec![7]]);

    // Edges inside each subgraph are kept.
    let g = subgraphs.iter().find(|g| g.is_active(1)).unwrap();
    assert!(g.are_connected(1,2));
    assert!(g.are_connected(2,3));
    assert!(!g.are_connected(1,3));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_find_isolated_subgraphs(){
    let graph = get_graph_for_tests();

    let subgraphs = find_isolated_subgraphs(&graph).unwrap();
    assert_eq!(get_sorted_vertex_sets(&subgraphs),
        vec![vec![0,1,2,3,4,5,6], vec![7]]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_find_nearly_isolated_components(){
    let graph = get_graph_for_tests();

    assert_eq!(find_nearly_isolated_components(&graph, Some(0)),
        vec![vec![1,2,3], vec![4,5,6], vec![7]]);
    assert_eq!(find_nearly_isolated_components(&graph, None),
        vec![vec![0,1,2,3,4,5,6], vec![7]]);

    // A later vertex joins two earlier vertices: 0-2-1 is one component.
    let mut graph = AdjacencyList::with_capacity(3);
    graph.connect(0,2);
    graph.connect(1,2);
    assert_eq!(find_nearly_isolated_components(&graph, None), 
        vec![vec![0,1,2]]);

    // Vertices without an entry form single-vertex components.
    let mut graph = AdjacencyList::with_capacity(5);
    graph.connect(1,3);
    assert_eq!(find_nearly_isolated_components(&graph, Some(4)),
        vec![vec![0], vec![1,3], vec![2]]);

    let graph = AdjacencyList::with_capacity(0);
    assert!(find_nearly_isolated_components(&graph, Some(0)).is_empty());
  }
  //----------------------------------------------------------------------------
  // This function returns the component label of each vertex of G\{root},
  // using union-find, as an independent reference.
  fn reference_component_labels(graph: &AdjacencyList, root: Option<usize>)
    -> Vec::<Option<usize>>
  {
    fn find(parent: &mut Vec::<usize>, v: usize) -> usize{
      let mut r = v;
      while parent[r] != r{ r = parent[r]; }
      parent[v] = r;
      r
    }
    let n = graph.len();
    let mut parent = (0..n).collect::<Vec::<usize>>();
    for v in 0..n{
      if Some(v) == root{ continue; }
      let Some(neighbors) = graph.get_neighbors(v) else{ continue; };
      for &w in neighbors.iter(){
        if Some(w) == root{ continue; }
        let (rv, rw) = (find(&mut parent, v), find(&mut parent, w));
        if rv != rw{ parent[rv.max(rw)] = rv.min(rw); }
      }
    }
    (0..n).map(|v| if Some(v) == root { None } 
        else { Some(find(&mut parent, v)) }).collect()
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_components_and_subgraphs_match_union_find_on_random_graphs(){
    use rand::Rng;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    let mut rng = ChaCha20Rng::seed_from_u64(0);
    for _trial in 0..200{
      let n = rng.gen_range(1..40);
      let p_edge = rng.gen_range(0.0..0.15);
      let mut graph = AdjacencyList::with_capacity(n);
      for v in 0..n{
        if rng.gen_bool(0.9){ graph.connect(v,v); } // Some vertices unset.
        for w in v+1..n{
          if rng.gen_bool(p_edge){ graph.connect(v,w); }
        }
      }
      let root = if n > 1 && rng.gen_bool(0.7) { Some(0) } else { None };

      let labels = reference_component_labels(&graph, root);
      let components = find_nearly_isolated_components(&graph, root);

      // Every non-root vertex is in exactly one component, components are
      // sorted and ordered by their smallest vertex, and two vertices share
      // a component iff union-find puts them together.
      let mut seen = vec![false; n];
      let mut previous_min = None;
      for component in components.iter(){
        assert!(component.windows(2).all(|w| w[0] < w[1]));
        assert!(previous_min < Some(component[0]));
        previous_min = Some(component[0]);
        for &v in component.iter(){
          assert!(!seen[v]);
          seen[v] = true;
          assert_eq!(labels[v], labels[component[0]]);
        }
      }
      for v in 0..n{
        assert_eq!(seen[v], Some(v) != root);
      }
      let n_labels = { let mut l = labels.iter().flatten().copied()
          .collect::<Vec::<usize>>(); l.sort(); l.dedup(); l.len() };
      assert_eq!(components.len(), n_labels);

      // The subgraphs have the same vertices and every edge of G\{root}.
      let subgraphs = find_nearly_isolated_subgraphs(&graph, root).unwrap();
      assert_eq!(subgraphs.len(), components.len());
      for (subgraph, component) in subgraphs.iter().zip(components.iter()){
        assert_eq!(&subgraph.get_active_vertices(), component);
        for &v in component.iter(){
          for &w in component.iter(){
            if v != w{
              assert_eq!(subgraph.are_connected(v,w), 
                  graph.are_connected(v,w));
            }
          }
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_find_nearly_isolated_components_scales_to_large_baths(){
    // 20000 isolated spins and a root connected to all of them.  The old 
    // search needed one 20001-entry AdjacencyList per spin (about 9.6 GB).
    let n = 20_000;
    let mut graph = AdjacencyList::with_capacity(n + 1);
    for v in 1..=n{ graph.connect(0,v); }

    let start = std::time::Instant::now();
    let components = find_nearly_isolated_components(&graph, Some(0));
    println!("{} components in {:?}", components.len(), start.elapsed());

    assert_eq!(components.len(), n);
    assert!(components.iter().enumerate()
        .all(|(ii,component)| *component == vec![ii + 1]));
  }
  //----------------------------------------------------------------------------
}
