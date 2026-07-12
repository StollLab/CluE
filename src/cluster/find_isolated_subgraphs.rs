use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;

use std::collections::HashMap;
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

  let mut subgraphs = Vec::<AdjacencyList>::new();

  for vertex in 0..graph.len(){
    if Some(vertex) == root_node{
      continue;
    }

    match get_subgraph_index(&subgraphs,vertex){
      None => subgraphs.push( 
          grow_subgraph_from_vertex(graph,vertex,root_node)? ),
      Some(_) => continue,
    }

  }

  Ok(subgraphs)
}
//------------------------------------------------------------------------------
fn get_subgraph_index(subgraphs: &[AdjacencyList], vertex: usize) 
  -> Option<usize>
{
  for (idx,g) in subgraphs.iter().enumerate(){
    if g.is_active(vertex){
      return Some(idx);
    }
  }
  None
}
//------------------------------------------------------------------------------
fn grow_subgraph_from_vertex(
    graph: &AdjacencyList,
    vertex: usize,
    root_node: Option<usize>,
    )
  -> Result<AdjacencyList,CluEError>
{

  let mut subgraph = AdjacencyList::with_capacity(graph.len());
  subgraph.activate(vertex);

  let mut neighbors = HashMap::<usize,bool>::with_capacity(graph.len()); 
  neighbors.insert(vertex,false);

  loop{

    let Some(v) = get_next_vertex(&neighbors) else{
      break;
    };

    if let Some(new_vertices) = graph.get_neighbors(v){
      insert_neighbors(&mut neighbors, new_vertices, root_node);
    }
  }

  for (&v,_) in neighbors.iter(){
    for (&w,_) in neighbors.iter(){
      if w <= v{
        continue;
      }

      if graph.are_connected(v,w){
        subgraph.connect(w,v)
      }
    }
  }

  Ok(subgraph)
}
//------------------------------------------------------------------------------
fn get_next_vertex(neighbors: &HashMap::<usize,bool>) -> Option::<usize>
{
  for (key,value) in neighbors.iter(){
    if !(*value){
      return Some(*key);
    }
  }

  None
}
//------------------------------------------------------------------------------
fn insert_neighbors(
    neighbors: &mut HashMap::<usize,bool>, 
    new_vertices: &[usize],
    root_node: Option<usize>,
    )
{
  for v in new_vertices.iter(){
    if Some(*v) == root_node{
      continue;
    }
    if neighbors.contains_key(v){
      continue;
    }
    neighbors.insert(*v,false);
  }
}

