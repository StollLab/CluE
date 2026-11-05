/// This struct is holds the vertex adjacency lists for a graph.
/// Each neighbor list is kept sorted and without duplicates.
#[derive(Debug,Clone)]
pub struct AdjacencyList{
  list: Vec::< Option<Vec::<usize>> >,
}

impl AdjacencyList{

  //----------------------------------------------------------------------------
  /// This method returns the number of vertices in the list.
  pub fn len(&self) -> usize{ self.list.len() }
  //----------------------------------------------------------------------------
  /// This function `true` iff there are no vertices in the list.
  pub fn is_empty(&self) -> bool{ self.list.is_empty() }
  //----------------------------------------------------------------------------
  /// This function returns an `AdjacencyList` with `n` active but empty 
  /// vertices.
  pub fn active_with_capacity(n: usize) -> Self{
    let mut list = Vec::< Option<Vec::<usize>> >::with_capacity(n);
    for _ii in 0..n{
      list.push(Some(Vec::new()));
    }

    AdjacencyList{list}
  }
  //----------------------------------------------------------------------------
  /// This function returns an `AdjacencyList` with `n` empty vertices.
  pub fn with_capacity(n: usize) -> Self{
    let mut list = Vec::< Option<Vec::<usize>> >::with_capacity(n);
    for _ii in 0..n{
      list.push(None);
    }

    AdjacencyList{list}
  }
  //----------------------------------------------------------------------------
  /// This method returns `true` iff vertices `m` and `n` have an edge between
  /// them.
  pub fn are_connected(&self, m: usize, n: usize) -> bool{
    self.contains(m,n)
  }

  //----------------------------------------------------------------------------
  // This method returns `true` iff vertices `m` and `n` have an edge between
  // them.
  fn contains(&self,m: usize, n: usize) -> bool{

    // Neighbor lists are sorted.
    if let Some(neighbors) = &self.list[m]{
      return neighbors.binary_search(&n).is_ok();
    }

    false
  }
  //----------------------------------------------------------------------------
  /// This function places an edge between vertices `m` and `n`.
  /// Connecting a vertex to itself activates it.
  pub fn connect(&mut self, m: usize, n: usize){
    if m==n {
      self.activate(m);
      return
    }

    self.insert(m,n);
    self.insert(n,m);
  }

  //----------------------------------------------------------------------------
  pub fn is_active(&self, n: usize) -> bool {
    self.list[n].is_some()
  }
  //----------------------------------------------------------------------------
  /// This method activates vertex `n`.
  pub fn activate(&mut self, n: usize){
    if self.list[n].is_none(){
      self.list[n] = Some(Vec::new());
    }
  }
  //----------------------------------------------------------------------------
  // This function places `n` onto the neighbor list of `m`, keeping the 
  // list sorted and without duplicates.
  fn insert(&mut self, m: usize, n: usize){

    if let Some(neighbors) = &mut self.list[m]{
      if let Err(position) = neighbors.binary_search(&n){
        neighbors.insert(position, n);
      }
    }else{
      self.list[m] = Some(Vec::from([n]));
    }
  }
  //----------------------------------------------------------------------------
  /// This function returns the neighbor list of vertex `n`.
  pub fn get_neighbors(& self, n: usize) -> Option<& Vec::<usize>>{
    if let Some(neighbors) = &self.list[n]{
      return Some(neighbors);
    }
    None
  }
  //----------------------------------------------------------------------------
  /// This function returns the number of active vertices.
  pub fn n_active_vertices(&self) -> usize{
    let mut counter = 0;
    for ii in 0..self.len(){
      if self.list[ii].is_some(){
        counter += 1;
      }
    }
    counter
  }
  //----------------------------------------------------------------------------
  /// This method returns a list of the active vertices.
  pub fn get_active_vertices(&self) -> Vec::<usize>{
    let counter = self.n_active_vertices();

    let mut vertices = Vec::<usize>::with_capacity(counter);

    for ii in 0..self.len(){
      if self.list[ii].is_some(){
        vertices.push(ii);
      }
    }

    vertices
  }
  //----------------------------------------------------------------------------

}

#[cfg(test)]
mod tests {
  use super::*;

  #[allow(non_snake_case)]
  #[test]
  fn test_AdjacencyList(){
    let mut adjacency_list = AdjacencyList::with_capacity(6);
    
    assert!(!adjacency_list.are_connected(0,1));
    
    let neighbors = adjacency_list.get_neighbors(0);
    assert_eq!(neighbors,None);

    adjacency_list.connect(0,1);
    assert!(adjacency_list.are_connected(0,1));
    assert!(adjacency_list.are_connected(1,0));

    assert!(!adjacency_list.are_connected(0,2));

    let neighbors = adjacency_list.get_neighbors(0).unwrap();
    assert_eq!(*neighbors,vec![1]);

    adjacency_list.connect(0,2);
    let neighbors = adjacency_list.get_neighbors(0).unwrap();
    assert_eq!(*neighbors,vec![1,2]);

    adjacency_list.connect(0,2);
    let neighbors = adjacency_list.get_neighbors(0).unwrap();
    assert_eq!(*neighbors,vec![1,2]);

    let neighbors = adjacency_list.get_neighbors(1).unwrap();
    assert_eq!(*neighbors,vec![0]);

    adjacency_list.connect(4,2);
    adjacency_list.connect(4,0);
    let neighbors = adjacency_list.get_neighbors(4).unwrap();
    assert_eq!(*neighbors,vec![0,2]);

    let vertices = adjacency_list.get_active_vertices();
    assert_eq!(vertices,vec![0,1,2,4]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_AdjacencyList_matches_sorted_set_reference(){
    use rand::Rng;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};
    use std::collections::BTreeSet;

    let mut rng = ChaCha20Rng::seed_from_u64(0);
    for _trial in 0..100{
      let n = rng.gen_range(1..30);
      let mut graph = AdjacencyList::with_capacity(n);
      let mut reference = vec![None::<BTreeSet<usize>>; n];

      // Random edges, with repeats, both orders, and self-links.
      for _ in 0..rng.gen_range(0..200){
        let (a, b) = (rng.gen_range(0..n), rng.gen_range(0..n));
        graph.connect(a,b);
        reference[a].get_or_insert_with(BTreeSet::new);
        if a != b{
          reference[a].as_mut().unwrap().insert(b);
          reference[b].get_or_insert_with(BTreeSet::new).insert(a);
        }
      }

      for v in 0..n{
        let expected = reference[v].as_ref()
          .map(|set| set.iter().copied().collect::<Vec::<usize>>());
        assert_eq!(graph.get_neighbors(v).cloned(), expected);
        for w in 0..n{
          let connected = reference[v].as_ref()
            .map_or(false, |set| set.contains(&w));
          assert_eq!(graph.are_connected(v,w), connected);
        }
      }
    }
  }
}
