use crate::CluEError;
use crate::space_3d::Vector3D;

use rand::seq::index::sample;
use rand_chacha::ChaCha20Rng;

use mcmf::{GraphBuilder, Vertex, Cost, Capacity};
use math::round;

// Largest integer cost of a datum-to-kluster edge in restricted k-means.
const MAX_RESTRICTED_KMEANS_COST: i64 = 1_000_000;

// The mcmf crate solves the flow problem with LEMON's NetworkSimplex with 
// 32-bit int costs.  LEMON adds artificial arcs with cost 
// (max_cost + 1)*(number of nodes) (lemon/network_simplex.h, ART_COST), 
// which must not overflow; this bound leaves a factor of two of headroom.
const LEMON_COST_BOUND: i64 = 1 << 30;

/// This function partitions `data` into `k` k-means clusters and returns the
/// cluster label of each datum.  If `do_restricted_kmeans` is true, every
/// cluster has at most ceil(n/k) members.  At most `iter_cap` 
/// assignment/update iterations are done after the first assignment; if the
/// assignments have not converged by then, a warning is printed and the last
/// assignment is returned.
pub fn get_kmeans_partition(
    rng: &mut ChaCha20Rng,
    k: usize,
    data: &[Vector3D],
    do_restricted_kmeans: bool,
    iter_cap: usize,
    ) -> Result<Vec::<usize>,CluEError>
{
  // Each kluster starts at a distinct data point, so k must be in 1..=n.
  if k == 0 || k > data.len(){
    return Err(CluEError::InvalidNumberOfKMeansClusters(k,data.len()));
  }
  if iter_cap == 0{
    return Err(CluEError::InvalidKMeansIterCap);
  }

  let mut klusters = initialize_klusters(rng, k, data);

  if do_restricted_kmeans{
    restricted_assign_data(&mut klusters, data)?;
  }else{  
    assign_data(&mut klusters, data)?;
  }

  move_kluster_centers(&mut klusters, data);

  let mut is_converged = false;
  for _iter in 0..iter_cap{

    is_converged = if do_restricted_kmeans{
      restricted_assign_data(&mut klusters, data)?
    }else{  
      assign_data(&mut klusters, data)?
    };

    if is_converged{ break; }

    move_kluster_centers(&mut klusters, data);
    
  }

  // Every assignment is a valid partition, so the last one is used.
  if !is_converged{
    println!("Warning: k-means did not converge within kmeans_iter_cap = {} \
iterations; the last assignment is used.", iter_cap);
  }

  let mut element_to_block = (0..data.len()).map(|_| 0)
    .collect::<Vec::<usize>>();

  for (hh,kluster) in klusters.iter().enumerate(){
    for idx in kluster.elements.iter(){
      element_to_block[*idx] = hh;
    }
  }

  Ok(element_to_block)
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// K-Means Clustering Cluster
// Since `Cluster` is already an important data structure in CluE Oxide,
// the k-means cluster is called `Kluster`.
#[derive(Debug,Clone)]
struct Kluster{
  center: Vector3D,
  elements: Vec::<usize>,
}

impl Kluster{
  //----------------------------------------------------------------------------
  // This function returns a `Kluster` centered at `center` with no elements.
  fn at(center: Vector3D) -> Self
  {
    Self{
      center,
      elements: Vec::<usize>::new(),
    }
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq,Eq,Ord,PartialOrd,Hash)]
enum VertexClass{
  Data,
  ClusterCenter,
}

#[derive(Debug,Clone,PartialEq,Eq,Ord,PartialOrd,Hash)]
struct KMNode{
  class: VertexClass,
  index: usize,
}

fn restricted_assign_data(
    klusters: &mut[Kluster],
    data: &[Vector3D], 
    ) 
  -> Result<bool,CluEError>
{ 
  // Clear old assignments.
  let mut old_elements = Vec::<Vec::<usize>>::with_capacity(klusters.len());
  for kluster in klusters.iter_mut(){
    old_elements.push(kluster.elements.clone());
    kluster.elements = Vec::<usize>::new();
  }

  let n_klusters = klusters.len();

  let mut graph = GraphBuilder::<KMNode>::new();

  let src_cap = round::ceil((data.len() as f64)/(n_klusters as f64),0) as i32;

  // Add edge from each kluster to the sink vertex with capacity `src_cap`. 
  for jj in 0..klusters.len(){
    let vertex_c = KMNode{class: VertexClass::ClusterCenter, index: jj};
    graph.add_edge(vertex_c.clone(),Vertex::Sink,Capacity(src_cap),Cost(0));
  }

  let costs = get_restricted_kmeans_costs(klusters, data);

  // Add edge from the source vertex each data point. 
  for (ii,costs_of_datum) in costs.iter().enumerate(){

    let vertex_x = KMNode{class: VertexClass::Data, index: ii};
    graph.add_edge(Vertex::Source,vertex_x.clone(),Capacity(1),Cost(0));
    
    // Add edge from each datum to each kluster.
    for (jj,&cost) in costs_of_datum.iter().enumerate(){
      let vertex_c = KMNode{class: VertexClass::ClusterCenter, index: jj};

      graph.add_edge(vertex_x.clone(),vertex_c.clone(),Capacity(1),Cost(cost));
    }
  }

  // Run min-cost max-flow.
  let (_cost, paths) = graph.mcmf();

  for p in paths.iter(){
    let mcmf::Vertex::Node(kluster) = p.vertices()[2] else{
      return Err(CluEError::FailedKMeans);
    };
    if kluster.class != VertexClass::ClusterCenter{
      return Err(CluEError::FailedKMeans);
    }

    let mcmf::Vertex::Node(datum) = p.vertices()[1] else{
      return Err(CluEError::FailedKMeans);
    };
    if datum.class != VertexClass::Data{
      return Err(CluEError::FailedKMeans);
    }
    
    klusters[kluster.index].elements.push(datum.index);
  }

  for kluster in klusters.iter_mut(){
    kluster.elements.sort();
  }

  for (ii,kluster) in klusters.iter().enumerate(){
    if kluster.elements != old_elements[ii] {
      return Ok(false);
    }
  }
  Ok(true)
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// This function returns the integer cost, costs[i][j], of assigning datum i
// to kluster j in restricted k-means.  The cost is the squared distance
// from the datum to the kluster center, as in constrained k-means
// (P. S. Bradley, K. P. Bennett, and A. Demiriz, "Constrained K-Means 
// Clustering", Microsoft Research Technical Report MSR-TR-2000-65 (2000)),
// so that the assignment and the center update (the mean) minimize the same
// sum of squared distances.  The flow solver needs integer costs, so the 
// squared distances are divided by the largest one and scaled to 
// 0..=max_cost; this makes the costs independent of the unit of length.
// max_cost is at most MAX_RESTRICTED_KMEANS_COST and is lowered for large
// problems to keep LEMON's artificial-arc cost within LEMON_COST_BOUND.
fn get_restricted_kmeans_costs(klusters: &[Kluster], data: &[Vector3D])
  -> Vec::<Vec::<i32>>
{
  let squared_distances = data.iter()
    .map(|x| klusters.iter()
        .map(|kluster| (x - &kluster.center).norm_squared())
        .collect::<Vec::<f64>>())
    .collect::<Vec::<Vec::<f64>>>();

  let max_squared_distance = squared_distances.iter().flatten()
    .fold(0.0_f64, |acc, &d2| f64::max(acc, d2));

  // Nodes: data, klusters, source, sink, and LEMON's root node.
  let n_nodes = (data.len() + klusters.len() + 3) as i64;
  let max_cost = i64::min(MAX_RESTRICTED_KMEANS_COST,
      LEMON_COST_BOUND/n_nodes - 1).max(0) as f64;

  squared_distances.iter()
    .map(|row| row.iter()
        .map(|&d2| {
          if max_squared_distance > 0.0 {
            (d2/max_squared_distance*max_cost).round() as i32
          }else{
            0
          }
        })
        .collect::<Vec::<i32>>())
    .collect::<Vec::<Vec::<i32>>>()
}
//------------------------------------------------------------------------------
fn move_kluster_centers(klusters: &mut [Kluster], data: &[Vector3D])
{
  for kluster in klusters.iter_mut(){
    if kluster.elements.is_empty(){
      continue;
    }
    let mut center = Vector3D::zeros();
    for idx in kluster.elements.iter(){
      center = &center + &data[*idx];
    }
    center = center.scale(1.0/(kluster.elements.len() as f64));

    kluster.center = center;
  }
}
//------------------------------------------------------------------------------
fn assign_data(klusters: &mut [Kluster], data: &[Vector3D])
  -> Result<bool,CluEError>
{

  let mut old_elements = Vec::<Vec::<usize>>::with_capacity(klusters.len());
  for kluster in klusters.iter_mut(){
    old_elements.push(kluster.elements.clone());
    kluster.elements = Vec::<usize>::new();
  }

  for (ii,x) in data.iter().enumerate(){
    let mut min_index = 0;
    let mut min_distance_2 = f64::INFINITY;

    for (hh,kluster) in klusters.iter().enumerate(){
      let distance_2 = (x - &kluster.center).norm_squared();

      if distance_2 < min_distance_2{
        min_distance_2 = distance_2;
        min_index = hh;
      }
    }
    klusters[min_index].elements.push(ii);
  }

  for (ii,kluster) in klusters.iter().enumerate(){
    if kluster.elements != old_elements[ii] {
      return Ok(false);
    }
  }
  Ok(true)
}
//------------------------------------------------------------------------------
// This function places the k initial kluster centers at k distinct data 
// points, chosen uniformly at random (Forgy initialization, E. W. Forgy, 
// Biometrics 21, 768 (1965)). Unlike centers drawn uniformly in the bounding 
// box, every kluster then starts with at least its own point (for distinct
// points), which avoids klusters that stay empty.
// Requires 1 <= k <= data.len().
fn initialize_klusters(
    rng: &mut ChaCha20Rng,
    k: usize,
    data: &[Vector3D],
    ) -> Vec::<Kluster>
{
  sample(rng, data.len(), k).into_iter()
    .map(|idx| Kluster::at(data[idx].clone()))
    .collect::<Vec::<Kluster>>()
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[cfg(test)]
mod tests{
  use super::*;
  use rand::SeedableRng;

  //----------------------------------------------------------------------------
  // TODO: Improve test so that it passes for all seed values.
  // Since k-means can settle into one of multiple minima, only some 
  // seeds give the expected result. 
  #[test]
  fn test_get_kmeans_partition(){
    let k = 4;
    let mut rng = ChaCha20Rng::seed_from_u64(2);
    let (data,_expected_klusters) = get_test_kluster_data();

    let elements_to_block = get_kmeans_partition(&mut rng,k,&data,false, 1000).unwrap();
    assert_eq!(elements_to_block.len(),12); 

    assert_eq!(elements_to_block[0],elements_to_block[1]);
    assert_eq!(elements_to_block[0],elements_to_block[2]);

    assert_eq!(elements_to_block[3],elements_to_block[4]);
    assert_eq!(elements_to_block[3],elements_to_block[5]);

    assert_eq!(elements_to_block[6],elements_to_block[7]);
    assert_eq!(elements_to_block[6],elements_to_block[8]);

    assert_eq!(elements_to_block[9],elements_to_block[10]);
    assert_eq!(elements_to_block[9],elements_to_block[11]);

    assert_ne!(elements_to_block[0],elements_to_block[3]);
    assert_ne!(elements_to_block[0],elements_to_block[6]);
    assert_ne!(elements_to_block[0],elements_to_block[9]);

    assert_ne!(elements_to_block[3],elements_to_block[6]);
    assert_ne!(elements_to_block[3],elements_to_block[9]);

    assert_ne!(elements_to_block[6],elements_to_block[9]);
  
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_restricted_kmeans_partition(){
    let k = 4;
    let mut rng = ChaCha20Rng::seed_from_u64(2);
    let (data,_expected_klusters) = get_test_kluster_data();

    let elements_to_block = get_kmeans_partition(&mut rng,k,&data,true, 1000).unwrap();
    assert_eq!(elements_to_block.len(),12); 

    assert_eq!(elements_to_block[0],elements_to_block[1]);
    assert_eq!(elements_to_block[0],elements_to_block[2]);

    assert_eq!(elements_to_block[3],elements_to_block[4]);
    assert_eq!(elements_to_block[3],elements_to_block[5]);

    assert_eq!(elements_to_block[6],elements_to_block[7]);
    assert_eq!(elements_to_block[6],elements_to_block[8]);

    assert_eq!(elements_to_block[9],elements_to_block[10]);
    assert_eq!(elements_to_block[9],elements_to_block[11]);

    assert_ne!(elements_to_block[0],elements_to_block[3]);
    assert_ne!(elements_to_block[0],elements_to_block[6]);
    assert_ne!(elements_to_block[0],elements_to_block[9]);

    assert_ne!(elements_to_block[3],elements_to_block[6]);
    assert_ne!(elements_to_block[3],elements_to_block[9]);

    assert_ne!(elements_to_block[6],elements_to_block[9]);

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_move_kluster_centers(){

    let (data,expected_klusters) = get_test_kluster_data();

    let mut klusters = vec![
      Kluster{
        center: Vector3D::zeros(),
        elements: vec![0,1,2],
      },
      Kluster{
        center: Vector3D::zeros(),
        elements: vec![3,4,5],
      },
      Kluster{
        center: Vector3D::zeros(),
        elements: vec![6,7,8],
      },
      Kluster{
        center: Vector3D::zeros(),
        elements: vec![9,10,11],
      },
    ];

    move_kluster_centers(&mut klusters, &data);

    for (ii,kluster) in klusters.iter().enumerate(){
      let err = (&kluster.center - &expected_klusters[ii].center).norm();
      assert!( err < 1e-12 );
    }

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_initialize_klusters(){
    let k = 4;
    let mut rng = ChaCha20Rng::from_entropy();
    let (data,_) = get_test_kluster_data();
    let klusters = initialize_klusters(&mut rng, k, &data);
    assert_eq!(klusters.len(), k);

    for kluster in klusters.iter(){
      let center = &kluster.center;
      // Centers are data points, so they can lie on the bounding box.
      assert!(center.x() >= 34.530);
      assert!(center.x() <= 38.970);
      assert!(center.y() >= 34.140);
      assert!(center.y() <= 40.260);
      assert!(center.z() >= 35.140);
      assert!(center.z() <= 39.500);
      assert!(kluster.elements.is_empty());
    }

    // Each center is a distinct data point.
    let mut center_indices = klusters.iter()
      .map(|kluster| data.iter().position(|x| *x == kluster.center).unwrap())
      .collect::<Vec::<usize>>();
    center_indices.sort();
    center_indices.dedup();
    assert_eq!(center_indices.len(), k);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_kmeans_partition_rejects_invalid_k(){
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    let (data,_) = get_test_kluster_data();
    for do_restricted_kmeans in [false, true]{
      for k in [0, data.len()+1]{
        assert_eq!(
            get_kmeans_partition(&mut rng,k,&data,do_restricted_kmeans, 1000).err(),
            Some(CluEError::InvalidNumberOfKMeansClusters(k,data.len())));
      }
      // k = n is allowed: every point is its own kluster.
      let labels = get_kmeans_partition(&mut rng, data.len(), &data,
          do_restricted_kmeans, 1000).unwrap();
      let mut sorted = labels.clone();
      sorted.sort();
      assert_eq!(sorted, (0..data.len()).collect::<Vec::<usize>>());
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_unrestricted_kmeans_has_no_empty_klusters_for_distinct_points(){
    // Before Forgy initialization, random centers in the bounding box
    // left klusters empty for most seeds.
    let (data,_) = get_test_kluster_data();
    for seed in 0..50{
      let mut rng = ChaCha20Rng::seed_from_u64(seed);
      for k in 1..=data.len(){
        let labels = get_kmeans_partition(&mut rng,k,&data,false, 1000).unwrap();
        let mut used = labels.clone();
        used.sort();
        used.dedup();
        assert_eq!(used.len(), k, "seed {}, k {}", seed, k);
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_kmeans_with_duplicate_points(){
    // Duplicate points can share a center, so a kluster may stay empty; 
    // the result must still label every point with a label below k.
    let p = Vector3D::from([1.0,2.0,3.0]);
    let q = Vector3D::from([5.0,2.0,3.0]);
    let data = vec![p.clone(), p.clone(), p, q];
    for seed in 0..20{
      for do_restricted_kmeans in [false, true]{
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let labels = get_kmeans_partition(&mut rng,3,&data,
            do_restricted_kmeans, 1000).unwrap();
        assert_eq!(labels.len(), 4);
        assert!(labels.iter().all(|&label| label < 3));
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_restricted_kmeans_in_si_units(){
    // CluE positions are in metres.  The old cost, round(distance), was 0 
    // for every pair, so the four groups were never found (N5).
    let (data,_) = get_test_kluster_data();
    let data = data.iter().map(|x| x.scale(1e-10))
      .collect::<Vec::<Vector3D>>();

    let mut rng = ChaCha20Rng::seed_from_u64(2);
    let labels = get_kmeans_partition(&mut rng,4,&data,true,1000).unwrap();

    for g in 0..4{
      assert_eq!(labels[3*g],labels[3*g+1]);
      assert_eq!(labels[3*g],labels[3*g+2]);
    }
    let mut group_labels = (0..4).map(|g| labels[3*g])
      .collect::<Vec::<usize>>();
    group_labels.sort();
    assert_eq!(group_labels, vec![0,1,2,3]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_restricted_kmeans_is_independent_of_length_unit(){
    let (data,_) = get_test_kluster_data();
    let data_si = data.iter().map(|x| x.scale(1e-10))
      .collect::<Vec::<Vector3D>>();

    for seed in 0..20{
      for k in [2,3,4,6]{
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let labels = get_kmeans_partition(&mut rng,k,&data,true,1000)
          .unwrap();
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let labels_si = get_kmeans_partition(&mut rng,k,&data_si,true,1000)
          .unwrap();
        assert_eq!(labels, labels_si, "seed {}, k {}", seed, k);
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_restricted_kmeans_respects_capacity(){
    let (data,_) = get_test_kluster_data();
    let n = data.len();
    for seed in 0..20{
      for k in 1..=n{
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let labels = get_kmeans_partition(&mut rng,k,&data,true,1000)
          .unwrap();
        let capacity = (n + k - 1)/k;
        for label in 0..k{
          let size = labels.iter().filter(|&&l| l == label).count();
          assert!(size <= capacity, "seed {}, k {}, label {}", seed, k, label);
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_restricted_kmeans_costs(){
    let data = vec![
      Vector3D::from([0.0,0.0,0.0]),
      Vector3D::from([1e-10,0.0,0.0]),
      Vector3D::from([2e-10,0.0,0.0]),
    ];
    let klusters = vec![Kluster::at(data[0].clone()), 
        Kluster::at(data[2].clone())];
    let costs = get_restricted_kmeans_costs(&klusters, &data);
    // Squared distances 0, 1, 4 (in Å^2), scaled so that 4 -> 10^6.
    assert_eq!(costs, vec![vec![0, 1_000_000], vec![250_000, 250_000],
        vec![1_000_000, 0]]);

    // All data at one point: all costs are 0.
    let data = vec![Vector3D::from([1.0,2.0,3.0]); 3];
    let klusters = vec![Kluster::at(data[0].clone())];
    assert_eq!(get_restricted_kmeans_costs(&klusters, &data), 
        vec![vec![0]; 3]);

    // Large problems lower the maximum cost to keep LEMON's artificial-arc
    // cost, (max_cost + 1)*(number of nodes), within LEMON_COST_BOUND.
    let n = 3000;
    let data = (0..n).map(|ii| Vector3D::from([ii as f64,0.0,0.0]))
      .collect::<Vec::<Vector3D>>();
    let klusters = vec![Kluster::at(data[0].clone())];
    let costs = get_restricted_kmeans_costs(&klusters, &data);
    let max_cost = costs.iter().flatten().max().unwrap();
    let n_nodes = (n + 1 + 3) as i64;
    assert_eq!(*max_cost as i64, LEMON_COST_BOUND/n_nodes - 1);
    assert!((*max_cost as i64 + 1)*n_nodes <= LEMON_COST_BOUND);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_kmeans_iter_cap(){
    let (data,_) = get_test_kluster_data();
    for do_restricted_kmeans in [false, true]{
      // A cap of 1 still returns a valid labelling.
      let mut rng = ChaCha20Rng::seed_from_u64(0);
      let labels = get_kmeans_partition(&mut rng,4,&data,
          do_restricted_kmeans,1).unwrap();
      assert_eq!(labels.len(), data.len());
      assert!(labels.iter().all(|&label| label < 4));

      let mut rng = ChaCha20Rng::seed_from_u64(0);
      assert_eq!(get_kmeans_partition(&mut rng,4,&data,
          do_restricted_kmeans,0).err(), 
          Some(CluEError::InvalidKMeansIterCap));
    }
  }
  //----------------------------------------------------------------------------
  fn get_test_kluster_data() -> (Vec::<Vector3D>,Vec::<Kluster>){  
    let data = vec![
      Vector3D::from([38.970,  35.510,  38.800]),
      Vector3D::from([37.380,  36.130,  39.500]),
      Vector3D::from([38.800,  37.180,  39.010]),

      Vector3D::from([37.110,  34.140,  37.860]),
      Vector3D::from([38.560,  34.210,  37.060]),
      Vector3D::from([37.060,  34.520,  36.030]),

      Vector3D::from([37.570,  38.870,  38.780]),
      Vector3D::from([35.840,  38.870,  38.940]),
      Vector3D::from([36.850,  40.260,  38.120]),

      Vector3D::from([35.260,  38.370,  35.140]),
      Vector3D::from([35.270,  39.920,  36.060]),
      Vector3D::from([34.530,  38.470,  36.820]),
    ];

    let klusters = vec![
      Kluster{
        center: (&data[0] + &(&data[1]+&data[2])).scale(1.0/3.0),
        elements: vec![0,1,2],
      },
      Kluster{
        center: (&data[3] + &(&data[4]+&data[5])).scale(1.0/3.0),
        elements: vec![3,4,5],
      },
      Kluster{
        center: (&data[6] + &(&data[7]+&data[8])).scale(1.0/3.0),
        elements: vec![6,7,8],
      },
      Kluster{
        center: (&data[9] + &(&data[10]+&data[11])).scale(1.0/3.0),
        elements: vec![9,10,11],
      },
    ];

    (data,klusters)
  }
  //----------------------------------------------------------------------------
}
