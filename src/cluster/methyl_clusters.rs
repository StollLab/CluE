use crate::clue_errors::CluEError;
use crate::cluster::Cluster;
use crate::cluster::cluster_set::ClusterSet;
use crate::math;
use crate::structure::{Structure,exchange_groups::ExchangeGroupManager};

use std::collections::HashMap;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This function takes in a `ClusterSet` and sorts the clusters by their
/// `ExchangeGroup` content labeled by "bath" for no groups, 
/// or a string that lists all the groups in the category.
/// For example in the category labeled by "methyl_3_4_5", 
/// every cluster will contains hydrogens 3, 4, and 5, and possibly 
/// hydrogens that are not part of an exchange group, but no cluster
/// will contain a hydrogen from another exchange group.  
/// Clusters that contain the hydrogens from "methyl_3_4_5" and the
/// hydrogens from "methyl_7_8_9" and no other exchange group hydrogens will
/// fall under the category labeled "methyl_3_4_5_methyl_7_8_9".
pub fn partition_cluster_set_by_exchange_groups(cluster_set: ClusterSet,
    exchange_group_manager: &ExchangeGroupManager,
    structure: &Structure)
  -> Result<HashMap::<String,ClusterSet>,CluEError>
{

  let max_size = cluster_set.clusters.len();


  // Count the number of partitions. 
  let n_groups = exchange_group_manager.exchange_groups.len() + 1;

  let mut n_cluster_partitions 
    = HashMap::<String,Vec::<usize>>::with_capacity(n_groups);

  // Loop over cluster sizes.
  for clu_size in 0..max_size{

    // Loop over clusters of a given size.
    for cluster in cluster_set.clusters[clu_size].iter(){

      let key = get_cluster_partition_key(cluster,
          exchange_group_manager,structure)?;

      let mut counts: Vec::<usize> = vec![0;max_size];

      if let Some(n) = n_cluster_partitions.get(&key){
        counts = n.clone();
      }

      counts[clu_size] += 1;

      n_cluster_partitions.insert(key,counts);
    }
  }


  // Collect the partitions.
  let n_groups = n_cluster_partitions.len();
  let mut cluster_partitions 
    = HashMap::<String,Vec::<Vec<Cluster>>>::with_capacity(n_groups);

    // Loop over cluster sizes.
    for clu_size in 0..cluster_set.clusters.len(){

      // Loop over clusters of a given size.
      for cluster in cluster_set.clusters[clu_size].iter(){

        let key = get_cluster_partition_key(cluster,
            exchange_group_manager, structure)?;

        let Some(counts) = n_cluster_partitions.get(&key) else
        {
          return Err(CluEError::InvalidClusterPartitionKey);
        };

        if let Some(clusters) = cluster_partitions.get_mut(&key){
          clusters[clu_size].push(cluster.clone());
        }else{
          let mut clusters = (0..max_size)
            .map(|size| Vec::<Cluster>::with_capacity(counts[size]))
            .collect::<Vec::<Vec::<Cluster>>>();

          clusters[clu_size].push(cluster.clone());
          cluster_partitions.insert(key.clone(), clusters);
        }
      }
    }


  // Transformation: Cluster -> ClusterSet.
  let mut cluster_set_partitions 
    = HashMap::<String,ClusterSet>::with_capacity(n_groups);
  
  for (key,clusters) in cluster_partitions{
    let cluster_set = ClusterSet::from(clusters);
    cluster_set_partitions.insert(key,cluster_set);
  }


  Ok(cluster_set_partitions)
}
//------------------------------------------------------------------------------
// This function takes a Cluster (in tensor indices) 
// and an ExchangeGroupManager (in structure indices) and returns a string
// identifyng the partition (in reference indices).
fn get_cluster_partition_key(cluster: &Cluster,
    exchange_group_manager: &ExchangeGroupManager, structure: &Structure) 
  -> Result<String,CluEError>
{
  let mut key_ids = Vec::<usize>::with_capacity(cluster.len());

  for &vertex in cluster.vertices().iter(){
    
    let structure_idx = structure.get_bath_index_of_nth_active(vertex)?;

    if let Some(id) = exchange_group_manager.exchange_group_ids[structure_idx]{
      key_ids.push(id);
    } 
  }

  key_ids = math::unique(key_ids);

  if key_ids.is_empty(){
    return Ok(String::from("bath"));
  }

  let mut key = exchange_group_manager.exchange_groups[key_ids[0]].to_string();
  for &id in key_ids.iter().skip(1){
    key = format!("{}_{}",key,
        exchange_group_manager.exchange_groups[id]);
  }

  Ok(key)

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// TODO: generalize to other exchange groups
/*
impl ClusterSet{
  /// This function removes all clusters from a cluster set that contain partial
  /// methyl group, 
  /// leaving only clusters that contain either no methyl hydrogens 
  /// or all three hydrogens from any methyl.
  pub fn remove_partial_methyls(&mut self,structure: &Structure) 
    -> Result<(),CluEError>
  {

    let Some(exchange_group_manager) = &structure.exchange_groups else{
      return Ok(());
    };

    for clu_size in 0..self.clusters.len(){

      let n_clusters = self.clusters[clu_size].len();
      let mut to_keep = Vec::<usize>::with_capacity(n_clusters);
      
      for (idx,cluster) in self.clusters[clu_size].iter().enumerate(){
        if cluster.contains_partial_methyl(exchange_group_manager, structure)?{ 
          continue; 
        }
        to_keep.push(idx);
      }

      let mut kept_clusters = Vec::<Cluster>::with_capacity(to_keep.len());
      let mut kept_cluster_indices = HashMap::with_capacity(to_keep.len());
      for (new_idx, old_idx) in to_keep.iter().enumerate(){
        let cluster = self.clusters[clu_size][*old_idx].clone();
        kept_cluster_indices.insert(cluster.vertices.clone(),new_idx);
        kept_clusters.push(cluster);
      }

      self.clusters[clu_size] = kept_clusters;
      self.cluster_indices[clu_size] = kept_cluster_indices;
    }

    Ok(())
  }
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// TODO: generalize to other exchange groups
/*    
// This function determines if the cluster contains a partial methyl.
impl Cluster{
pub fn contains_partial_methyl(&self,
    exchange_group_manager: &ExchangeGroupManager,
    structure: &Structure) -> Result<bool,CluEError>
{
  let mut count_map = HashMap::<usize,usize>::with_capacity(self.len());
  for &vertex in self.vertices().iter(){
    let idx = structure.get_bath_index_of_nth_active(vertex)?;

    if let Some(key) = exchange_group_manager.exchange_group_ids[idx]{
      
        let mut count: usize = 0;
        if let Some(n) = count_map.get(&key){
          count = *n;
        }
        count += 1;

        count_map.insert(key,count);
    }
  }

  for (_key, count) in count_map.iter(){
    if count % 3 != 0 {
      return Ok(true);
    }
  }
  Ok(false)
}
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  use crate::space_3d::Vector3D;
  use crate::structure::exchange_groups::{C3Rotor,ExchangeGroup};
  use crate::structure::particle::Particle;
  use crate::cluster::adjacency::AdjacencyList;
  use crate::elements::Element;
  
  //----------------------------------------------------------------------------
  // This test checks the basics of partition_cluster_set_by_exchange_groups.
  #[test]
  fn test_partition_cluster_set_by_exchange_groups_1(){
    let exchange_group_manager =  generate_test_exchange_group_manager();
    let structure = get_test_structure();

    let clusters = vec![
      vec![Cluster::from(vec![1]),Cluster::from(vec![5])],
      vec![Cluster::from(vec![1,2]),Cluster::from(vec![1,5]),
        Cluster::from(vec![5,9])],
      vec![Cluster::from(vec![1,2,3]),Cluster::from(vec![1,4,5]),
        Cluster::from(vec![4,7,8])],
      vec![Cluster::from(vec![1,2,3,4])],
      vec![Cluster::from(vec![1,2,3,4,5])],
      vec![Cluster::from(vec![1,2,3,4,5,6])],
      vec![Cluster::from(vec![1,2,3,4,5,6,7])],
      vec![Cluster::from(vec![1,2,3,4,5,6,7,8])],
      vec![Cluster::from(vec![1,2,3,4,5,6,7,8,9])],
    ];

    // ref_idx: 1 2 3 4 5 6 7 8 9
    // CH3_id#: 0 0 0 1 - - 1 1 - 
    let mut ref_partitions 
      = HashMap::<String,ClusterSet>::with_capacity(4);


    ref_partitions.insert("bath".to_string(), ClusterSet::from(vec![
      vec![Cluster::from(vec![5])],
      vec![Cluster::from(vec![5,9])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
    ]));

    ref_partitions.insert("methyl_1_2_3".to_string(), ClusterSet::from(vec![
      vec![Cluster::from(vec![1])],
      vec![Cluster::from(vec![1,2]),Cluster::from(vec![1,5])],
      vec![Cluster::from(vec![1,2,3])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
    ]));

    ref_partitions.insert("methyl_7_4_8".to_string(), ClusterSet::from(vec![
      vec![],
      vec![],
      vec![Cluster::from(vec![4,7,8])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
    ]));

    ref_partitions.insert("methyl_1_2_3_methyl_7_4_8".to_string(), 
     ClusterSet::from(vec![
      vec![],
      vec![],
      vec![Cluster::from(vec![1,4,5])],
      vec![Cluster::from(vec![1,2,3,4])],
      vec![Cluster::from(vec![1,2,3,4,5])],
      vec![Cluster::from(vec![1,2,3,4,5,6])],
      vec![Cluster::from(vec![1,2,3,4,5,6,7])],
      vec![Cluster::from(vec![1,2,3,4,5,6,7,8])],
      vec![Cluster::from(vec![1,2,3,4,5,6,7,8,9])],
    ]));

    let ref_counts = clusters.iter().map(|v| v.len())
      .collect::<Vec::<usize>>();

    let cluster_set = ClusterSet::from(clusters);

    let cluster_partitions: HashMap::<String,ClusterSet> = 
      partition_cluster_set_by_exchange_groups(cluster_set,
          &exchange_group_manager, &structure).unwrap();


    let mut counts = ref_counts.iter().map(|_x| 0)
      .collect::<Vec::<usize>>();

    for (_key, value_cluster_set) in cluster_partitions.iter(){
      for (clu_size,clusters) in value_cluster_set.clusters.iter().enumerate(){
        counts[clu_size] += clusters.len();
      }
    }
    assert_eq!(counts,ref_counts);


    for (ref_key,ref_cluster_set) in ref_partitions.iter(){
      assert!(cluster_partitions.contains_key(ref_key));
      let cluster_set = &cluster_partitions[ref_key];
      assert_eq!(*cluster_set, *ref_cluster_set);
    }
    for (clu_key,cluster_set) in cluster_partitions.iter(){
      assert!(ref_partitions.contains_key(clu_key));
      let ref_cluster_set = &ref_partitions[clu_key];
      assert_eq!(*cluster_set, *ref_cluster_set);
    }

  }
  //----------------------------------------------------------------------------
  // This test checks partition_cluster_set_by_exchange_groups when 
  // methyl-methyl clusters are not included.
  #[test]
  fn test_partition_cluster_set_by_exchange_groups_2(){
    let exchange_group_manager =  generate_test_exchange_group_manager();
    let structure = get_test_structure();

    // methyl 1: 1,2,3
    // methyl 2: 4,7,8
    let clusters = vec![
      vec![Cluster::from(vec![1]),Cluster::from(vec![5])],
      vec![Cluster::from(vec![1,2]),Cluster::from(vec![1,5]),
        Cluster::from(vec![5,9])],
      vec![Cluster::from(vec![1,2,3]),Cluster::from(vec![4,7,8])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
    ];

    // ref_idx: 1 2 3 4 5 6 7 8 9
    // CH3_id#: 0 0 0 1 - - 1 1 - 
    let mut ref_partitions 
      = HashMap::<String,ClusterSet>::with_capacity(4);


    ref_partitions.insert("bath".to_string(), ClusterSet::from(vec![
      vec![Cluster::from(vec![5])],
      vec![Cluster::from(vec![5,9])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
    ]));

    ref_partitions.insert("methyl_1_2_3".to_string(), ClusterSet::from(vec![
      vec![Cluster::from(vec![1])],
      vec![Cluster::from(vec![1,2]),Cluster::from(vec![1,5])],
      vec![Cluster::from(vec![1,2,3])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
    ]));

    ref_partitions.insert("methyl_7_4_8".to_string(), ClusterSet::from(vec![
      vec![],
      vec![],
      vec![Cluster::from(vec![4,7,8])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![],
    ]));

    let ref_counts = clusters.iter().map(|v| v.len())
      .collect::<Vec::<usize>>();

    let cluster_set = ClusterSet::from(clusters);

    let cluster_partitions: HashMap::<String,ClusterSet> = 
      partition_cluster_set_by_exchange_groups(cluster_set,
          &exchange_group_manager, &structure).unwrap();


    let mut counts = ref_counts.iter().map(|_x| 0)
      .collect::<Vec::<usize>>();

    for (_key, value_cluster_set) in cluster_partitions.iter(){
      for (clu_size,clusters) in value_cluster_set.clusters.iter().enumerate(){
        counts[clu_size] += clusters.len();
      }
    }
    assert_eq!(counts,ref_counts);


    for (ref_key,ref_cluster_set) in ref_partitions.iter(){
      assert!(cluster_partitions.contains_key(ref_key));
      let cluster_set = &cluster_partitions[ref_key];
      assert_eq!(*cluster_set, *ref_cluster_set);
    }
    for (clu_key,cluster_set) in cluster_partitions.iter(){
      assert!(ref_partitions.contains_key(clu_key));
      let ref_cluster_set = &ref_partitions[clu_key];
      assert_eq!(*cluster_set, *ref_cluster_set);
    }

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_cluster_partition_key(){
    let exchange_group_manager =  generate_test_exchange_group_manager();
    let structure = get_test_structure();

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![1,2,3]),
        &exchange_group_manager, &structure), Ok("methyl_1_2_3".to_string()));

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![1,2,3,5]),
        &exchange_group_manager, &structure), Ok("methyl_1_2_3".to_string()));

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![5]),
        &exchange_group_manager, &structure), Ok("bath".to_string()));

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![2]),
        &exchange_group_manager, &structure), Ok("methyl_1_2_3".to_string()));

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![1,2,3,4]),
        &exchange_group_manager, &structure), 
        Ok("methyl_1_2_3_methyl_7_4_8".to_string()));

    assert_eq!(get_cluster_partition_key(&Cluster::from(vec![4,3]),
        &exchange_group_manager, &structure), 
        Ok("methyl_1_2_3_methyl_7_4_8".to_string()));

  }
  //----------------------------------------------------------------------------
  /*
  #[test]
  fn test_remove_partial_methyls(){

    let clusters = vec![
      vec![Cluster::from(vec![1]),Cluster::from(vec![5])],
      vec![Cluster::from(vec![1,2]),Cluster::from(vec![1,5]),
        Cluster::from(vec![5,9])],
      vec![Cluster::from(vec![1,2,3]),Cluster::from(vec![1,4,5]),
        Cluster::from(vec![4,7,8])],
      vec![Cluster::from(vec![1,2,3,4])],
      vec![Cluster::from(vec![1,2,3,4,5])],
      vec![Cluster::from(vec![1,2,3,4,5,5])],
      vec![Cluster::from(vec![1,2,3,4,5,6,7])],
      vec![Cluster::from(vec![1,2,3,4,5,6,7,8])],
      vec![Cluster::from(vec![1,2,3,4,5,6,7,8,9])],
    ];

    let mut cluster_set = ClusterSet::from(clusters);

    let clusters = vec![
      vec![Cluster::from(vec![5])],
      vec![Cluster::from(vec![5,9])],
      vec![Cluster::from(vec![1,2,3]),Cluster::from(vec![4,7,8])],
      vec![],
      vec![],
      vec![],
      vec![],
      vec![Cluster::from(vec![1,2,3,4,5,6,7,8])],
      vec![Cluster::from(vec![1,2,3,4,5,6,7,8,9])],
    ];
    let ref_cluster_set = ClusterSet::from(clusters);

    let structure = get_test_structure();

    cluster_set.remove_partial_methyls(&structure).unwrap();
    assert_eq!(cluster_set,ref_cluster_set);
  }
  */
  //----------------------------------------------------------------------------
  fn get_test_structure() -> Structure{
    let bath_particles = (0..10).map(|_| 
        {
          let particle = Particle::new(Element::Hydrogen, 0.0, 0.0, 0.0);
          particle
        }).collect::<Vec::<Particle>>();
    let connections = AdjacencyList::with_capacity(1);
    let cell_offsets = Vec::<Vector3D>::new();
    let mut structure = Structure::new(
        bath_particles,
        connections,
        cell_offsets
        );
    
    structure.map_nth_active_to_reference_indices();

    let exchange_group_manager =  generate_test_exchange_group_manager();
    structure.exchange_groups = Some(exchange_group_manager);

    structure
  }
  //----------------------------------------------------------------------------
  /*
  #[test]
  fn test_contains_partial_methyl(){

    let exchange_group_manager =  generate_test_exchange_group_manager();
    let structure = get_test_structure();

    let cluster = Cluster::from(vec![5]);
    assert!(!cluster.contains_partial_methyl(&exchange_group_manager,
          &structure).unwrap());

    let cluster = Cluster::from(vec![1,2,3]);
    assert!(!cluster.contains_partial_methyl(&exchange_group_manager,
          &structure).unwrap());

    let cluster = Cluster::from(vec![1]);
    assert!(cluster.contains_partial_methyl(&exchange_group_manager,
          &structure).unwrap());

    let cluster = Cluster::from(vec![2,3]);
    assert!(cluster.contains_partial_methyl(&exchange_group_manager,
          &structure).unwrap());

    let cluster = Cluster::from(vec![1,2,4]);
    assert!(cluster.contains_partial_methyl(&exchange_group_manager,
          &structure).unwrap());
   
    let cluster = Cluster::from(vec![1,2,3,4,5,6,7,8,9]);
    assert!(!cluster.contains_partial_methyl(&exchange_group_manager,
          &structure).unwrap());
  }
  */
  //----------------------------------------------------------------------------
  fn generate_test_exchange_group_manager() -> ExchangeGroupManager{
    
    let exchange_groups = vec![
      ExchangeGroup::Methyl(C3Rotor{
        center: Vector3D::from([0.0,0.0,0.0]),
        normal: Vector3D::from([0.0,0.0,1.0]),
        indices: [0,1,2]
      }),
      ExchangeGroup::Methyl(C3Rotor{
        center: Vector3D::from([1.0,1.0,1.0]),
        normal: Vector3D::from([1.0,0.0,0.0]),
        indices: [6,3,7]
      })
    ];

    let exchange_group_ids = vec![
      Some(0),Some(0),Some(0),
      Some(1),None, None, 
      Some(1),Some(1),None];

    let exchange_couplings = exchange_group_ids.iter()
      .map(|_id| 0.0).collect::<Vec::<f64>>();

    ExchangeGroupManager{
      exchange_groups,
      exchange_group_ids,
      exchange_couplings
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
