use crate::cluster::{
  adjacency::AdjacencyList,
  Cluster,
  cluster_set::ClusterSet,
};
use crate::clue_errors::CluEError;
use crate::config::Config;
use crate::config::SAVE_FILE_PARTITION_TABLE;
use crate::math;
use crate::quantum::tensors::HamiltonianTensors;
use crate::structure::{
  exchange_groups::GetIndices,
  Structure,
};

use crate::cluster::find_isolated_subgraphs::find_nearly_isolated_components;
use crate::kmeans::get_kmeans_partition;
use crate::math::unique;
use crate::space_3d::Vector3D;
use crate::io::write_data;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use rand_chacha::ChaCha20Rng;


// TODO: Change this to a value in Config.
const DROP_ALL_SPINS_FROM_METHYLS_THAT_ARE_NOT_CLUSTERS: bool = false;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq)]
pub enum PartitioningMethod{
  Particles, 
  ExchangeGroupsAndParticles,
  KMeans(usize),
  RestrictedKMeans(usize),
}
/*
impl PartitioningMethod{
  pub fn from_toml_value(value: toml::Value) -> Result<Self,CluEError> {

    match value{
      toml::Value::String(s) => {
        if s == KEY_PARTITION_PARTICLE.to_string({
            Ok(Self::Particles)
        }else if s == KEY_PARTITION_EX_GROUPS.to_string(){
          Ok(Self::ExchangeGroupsAndParticles)
        }else{
          return Err(CluEError::CannotParsePartitioningMethod(0,s.to_string()));     
        }
      },
      toml::Value::Table(table) => from_toml_table(table),
      _ => Err(CluEError::CannotParsePartitioningMethod(0,value.to_string())),   
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_table(table: toml::Table,unit_of_distance: f64) 
      -> Result<Self,CluEError>
  {
    let mut out = Self::new();
    out.set_from_toml_table(table,unit_of_distance)?;
    Ok(out)
  }
  //----------------------------------------------------------------------------
  pub fn set_from_toml_table(&mut self, 
      table: toml::Table,
      unit_of_distance: f64)
      -> Result<(),CluEError>
  {
    if let Some(value) = table.get(KEY_PARTITIONING_METHOD){
      let Some(array) = value.as_array() else{
        return Err(CluEError::ExpectedTOMLArray(value.type_str().to_string()));
      };
      self.indices = array.iter().filter_map(|v| v.as_integer())
        .map(|n| n as usize).collect::<Vec::<usize>>();
    }

    if let Some(value) = table.get(KEY_SELE_NOT_INDICES){
      let Some(array) = value.as_array() else{
        return Err(CluEError::ExpectedTOMLArray(value.type_str().to_string()));
      }; 
      self.not_indices = array.iter().filter_map(|v| v.as_integer())
        .map(|n| n as usize).collect::<Vec::<usize>>();
    }
  }
  //----------------------------------------------------------------------------
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq)]
pub enum PartitionBlock{
  Unset,
  Detection,
  Bath(usize),
  Preset(usize),
  Inactive,
}
impl PartitionBlock{
  //----------------------------------------------------------------------------
  fn to_number(&self) -> Option<usize>{
    match self{
      PartitionBlock::Unset => None,
      PartitionBlock::Detection => Some(0),
      PartitionBlock::Bath(n) => Some(*n),
      PartitionBlock::Preset(n) => Some(*n),
      PartitionBlock::Inactive => None,
    }
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)]
pub struct PartitionTable{
  pub element_to_block: Vec::<PartitionBlock>,
  pub block_to_elements: Vec<Vec::<usize>>,
}

impl PartitionTable{
  //----------------------------------------------------------------------------
  pub fn number_blocks(&self) -> usize
  { 
    self.block_to_elements.len()
  }
  //----------------------------------------------------------------------------
  pub fn number_elements(&self) -> usize{
    self.element_to_block.len()
  }
  //----------------------------------------------------------------------------
  pub fn from_indices(element_to_block_indices: Vec::<usize>) 
    -> Result<Self,CluEError>
  {
    let element_to_block: Vec::<PartitionBlock> 
      = element_to_block_indices.iter().map(|idx|
        if *idx == 0{
          PartitionBlock::Detection
        }else{
          PartitionBlock::Bath(*idx)
        }
    ).collect();

    Self::from(element_to_block)
  }
  //----------------------------------------------------------------------------
  pub fn from(element_to_block: Vec::<PartitionBlock>) -> Result<Self,CluEError>{
  
  let element_to_num_block = element_to_block.iter()
    .filter_map(|p_blk| p_blk.to_number())
    .collect::<Vec::<usize>>();

  let n_blocks = {
    let partitions = math::unique(element_to_num_block.clone());
    for (ii,p) in partitions.iter().enumerate(){
      if ii != *p{
        return Err(CluEError::PartitionIsIncomplete);
      }
    }
    partitions.len()
  };


  //let mut block_to_elements = Vec<Vec::<usize>>::with_capacity(n_blocks);
  let mut block_to_elements = (0..n_blocks).map(|_| Vec::<usize>::new())
    .collect::< Vec<Vec::<usize>> >();

  for (element_idx, part_block) in element_to_block.iter().enumerate(){
    if *part_block == PartitionBlock::Unset{
      return Err(CluEError::PartitionIsIncomplete);
    }

    let Some(block_idx) = part_block.to_number() else{
      continue;
    };
    block_to_elements[block_idx].push(element_idx);
  }

    Ok(Self{
        element_to_block,
        block_to_elements
    })

  }
  //----------------------------------------------------------------------------
  pub fn to_csv(&self,save_path: &str,structure: &Structure)
    -> Result<(),CluEError> 
  {
    let n = self.element_to_block.len();
    let mut indices = Vec::<usize>::with_capacity(n);
    let mut blocks = Vec::<usize>::with_capacity(n);

    for (active_idx,block) in self.element_to_block.iter().enumerate(){
      if active_idx == 0{continue;}
      let Some(b) = block.to_number() else { continue; };
      let idx = structure.nth_active_to_reference_index[active_idx];
      indices.push(idx);
      blocks.push(b);
    }

    let headers = vec!["index".to_string(),"block".to_string()];
    let filename = format!("{}/{}.csv",save_path, 
        SAVE_FILE_PARTITION_TABLE);
    write_data(&[indices,blocks], &filename, headers)?;


    let filename = &format!("{}/{}.toml",save_path, 
        SAVE_FILE_PARTITION_TABLE);

    let file = match File::create(filename){
      Ok(file) => file,
      Err(err) => return Err(CluEError::CannotCreateFile(
          filename.to_string(), err.to_string())),
    };

    let mut toml_string = "".to_string();

    for (blk,elements) in self.block_to_elements.iter().enumerate(){
      if elements[0] == 0{ continue;}
      let idx = structure.nth_active_to_reference_index[elements[0]];
      let mut elements_string = format!("{}",idx);
      for e in elements.iter().skip(1){
        if *e == 0{ continue;}
        let idx = structure.nth_active_to_reference_index[*e];
        elements_string = format!("{},{}",elements_string,idx);
      }

      toml_string = format!("{}block_{} = [{}]\n",toml_string,blk,elements_string);
    }
    let mut stream = BufWriter::with_capacity(toml_string.len(),file);
    if let Err(err) = stream.write_all(toml_string.as_bytes()){
      return Err(CluEError::CannotWriteFile(filename.to_string(),
            err.to_string()));
    };
    // Dropping a BufWriter ignores write errors, so flush explicitly.
    if let Err(err) = stream.flush(){
      return Err(CluEError::CannotWriteFile(filename.to_string(),
            err.to_string()));
    }

    Ok(())
  }

  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//------------------------------------------------------------------------------
/// This function builds the partition table for use in pCCE.
pub fn get_partition_table(rng: &mut ChaCha20Rng,
    spin_adjacency_list: &AdjacencyList,
    tensors: &HamiltonianTensors,
    structure: &Structure, config: &Config) 
    -> Result<PartitionTable,CluEError>
{

  let Some(partitioning) = &config.partitioning else{
    return Err(CluEError::NoPartitioningMethod)
  };

  // element_to_block: Vec::<PartitionBlock>, has an entry for each particle
  // specified in tensors.  
  // max_block_index specifies the highest block id used.
  let (mut element_to_block, max_block_index) 
      = get_partition_pretable(tensors, structure, config)?;

  match partitioning{
    PartitioningMethod::Particles => set_element_to_block_singles(
        &mut element_to_block, max_block_index),
    PartitioningMethod::ExchangeGroupsAndParticles 
      => set_element_to_block_exchange_groups(
          &mut element_to_block, max_block_index,
          spin_adjacency_list,structure)?,
    PartitioningMethod::KMeans(kluster_size) 
      => set_element_to_block_kmeans(
          rng, 
          &mut element_to_block, max_block_index,
          *kluster_size,false,get_kmeans_iter_cap(config)?,
          spin_adjacency_list,structure)?,
    PartitioningMethod::RestrictedKMeans(kluster_size) 
      => set_element_to_block_kmeans(
          rng, 
          &mut element_to_block, max_block_index,
          *kluster_size,true,get_kmeans_iter_cap(config)?,
          spin_adjacency_list,structure)?,
  };

  PartitionTable::from(element_to_block)
}
//------------------------------------------------------------------------------
// This function extracts and set user specified blocks.
fn get_partition_pretable(tensors: &HamiltonianTensors, 
    structure: &Structure, config: &Config)
    ->  Result<(Vec::<PartitionBlock>, usize),CluEError>{

  
  // Initialize unset partition.  
  let mut element_to_block : Vec::<PartitionBlock>
      = (0..tensors.len()).map(|_| PartitionBlock::Unset).collect();

  // Set Block 0 as the detection block.
  element_to_block[0] = PartitionBlock::Detection;

  let mut max_block_index = 0;

  // partition_table: Vec::<serial: usize, block_id: usize>
  let partition_table = config.get_partition_table()?;

  // The partition_table is user specifed in terms of PDB serial ids,
  // but CluE needs the partition table in terms of the indices used
  // in tensors: &HamiltonianTensors.
  for (serial,blk) in partition_table.iter(){

    // Block 0 is the detection block.
    if *blk == 0{
      return Err(CluEError::InvalidPartitionTableEntry(*serial,*blk,
            "block 0 is reserved for the detected spin; bath blocks start \
at 1".to_string()));
    }

    // PDB serial ids are u32.
    let Ok(serial_u32) = u32::try_from(*serial) else{
      return Err(CluEError::InvalidPartitionTableEntry(*serial,*blk,
            "serial numbers cannot be larger than 4294967295".to_string()));
    };

    // Converert PDB serial id to bath_index used on structure.
    let Some(bath_idx) = structure.map_serial_id_to_bath_index
        .get(&(serial_u32,0_usize))
    else{
      return Err(CluEError::InvalidPartitionTableEntry(*serial,*blk,
            "no bath particle in the primary cell has this serial number"
            .to_string()));
    };

    // Convert bath_index to active_idx, used in tensors.
    // An inactive particle has no tensors and so no block; it is skipped,
    // so that a table written for a larger bath can be reused.
    let Some(active_idx) = structure.bath_indices_to_active_indices[*bath_idx]
    else{
      continue;
    };
    element_to_block[active_idx] = PartitionBlock::Preset(*blk);

    // Update max_block_index.
    if *blk > max_block_index{
      max_block_index = *blk;
    }
  }

  Ok( (element_to_block, max_block_index) )
}

//------------------------------------------------------------------------------
// This function looks for all entries in element_to_block that are unset
// and gives them a unique block id.
fn set_element_to_block_singles(
    element_to_block: &mut [PartitionBlock], mut max_block_index: usize)
{
  for blk in element_to_block.iter_mut(){
    if *blk != PartitionBlock::Unset{
      continue;
    }
    max_block_index += 1;
    *blk = PartitionBlock::Bath(max_block_index);
  }
}
//------------------------------------------------------------------------------
fn set_element_to_block_exchange_groups(
    element_to_block: &mut [PartitionBlock], mut max_id: usize,
    spin_adjacency_list: &AdjacencyList,
    structure: &Structure,
    ) -> Result<(),CluEError>
{
  // Look for an exchange group manager.
  let Some(exchange_group_manager) = &structure.exchange_groups else {
    set_element_to_block_singles(element_to_block, max_id);
    return Ok(());
  };

  // Check that there are exchange groups. 
  if exchange_group_manager.exchange_groups.is_empty(){
    set_element_to_block_singles(element_to_block, max_id);
    return Ok(());
  }
 

  // Loop through all exchange groups.
  for exchange_group in exchange_group_manager.exchange_groups.iter()
  {

    // Get indices for exchange group.
    let indices = exchange_group.indices();
    
    // Check exchange group status.
    if indices.is_empty(){
      continue;
    }

    // Initialize list of spin indices. 
    let mut spins = Vec::<usize>::with_capacity(3);

    // Loop throught the indices of each methyl/ammonium group.
    for &bath_idx in indices.iter(){

      // Get tensor/active index to determine which entry of element_to_block
      // to access. 
      let Some(spin_idx)  = structure.bath_indices_to_active_indices[bath_idx]
          else { continue;};

      match element_to_block[spin_idx]{

        // If the spin has no block, save its index for later.
        PartitionBlock::Unset => spins.push(spin_idx),
        
        // If the spin is in the detection or bath blocks,
        // run some validation.  
        PartitionBlock::Detection | PartitionBlock::Bath(_) =>{
          let Some(block_n) = element_to_block[spin_idx].to_number() else{
            return Err(CluEError::SpinAlreadyPartitioned(bath_idx+1,-1));
          }; 
          return Err(
              CluEError::SpinAlreadyPartitioned(bath_idx+1, block_n as i32));
        },

        // If the spin's block is already set, do nothing.
        PartitionBlock::Preset(_) => (),
        
        // If the spin is inactive, do nothing.
        PartitionBlock::Inactive => (),
      }
    }

    // Get the number of edges in the adjacency subgraph containing the
    // spins.
    let mut num_edges = 0;
    for &spin_idx0 in spins.iter(){
      for &spin_idx1 in spins.iter(){
        if spin_idx0 >= spin_idx1 { continue; }
          num_edges += spin_adjacency_list
            .are_connected(spin_idx0,spin_idx1) as usize;
      }
    }

    // Check if the spins do not form a connected methyl group.
    if spins.len() != 3 || num_edges < 2{
      for spin_idx in spins.iter(){
        if DROP_ALL_SPINS_FROM_METHYLS_THAT_ARE_NOT_CLUSTERS{
          element_to_block[*spin_idx] = PartitionBlock::Inactive;
        }else{
          max_id += 1;
          element_to_block[*spin_idx] = PartitionBlock::Bath(max_id);
        }
      }
    }else{
      // Assign block ID.
      max_id += 1;

      // Write block ID to the partition table.
      for spin_idx in spins.iter(){
        element_to_block[*spin_idx] = PartitionBlock::Bath(max_id);
      }
    }

  }

  // Loop though partition table, 
  // and assign all remaining bath spins a non-zero ID.
  for el_to_blk in element_to_block.iter_mut().skip(1){
    if *el_to_blk != PartitionBlock::Unset{
      continue;
    }
    max_id += 1;
    *el_to_blk = PartitionBlock::Bath(max_id);
  }

  Ok(())    
}
//------------------------------------------------------------------------------
fn set_element_to_block_kmeans(
    rng: &mut ChaCha20Rng,
    element_to_block: &mut [PartitionBlock], mut max_id: usize,
    kluster_size: usize,do_restricted_kmeans: bool, iter_cap: usize,
    spin_adjacency_list: &AdjacencyList,
    structure: &Structure,
    ) -> Result<(),CluEError>
{
  // Separate the bath spins into connected subgraphs, excluding the
  // detected spin (vertex 0); only the vertex sets are needed.
  let components = find_nearly_isolated_components(spin_adjacency_list,
      Some(0));

  // Run k-means clustering on each subgraph and collect the results.
  for component in components.iter(){
    // Vertices of the spin adjacency list are active (tensor) indices.
    // Spins whose block is already set, e.g. by a user partition_table
    // (PartitionBlock::Preset), keep it; only unset spins are clustered.
    let indices = component.iter().copied()
      .filter(|&index| element_to_block[index] == PartitionBlock::Unset)
      .collect::<Vec::<usize>>();

    if indices.is_empty(){
      continue;
    }

    let positions = get_positions_of_active_indices(&indices, structure)?;

    // At least one kluster, so that a subgraph with fewer than about
    // kluster_size/2 spins becomes a single block.
    let number_klusters = usize::max(1,
        (indices.len() as f64 / kluster_size as f64).round() as usize);

    let klusters = get_kmeans_partition(rng, number_klusters, &positions,
        do_restricted_kmeans, iter_cap)?;

    // Klusters can be empty, so the labels need not be 0, 1, ..., k-1.
    // Renumber them so that the block IDs are unique and without gaps.
    let (labels, number_used) = renumber_labels_contiguously(&klusters);

    for (ii,index) in indices.iter().enumerate(){
      element_to_block[*index] 
          = PartitionBlock::Bath(labels[ii] + max_id + 1);
    }
    max_id += number_used;
  }

  Ok(())
}
//------------------------------------------------------------------------------
// This function returns `config.kmeans_iter_cap`, the maximum number of 
// k-means iterations.
fn get_kmeans_iter_cap(config: &Config) -> Result<usize,CluEError>
{
  match config.kmeans_iter_cap{
    Some(iter_cap) => Ok(iter_cap),
    None => Err(CluEError::NoKMeansIterCap),
  }
}
//------------------------------------------------------------------------------
// This function maps the distinct values of `labels`, in increasing order,
// to 0, 1, ..., m-1, where m is the number of distinct values.
// It returns the renumbered labels and m.
fn renumber_labels_contiguously(labels: &[usize]) -> (Vec::<usize>,usize)
{
  let used = unique(labels.to_vec());

  let rank: HashMap::<usize,usize> = used.iter().enumerate()
    .map(|(new_label,&old_label)| (old_label,new_label))
    .collect();

  let renumbered = labels.iter().map(|label| rank[label])
    .collect::<Vec::<usize>>();

  (renumbered, used.len())
}
//------------------------------------------------------------------------------
// This function returns the coordinates of the particles with the input
// active (tensor) indices, where index n >= 1 is the nth active bath
// particle.  Index 0 is the detected spin, which is not a bath particle, and
// is an error.  The active indices are converted to bath indices before
// `bath_particles` is accessed.
fn get_positions_of_active_indices(active_indices: &[usize],
    structure: &Structure) -> Result<Vec::<Vector3D>,CluEError>
{
  let bath_indices = active_indices.iter()
    .map(|&n| {
      if n == 0 {
        return Err(CluEError::CannotFindRefIndexFromNthActive(n));
      }
      structure.get_bath_index_of_nth_active(n)
    })
    .collect::<Result<Vec::<usize>,CluEError>>()?;

  Ok(structure.get_positions(&bath_indices))
}
//------------------------------------------------------------------------------
/// This function takes a graph (as a `&AdjacencyList`), and a 
/// `&PartitionTable`, and builds a coarse-grained graph, 
/// another `AdjacencyList`, where the vertices of the new graph are the
/// blocks of the partition of the old graph.  
/// Edges between vertices A and B are placed if any edges were present in the 
/// original graph that connected an element of block A to an element of block B.
pub fn partition_system(
    neighbor_list: &AdjacencyList,
    partition_table: &PartitionTable) -> Result<AdjacencyList,CluEError>
{

  if neighbor_list.len() != partition_table.number_elements(){
    return Err(CluEError::NeighborListAndPartitionTableAreNotCompatible);
  }

  let n_blocks = partition_table.number_blocks();

  let mut block_list = AdjacencyList::with_capacity(n_blocks);
  for v in neighbor_list.get_active_vertices(){
    let Some(idx) = partition_table.element_to_block[v].to_number() else{
      continue;
    };
    block_list.connect(idx,idx);
  }

  for (particle_idx,partition_blk) in partition_table.element_to_block
      .iter().enumerate(){

    let Some(partition_idx) = partition_blk.to_number() else{
      continue;
    };

    let Some(neighbors) = neighbor_list.get_neighbors(particle_idx) else{
      continue;
    };

    for neighbor in neighbors.iter(){
      let Some(neighbor_partition_idx) = partition_table.element_to_block[*neighbor]
       .to_number() else{
         continue;
       }; 

      block_list.connect( partition_idx, neighbor_partition_idx);
    }
  }

  Ok(block_list)

}
//------------------------------------------------------------------------------
/// Given a set that is partitioned into a number of blocks, 
/// as defined by a `PartitionTable` that maps individual elements to blocks,
/// this function takes a `ClusterSet` of block clusters and a `&PartitionTable`,
/// and retturns `ClusterSet` of element clusters.
//------------------------------------------------------------------------------
pub fn expand_block_clusters(
    mut block_cluster_set: ClusterSet, 
    partition_table: &PartitionTable)
  -> Result<ClusterSet,CluEError>
{
  let n_clusters = count_expanded_clusters(&block_cluster_set,
      partition_table);

  let mut cluster_indices 
      = Vec::<HashMap::<Vec::<usize>,usize>>::with_capacity(n_clusters.len());

  // Loop through cluster sizes.
  for (ii, block_clusters) in block_cluster_set.clusters.iter_mut()
      .enumerate()
  {

    cluster_indices.push(HashMap::<Vec::<usize>,usize>::with_capacity(
          block_cluster_set.cluster_indices[ii].len()));

    // Loop through all cluster of the given size.
    for block_cluster in block_clusters.iter_mut(){

      let Some(index) = block_cluster_set.cluster_indices[ii]
          .get(&block_cluster.vertices) 
      else{
        return Err(CluEError::CannotExpandBlockClusters);    
      };

      *block_cluster 
          = expand_block_cluster(block_cluster.clone(), partition_table); 

      cluster_indices[ii].insert(block_cluster.vertices.clone(),*index);
    }
  }

  block_cluster_set.cluster_indices = cluster_indices;

  Ok(block_cluster_set) 
}
//------------------------------------------------------------------------------
// This function counts the number of clusters of each size that a `ClusterSet`
// of block clusters will expand into when converted to element clusters.
// Note that the total number of clusters should remain constant, but the
// numbers of clusters of a particular size may change.
fn count_expanded_clusters(
    block_cluster_set: &ClusterSet, 
    partition_table: &PartitionTable) -> Vec::<usize>
{
  let mut n_clusters = Vec::<usize>::new();

  // Loop through cluster sizes.
  for clusters in block_cluster_set.clusters.iter(){

    // Loop through all cluster of the given size.
    for cluster in clusters.iter(){

      // Initialize size of the expanded cluster.
      let mut size = 0;

      // Loop though all block indices in the block cluster.
      for block_idx in cluster.vertices.iter(){
        // Accumulate the size of each expanded block.
        size += partition_table.block_to_elements[*block_idx].len();
      }
      // Ensure n_clusters can hold the extra count.
      loop{
        if n_clusters.len() > size{ break; }
        n_clusters.push(0);
      }

      // Count the expanded cluster as the correct size.
      n_clusters[size] += 1;
    }
  } 
  n_clusters
}
//------------------------------------------------------------------------------
// This function expands a block cluster into an element cluster.
fn expand_block_cluster(block_cluster: Cluster, partition_table: &PartitionTable) 
  -> Cluster{

  // Get the block cluster vertices.
  let block_vertices = &block_cluster.vertices;
  
  // Determine how many spins the union of the blocks contain,
  let mut n_elements = 0;
  for block_idx in block_vertices.iter(){
    n_elements += partition_table.block_to_elements[*block_idx].len();
  }

  // Collect the spins in the union of the blocks.
  let mut expanded_vertices = Vec::<usize>::with_capacity(n_elements);
  for block_idx in block_vertices.iter(){
    for element in partition_table.block_to_elements[*block_idx].iter(){
      expanded_vertices.push(*element);
    }
  }

  // Keep most of the block cluster: any signal(not yet calculated),
  // would be the same, and the subcluster indices will also remain the same.
  // Let B be a block cluster that maps to the spin cluster S.
  // For any subcluster, b, of B, the spin cluster, s, that b expands to will
  // be a subcluster of S, since by design the expansion is an injective map 
  // that preserves the cluster hierarchy; 
  // however the map is not surjective and so not all subclusters of S come
  // from a b.  
  // The missed subclusters are those that do not contain all or none of the
  // spins from a block, and are the subclusters that are intended to be
  // elided.
  let mut expanded_cluster = block_cluster;
  expanded_cluster.vertices = expanded_vertices;

  expanded_cluster
}
 
//------------------------------------------------------------------------------


#[cfg(test)]
mod tests{
  use super::*;
  use crate::find_clusters;
  use crate::build_adjacency_list;
  use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};
  use crate::io::FromTOMLString;

  //----------------------------------------------------------------------------
  //    block 0-----block 1     block 2---block 3---block 4
  //     1---0-------4---3        6--------7--8--------9
  //      \ /         \ /
  //       2           5
  //
  //
  fn get_neighbor_list_for_tests() -> (AdjacencyList,PartitionTable){
    let mut neighbor_list = AdjacencyList::with_capacity(10);
    neighbor_list.connect(0,1);
    neighbor_list.connect(0,2);
    neighbor_list.connect(1,2);

    neighbor_list.connect(0,4);

    neighbor_list.connect(3,4);
    neighbor_list.connect(4,5);
    neighbor_list.connect(3,5);

    neighbor_list.connect(6,7);
    neighbor_list.connect(8,9);

    let partition_table 
      = PartitionTable::from_indices(vec![0,0,0,1,1,1,2,3,3,4])
      .unwrap();

    (neighbor_list,partition_table)
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_partition_system(){
    let (neighbor_list, partition_table) = get_neighbor_list_for_tests();

    let block_list = partition_system(&neighbor_list, &partition_table).unwrap();

    assert_eq!(block_list.len(), 5);

    assert_eq!(*(block_list.get_neighbors(0).unwrap()), vec![1]);
    assert_eq!(*(block_list.get_neighbors(1).unwrap()), vec![0]);
    assert_eq!(*(block_list.get_neighbors(2).unwrap()), vec![3]);
    assert_eq!(*(block_list.get_neighbors(3).unwrap()), vec![2,4]);
    assert_eq!(*(block_list.get_neighbors(4).unwrap()), vec![3]);
  }
  //----------------------------------------------------------------------------
  //    block 0-----block 1     block 2---block 3---block 4
  //     1---0-------4---3        6--------7--8--------9
  //      \ /         \ /
  //       2           5
  //
  //
  #[test]
  fn test_expand_block_clusters(){
    let (neighbor_list, partition_table) = get_neighbor_list_for_tests();
    let block_list = partition_system(&neighbor_list, &partition_table).unwrap();
    let block_cluster_set = find_clusters(&block_list,2).unwrap();

    let cluster_set = expand_block_clusters(
        block_cluster_set, &partition_table).unwrap();

    assert_eq!(cluster_set.clusters[1][0].vertices, vec![0,1,2]);
    assert_eq!(cluster_set.clusters[1][1].vertices, vec![3,4,5]);
    assert_eq!(cluster_set.clusters[1][2].vertices, vec![6]);
    assert_eq!(cluster_set.clusters[1][3].vertices, vec![7,8]);
    assert_eq!(cluster_set.clusters[1][4].vertices, vec![9]);

    assert_eq!(cluster_set.clusters[2][0].vertices, vec![0,1,2,3,4,5]);
    assert_eq!(cluster_set.clusters[2][1].vertices, vec![6,7,8]);
    assert_eq!(cluster_set.clusters[2][2].vertices, vec![7,8,9]);


    assert_eq!(cluster_set.cluster_indices[1][&vec![0,1,2]],0);
    assert_eq!(cluster_set.cluster_indices[1][&vec![3,4,5]],1);
    assert_eq!(cluster_set.cluster_indices[1][&vec![6]],2);
    assert_eq!(cluster_set.cluster_indices[1][&vec![7,8]],3);
    assert_eq!(cluster_set.cluster_indices[1][&vec![9]],4);

    assert_eq!(cluster_set.cluster_indices[2][&vec![0,1,2,3,4,5]],0);
    assert_eq!(cluster_set.cluster_indices[2][&vec![6,7,8]],1);
    assert_eq!(cluster_set.cluster_indices[2][&vec![7,8,9]],2);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_count_expanded_clusters(){
    let (neighbor_list, partition_table) = get_neighbor_list_for_tests();
    let block_list = partition_system(&neighbor_list, &partition_table).unwrap();
    let block_cluster_set = find_clusters(&block_list,2).unwrap();

    let n_clusters = count_expanded_clusters(&block_cluster_set,
        &partition_table);

    assert_eq!(n_clusters.len(),7);
    assert_eq!(n_clusters,vec![1,2,1,4,0,0,1]);

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_expand_block_cluster(){
  
    let partition_table 
      = PartitionTable::from_indices(vec![0,0,0,1,1,1,2,3,3,4])
      .unwrap();

    let block_cluster = Cluster::from(vec![0]);
    let cluster = expand_block_cluster(block_cluster, &partition_table);

    assert_eq!(cluster.vertices, vec![0,1,2]);
    assert_eq!(cluster.signal, Ok(None));

    let block_cluster = Cluster::from(vec![1]);
    let cluster = expand_block_cluster(block_cluster, &partition_table);
    assert_eq!(cluster.vertices, vec![3,4,5]);

    let block_cluster = Cluster::from(vec![2]);
    let cluster = expand_block_cluster(block_cluster, &partition_table);
    assert_eq!(cluster.vertices, vec![6]);

    let block_cluster = Cluster::from(vec![3]);
    let cluster = expand_block_cluster(block_cluster, &partition_table);
    assert_eq!(cluster.vertices, vec![7,8]);

    let block_cluster = Cluster::from(vec![4]);
    let cluster = expand_block_cluster(block_cluster, &partition_table);
    assert_eq!(cluster.vertices, vec![9]);

    let block_cluster = Cluster::from(vec![0,1]);
    let cluster = expand_block_cluster(block_cluster, &partition_table);
    assert_eq!(cluster.vertices, vec![0,1,2,3,4,5]);

    let block_cluster = Cluster::from(vec![2,4]);
    let cluster = expand_block_cluster(block_cluster, &partition_table);
    assert_eq!(cluster.vertices, vec![6,9]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_set_element_to_block_exchange_groups(){
    let (adjacency_list,tensors,structure,config) = get_tempo();

    let (mut element_to_block, max_block_index) 
        = get_partition_pretable(&tensors, &structure, &config).unwrap();

    set_element_to_block_exchange_groups(
        &mut element_to_block, max_block_index,
        &adjacency_list, &structure).unwrap();

    let expected = vec![
        PartitionBlock::Detection,
        // C1,C3
        PartitionBlock::Preset(1),PartitionBlock::Preset(1),
        PartitionBlock::Preset(1), // H6,H4,H5
        // C2
        PartitionBlock::Preset(1),PartitionBlock::Preset(1),
        PartitionBlock::Preset(1), // H2,H3,H1
        // C4
        PartitionBlock::Bath(4),PartitionBlock::Bath(5), // H8,H7
        // C5
        PartitionBlock::Bath(6),PartitionBlock::Bath(7), // H10, H9
        // C6
        PartitionBlock::Bath(8),PartitionBlock::Bath(9), // H11,H12
        //C7,C9
        PartitionBlock::Bath(2),PartitionBlock::Bath(2),
        PartitionBlock::Bath(2), // H18,H16,H17 
        // C8
        PartitionBlock::Bath(3),PartitionBlock::Bath(3),
        PartitionBlock::Bath(3),// H14,H13,H15
        PartitionBlock::Bath(10), // N
        // O
    ];
    assert_eq!(element_to_block.len(),expected.len());
    assert_eq!(element_to_block, expected);
  }
  //----------------------------------------------------------------------------
  fn get_tempo() -> (AdjacencyList,HamiltonianTensors,Structure,Config)
  {
    let config = Config::from_toml_string(r##"
      input_structure_file = "./assets/TEMPO.pdb"
      detected_spin.position = [28,29]
      radius = 73.5676
      magnetic_field = 1.2
      replicate_unit_cell = false
      partitioning = "exchange_groups"
      partition_table = [
        3,1,
        4,1,
        5,1,
        7,1,
        8,1,
        9,1,
      ]
      
      pulse_sequence = "hahn"
      tau_increments = [1e-2]
      number_timepoints = [101]
      [pair_cutoffs]
      hahn_mod_depth = 3.23e-5
      hahn_taylor_4 = 6.416238909177711e-11
      
      [[groups]]
      name = "hydrogens"
      selection.elements = ["H"]
      1H.c3_tunnel_splitting = 75e-3
    "##).unwrap();

    let mut rng = ChaCha20Rng::seed_from_u64(0);
    let mut structure = Structure::build_structure(&mut rng, &config).unwrap();

    assert_eq!(structure.number_active(),19);
    structure.map_nth_active_to_reference_indices();

    assert_eq!(structure.bath_particles.len(),29);
    let tensors = HamiltonianTensors::generate(&mut rng, &structure, &config)
        .unwrap();

    let adjacency_list = build_adjacency_list(&tensors, &structure, &config)
        .unwrap();

    (adjacency_list,tensors,structure,config)

  }
  //----------------------------------------------------------------------------
  // Invalid partition_table entries are errors instead of panics or
  // wrapped serial numbers; inactive particles are skipped.
  #[test]
  fn test_get_partition_pretable_checks_entries(){
    let (_adjacency_list,tensors,structure,mut config) = get_tempo();
    let s = |x: &str| x.to_string();

    // Serial 3 (H6) is active; serial 1 (C1, 12C) is not.
    let active = structure.map_serial_id_to_bath_index.get(&(3,0)).unwrap();
    assert!(structure.bath_indices_to_active_indices[*active].is_some());
    let inactive = structure.map_serial_id_to_bath_index.get(&(1,0)).unwrap();
    assert!(structure.bath_indices_to_active_indices[*inactive].is_none());

    let cases = vec![
      (vec![(3,1),(999999,2)], CluEError::InvalidPartitionTableEntry(999999,2,
          s("no bath particle in the primary cell has this serial number"))),
      (vec![(3,0)], CluEError::InvalidPartitionTableEntry(3,0,
          s("block 0 is reserved for the detected spin; bath blocks start \
at 1"))),
      // 4294967299 = 2^32 + 3 used to be read as serial 3.
      (vec![(4294967299,1)], CluEError::InvalidPartitionTableEntry(
          4294967299,1,
          s("serial numbers cannot be larger than 4294967295"))),
    ];
    for (table, expected) in cases{
      config.partition_table
        = Some(crate::config::PartitionTableConfig::IndexBlockList(
              table.clone()));
      let result = get_partition_pretable(&tensors, &structure, &config);
      assert_eq!(result.err(), Some(expected), "{:?}", table);
    }

    // The inactive particle is skipped, and the active one is preset.
    config.partition_table
      = Some(crate::config::PartitionTableConfig::IndexBlockList(
            vec![(1,2),(3,1)]));
    let (element_to_block, max_block_index)
      = get_partition_pretable(&tensors, &structure, &config).unwrap();
    assert_eq!(max_block_index,1);
    let active_idx = structure.bath_indices_to_active_indices[*active].unwrap();
    assert_eq!(element_to_block[active_idx], PartitionBlock::Preset(1));
    assert_eq!(element_to_block.iter()
        .filter(|b| **b != PartitionBlock::Unset).count(), 2); // Detection.

    let err = CluEError::InvalidPartitionTableEntry(999999,2,
        s("no bath particle in the primary cell has this serial number"));
    assert_eq!(err.to_string(), "partition_table entry (serial 999999, \
block 2): no bath particle in the primary cell has this serial number");
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_positions_of_active_indices(){
    // TEMPO: 29 bath particles; the carbons and the oxygen are inactive, 
    // so active indices and bath indices differ.
    let (_adjacency_list,_tensors,structure,_config) = get_tempo();

    let active_indices = (1..=structure.number_active())
      .collect::<Vec::<usize>>();
    let positions 
      = get_positions_of_active_indices(&active_indices, &structure).unwrap();
    assert_eq!(positions.len(), active_indices.len());

    for (&n, r) in active_indices.iter().zip(positions.iter()){
      let bath_idx = structure.nth_active_to_reference_index[n] - 1;
      let particle = &structure.bath_particles[bath_idx];
      assert!(particle.active);
      assert_eq!(*r, particle.coordinates);
    }

    // Active index 1 is H6 (bath index 2), not C3 (bath index 1).
    assert_eq!(structure.nth_active_to_reference_index[1], 3);
    assert_eq!(positions[0], structure.bath_particles[2].coordinates);
    assert_ne!(positions[0], structure.bath_particles[1].coordinates);

    // The last active particle, N, is bath index 27.
    assert_eq!(*positions.last().unwrap(), 
        structure.bath_particles[27].coordinates);

    // Indices beyond the active particles are an error, not a panic.
    let n_out = structure.number_active() + 1;
    assert_eq!(get_positions_of_active_indices(&[n_out], &structure).err(),
        Some(CluEError::CannotFindRefIndexFromNthActive(n_out)));

    // The detected spin (index 0) is not a bath particle.
    assert_eq!(get_positions_of_active_indices(&[0], &structure).err(),
        Some(CluEError::CannotFindRefIndexFromNthActive(0)));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_set_element_to_block_kmeans_terminates_and_assigns_all_spins(){
    let (adjacency_list,tensors,structure,_config) = get_tempo();
    let mut rng = ChaCha20Rng::seed_from_u64(0);

    for do_restricted_kmeans in [false, true]{
      let mut element_to_block = (0..tensors.len())
        .map(|_| PartitionBlock::Unset).collect::<Vec::<PartitionBlock>>();
      element_to_block[0] = PartitionBlock::Detection;

      set_element_to_block_kmeans(&mut rng, &mut element_to_block, 0, 2,
          do_restricted_kmeans, 1000, &adjacency_list, &structure).unwrap();

      assert_eq!(element_to_block[0], PartitionBlock::Detection);
      for blk in element_to_block.iter().skip(1){
        assert!(matches!(blk, PartitionBlock::Bath(id) if *id > 0), 
            "{:?}", blk);
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_renumber_labels_contiguously(){
    // Labels 1 and 3 are unused (empty klusters).
    let (labels, n) = renumber_labels_contiguously(&[4,0,2,2,4,0]);
    assert_eq!(labels, vec![2,0,1,1,2,0]);
    assert_eq!(n, 3);

    let (labels, n) = renumber_labels_contiguously(&[5]);
    assert_eq!(labels, vec![0]);
    assert_eq!(n, 1);

    let (labels, n) = renumber_labels_contiguously(&[0,1,2]);
    assert_eq!(labels, vec![0,1,2]);
    assert_eq!(n, 3);
  }
  //----------------------------------------------------------------------------
  // This function checks that a k-means partition of TEMPO is a valid 
  // partition table: every active bath spin is in a bath block, and the
  // block IDs are 0, 1, ..., n-1 without gaps (`PartitionTable::from`).
  // It returns the number of blocks.
  fn check_kmeans_partition_of_tempo(rng: &mut ChaCha20Rng, 
      kluster_size: usize, do_restricted_kmeans: bool,
      adjacency_list: &AdjacencyList, n_spins: usize, structure: &Structure)
    -> usize
  {
    let mut element_to_block = (0..n_spins)
      .map(|_| PartitionBlock::Unset).collect::<Vec::<PartitionBlock>>();
    element_to_block[0] = PartitionBlock::Detection;

    set_element_to_block_kmeans(rng, &mut element_to_block, 0, kluster_size,
        do_restricted_kmeans, 1000, adjacency_list, structure).unwrap();

    for blk in element_to_block.iter().skip(1){
      assert!(matches!(blk, PartitionBlock::Bath(id) if *id > 0), "{:?}", blk);
    }

    let partition_table = PartitionTable::from(element_to_block).unwrap();
    partition_table.number_blocks()
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_set_element_to_block_kmeans_gives_valid_partition(){
    let (adjacency_list,tensors,structure,_config) = get_tempo();
    let n_spins = tensors.len();
    assert_eq!(n_spins, 20);

    // A graph on the same spins with three isolated subgraphs, 
    // so that block IDs from different subgraphs must not collide.
    let mut three_subgraphs = AdjacencyList::with_capacity(n_spins);
    for ii in 2..=9{ three_subgraphs.connect(1,ii); }
    for ii in 11..=18{ three_subgraphs.connect(10,ii); }
    three_subgraphs.connect(19,19);

    // A 16-spin subgraph followed by a 3-spin subgraph. With 
    // kluster_size = 3, restricted k-means uses k = 5 klusters with 
    // capacity 4 for the 16 spins, so a kluster can stay empty; the 
    // block IDs of the next subgraph must not collide or leave a gap.
    let mut sixteen_and_three = AdjacencyList::with_capacity(n_spins);
    for ii in 2..=16{ sixteen_and_three.connect(1,ii); }
    sixteen_and_three.connect(17,18);
    sixteen_and_three.connect(18,19);

    for seed in 0..20{
      let mut rng = ChaCha20Rng::seed_from_u64(seed);
      for do_restricted_kmeans in [false, true]{
        for graph in [&adjacency_list, &three_subgraphs, &sixteen_and_three]{
          // kluster_size = 100 gives round(n/100) = 0 for every subgraph,
          // which must give one block per subgraph instead of a panic.
          for kluster_size in [1,2,3,4,100]{
            let n_blocks = check_kmeans_partition_of_tempo(&mut rng, 
                kluster_size, do_restricted_kmeans, graph, n_spins, 
                &structure);
            // Detection block plus at least one block per subgraph.
            let n_subgraphs 
              = find_nearly_isolated_components(graph, Some(0)).len();
            assert!(n_blocks > n_subgraphs);
            if kluster_size == 1{
              assert_eq!(n_blocks, n_spins);
            }
            if kluster_size == 100{
              assert_eq!(n_blocks, 1 + n_subgraphs);
            }
          }
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_kmeans_keeps_preset_blocks(){
    // get_tempo() presets the spins with serials 3, 4, 5, 7, 8, 9 (the two
    // methyl groups on C1 and C2) to block 1.
    let (adjacency_list,tensors,structure,mut config) = get_tempo();

    let preset = [3_u32,4,5,7,8,9].iter().map(|&serial| {
      let bath_idx = structure.map_serial_id_to_bath_index[&(serial,0)];
      structure.bath_indices_to_active_indices[bath_idx].unwrap()
    }).collect::<Vec::<usize>>();

    for partitioning in [PartitioningMethod::KMeans(3), 
        PartitioningMethod::RestrictedKMeans(3)]
    {
      config.partitioning = Some(partitioning.clone());
      for seed in 0..10{
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let partition_table = get_partition_table(&mut rng, &adjacency_list,
            &tensors, &structure, &config).unwrap();

        for &idx in preset.iter(){
          assert_eq!(partition_table.element_to_block[idx], 
              PartitionBlock::Preset(1), "{:?}", partitioning);
        }
        let mut block_1 = partition_table.block_to_elements[1].clone();
        block_1.sort();
        let mut expected = preset.clone();
        expected.sort();
        assert_eq!(block_1, expected);

        // Every other bath spin is in a k-means block after the presets.
        for blk in partition_table.element_to_block.iter().skip(1){
          match blk{
            PartitionBlock::Preset(1) => (),
            PartitionBlock::Bath(id) => assert!(*id > 1),
            _ => panic!("{:?}", blk),
          }
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_kmeans_skips_fully_preset_subgraph(){
    let (_adjacency_list,tensors,structure,_config) = get_tempo();
    let n_spins = tensors.len();

    // Subgraph {1,2,3} is fully preset; subgraph {4,...,19} is not.
    let mut graph = AdjacencyList::with_capacity(n_spins);
    graph.connect(1,2);
    graph.connect(2,3);
    for ii in 5..n_spins{ graph.connect(4,ii); }

    let mut element_to_block = (0..n_spins)
      .map(|_| PartitionBlock::Unset).collect::<Vec::<PartitionBlock>>();
    element_to_block[0] = PartitionBlock::Detection;
    for ii in 1..=3{ element_to_block[ii] = PartitionBlock::Preset(1); }

    let mut rng = ChaCha20Rng::seed_from_u64(0);
    set_element_to_block_kmeans(&mut rng, &mut element_to_block, 1, 4,
        true, 1000, &graph, &structure).unwrap();

    for ii in 1..=3{
      assert_eq!(element_to_block[ii], PartitionBlock::Preset(1));
    }
    for blk in element_to_block.iter().skip(4){
      assert!(matches!(blk, PartitionBlock::Bath(id) if *id > 1), "{:?}", blk);
    }
    // 16 spins, kmeans_size 4: restricted k-means gives 4 blocks of 4.
    let table = PartitionTable::from(element_to_block).unwrap();
    assert_eq!(table.number_blocks(), 1 + 1 + 4);
  }
  //----------------------------------------------------------------------------
}
