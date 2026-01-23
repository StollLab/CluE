use crate::config::Config;
//use crate::config::{Config,NeighborCutoff};
use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::quantum::tensors::{
  get_perpendicular_dipole_dipole_frequency, HamiltonianTensors};
use crate::signal::calculate_analytic_restricted_2cluster_signals::{
  hahn_three_spin_modulation_depth,
  hahn_three_spin_fourth_order_coefficient,
};
use crate::structure::Structure;
use crate::structure::exchange_groups::GetIndices;

//------------------------------------------------------------------------------
/// This function builds the clusterization graph, where spins are vertices,
/// and edges are crerated based on criteria defined in `config`.
pub fn build_adjacency_list(tensors: &HamiltonianTensors, 
    structure: &Structure, config: &Config) 
  -> Result<AdjacencyList, CluEError>
{

  let n_spins = tensors.len();
  let mut adjacency_list = AdjacencyList::with_capacity(n_spins);

  for idx0 in 1..n_spins{
    adjacency_list.connect(idx0,idx0);


    for idx1 in idx0+1..n_spins{

      if are_spins_neighbors(idx0,idx1, tensors,structure,config)?{
        adjacency_list.connect(idx0,idx1);
      }

    }
  }

  if config.connect_exchange_groups == Some(true){
    connect_exchange_groups(&mut adjacency_list,structure)?;
  }
  Ok(adjacency_list)
}
//------------------------------------------------------------------------------
// This function determines if an edge should be placed between the vertices
// of spins with tensor indices `idx0` and `idx1`.
fn are_spins_neighbors(idx0: usize,idx1: usize,
    tensors: &HamiltonianTensors, structure: &Structure, config: &Config) 
  -> Result<bool, CluEError>
{

  // CluE can only handle clusters where all spin have the same 
  // spin multiplicity. 
  if tensors.spin_multiplicities[idx0] != tensors.spin_multiplicities[idx1]{
    return Ok(false);
  }

  let Some(zeeman0) = tensors.spin1_tensors.get(idx0) else {
    return Ok(false);
  };

  let Some(zeeman1) = tensors.spin1_tensors.get(idx1) else {
    return Ok(false);
  };

  let delta_zeeman = zeeman1 - zeeman0;
  let sum_zeeman = zeeman1 + zeeman0;

  // Enforce conservation of Zeeman energy.
  if 2.0*delta_zeeman.norm()/sum_zeeman.norm() > 1e-12 {
    return Ok(false);
  }

  let Some(dipdip) = tensors.spin2_tensors.get(idx0,idx1) else {
    return Ok(false);
  };

  let Some(hf0) = tensors.spin2_tensors.get(0,idx0) else {
    return Ok(false);
  };
  
  let Some(hf1) = tensors.spin2_tensors.get(0,idx1) else {
    return Ok(false);
  };

  let delta_hf = (hf0.zz() -hf1.zz()).abs();
  let b = ( -dipdip.xx() - dipdip.yy() ).abs();

  if let Some(cutoff) = &config.neighbor_cutoff_delta_hyperfine_zz{
    if delta_hf < *cutoff {return Ok(false);}
  }
  if let Some(cutoff) = &config.neighbor_cutoff_coupling_xx_yy{
    if b < *cutoff {return Ok(false);}
  }
  if let Some(cutoff) = &config.neighbor_cutoff_point_dipole_perpendicular{
    let structure_index0 = structure.get_bath_index_of_nth_active(idx0)?;
    let structure_index1 = structure.get_bath_index_of_nth_active(idx1)?;

    let r0 = &structure.bath_particles[structure_index0].coordinates;
    let r1 = &structure.bath_particles[structure_index1].coordinates;

    let delta_r = (r1 - r0).norm();

    let gamma_0 = structure.bath_particles[structure_index0].isotope
      .gyromagnetic_ratio();
    let gamma_1 = structure.bath_particles[structure_index1].isotope
      .gyromagnetic_ratio();

    let b_perp = get_perpendicular_dipole_dipole_frequency(
        gamma_0, gamma_1, delta_r);
    if b_perp < *cutoff {return Ok(false);}
  }

  if let Some(cutoff) = &config.neighbor_cutoff_3_spin_hahn_mod_depth{
    let k = hahn_three_spin_modulation_depth(delta_hf,b);
    if k < *cutoff {return Ok(false);}
  }
  if let Some(cutoff) = &config.neighbor_cutoff_3_spin_hahn_taylor_4{
    let kom4 = hahn_three_spin_fourth_order_coefficient(delta_hf,b);
    if kom4 < *cutoff {return Ok(false);}
  }

  if let Some(cutoff) = &config.neighbor_cutoff_distance{
    let structure_index0 = structure.get_bath_index_of_nth_active(idx0)?;
    let structure_index1 = structure.get_bath_index_of_nth_active(idx1)?;

    let r0 = &structure.bath_particles[structure_index0].coordinates;
    let r1 = &structure.bath_particles[structure_index1].coordinates;

    let delta_r = (r1 - r0).norm();

    if delta_r > *cutoff { return Ok(false); } 
  }

  Ok(true)
}
//------------------------------------------------------------------------------
fn connect_exchange_groups(
    adjacency_list: &mut AdjacencyList, structure: &Structure)
  -> Result<(),CluEError>
{
  let Some(manager) = &structure.exchange_groups else{
    return Ok(());
  };

  'ex_grp: for exchange_group in manager.exchange_groups.iter(){

    for idx in (*exchange_group).indices(){
      if !structure.bath_particles[idx].active{ continue 'ex_grp;}
    } 

    let indices = (*exchange_group).indices().iter()
        .filter_map(|idx| structure.bath_indices_to_active_indices[*idx])
        .collect::<Vec::<usize>>();

    for &idx0 in indices.iter(){
      for &idx1 in indices.iter(){
        if idx1 <= idx0{ continue; }

        adjacency_list.connect(idx0,idx1)
      }
    }
  }

  Ok(())
}  




