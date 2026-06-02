use crate::Config;
use crate::config::SAVE_DIR_AUXILIARY_SIGNALS;
use crate::cluster::cluster_set::ClusterSet;
use crate::signal::{Signal, load_batch_signals, write_batch_signals};
use crate::structure::Structure;
use crate::HamiltonianTensors;
use crate::CluEError;
use crate::physical_constants::{ONE,PI};
use crate::math;

use num_complex::Complex;
use rayon::prelude::*;
use std::path::Path;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This function calculates 2-CCE, using an analytic solution
/// to the 3-spin Hamiltonian from
/// W. Witzel and S. Das Sarma, “Quantum theory for electron spin decoherence 
/// induced by nuclear spin dynamics in semiconductor quantum computer 
/// architectures: Spectral diffusion of localized electron spins in the 
/// nuclear solid-state environment,” Phys. Rev. B 74, 035322 (2006).
pub fn calculate_analytic_restricted_2cluster_signals(
    cluster_set: &mut ClusterSet, 
    tensors: &HamiltonianTensors, config: &Config, 
    save_path_opt: &Option<String>,
    structure: &Structure,
    ) ->Result<Vec::<Signal>, CluEError>
{

  
  let clusters = &mut cluster_set.clusters;

  let cluster_size = 2;
  if clusters.len() != cluster_size {
    return Err(CluEError::WrongClusterSizeForAnalyticCCE(clusters.len()));
  }

  let Some(batch_size) = config.cluster_batch_size else{
    return Err(CluEError::NoClusterBatchSize);
  };
  
  
  let n_tot = config.get_total_number_timesteps();
  let mut signal = Signal::ones(n_tot);
  let mut order_n_signals = Vec::<Signal>::with_capacity(2);
  order_n_signals.push(signal.clone());

  let n_clusters = clusters[cluster_size-1].len();
  let n_batches = math::ceil( (n_clusters as f64)/(batch_size as f64)) as usize;

  for ibatch in 0..n_batches{
    let idx = ibatch*batch_size;

    if let Some(path) = &save_path_opt{
      if config.write_auxiliary_signals == Some(true){
        let load_dir = format!("{}/{}",path, SAVE_DIR_AUXILIARY_SIGNALS);
          let aux_filename = format!("{}/cluster_size_{}_batch_{}.csv",
              load_dir, cluster_size,ibatch);

          if Path::new(&aux_filename).exists(){
             load_batch_signals(&mut clusters[cluster_size - 1],
                 idx,batch_size,&aux_filename)?;

             continue;
          }
        }
    }

    clusters[1].par_iter_mut().skip(idx).take(batch_size).for_each(
        |cluster| {
          cluster.signal = analytic_restricted_2cluster_signal(
            cluster.vertices(), tensors, config);
        });

    let idx_end = if idx+batch_size < clusters[cluster_size-1].len() {
      idx + batch_size
    }else{
      clusters[1].len()
    };

    for ii in idx..idx_end{
      let Ok(Some(aux_signal)) = &clusters[1][ii].signal else{
        return Err(CluEError::ClusterHasNoSignal(
              clusters[cluster_size-1][ii].to_string()));
      };

      signal = &signal * aux_signal;
    }

    if let Some(path) = &save_path_opt{
        if config.write_auxiliary_signals == Some(true){
          let save_dir = format!("{}/{}",path, SAVE_DIR_AUXILIARY_SIGNALS);
          match std::fs::create_dir_all(save_dir.clone()){
            Ok(_) => (),
            Err(_) => return Err(CluEError::CannotCreateDir(save_dir)),
          }
          let aux_filename = format!("{}/cluster_size_{}_batch_{}.csv",
              save_dir, cluster_size,ibatch);

          write_batch_signals(&clusters[cluster_size-1], n_tot, idx, batch_size,
              &aux_filename, structure)?;
        }
    }
  }

  order_n_signals.push(signal);

  Ok(order_n_signals)
}
//------------------------------------------------------------------------------
/// This function applies Witzel Das Sarma's analytic solution to calculate
/// the 2-CCE auxiliary signal the specified cluster.
pub fn analytic_restricted_2cluster_signal(vertices: &[usize],
    tensors: &HamiltonianTensors, config: &Config) 
  -> Result<Option<Signal>,CluEError>
{
  if vertices.len() != 2 {
    return Err(CluEError::WrongClusterSizeForAnalyticCCE(vertices.len()));
  }
  let idx0 = vertices[0];
  let idx1 = vertices[1];

  let Some(dipdip) = &tensors.spin2_tensors.get(idx0,idx1) else {
    return Ok(Some(Signal::new()));
  };
  let Some(hf0) =  &tensors.spin2_tensors.get(0,idx0) else {
    return Ok(Some(Signal::new()));
  }; 
  let Some(hf1) =  &tensors.spin2_tensors.get(0,idx1) else {
    return Ok(Some(Signal::new()));
  }; 

  // For a point dipole b_zz = -(b_xx + b_yy), but
  // -(b_xx + b_yy) account for isotropic coupling as well.
  let b = -(dipdip.xx() + dipdip.yy() );
  let delta_hf = (hf0.zz() - hf1.zz()).abs();
  
  //let omega = (delta_hf*delta_hf + b*b).sqrt();
  //let k = (2.0*b*delta_hf/( omega*omega )).powi(2);
  let omega = 2.0*PI*hahn_three_spin_modulation_frequency(delta_hf,b);
  let k = hahn_three_spin_modulation_depth(delta_hf,b);

  let tau_axis = config.get_tau_axis_as_ref()?;
  let nt = tau_axis.len();
  let mut data = Vec::<Complex<f64>>::with_capacity(nt);

  for t in tau_axis.iter(){
    let twotau = 2.0*(*t);
    let s4 = (omega*twotau).sin().powi(4);
    data.push(ONE - k*s4 );
  }

  Ok(Some(Signal{data}))

}
//------------------------------------------------------------------------------
/// This function calculate the modulation depth for a single 2-cluster.
pub fn hahn_three_spin_modulation_depth(delta_hf: f64,b:f64) -> f64{
  (2.0*b*delta_hf/( delta_hf*delta_hf + b*b) ).powi(2)
}
//------------------------------------------------------------------------------
/// This function return the modulation frequency to the three spin Hahn echo
/// in unit matching the inputs, which are assumed to have the same units as
/// each other.
pub fn hahn_three_spin_modulation_frequency(delta_hf: f64,b:f64) -> f64{
  (delta_hf*delta_hf + b*b).sqrt()/8.0
}
//------------------------------------------------------------------------------
/// This function calculates the 4th order Taylor series coefficient 
/// for a single 2-cluster.
pub fn hahn_three_spin_fourth_order_coefficient(delta_hf: f64,b:f64) -> f64{
  let k = hahn_three_spin_modulation_depth(delta_hf,b);
  let omega = 2.0*PI*hahn_three_spin_modulation_frequency(delta_hf,b);
  k*(omega.powi(4))
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#[cfg(test)]
mod tests{
  use super::*;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  #[test]
  fn test_hahn_three_spin_modulation_depth(){

    let delta_hfs = [1e-2, 1e-1, 1.0, 1e1, 1e2];
    let bs = [1e-2, 1e-1, 1.0, 1e1, 1e2];

    for &delta_hf in  delta_hfs.iter(){
      assert_eq!(hahn_three_spin_modulation_depth(delta_hf,delta_hf), 1.0);

      for &b in bs.iter(){
        assert!(hahn_three_spin_modulation_depth(delta_hf,b) <= 1.0);
      }
    }
  }
  //----------------------------------------------------------------------------
}







