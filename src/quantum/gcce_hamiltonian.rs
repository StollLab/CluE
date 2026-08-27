use crate::CluEError;
use crate::config::{
  Config,
};
use crate::config::pulse_sequence::PulseSequence;
use crate::math::cxmat_pow_n;
use crate::physical_constants::{BOLTZMANN,HBAR,I,ZERO};
use crate::quantum::spin_hamiltonian::{
  get_propagators_from_eig,
  get_propagators_complex_time_from_eig,
};
use crate::quantum::pulse_sequences::{
  PulseStep,  
};
use crate::quantum::cluster_operators::{
  ClusterSpinOperators,
  SpinOp,
};
use crate::signal::Signal;
use crate::quantum::tensors::HamiltonianTensors;

use num_complex::Complex;
use ndarray::{Array1,Array2};
use ndarray_linalg::{Eigh,Trace,UPLO};

type CxMat = Array2::<Complex<f64>>;

pub fn propagate_custom_pulse_sequence(
    pulse_sequence: &[PulseStep],
    density_matrix: &CxMat, //hamiltonian: &CxMat, 
    h_eigvals: &Array1::<f64>, h_eigvecs: &CxMat,
    config: &Config)
  -> Result<Signal,CluEError>
{

  let tau_increments = &config.tau_increments;
  if tau_increments.is_empty(){
    return Err(CluEError::NoTimeIncrements);
  }

  let number_timepoints = &config.number_timepoints;
  if number_timepoints.is_empty(){
    return Err(CluEError::NoTimepoints);
  }

  let tau2_increments = &config.tau2_increments;
  if tau2_increments.is_empty(){
    return Err(CluEError::NoTimeIncrements2);
  }

  let number_timepoints2 = &config.number_timepoints2;
  if number_timepoints2.is_empty(){
    return Err(CluEError::NoTimepoints2);
  }


  let n_tot = config.get_total_number_timesteps();

  let mut signal = Vec::<Complex<f64>>::with_capacity(n_tot);

  let dus = get_propagators_from_eig(h_eigvals,h_eigvecs,tau_increments)?;

  let du2s = get_propagators_from_eig(h_eigvals,h_eigvecs,tau2_increments)?;

  let mut u_of_tau = CxMat::eye(dus[0].dim().0);
  for (idt,_dt) in tau_increments.iter().enumerate(){
    let n_timepoints = number_timepoints[idt];
    for _inumt in 0..n_timepoints{

    let mut u_of_tau2 = CxMat::eye(du2s[0].dim().0);
      for (idt2,_dt2) in tau2_increments.iter().enumerate(){
        let n_timepoints2 = number_timepoints2[idt2];
        
        for _inumt2 in 0..n_timepoints2{
    
          let mut v = ZERO;
          let mut u_sequence = CxMat::eye(dus[0].dim().0);

          for step in pulse_sequence.iter(){
            match step{
              PulseStep::Pulse(p) => u_sequence = p.dot(&u_sequence),  
              PulseStep::FixedDelay(number, index_opt) =>{
                let delay_idx = match index_opt{ 
                  Some(i) => *i,
                  None => idt,
                };

                let u_step = cxmat_pow_n(&dus[delay_idx],*number);
                u_sequence = u_step.dot(&u_sequence); 

              }, 
              PulseStep::InvFixedDelay(number,index_opt) =>{
                let delay_idx = match index_opt{ 
                  Some(i) => *i,
                  None => idt,
                };

                let u = dus[delay_idx].t().map(|u_ij| u_ij.conj() );
                let u_step = cxmat_pow_n(&u,*number);
                u_sequence = u_step.dot(&u_sequence); 
              }, 
              PulseStep::FixedDelay2(number, index_opt) =>{
                let delay_idx = match index_opt{ 
                  Some(i) => *i,
                  None => idt2,
                };

                let u_step = cxmat_pow_n(&du2s[delay_idx],*number);
                u_sequence = u_step.dot(&u_sequence); 

              }, 
              PulseStep::InvFixedDelay2(number,index_opt) =>{
                let delay_idx = match index_opt{ 
                  Some(i) => *i,
                  None => idt2,
                };

                let u = du2s[delay_idx].t().map(|u_ij| u_ij.conj() );
                let u_step = cxmat_pow_n(&u,*number);
                u_sequence = u_step.dot(&u_sequence); 
              }, 
              PulseStep::TauDelay => {
                  u_sequence = u_of_tau.dot(&u_sequence); 
              },
              PulseStep::InvTauDelay => {
                  let u_of_tau_dag = u_of_tau.t().map(|u_ij| u_ij.conj() );
                  u_sequence = u_of_tau_dag.dot(&u_sequence); 
              },
              PulseStep::Tau2Delay => {
                  u_sequence = u_of_tau2.dot(&u_sequence); 
              },
              PulseStep::InvTau2Delay => {
                  let u_of_tau2_dag = u_of_tau2.t().map(|u_ij| u_ij.conj() );
                  u_sequence = u_of_tau2_dag.dot(&u_sequence); 
              },
              PulseStep::Detect(detection_op0) => {
                let u_sequence_dag = u_sequence.t().map(|u_ij| u_ij.conj() );
                let detection_op = u_sequence_dag.dot(&detection_op0.dot(&u_sequence));
                let it = std::iter::zip(density_matrix,&detection_op);
                v += it.map(|(rho_ij,u_ij)| rho_ij*u_ij).sum::<Complex<f64>>();
              },
              /*
              PulseStep::Integrate(delay,detection_op0) =>{
                let delay_idx = match delay.index{ 
                  Some(i) => i,
                  None => idt,
                };
    
                for _ii in 0..delay.number{
                  let u_sequence_dag = u_sequence.t().map(|u_ij| u_ij.conj() );
                  let detection_op = u_sequence_dag.dot(&detection_op0.dot(&u_sequence));
                  let it = std::iter::zip(density_matrix,&detection_op);
                  v += it.map(|(rho_ij,u_ij)| rho_ij*u_ij).sum::<Complex<f64>>();
                  u_sequence = dus[delay_idx].dot(&u_sequence); 
                }
              },
             */ 
              
            }
          }

          signal.push(v);

          u_of_tau2 = du2s[idt2].dot(&u_of_tau2);
        }
      }
      u_of_tau = dus[idt].dot(&u_of_tau);

    }
  }

  let v0 = signal[0];
  for v in signal.iter_mut(){
    *v /= v0;
  }

  Ok(Signal{data: signal})
}
/// This function is still under construction.
/// Please be patient.
pub fn propagate_pulse_sequence(
    pulses: &[PulseStep],
    density_matrix: &CxMat, //hamiltonian: &CxMat, 
    h_eigvals: &Array1::<f64>, h_eigvecs: &CxMat,
    config: &Config)
  -> Result<Signal,CluEError>
{  
  let Some(pulse_sequence_name) = &config.pulse_sequence else{
    return Err(CluEError::NoPulseSequence);
  };

  let PulseStep::Pulse(u_half_pi) =  pulses[0] else{
    return Err(CluEError::InvalidPulse("pi/2 pulse".to_string()))
  };
  let PulseStep::Pulse(u_pi) =  pulses[1] else{
    return Err(CluEError::InvalidPulse("pi pulse".to_string()))
  };
  let PulseStep::Detect(detection_op0) =  pulses[2] else{
    return Err(CluEError::InvalidPulse("detection operator".to_string()))
  };


  let (_dus,u_of_taus) = get_free_evolutions_propagators(
      h_eigvals, h_eigvecs, config)?; 

  if u_of_taus[0].shape() != density_matrix.shape(){
    return Err(CluEError::PropagatorAndDensityNotSameDimension(
          u_of_taus[0].shape()[0],density_matrix.shape()[0]));
  }

  let (_du2s,u_of_tau2s) = match pulse_sequence_name{
    PulseSequence::RefocusedHahnEcho =>
      get_2nd_free_evolutions_propagators(h_eigvals, h_eigvecs, config)?,

    _ => (vec![CxMat::eye(0)],vec![CxMat::eye(0)]),
  };  

  let mut signal = Vec::<Complex<f64>>::with_capacity(u_of_taus.len());

  for u_of_tau in u_of_taus.iter(){
    
    let mut u_sequence = u_half_pi.clone();

    match pulse_sequence_name{
      // U(𝜏) = U0(𝜏)  
      PulseSequence::CarrPurcell(0) => // FID
        u_sequence = u_of_tau.dot(&u_sequence),

      // U(2𝜏) = U0(𝜏) Upi U0(𝜏) 
      PulseSequence::CarrPurcell(1) => { // Hahn echo
        let u_seg = u_of_tau.dot(&u_pi.dot(u_of_tau));
        u_sequence = u_seg.dot(&u_sequence);
      },

      // U(4𝜏) = U0(𝜏) Upi U0(2𝜏) Upi U0(𝜏) = ( U0(𝜏) Upi U0(𝜏) )^2
      PulseSequence::CarrPurcell(2) => {// CP-2  
        let mut u_seg = u_of_tau.dot(&u_pi.dot(u_of_tau));
        u_seg = u_seg.dot(&u_seg);     
        u_sequence = u_seg.dot(&u_sequence);
        },

      // U(2n𝜏) = ( U0(𝜏) Upi U0(𝜏) )^n
      PulseSequence::CarrPurcell(n_pi) => { // CP-n
        let u_seg = u_of_tau.dot(&u_pi.dot(u_of_tau));
        u_sequence = cxmat_pow_n(&u_seg,*n_pi);
      }

      PulseSequence::RefocusedHahnEcho =>{
        let u_seg1 = u_of_tau.dot(&u_pi.dot(u_of_tau));
        for u_of_tau2 in u_of_tau2s.iter(){
          let u_seg2 = u_of_tau2.dot(&u_pi.dot(u_of_tau2));
          u_sequence = u_seg2.dot(&u_seg1.dot(&u_sequence));
          
          let u_sequence_dag = u_sequence.t().map(|u_ij| u_ij.conj() );
    
          // Find time evolved detection operator:
          // D(T) = U(T)^* D U(T). 
          let detection_op = u_sequence_dag.dot(&detection_op0.dot(&u_sequence));

          let it = std::iter::zip(density_matrix,&detection_op);
          let v = it.map(|(rho_ij,u_ij)| rho_ij*u_ij).sum::<Complex<f64>>();
          signal.push(v);
        }
      },
      PulseSequence::FreeEvolution => u_sequence = u_of_tau.clone(),
        
      PulseSequence::Custom(_) => return Err(
          CluEError::PulseSequenceNotSupported(
            "propagate_pulse_sequence".to_string(),
            "custom pulse sequences".to_string()
            )),
    }

    if *pulse_sequence_name == PulseSequence::RefocusedHahnEcho { continue; }

    let u_sequence_dag = u_sequence.t().map(|u_ij| u_ij.conj() );
    
    // Find time evolved detection operator:
    // D(T) = U(T)^* D U(T). 
    let detection_op = u_sequence_dag.dot(&detection_op0.dot(&u_sequence));

    let it = std::iter::zip(density_matrix,&detection_op);
    let v = it.map(|(rho_ij,u_ij)| rho_ij*u_ij).sum::<Complex<f64>>();
    signal.push(v);
  }

  let v0 = signal[0];
  for v in signal.iter_mut(){
    *v /= v0;
  }

  Ok(Signal{data: signal})
}
//------------------------------------------------------------------------------
pub fn get_free_evolutions_propagators(//hamiltonian: &CxMat,
    eigvals: &Array1::<f64>, eigvecs: &CxMat,
    config: &Config) -> Result<(Vec::<CxMat>,Vec::<CxMat>),CluEError>
{
  let tau_increments = &config.tau_increments;
  if tau_increments.is_empty(){
    return Err(CluEError::NoTimeIncrements);
  }

  let number_timepoints = &config.number_timepoints;
  if number_timepoints.is_empty(){
    return Err(CluEError::NoTimepoints);
  }
  let n_tot = config.get_total_number_timesteps();

  let dus = get_propagators_from_eig(eigvals,eigvecs,tau_increments)?;

  let mut u_of_taus = Vec::<CxMat>::with_capacity(n_tot);
  u_of_taus.push( CxMat::eye(dus[0].dim().0) );

  for (idt,_dt) in tau_increments.iter().enumerate(){
    let n_timepoints = number_timepoints[idt];
    for ii in 1..n_timepoints{
      u_of_taus.push( dus[idt].dot(&u_of_taus[ii-1]));
    }
  }

  Ok((dus,u_of_taus))
}
//------------------------------------------------------------------------------
pub fn get_2nd_free_evolutions_propagators(//hamiltonian: &CxMat,
    eigvals: &Array1::<f64>, eigvecs: &CxMat,
    config: &Config) -> Result<(Vec::<CxMat>,Vec::<CxMat>),CluEError>
{
  let tau_increments = &config.tau2_increments;
  if tau_increments.is_empty(){
    return Err(CluEError::NoTimeIncrements);
  }

  let number_timepoints = &config.number_timepoints2;
  if number_timepoints.is_empty(){
    return Err(CluEError::NoTimepoints);
  }
  let n_tot = config.get_total_number_timesteps();

  let dus = get_propagators_from_eig(eigvals,eigvecs,tau_increments)?;

  let mut u_of_taus = Vec::<CxMat>::with_capacity(n_tot);
  u_of_taus.push( CxMat::eye(dus[0].dim().0) );

  for (idt,_dt) in tau_increments.iter().enumerate(){
    let n_timepoints = number_timepoints[idt];
    for ii in 1..n_timepoints{
      u_of_taus.push( dus[idt].dot(&u_of_taus[ii-1]));
    }
  }

  Ok((dus,u_of_taus))
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/*
fn build_pump_pulse(spin_indices: &[usize], pumped_spins: &[usize],
    spin_ops: &ClusterSpinOperators, config: &Config)
  -> Option<CxMat>
{

  // Loop through spins and build up Hamiltonian,
  for (sop_idx0, &ten_idx0) in spin_indices.iter().enumerate(){
    if !pumped_spins.contains(ten_ind0){ continue} 
  }

}
*/
//------------------------------------------------------------------------------
/// This function builds the cluster spin Hamiltonian without assuming
/// < mS | H | mS' > = 0, for mS != mS',
pub fn build_spin_hamiltonian(spin_indices: &[usize],
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors)
  -> Result<(Array1::<f64>,CxMat),CluEError>
{

  // Get list of spin multiplicities.
  let spin_multiplicities: Vec::<usize> =
    spin_indices.iter().map(|idx| tensors.spin_multiplicities[*idx]).collect();

  // Find Hilbert space dimensionality.
  let mut dim: usize = 1;
  spin_multiplicities.iter().for_each(|spin_mul| dim *= spin_mul);

  // Initialize Hamiltonian.
  let mut ham = CxMat::zeros((dim,dim));

  let cluster_size = spin_indices.len();

  // Loop through spins and build up Hamiltonian,
  for (sop_idx0, &ten_idx0) in spin_indices.iter().enumerate(){


    let spin_mult0 = if spin_multiplicities.len() == 1{
      spin_multiplicities[0]
    }else{
      spin_multiplicities[1]
    };
    //let spin_mult0 = tensors.spin_multiplicities[ten_idx0];
    let sx0 = spin_ops.get(&SpinOp::Sx,spin_mult0,cluster_size,sop_idx0)?;
    let sy0 = spin_ops.get(&SpinOp::Sy,spin_mult0,cluster_size,sop_idx0)?;
    let sz0 = spin_ops.get(&SpinOp::Sz,spin_mult0,cluster_size,sop_idx0)?;

    // Zeeman
    if let Some(vec) = tensors.spin1_tensors.get(ten_idx0){
      ham = ham + sx0*vec.x();
      ham = ham + sy0*vec.y();
      ham = ham + sz0*vec.z();
    }

    if let Some(vec) = tensors.get_mean_field_couplings(ten_idx0,spin_indices){
      ham = ham + sx0*vec.x();
      ham = ham + sy0*vec.y();
      ham = ham + sz0*vec.z();
    }

    for (sop_idx1, &ten_idx1) in spin_indices.iter().enumerate().skip(sop_idx0){

      let sx1 = spin_ops.get(&SpinOp::Sx,spin_mult0,cluster_size,sop_idx1)?;
      let sy1 = spin_ops.get(&SpinOp::Sy,spin_mult0,cluster_size,sop_idx1)?;
      let sz1 = spin_ops.get(&SpinOp::Sz,spin_mult0,cluster_size,sop_idx1)?;

      // hyperfine, dipole-dipole, electric quadrupole, and zero-field
      if let Some(ten) = tensors.spin2_tensors.get(ten_idx0,ten_idx1){
        ham = ham + sx0.dot(sx1)*ten.xx();
        ham = ham + sx0.dot(sy1)*ten.xy();
        ham = ham + sx0.dot(sz1)*ten.xz();

        ham = ham + sy0.dot(sx1)*ten.yx();
        ham = ham + sy0.dot(sy1)*ten.yy();
        ham = ham + sy0.dot(sz1)*ten.yz();

        ham = ham + sz0.dot(sx1)*ten.zx();
        ham = ham + sz0.dot(sy1)*ten.zy();
        ham = ham + sz0.dot(sz1)*ten.zz();
      }
    }
  }

  let Ok((eigvals, eigvecs)) = ham.eigh(UPLO::Lower) else{
    return Err(
        CluEError::CannotDiagonalizeHamiltonian(ham.to_string()));
  };

  Ok((eigvals, eigvecs))
}
//------------------------------------------------------------------------------
pub fn get_electron_cluster_thermal_density_matrix(
    h_eigvals: &Array1::<f64>, h_eigvecs: &CxMat, config: &Config)
  -> Result<CxMat,CluEError>
{

  let Some(temperature) = config.temperature else {
    return Err(CluEError::NoTemperature);
  };
  let beta = I*HBAR/(temperature*BOLTZMANN);

  let mut rho_list = get_propagators_complex_time_from_eig(
      h_eigvals, h_eigvecs, &[-beta])?;

  let mut density_matrix = rho_list.swap_remove(0);

  let Ok(z) = density_matrix.trace() else{
    return Err(CluEError::CannotTakeTrace(format!("{}",density_matrix)));
  };
  if z.norm() < 1e-12{
    return Err(CluEError::CannotTakeTrace(format!("{}",density_matrix)));
  }
  density_matrix /= z;
  Ok(density_matrix)


}  
//------------------------------------------------------------------------------


//==============================================================================
#[cfg(test)]
mod tests{
  use super::*;
  use crate::io::FromTOMLString;
  use crate::physical_constants::{SQRT2_INV};
  use crate::space_3d::{Vector3D,SymmetricTensor3D};
  use crate::quantum::tensors::*;
  use crate::quantum::spin_hamiltonian::{
    build_block_diag_hamiltonian,
  };
  use crate::cluster_methods::appa::appa_hahn;


  //----------------------------------------------------------------------------
  #[test]
  fn test_propagate_custom_pulse_sequence(){
    let mut config= Config::from_toml_string(r##"
       cluster_method = "cce"
       magnetic_field = 1.2
       pulse_sequence = "hahn"
       number_timepoints = [101]
       tau_increments = [5e-2]
       populations = "thermal"
       temperature = 20
       [detected_spin]
         transition = [0,1]
     "##).unwrap();
    config.set_tau_axis().unwrap();

    let z0 = 33.0e9;
    let z1 = 80.0e6;
    let a1 = 10.0e6;
    let a2 = -10.0e6;
    let b = 10.0e3;

    let tensors = build_restricted_three_spin_tensors(z0, z1, a1, a2, b);

    let spin_indices = vec![0,1,2];
    let spin_ops = ClusterSpinOperators::new(1,&vec![2],3,&config).unwrap();

    let (h_eigvals, h_eigvecs) = build_spin_hamiltonian(
        &spin_indices,&spin_ops, &tensors,
        ).unwrap();
    
    let dim = h_eigvecs.dim().0;

    let density_matrix = get_electron_cluster_thermal_density_matrix(
        &h_eigvals, &h_eigvecs, &config).unwrap();

    let u_pi = (-2.0*I)*spin_ops.get(&SpinOp::Sy,2,3,0).unwrap().clone();

    let u_half_pi = Complex::<f64>{re:SQRT2_INV,im:0.0}*(CxMat::eye(8) 
        + u_pi.clone());
    
    let detection_op0 = spin_ops.get(&SpinOp::Sp,2,3,0).unwrap().clone();

    let pulse_sequence = vec![
      PulseStep::Pulse(&u_half_pi),
      PulseStep::TauDelay,
      PulseStep::Pulse(&u_pi),
      PulseStep::TauDelay,
      PulseStep::Detect(&detection_op0),
    ];

    let signal = propagate_custom_pulse_sequence(
        &pulse_sequence, &density_matrix, &h_eigvals,&h_eigvecs, &config)
        .unwrap();

    let spin_indices = vec![1,2];

    let ref_signal = appa_hahn(
        &spin_indices,&tensors,&config).unwrap().unwrap();

    for (ii,v) in signal.data.iter().enumerate(){
      let v0 = ref_signal.data[ii];
      let err: f64;
      if (v+v0).norm() < 1e12{
        err = (v-v0).norm();
      }else{
        err = (2.0+(v-v0)/(v+v0)).norm();
      }
      assert!(err < 1e-9);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_propagate_pulse_sequence(){
    let mut config= Config::from_toml_string(r##"
       cluster_method = "cce"
       magnetic_field = 1.2
       pulse_sequence = "hahn"
       number_timepoints = [101]
       tau_increments = [5e-2]
       populations = "thermal"
       temperature = 20
       [detected_spin]
         transition = [0,1]
     "##).unwrap();
    config.set_tau_axis().unwrap();

    let z0 = 33.0e9;
    let z1 = 80.0e6;
    let a1 = 10.0e6;
    let a2 = -10.0e6;
    let b = 10.0e3;

    let tensors = build_restricted_three_spin_tensors(z0, z1, a1, a2, b);
    
    let spin_indices = vec![0,1,2];
    let spin_ops = ClusterSpinOperators::new(1,&vec![2],3, &config).unwrap();

    let (h_eigvals, h_eigvecs) = build_spin_hamiltonian(
        &spin_indices,&spin_ops, &tensors).unwrap();

    let dim = h_eigvecs.dim().0;

    let density_matrix = get_electron_cluster_thermal_density_matrix(
        &h_eigvals, &h_eigvecs, &config).unwrap();

    let u_pi = (-2.0*I)*spin_ops.get(&SpinOp::Sy,2,3,0).unwrap().clone();

    let u_half_pi = Complex::<f64>{re:SQRT2_INV,im:0.0}*(CxMat::eye(8) 
        + u_pi.clone());
    let detection_op0 = spin_ops.get(&SpinOp::Sp,2,3,0).unwrap().clone();

    let pulses = vec![
      PulseStep::Pulse(&u_half_pi),
      PulseStep::Pulse(&u_pi),
      PulseStep::Detect(&detection_op0),
    ];

    let signal = propagate_pulse_sequence(
        &pulses, &density_matrix, &h_eigvals, &h_eigvecs, &config)
        .unwrap();

    let spin_indices = vec![1,2];

    let ref_signal = appa_hahn(
        &spin_indices,&tensors,&config).unwrap().unwrap();

    for (ii,v) in signal.data.iter().enumerate(){
      let v0 = ref_signal.data[ii];
      let err: f64;
      if (v+v0).norm() < 1e12{
        err = (v-v0).norm();
      }else{
        err = (2.0+(v-v0)/(v+v0)).norm();
      }
      assert!(err < 1e-9);
    }

  }
  
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_spin_hamiltonian(){

    let z0 = 1000.0;
    let z1 = 100.0;
    let a1 = 2.0;
    let a2 = -1.0;
    let b = 0.1;

    let tensors = build_restricted_three_spin_tensors(z0, z1, a1, a2, b);
  
    let mut config = Config::new();
    config.set_defaults().unwrap();

    let spin_indices = vec![0,1,2];
    let spin_ops = ClusterSpinOperators::new(1,&vec![2],3,&config).unwrap();

    let (h_eigvals, h_eigvecs) = build_spin_hamiltonian(
        &spin_indices,&spin_ops, &tensors).unwrap();

    let h_inv_eigvecs = h_eigvecs.t().map(|v| v.conj());
    let eigval_matrix = CxMat::from_diag(&
        h_eigvals.map(|v| Complex::<f64>{re: *v, im: 0.0})); 
    let hamiltonian = h_eigvecs.dot(&eigval_matrix.dot(&h_inv_eigvecs));
      
    let mut config = Config::new();
    config.set_defaults().unwrap();

    let spin_indices = vec![1,2];

    let block_hamiltonian = build_block_diag_hamiltonian(&spin_indices,&spin_ops, &tensors,
        &config).unwrap();


    assert_eq!(hamiltonian.len(), 4*block_hamiltonian.beta_eigvecs.len());

    let beta = hamiltonian.slice(ndarray::s![4..,4..]).to_owned();

    let alpha = hamiltonian.slice(ndarray::s![0..4,0..4]).to_owned();

    let halpha = block_hamiltonian.alpha();
    assert!(approx_eq(&halpha, &alpha, 1e-12));

    let hbeta = block_hamiltonian.beta();
    assert!(approx_eq(&hbeta, &beta, 1e-12));


  }
  //----------------------------------------------------------------------------
  fn build_restricted_three_spin_tensors(z0: f64, z1: f64, a1: f64, a2: f64,
      b: f64) -> HamiltonianTensors{
    let spin_multiplicities = vec![2,2,2];

    let mut spin1_tensors = Spin1Tensors::new(3);
    let zeeman0 = Vector3D::from([0.0, 0.0, z0]);
    let zeeman1 = Vector3D::from([0.0, 0.0, z1]);
    spin1_tensors.set(0,zeeman0);
    spin1_tensors.set(1,zeeman1.clone());
    spin1_tensors.set(2,zeeman1);

    let mut spin2_tensors = Spin2Tensors::new(3);
    let hf1 = SymmetricTensor3D::from([ 0.0, 0.0, 0.0,
                                             0.0, 0.0,
                                                    a1]);
    let hf2 = SymmetricTensor3D::from([ 0.0, 0.0, 0.0,
                                             0.0, 0.0,
                                                   a2]);

    let dip = SymmetricTensor3D::from([ -b/2.0,    0.0, 0.0,
                                                -b/2.0, 0.0,
                                                         b]);

    spin2_tensors.set(0,1,hf1);
    spin2_tensors.set(0,2,hf2);
    spin2_tensors.set(1,2,dip);

    let ge = -1.7609e+11;
    HamiltonianTensors{
      spin_multiplicities,
      spin1_tensors,
      spin2_tensors,
      detected_gamma_matrix: SymmetricTensor3D::from([ge, 0.0, 0.0,
                                                           ge, 0.0,
                                                               ge]),
      magnetic_field: Vector3D::from([0.0,0.0,1.2]),
      mean_field_couplings: None,
      }
  }
  //----------------------------------------------------------------------------
  fn approx_eq(m0: &CxMat, m1: &CxMat, tol: f64) -> bool {

    assert!(tol > 0.0);

    if m0.ncols() != m1.ncols() || m0.nrows() != m1.nrows() {
      return false;
    }

    for irow in 0..m0.nrows() {
      for icol in 0..m0.ncols(){
        let err = ( m0[[irow,icol]] - m1[[irow,icol]] ).norm();
        if err >= tol {
          return false;
        }
      }
    }
    true
  }
  //----------------------------------------------------------------------------
}
