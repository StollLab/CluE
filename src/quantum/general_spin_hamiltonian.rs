use crate::CluEError;
use crate::config::{
  Config,
};
use crate::config::pulse_sequence::PulseSequence;
use crate::math::cxmat_pow_n;
use crate::physical_constants::{HALF,ZERO};
use crate::quantum::block_diagonal_spin_hamiltonian::{
  get_propagators_from_eig,
  thermal_density_matrix,
  thermal_density_matrix_from_eig,
};
use crate::quantum::pulse_sequences::{
  PulseStep,  
};
use crate::quantum::cluster_operators::{
  ClusterSpinOperators,
  SpinOp,
  spin_x,
  spin_y,
  spin_z,
  spin_ist,
};
use crate::signal::Signal;
use crate::quantum::tensors::HamiltonianTensors;

use num_complex::Complex;
use ndarray::{Array1,Array2};
use ndarray::linalg::kron;
use ndarray_linalg::UPLO;
use crate::ndarray_linalg_eigh_workaround::eigh_hermitian;

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
                // Tr[ρ D] = Σ_ij ρ_ji D_ij.
                let it = std::iter::zip(density_matrix.t(),&detection_op);
                v += it.map(|(rho_ji,u_ij)| rho_ji*u_ij).sum::<Complex<f64>>();
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
                  let it = std::iter::zip(density_matrix.t(),&detection_op);
                  v += it.map(|(rho_ji,u_ij)| rho_ji*u_ij).sum::<Complex<f64>>();
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
    return Err(CluEError::UnexpectedPulseStep(0,"a pi/2 pulse".to_string()))
  };
  let PulseStep::Pulse(u_pi) =  pulses[1] else{
    return Err(CluEError::UnexpectedPulseStep(1,"a pi pulse".to_string()))
  };
  let PulseStep::Detect(detection_op0) =  pulses[2] else{
    return Err(CluEError::UnexpectedPulseStep(2,
          "a detection operator".to_string()))
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
        let mut u_seg = u_of_tau.dot(&u_pi.dot(u_of_tau));
        u_seg = cxmat_pow_n(&u_seg,*n_pi);
        u_sequence = u_seg.dot(&u_sequence);
      }

      PulseSequence::RefocusedHahnEcho =>{
        let u_seg1 = u_of_tau.dot(&u_pi.dot(u_of_tau));
        for u_of_tau2 in u_of_tau2s.iter(){
          let mut u_sequence = u_half_pi.clone();
          let u_seg2 = u_of_tau2.dot(&u_pi.dot(u_of_tau2));
          u_sequence = u_seg2.dot(&u_seg1.dot(&u_sequence));
          
          let u_sequence_dag = u_sequence.t().map(|u_ij| u_ij.conj() );
    
          // Find time evolved detection operator:
          // D(T) = U(T)^* D U(T). 
          let detection_op = u_sequence_dag.dot(&detection_op0.dot(&u_sequence));

          let it = std::iter::zip(density_matrix.t(),&detection_op);
          let v = it.map(|(rho_ji,u_ij)| rho_ji*u_ij).sum::<Complex<f64>>();
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

    // Tr[ρ D(T)] = Σ_ij ρ_ji D(T)_ij.
    let it = std::iter::zip(density_matrix.t(),&detection_op);
    let v = it.map(|(rho_ji,u_ij)| rho_ji*u_ij).sum::<Complex<f64>>();
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
    for _ in 1..n_timepoints{
      let idx_u = u_of_taus.len() - 1;
      u_of_taus.push( dus[idt].dot(&u_of_taus[idx_u]));
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
    return Err(CluEError::NoTimeIncrements2);
  }

  let number_timepoints = &config.number_timepoints2;
  if number_timepoints.is_empty(){
    return Err(CluEError::NoTimepoints2);
  }
  let n_tot = config.get_total_number_timesteps();

  let dus = get_propagators_from_eig(eigvals,eigvecs,tau_increments)?;

  let mut u_of_taus = Vec::<CxMat>::with_capacity(n_tot);
  u_of_taus.push( CxMat::eye(dus[0].dim().0) );

  for (idt,_dt) in tau_increments.iter().enumerate(){
    let n_timepoints = number_timepoints[idt];
    for _ in 1..n_timepoints{
      let idx_u = u_of_taus.len() - 1;
      u_of_taus.push( dus[idt].dot(&u_of_taus[idx_u]));
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
  let ham = build_spin_hamiltonian_matrix(spin_indices, spin_ops, tensors,
      true, true)?;

  let Ok((eigvals, eigvecs)) = eigh_hermitian(&ham, UPLO::Lower) else{
    return Err(
        CluEError::CannotDiagonalizeHamiltonian(
          crate::math::cx_matrix_summary(&ham)));
  };

  Ok((eigvals, eigvecs))
}
//------------------------------------------------------------------------------
/// This function builds 1_S ⊗ H_C⁰ (Hz) for a gCCE cluster `spin_indices`
/// (position 0 is the detected spin): the gCCE cluster Hamiltonian without
/// the detected-spin terms, without the hyperfine couplings, and without
/// mean fields. H_C⁰ holds the one-spin terms of the bath spins (Zeeman,
/// quadrupole/zero-field, spherical terms up to `max_spherical_tensor_rank`)
/// and the couplings within the cluster. It sets the bath thermal state of
/// ensemble gCCE with a given detected-spin state
/// (`get_gcce_product_density_matrix`).
pub fn build_bath_spin_hamiltonian(spin_indices: &[usize],
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors)
  -> Result<CxMat,CluEError>
{
  build_spin_hamiltonian_matrix(spin_indices, spin_ops, tensors, false, false)
}
//------------------------------------------------------------------------------
// This function builds the gCCE cluster Hamiltonian matrix (Hz), made 
// Hermitian. With `include_detected_spin = false`, the one-spin terms and
// couplings of position 0 (the detected spin) are left out; with 
// `include_mean_field = false`, the mean fields are left out.
fn build_spin_hamiltonian_matrix(spin_indices: &[usize],
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors,
    include_detected_spin: bool, include_mean_field: bool)
  -> Result<CxMat,CluEError>
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

    // Position 0 is the detected spin; skipping it leaves out its one-spin
    // terms and, since pairs start at sop_idx1 = sop_idx0, its couplings.
    if sop_idx0 == 0 && !include_detected_spin{
      continue;
    }


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

    if include_mean_field{
      if let Some(vec) = tensors.get_mean_field_couplings(ten_idx0,
          spin_indices)
      {
        ham = ham + sx0*vec.x();
        ham = ham + sy0*vec.y();
        ham = ham + sz0*vec.z();
      }
    }

    if let Some(tens) = tensors.spin_spherical_tensors
        .get(&(ten_idx0,ten_idx0)){

      for coefs in tens.iter(){
        let il = coefs.dim().0;
        let l = (il as i32 -1)/2;
        if l > tensors.max_spherical_tensor_rank{
          continue;
        }
        for (im,&c) in coefs.iter().enumerate(){
          let m = im as i32 -l;
          let tlm = spin_ops.get(&SpinOp::T(l,m),spin_mult0,cluster_size,
              sop_idx0)?;
          ham = ham + tlm*c;
        }
      }
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

  // Remove any floating point errors that break Hermiticity.
  let ham_dag = ham.t().map(|u_ij| u_ij.conj() );
  ham = HALF*(ham + ham_dag);

  Ok(ham)
}
//------------------------------------------------------------------------------
/// This function builds the one-spin Hamiltonian (Hz) of spin `index`
/// in the Zeeman basis (index i is m = S - i): the Zeeman term
/// (`spin1_tensors`), the Cartesian quadrupole or zero-field term
/// (`spin2_tensors(index,index)`), and, if `include_spherical`, the
/// spherical-tensor terms up to `max_spherical_tensor_rank`.
/// It has no mean field and no coupling to other spins.
/// Bath spin states are sampled from the eigenstates of this Hamiltonian
/// (`SpinStates::generate`); gCCE includes the spherical terms because
/// `build_spin_hamiltonian` propagates them, while CCE does not, because
/// `build_block_diag_hamiltonian` does not.
pub fn build_single_spin_hamiltonian(tensors: &HamiltonianTensors,
    index: usize, include_spherical: bool)
  -> Result<CxMat,CluEError>
{
  let Some(&spin_multiplicity) = tensors.spin_multiplicities.get(index) else{
    return Err(CluEError::NoSpinWithIndex(index,
          tensors.spin_multiplicities.len()));
  };

  // Initialize Hamiltonian.
  let mut ham = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  let sx = spin_x(spin_multiplicity);
  let sy = spin_y(spin_multiplicity);
  let sz = spin_z(spin_multiplicity);

  // Zeeman
  if let Some(vec) = tensors.spin1_tensors.get(index){
    ham = ham + &sx*vec.x();
    ham = ham + &sy*vec.y();
    ham = ham + &sz*vec.z();
  }

  if include_spherical{
    if let Some(tens) = tensors.spin_spherical_tensors
        .get(&(index,index)){

      for coefs in tens.iter(){
        let il = coefs.dim().0;
        let l = (il as i32 -1)/2;
        if l > tensors.max_spherical_tensor_rank{
          continue;
        }
        for (im,&c) in coefs.iter().enumerate(){
          let m = im as i32 -l;
          let tlm = spin_ist(spin_multiplicity,l,m);
          ham = ham + tlm*c;
        }
      }
    }
  }

  // quadrupole or zero-field
  if let Some(ten) = tensors.spin2_tensors.get(index,index){
    ham = ham + sx.dot(&sx)*ten.xx();
    ham = ham + sx.dot(&sy)*ten.xy();
    ham = ham + sx.dot(&sz)*ten.xz();

    ham = ham + sy.dot(&sx)*ten.yx();
    ham = ham + sy.dot(&sy)*ten.yy();
    ham = ham + sy.dot(&sz)*ten.yz();

    ham = ham + sz.dot(&sx)*ten.zx();
    ham = ham + sz.dot(&sy)*ten.zy();
    ham = ham + sz.dot(&sz)*ten.zz();
  }

  // Remove any floating point errors that break Hermiticity.
  let ham_dag = ham.t().map(|u_ij| u_ij.conj() );
  ham = HALF*(ham + ham_dag);

  Ok(ham)
}
//------------------------------------------------------------------------------
/// This function builds and diagonalizes the Hamiltonian of the detected
/// spin (index 0): `build_single_spin_hamiltonian` with the spherical terms,
/// plus the mean field of the bath, if mean-field couplings are set.
/// The eigenvalues are in ascending order.
pub fn build_detected_spin_hamiltonian(tensors: &HamiltonianTensors)
  -> Result<(Array1::<f64>,CxMat),CluEError>
{

  let ten_idx0 = 0;

  let spin_multiplicity: usize = tensors.spin_multiplicities[ten_idx0];

  let mut ham = build_single_spin_hamiltonian(tensors, ten_idx0, true)?;

  if let Some(vec) = tensors.get_mean_field_couplings(ten_idx0,&[ten_idx0]){
    ham = ham + &spin_x(spin_multiplicity)*vec.x();
    ham = ham + &spin_y(spin_multiplicity)*vec.y();
    ham = ham + &spin_z(spin_multiplicity)*vec.z();
  }

  let Ok((eigvals, eigvecs)) = eigh_hermitian(&ham, UPLO::Lower) else{
    return Err(
        CluEError::CannotDiagonalizeHamiltonian(
          crate::math::cx_matrix_summary(&ham)));
  };

  Ok((eigvals, eigvecs))
}
//------------------------------------------------------------------------------
/// This function returns the thermal density matrix exp(-βH)/Z of the
/// Hamiltonian with eigenvalues `h_eigvals` (Hz) and eigenvectors
/// `h_eigvecs`, with β = h/(k_B T) and T = `conditions.temperature`.
/// The energies are shifted by their minimum before exponentiating.
pub fn get_electron_cluster_thermal_density_matrix(
    h_eigvals: &Array1::<f64>, h_eigvecs: &CxMat, config: &Config)
  -> Result<CxMat,CluEError>
{

  let Some(temperature) = config.temperature else {
    return Err(CluEError::NoTemperature);
  };
  thermal_density_matrix_from_eig(h_eigvals, h_eigvecs, temperature)
}
//------------------------------------------------------------------------------
/// This function returns the initial state of an ensemble-gCCE cluster when
/// the detected-spin state ρ_S (`rho_s`, d_S × d_S) is given:
/// ρ_S ⊗ exp(-βH_C⁰)/Tr with `temperature`, or ρ_S ⊗ 1/d_C without.
/// `bath_ham` is 1_S ⊗ H_C⁰ (`build_bath_spin_hamiltonian`). The thermal 
/// state of 1_S ⊗ H_C⁰ is G = (1_S/d_S) ⊗ g_C with g_C = exp(-βH_C⁰)/Tr,
/// and ρ_S ⊗ 1_C commutes with it, so d_S (ρ_S ⊗ 1_C) G = ρ_S ⊗ g_C.
pub fn get_gcce_product_density_matrix(rho_s: &CxMat, bath_ham: &CxMat,
    temperature: Option<f64>)
  -> Result<CxMat,CluEError>
{
  let d_s = rho_s.dim().0;
  let dim = bath_ham.dim().0;
  if d_s == 0 || dim % d_s != 0 || rho_s.dim().1 != d_s 
    || bath_ham.dim().1 != dim
  {
    return Err(CluEError::DetectedSpinStateDoesNotFitCluster(dim, d_s));
  }
  let d_c = dim/d_s;

  match temperature{
    None => Ok(kron(rho_s, &(CxMat::eye(d_c)/(d_c as f64)))),
    Some(temperature) => {
      let g = thermal_density_matrix(bath_ham, temperature)?;
      let rho_s_bath = kron(rho_s, &CxMat::eye(d_c));
      Ok(rho_s_bath.dot(&g)*Complex::<f64>{re: d_s as f64, im: 0.0})
    },
  }
}
//------------------------------------------------------------------------------


//==============================================================================
#[cfg(test)]
mod tests{
  use super::*;
  use crate::io::FromTOMLString;
  use crate::physical_constants::{I,SQRT2_INV};
  use crate::space_3d::{Vector3D,SymmetricTensor3D};
  use crate::quantum::tensors::*;
  use crate::quantum::block_diagonal_spin_hamiltonian::{
    build_block_diag_hamiltonian,
  };
  use crate::cluster_methods::appa::appa_hahn;


  //----------------------------------------------------------------------------
  #[test]
  fn test_propagate_custom_pulse_sequence(){
    let mut config= Config::from_toml_string(r##"
       cluster_expansion.method = "cce"
       conditions.magnetic_field = 1.2
       pulse_sequence.sequence = "hahn"
       pulse_sequence.number_timepoints = [101]
       pulse_sequence.tau_increments = [5e-2]
       conditions.temperature = 20
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
    let det_h = (Array1::<f64>::zeros(1) ,CxMat::eye(1));
    let spin_ops = ClusterSpinOperators::new(det_h,&[2],3,&config).unwrap();

    let (h_eigvals, h_eigvecs) = build_spin_hamiltonian(
        &spin_indices,&spin_ops, &tensors,
        ).unwrap();
    
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
      let err: f64 = if (v+v0).norm() < 1e12{
        (v-v0).norm()
      }else{
        (2.0+(v-v0)/(v+v0)).norm()
      };
      assert!(err < 1e-9);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_propagate_pulse_sequence(){
    let mut config= Config::from_toml_string(r##"
       cluster_expansion.method = "cce"
       conditions.magnetic_field = 1.2
       pulse_sequence.sequence = "hahn"
       pulse_sequence.number_timepoints = [101]
       pulse_sequence.tau_increments = [5e-2]
       conditions.temperature = 20
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
    let det_h = (Array1::<f64>::zeros(1) ,CxMat::eye(1));
    let spin_ops = ClusterSpinOperators::new(det_h,&[2],3, &config).unwrap();

    let (h_eigvals, h_eigvecs) = build_spin_hamiltonian(
        &spin_indices,&spin_ops, &tensors).unwrap();

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
      let err: f64 = if (v+v0).norm() < 1e12{
        (v-v0).norm()
      }else{
        (2.0+(v-v0)/(v+v0)).norm()
      };
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
    let det_h = (Array1::<f64>::zeros(1) ,CxMat::eye(1));
    let spin_ops = ClusterSpinOperators::new(det_h,&[2],3,&config).unwrap();

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
      // No bath gamma matrices: every bath Zeeman vector is fixed.
      bath_gamma_matrices: vec![None; spin1_tensors.len()],
      spin1_tensors,
      spin2_tensors,
      detected_gamma_matrix: SymmetricTensor3D::from([ge, 0.0, 0.0,
                                                           ge, 0.0,
                                                               ge]),
      magnetic_field: Vector3D::from([0.0,0.0,1.2]),
      mean_field_couplings: None,
      spin_spherical_tensors: SpinSphericalTensors::new(),
      max_spherical_tensor_rank: 0,
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
  //----------------------------------------------------------------------------
  #[test]
  fn test_2nd_free_evolution_errors(){
    let eigvals = Array1::<f64>::zeros(2);
    let eigvecs = CxMat::eye(2);
    let mut config = Config::new();
    let result = get_2nd_free_evolutions_propagators(&eigvals, &eigvecs, 
        &config);
    assert_eq!(result.err(), Some(CluEError::NoTimeIncrements2));

    config.tau2_increments = vec![1e-8];
    let result = get_2nd_free_evolutions_propagators(&eigvals, &eigvecs, 
        &config);
    assert_eq!(result.err(), Some(CluEError::NoTimepoints2));
  }
  //----------------------------------------------------------------------------
  // The signal must be Tr[ρ U^† D U], not Tr[ρ^T U^† D U]; the two differ
  // when ρ is not real-symmetric.
  // One spin-1/2 with H = ν Sz, ρ = Sx + Sy, D = Sx, and identity pulses:
  // U^† Sx U = Sx cos(φ) - Sy sin(φ), with φ = 2πνt, so
  // Tr[ρ D(t)]/Tr[ρ D(0)] = cos(φ) - sin(φ).
  // Tr[ρ^T D(t)] would give cos(φ) + sin(φ).
  #[test]
  fn test_signal_is_trace_of_rho_times_detection_operator(){
    let mut config = Config::from_toml_string(r##"
      cluster_expansion.method = "gcce"
      conditions.magnetic_field = 1.2
      pulse_sequence.sequence = "fid"
      pulse_sequence.number_timepoints = [8]
      pulse_sequence.tau_increments = [0.125]
    "##).unwrap();
    config.set_tau_axis().unwrap();

    let nu = 1.0e6;
    // H = ν Sz = diag(ν/2, -ν/2); eigh returns ascending eigenvalues.
    let eigvals = Array1::<f64>::from_vec(vec![-0.5*nu, 0.5*nu]);
    let one = Complex::<f64>::new(1.0,0.0);
    let zero = Complex::<f64>::new(0.0,0.0);
    let eigvecs = ndarray::array![[zero, one],[one, zero]];

    let e = CxMat::eye(2);
    let sx = spin_x(2);
    let rho = &sx + &spin_y(2);

    let expected: Vec::<f64> = config.get_tau_axis().unwrap().iter()
      .map(|t| {
        let phi = 2.0*crate::physical_constants::PI*nu*t;
        phi.cos() - phi.sin()
      }).collect();

    // Standard pulse sequence (FID).
    let pulses = vec![PulseStep::Pulse(&e),PulseStep::Pulse(&e),
        PulseStep::Detect(&sx)];
    let signal = propagate_pulse_sequence(&pulses, &rho, &eigvals, &eigvecs,
        &config).unwrap();
    assert_eq!(signal.len(), expected.len());
    for (v,v0) in signal.data.iter().zip(expected.iter()){
      assert!((v - v0).norm() < 1e-12);
    }

    // Custom pulse sequence.
    let pulse_sequence = vec![PulseStep::TauDelay, PulseStep::Detect(&sx)];
    let signal = propagate_custom_pulse_sequence(&pulse_sequence, &rho,
        &eigvals, &eigvecs, &config).unwrap();
    assert_eq!(signal.len(), expected.len());
    for (v,v0) in signal.data.iter().zip(expected.iter()){
      assert!((v - v0).norm() < 1e-12);
    }
  }
  //----------------------------------------------------------------------------
  // This function builds tensors for one spin-1 (index 0) and one spin-1/2
  // (index 1). Spin 0 has a tilted Zeeman vector, a tilted axial Q
  // (Q = d(n n^T - 1/3), n at polar angle theta in the xz-plane), and a
  // rank-4 spherical term; spin 1 has only a Zeeman vector.
  fn build_single_spin_test_tensors(nu: Vector3D, d: f64, theta: f64,
      c4: f64, max_rank: i32) -> HamiltonianTensors
  {
    let mut spin1_tensors = Spin1Tensors::new(2);
    spin1_tensors.set(0, nu);
    spin1_tensors.set(1, Vector3D::from([0.0, 0.0, 5.0e6]));

    let (st, ct) = (theta.sin(), theta.cos());
    let mut spin2_tensors = Spin2Tensors::new(2);
    spin2_tensors.set(0, 0, SymmetricTensor3D::from([
        d*(st*st - 1.0/3.0), 0.0, d*st*ct,
                        -d/3.0, 0.0,
                                d*(ct*ct - 1.0/3.0)]));

    let mut spin_spherical_tensors = SpinSphericalTensors::new();
    let mut coefs = CxMat::zeros((9,1));
    // Index l + m: m = 4, m = -4, and m = 0.
    coefs[[8,0]] = Complex::<f64>{re: c4, im: 0.0};
    coefs[[0,0]] = Complex::<f64>{re: c4, im: 0.0};
    coefs[[4,0]] = Complex::<f64>{re: 0.5*c4, im: 0.0};
    spin_spherical_tensors.insert((0,0), vec![coefs]);

    HamiltonianTensors{
      spin_multiplicities: vec![3,2],
      // No bath gamma matrices: every bath Zeeman vector is fixed.
      bath_gamma_matrices: vec![None; spin1_tensors.len()],
      spin1_tensors,
      spin2_tensors,
      spin_spherical_tensors,
      max_spherical_tensor_rank: max_rank,
      detected_gamma_matrix: SymmetricTensor3D::eye(),
      magnetic_field: Vector3D::from([0.0,0.0,1.2]),
      mean_field_couplings: None,
    }
  }
  //----------------------------------------------------------------------------
  // The detected-spin Hamiltonian as it was built before
  // build_single_spin_hamiltonian existed (commit cecaa49), kept here as the
  // reference for the regression test.
  fn reference_detected_spin_hamiltonian(tensors: &HamiltonianTensors)
    -> CxMat
  {
    let m = tensors.spin_multiplicities[0];
    let mut ham = CxMat::zeros((m,m));
    let sx = spin_x(m);
    let sy = spin_y(m);
    let sz = spin_z(m);
    if let Some(vec) = tensors.spin1_tensors.get(0){
      ham = ham + &sx*vec.x();
      ham = ham + &sy*vec.y();
      ham = ham + &sz*vec.z();
    }
    if let Some(vec) = tensors.get_mean_field_couplings(0,&[0]){
      ham = ham + &sx*vec.x();
      ham = ham + &sy*vec.y();
      ham = ham + &sz*vec.z();
    }
    if let Some(tens) = tensors.spin_spherical_tensors.get(&(0,0)){
      for coefs in tens.iter(){
        let il = coefs.dim().0;
        let l = (il as i32 -1)/2;
        if l > tensors.max_spherical_tensor_rank{ continue; }
        for (im,&c) in coefs.iter().enumerate(){
          ham = ham + spin_ist(m,l,im as i32 - l)*c;
        }
      }
    }
    if let Some(ten) = tensors.spin2_tensors.get(0,0){
      ham = ham + sx.dot(&sx)*ten.xx();
      ham = ham + sx.dot(&sy)*ten.xy();
      ham = ham + sx.dot(&sz)*ten.xz();
      ham = ham + sy.dot(&sx)*ten.yx();
      ham = ham + sy.dot(&sy)*ten.yy();
      ham = ham + sy.dot(&sz)*ten.yz();
      ham = ham + sz.dot(&sx)*ten.zx();
      ham = ham + sz.dot(&sy)*ten.zy();
      ham = ham + sz.dot(&sz)*ten.zz();
    }
    let ham_dag = ham.t().map(|u_ij| u_ij.conj() );
    HALF*(ham + ham_dag)
  }
  //----------------------------------------------------------------------------
  fn sorted_eigenvalues(ham: &CxMat) -> Vec::<f64>{
    let (e, _v) = eigh_hermitian(ham, UPLO::Lower).unwrap();
    e.to_vec()
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_detected_spin_hamiltonian_unchanged(){
    use ndarray_linalg::Norm;
    let nu = Vector3D::from([1.0e8, -2.0e7, 9.5e9]);
    for max_rank in [0, 4]{
      let tensors = build_single_spin_test_tensors(nu.clone(), 2.0e9, 0.7,
          3.0e7, max_rank);
      let reference = reference_detected_spin_hamiltonian(&tensors);

      // Without mean fields, the detected-spin Hamiltonian is the one-spin
      // Hamiltonian with the spherical terms.
      let ham = build_single_spin_hamiltonian(&tensors, 0, true).unwrap();
      assert!((&ham - &reference).norm() <= 1e-12*reference.norm());

      // The eigendecomposition equals that of the reference matrix, done
      // the same way (eigh, UPLO::Lower).
      let (e, v) = build_detected_spin_hamiltonian(&tensors).unwrap();
      let (e_ref, v_ref) = eigh_hermitian(&reference, UPLO::Lower).unwrap();
      let n = e.len();
      let rebuild = |e: &Array1::<f64>, v: &CxMat| -> CxMat {
        let mut d = CxMat::zeros((n,n));
        for k in 0..n{ d[[k,k]] = Complex::<f64>{re: e[k], im: 0.0}; }
        v.dot(&d).dot(&v.t().map(|x| x.conj()))
      };
      for k in 0..n{
        assert!((e[k] - e_ref[k]).abs() <= 1e-12*reference.norm());
      }
      assert!((rebuild(&e,&v) - rebuild(&e_ref,&v_ref)).norm() 
          <= 1e-12*reference.norm());
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_single_spin_hamiltonian(){
    use ndarray_linalg::Norm;
    let d = 2.0e9;
    let nu = 9.5e9;

    // Zeeman along z and Q along z: E(m) = nu*m + d*(m^2 - 2/3).
    let tensors = build_single_spin_test_tensors(
        Vector3D::from([0.0,0.0,nu]), d, 0.0, 0.0, 0);
    let ham = build_single_spin_hamiltonian(&tensors, 0, false).unwrap();
    let mut expected = [-2.0*d/3.0, -nu + d/3.0, nu + d/3.0];
    expected.sort_by(|a,b| a.partial_cmp(b).unwrap());
    for (a,b) in sorted_eigenvalues(&ham).iter().zip(expected.iter()){
      assert!((a - b).abs() <= 1e-12*nu);
    }
    // Index i is m = S - i: the diagonal is (nu + d/3, -2d/3, -nu + d/3).
    assert!((ham[[0,0]].re - (nu + d/3.0)).abs() <= 1e-12*nu);
    assert!((ham[[1,1]].re + 2.0*d/3.0).abs() <= 1e-12*nu);
    assert!((ham[[2,2]].re - (-nu + d/3.0)).abs() <= 1e-12*nu);

    // No field, Q tilted by 0.7 rad: the levels do not depend on the tilt.
    let tensors = build_single_spin_test_tensors(
        Vector3D::from([0.0,0.0,0.0]), d, 0.7, 0.0, 0);
    let ham = build_single_spin_hamiltonian(&tensors, 0, false).unwrap();
    let expected = [-2.0*d/3.0, d/3.0, d/3.0];
    for (a,b) in sorted_eigenvalues(&ham).iter().zip(expected.iter()){
      assert!((a - b).abs() <= 1e-12*d);
    }

    // The spherical terms: with rank 4 the matrices with and without them
    // differ by exactly sum_m c_4m T_4m; with rank 0 they are identical.
    let c4 = 3.0e7;
    let tensors = build_single_spin_test_tensors(
        Vector3D::from([1.0e8,-2.0e7,nu]), d, 0.7, c4, 4);
    let with = build_single_spin_hamiltonian(&tensors, 0, true).unwrap();
    let without = build_single_spin_hamiltonian(&tensors, 0, false).unwrap();
    let c = Complex::<f64>{re: c4, im: 0.0};
    let term = spin_ist(3,4,4)*c + spin_ist(3,4,-4)*c + spin_ist(3,4,0)*(0.5*c);
    let term = HALF*(&term + &term.t().map(|x| x.conj()));
    assert!((&with - &without - &term).norm() <= 1e-12*with.norm());

    let tensors = build_single_spin_test_tensors(
        Vector3D::from([1.0e8,-2.0e7,nu]), d, 0.7, c4, 0);
    let with = build_single_spin_hamiltonian(&tensors, 0, true).unwrap();
    let without = build_single_spin_hamiltonian(&tensors, 0, false).unwrap();
    assert_eq!(with, without);

    // The spin-1/2 at index 1 has only its Zeeman term.
    let ham = build_single_spin_hamiltonian(&tensors, 1, true).unwrap();
    assert!((&ham - &(spin_z(2)*Complex::<f64>{re: 5.0e6, im: 0.0})).norm() 
        <= 1e-12*5.0e6);

    assert_eq!(build_single_spin_hamiltonian(&tensors, 2, true),
        Err(CluEError::NoSpinWithIndex(2,2)));
  }
  //----------------------------------------------------------------------------
  // Partial traces of a (d_s·d_c)×(d_s·d_c) matrix over the bath (Tr_C) and
  // over the detected spin (Tr_S); the detected spin is the first factor.
  fn trace_bath(rho: &CxMat, d_s: usize) -> CxMat{
    let d_c = rho.dim().0/d_s;
    let mut out = CxMat::zeros((d_s,d_s));
    for i in 0..d_s{ for j in 0..d_s{ for k in 0..d_c{
      out[[i,j]] += rho[[i*d_c + k, j*d_c + k]];
    }}}
    out
  }
  fn trace_detected(rho: &CxMat, d_s: usize) -> CxMat{
    let d_c = rho.dim().0/d_s;
    let mut out = CxMat::zeros((d_c,d_c));
    for k in 0..d_c{ for l in 0..d_c{ for i in 0..d_s{
      out[[k,l]] += rho[[i*d_c + k, i*d_c + l]];
    }}}
    out
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_gcce_product_density_matrix(){
    use crate::physical_constants::{BOLTZMANN,JOULES_TO_HERTZ};
    use crate::quantum::block_diagonal_spin_hamiltonian::{
      get_bath_thermal_density_matrix,
    };
    use ndarray::linalg::kron;
    use ndarray_linalg::Norm;

    let (z0, z1, b) = (1000.0, 100.0, 30.0);
    let beta = 0.01; // 1/Hz
    let temperature = 1.0/(JOULES_TO_HERTZ*BOLTZMANN*beta);
    let mut config = Config::new();
    config.set_defaults().unwrap();
    config.temperature = Some(temperature);

    let det_h = (Array1::<f64>::zeros(1) ,CxMat::eye(1));
    let spin_ops = ClusterSpinOperators::new(det_h,&[2],3,&config).unwrap();
    let spin_indices = vec![0,1,2];

    // (a) No hyperfine coupling: ρ_S ⊗ exp(-βH_C⁰)/Z equals the joint
    // thermal state exp(-βH_{S+C})/Z, with ρ_S = exp(-βH_S)/Z.
    let tensors = build_restricted_three_spin_tensors(z0, z1, 0.0, 0.0, b);
    let h_s = build_single_spin_hamiltonian(&tensors, 0, true).unwrap();
    let rho_s = thermal_density_matrix(&h_s, temperature).unwrap();
    let bath_ham = build_bath_spin_hamiltonian(&spin_indices, &spin_ops,
        &tensors).unwrap();
    let product = get_gcce_product_density_matrix(&rho_s, &bath_ham,
        Some(temperature)).unwrap();
    let (e, v) = build_spin_hamiltonian(&spin_indices, &spin_ops, &tensors)
      .unwrap();
    let joint = get_electron_cluster_thermal_density_matrix(&e, &v, &config)
      .unwrap();
    assert!((&product - &joint).norm() < 1e-12);

    // (b) With hyperfine coupling: Tr_C ρ = ρ_S and Tr_S ρ = exp(-βH_C⁰)/Z,
    // where the latter is computed by the CCE code (block-diagonal ham0).
    let tensors = build_restricted_three_spin_tensors(z0, z1, 40.0, -25.0, b);
    let bath_ham = build_bath_spin_hamiltonian(&spin_indices, &spin_ops,
        &tensors).unwrap();
    let rho = get_gcce_product_density_matrix(&rho_s, &bath_ham,
        Some(temperature)).unwrap();
    assert!((trace_bath(&rho, 2) - &rho_s).norm() < 1e-12);
    let cce_ops = ClusterSpinOperators::new((Array1::<f64>::zeros(1),
        CxMat::eye(1)), &[2], 2, &config).unwrap();
    let rho_c = get_bath_thermal_density_matrix(&[1,2], &cce_ops, &tensors,
        temperature).unwrap();
    assert!((trace_detected(&rho, 2) - &rho_c).norm() < 1e-12);
    assert!((&rho - &kron(&rho_s, &rho_c)).norm() < 1e-12);

    // (c) No temperature: ρ_S ⊗ 1/d_C.
    let rho = get_gcce_product_density_matrix(&rho_s, &bath_ham, None)
      .unwrap();
    assert_eq!(rho, kron(&rho_s, &(CxMat::eye(4)/4.0)));

    // The empty cluster gives ρ_S.
    let bath_ham0 = build_bath_spin_hamiltonian(&[0], &spin_ops, &tensors)
      .unwrap();
    for t in [Some(temperature), None]{
      let rho = get_gcce_product_density_matrix(&rho_s, &bath_ham0, t)
        .unwrap();
      assert!((&rho - &rho_s).norm() < 1e-12);
    }

    // (d) The bath Hamiltonian is the full gCCE Hamiltonian minus the
    // detected-spin terms, z0·S_z + a1·S_z·I_z1 + a2·S_z·I_z2.
    let (a1, a2) = (40.0, -25.0);
    let (e, v) = build_spin_hamiltonian(&spin_indices, &spin_ops, &tensors)
      .unwrap();
    let mut d = CxMat::zeros((8,8));
    for k in 0..8{ d[[k,k]] = Complex::<f64>{re: e[k], im: 0.0}; }
    let full = v.dot(&d).dot(&v.t().map(|x| x.conj()));
    let sz0 = spin_ops.get(&SpinOp::Sz,2,3,0).unwrap();
    let sz1 = spin_ops.get(&SpinOp::Sz,2,3,1).unwrap();
    let sz2 = spin_ops.get(&SpinOp::Sz,2,3,2).unwrap();
    let detected = sz0*z0 + sz0.dot(sz1)*a1 + sz0.dot(sz2)*a2;
    assert!((&full - &detected - &bath_ham).norm() < 1e-12*full.norm());

    // Dimension error.
    assert_eq!(get_gcce_product_density_matrix(&CxMat::eye(3), &bath_ham,
          None), Err(CluEError::DetectedSpinStateDoesNotFitCluster(8,3)));
  }
}
