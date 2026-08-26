
use crate::quantum::cluster_operators::*;

use crate::physical_constants::*;
use crate::clue_errors::*;
use crate::config::Config;
use crate::config::pulse_sequence::PulseSequence;
use crate::HamiltonianTensors;
use crate::signal::Signal;
use crate::quantum::gcce_hamiltonian::{
  get_free_evolutions_propagators,
  get_2nd_free_evolutions_propagators,
};

use ndarray::{Array1,Array2};
use ndarray_linalg::{Eigh, UPLO, Trace};
use num_complex::Complex;

type CxMat = Array2::<Complex<f64>>;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This function evolves the quantum system according to the specified
/// pulse sequence and returns the time dependent expectation value of the
/// measurement.
pub fn propagate_pulse_sequence_block_diag(
    density_matrix: &CxMat, hamiltonian: &BlockDiagSpinHamiltonian, config: &Config)
  -> Result<Signal,CluEError>
{

  let Some(pulse_sequence) = &config.pulse_sequence else{
    return Err(CluEError::NoPulseSequence);
  };

  let tau_increments = &config.tau_increments;
  if tau_increments.is_empty(){
    return Err(CluEError::NoTimeIncrements);
  }

  let number_timepoints = &config.number_timepoints;
  if number_timepoints.is_empty(){
    return Err(CluEError::NoTimepoints);
  }
  let n_tot = config.get_total_number_timesteps();

  let mut signal = Vec::<Complex<f64>>::with_capacity(n_tot);

  let (_du_betas,u_betas) = get_free_evolutions_propagators(
      &hamiltonian.beta_eigvals,&hamiltonian.beta_eigvecs,
      config)?;

  let (_du_alphas,u_alphas) = get_free_evolutions_propagators(
      &hamiltonian.alpha_eigvals,&hamiltonian.alpha_eigvecs,
      config)?;

  let (_du2_betas,u2_betas) = match pulse_sequence{
    PulseSequence::RefocusedHahnEcho =>
      get_2nd_free_evolutions_propagators(
          &hamiltonian.beta_eigvals,&hamiltonian.beta_eigvecs, config)?,

    _ => (vec![CxMat::eye(0)],vec![CxMat::eye(0)]),
  }; 

  let (_du2_alphas,u2_alphas) = match pulse_sequence{
    PulseSequence::RefocusedHahnEcho =>
      get_2nd_free_evolutions_propagators(
          &hamiltonian.alpha_eigvals,&hamiltonian.alpha_eigvecs, config)?,

    _ => (vec![CxMat::eye(0)],vec![CxMat::eye(0)]),
  }; 


  for (idt,u_beta) in u_betas.iter().enumerate(){
    let u_alpha = &u_alphas[idt];
    let u_beta_dag = u_beta.t().map(|u_ij| u_ij.conj() );
    let u_alpha_dag = u_alpha.t().map(|u_ij| u_ij.conj() );
    
    let u: CxMat;
    match pulse_sequence{
      PulseSequence::CarrPurcell(0) => // FID
        u = u_alpha_dag.dot(u_beta),

      PulseSequence::CarrPurcell(1) => // Hahn echo
        u = u_alpha_dag.dot(&u_beta_dag.dot(&u_alpha.dot(u_beta))),  
      
      PulseSequence::CarrPurcell(2) => // CP-2
        u = u_beta_dag.dot(&u_alpha_dag.dot(&u_alpha_dag.dot(&u_beta_dag
              .dot(&u_alpha.dot(&u_beta.dot(&u_alpha.dot(u_beta))))))),
      
      PulseSequence::CarrPurcell(n_pi) => { // CP-n
        let u_aa = u_alpha.dot(u_alpha);
        let u_bb = u_beta.dot(u_beta);
        let exponent = ((*n_pi as f64 - 1.0)/2.0) as usize;
        let aabb1 = u_aa.dot(&u_bb);
        let bbaa1 = u_bb.dot(&u_aa);
        let mut aabb = aabb1.clone();
        let mut bbaa = bbaa1.clone();
        for _ii in 1..exponent{
          aabb = aabb1.dot(&aabb);
          bbaa = bbaa1.dot(&bbaa);
        }
        if n_pi%2 == 0{
          u = ((u_alpha.dot(&bbaa.dot(&u_bb.dot(u_alpha))))
              .t().map(|v| v.conj()))
            .dot( &(u_beta.dot(&aabb.dot(&u_aa.dot(u_beta)))));
        }else{
          u = ((u_beta.dot(&aabb.dot(u_alpha))).t().map(|v| v.conj()))
            .dot( &(u_alpha.dot(&bbaa.dot(u_beta))) );
        }
      },

      PulseSequence::RefocusedHahnEcho => {
        // Initialize u to satisfy the compiler.
        u = CxMat::eye(0);

        for (idt2,u2_beta) in u2_betas.iter().enumerate(){
          let u2_alpha = &u2_alphas[idt2];
          let u2_beta_dag = u2_beta.t().map(|u_ij| u_ij.conj() );
          let u2_alpha_dag = u2_alpha.t().map(|u_ij| u_ij.conj() );

          // In the full spin Hamiltonian,
          // U(2τ1 + 2τ2) = U0(τ2)Up(π)U0(τ2)U0(τ1)Up(π)U0(τ1)Up(π/2).
          // In the reduced spin-space where only spin Hamiltonian
          // is block diagonal is the electron spin's ms, U reduces to
          // U(2τ1 + 2τ2,-) = U0(τ2,-)U0(τ2,+)U0(τ1,+)U0(τ1,-).
          // or
          // U(2τ1 + 2τ2,+) = U0(τ2,+)U0(τ2,-)U0(τ1,-)U0(τ1,+)
          // with
          // U(2τ1 + 2τ2,+)^† = U0(τ1,+)^†U0(τ1,-)^†U0(τ2,-)^†U0(τ2,+)^†
          // With S+ as the detection operator,
          // <S+(2τ1 + 2τ2)> = <U(2τ1 + 2τ2,+)^† U(2τ1 + 2τ2,-)>.

          // U0(τ2,-)U0(τ2,+)U0(τ1,+)U0(τ1,-)
          let u_baab = u2_beta.dot(&u2_alpha.dot(&u_alpha.dot(u_beta)));

          // U0(τ1,+)^† U0(τ1,-)^† U0(τ2,-)^† U0(τ2,+)^†
          let u_abba_dag 
            = u_alpha_dag.dot(&u_beta_dag.dot(&u2_beta_dag.dot(&u2_alpha_dag)));
  
          let u_re = u_abba_dag.dot(&u_baab);
          let it = std::iter::zip(density_matrix,&u_re);
          let v = it.map(|(rho_ij,u_ij)| rho_ij*u_ij).sum::<Complex<f64>>();
          signal.push(v);
        }
      },
      PulseSequence::FreeEvolution => return Err(
          CluEError::PulseSequenceNotSupported(
            "propagate_pulse_sequence_block_diag".to_string(),
            "free evolution".to_string(),
            )),
      PulseSequence::Custom(_) => return Err(
          CluEError::PulseSequenceNotSupported(
            "propagate_pulse_sequence_block_diag".to_string(),
            "custom pulse sequences".to_string(),
            )),
    }

    if *pulse_sequence == PulseSequence::RefocusedHahnEcho { continue; }
    let it = std::iter::zip(density_matrix,&u);
    let v = it.map(|(rho_ij,u_ij)| rho_ij*u_ij).sum::<Complex<f64>>();
    signal.push(v);
  }

  Ok(Signal{data: signal})
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//------------------------------------------------------------------------------
/// This function builds the density matrix.
pub fn get_cluster_thermal_density_matrix(
    hamiltonian: &BlockDiagSpinHamiltonian, config: &Config)
  -> Result<CxMat,CluEError>
{

  let Some(temperature) = config.temperature else {
    return Err(CluEError::NoTemperature);
  };

  let beta = I/(temperature*BOLTZMANN/HBAR);
  
  let rho_alpha = get_propagators_complex_time_from_eig(
      &hamiltonian.alpha_eigvals, &hamiltonian.alpha_eigvecs,
      &[-beta])?;

  let rho_beta = get_propagators_complex_time_from_eig(
      &hamiltonian.beta_eigvals, &hamiltonian.beta_eigvecs,
      &[-beta])?;

  let mut density_matrix = &rho_alpha[0] - &rho_beta[0];

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
/// This function takes a Hamiltonian (Hz) and a vector of times (s), 
/// and calculates the propagator for each time.
/// For each time _t_, the propagator _U_(_t_) = exp(-i2π_tH_), 
/// where _H_ is a Hamiltonian in frequency units.
pub fn get_propagators(hamiltonian: &CxMat, times: &[f64]  )
  -> Result<Vec::<CxMat>,CluEError> 
{
  let Ok((eigvals, eigvecs)) = hamiltonian.eigh(UPLO::Lower) else{
    return Err(
        CluEError::CannotDiagonalizeHamiltonian(hamiltonian.to_string()));
  };

  let mut propagators = Vec::<CxMat>::with_capacity(times.len());

  let inv_eigvecs = eigvecs.t().map(|v| v.conj());

  for &t in times.iter(){
    let u_eig = CxMat::from_diag(&eigvals.map(|nu| 
          { let i_phase: Complex<f64> = (-I*2.0*PI*nu)*t;
            i_phase.exp() 
          }  
          )
        );

    let u = eigvecs.dot( &u_eig.dot( &inv_eigvecs) );

    propagators.push(u);

  }

  Ok(propagators)
}
//------------------------------------------------------------------------------
/// This function takes a Hamiltonian (Hz) and a vector of times (s), 
/// and calculates the propagator for each time.
/// For each time _t_, the propagator _U_(_t_) = exp(-i2π_tH_), 
/// where _H_ is a Hamiltonian in frequency units.
pub fn get_propagators_from_eig(
    eigvals: &Array1::<f64>, eigvecs: &CxMat, 
    times: &[f64]  )
  -> Result<Vec::<CxMat>,CluEError> 
{

  let mut propagators = Vec::<CxMat>::with_capacity(times.len());

  let inv_eigvecs = eigvecs.t().map(|v| v.conj());

  for &t in times.iter(){
    let u_eig = CxMat::from_diag(&eigvals.map(|nu| 
          { let i_phase: Complex<f64> = (-I*2.0*PI*nu)*t;
            i_phase.exp() 
          }  
          )
        );

    let u = eigvecs.dot( &u_eig.dot( &inv_eigvecs) );

    propagators.push(u);

  }

  Ok(propagators)
}
//------------------------------------------------------------------------------
/// This function takes a Hamiltonian (Hz) and a vector of times (s), 
/// and calculates the propagator for each time.
/// For each time _t_, the propagator _U_(_t_) = exp(-i2π_tH_), 
/// where _H_ is a Hamiltonian in frequency units.
pub fn get_propagators_complex_time(hamiltonian: &CxMat, 
    times: &[Complex<f64>]  )
  -> Result<Vec::<CxMat>,CluEError> 
{
  let Ok((eigvals, eigvecs)) = hamiltonian.eigh(UPLO::Lower) else{
    return Err(
        CluEError::CannotDiagonalizeHamiltonian(hamiltonian.to_string()));
  };

  let mut propagators = Vec::<CxMat>::with_capacity(times.len());

  let inv_eigvecs = eigvecs.t().map(|v| v.conj());

  for &t in times.iter(){
    let u_eig = CxMat::from_diag(&eigvals.map(|nu| 
          { let i_phase: Complex<f64> = (-I*2.0*PI*nu)*t;
            i_phase.exp() 
          }  
          )
        );

    let u = eigvecs.dot( &u_eig.dot( &inv_eigvecs) );

    propagators.push(u);

  }

  Ok(propagators)
}
//------------------------------------------------------------------------------
/// This function takes a Hamiltonian (Hz) and a vector of times (s), 
/// and calculates the propagator for each time.
/// For each time _t_, the propagator _U_(_t_) = exp(-i2π_tH_), 
/// where _H_ is a Hamiltonian in frequency units.
pub fn get_propagators_complex_time_from_eig(
    eigvals: &Array1::<f64>, eigvecs: &CxMat,
    times: &[Complex<f64>]  )
  -> Result<Vec::<CxMat>,CluEError> 
{

  let mut propagators = Vec::<CxMat>::with_capacity(times.len());

  let inv_eigvecs = eigvecs.t().map(|v| v.conj());

  for &t in times.iter(){
    let u_eig = CxMat::from_diag(&eigvals.map(|nu| 
          { let i_phase: Complex<f64> = (-I*2.0*PI*nu)*t;
            i_phase.exp() 
          }  
          )
        );

    let u = eigvecs.dot( &u_eig.dot( &inv_eigvecs) );

    propagators.push(u);

  }

  Ok(propagators)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `BlockDiagSpinHamiltonian` defines a block diagonal spin Hamiltonian, where
/// the detected spin's eigenstates are beta and alpha.
pub struct BlockDiagSpinHamiltonian{
  pub beta_eigvals: Array1::<f64>,
  pub beta_eigvecs: CxMat,
  pub alpha_eigvals: Array1::<f64>,
  pub alpha_eigvecs: CxMat,
}
impl BlockDiagSpinHamiltonian{
  pub fn new(beta: &CxMat,alpha: &CxMat) -> Result<Self,CluEError>
  {

    let Ok((beta_eigvals, beta_eigvecs)) = beta.eigh(UPLO::Lower) else{
      return Err(
          CluEError::CannotDiagonalizeHamiltonian(beta.to_string()));
    };
    let Ok((alpha_eigvals, alpha_eigvecs)) = alpha.eigh(UPLO::Lower) else{
      return Err(
          CluEError::CannotDiagonalizeHamiltonian(alpha.to_string()));
    };
  
    Ok(BlockDiagSpinHamiltonian{
      beta_eigvals, beta_eigvecs,
      alpha_eigvals, alpha_eigvecs,
    })
  }
  //----------------------------------------------------------------------------
  pub fn beta(&self) -> CxMat{
    let eigvals = CxMat::from_diag(&self.beta_eigvals
        .map(|v| Complex::<f64>{re: *v, im: 0.0}));
    let inv_eigvecs = self.beta_eigvecs.t().map(|v| v.conj());  

    self.beta_eigvecs.dot(&eigvals.dot(&inv_eigvecs))
  }
  //----------------------------------------------------------------------------
  pub fn alpha(&self) -> CxMat{
    let eigvals = CxMat::from_diag(&self.alpha_eigvals
        .map(|v| Complex::<f64>{re: *v, im: 0.0}));
    let inv_eigvecs = self.alpha_eigvecs.t().map(|v| v.conj());  

    self.alpha_eigvecs.dot(&eigvals.dot(&inv_eigvecs))
  }
}
//------------------------------------------------------------------------------
/// This function builds the cluster spin Hamiltonian assuming
/// < mS | H | mS' > = 0, for mS != mS',
pub fn build_block_diag_hamiltonian(spin_indices: &[usize],
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors, 
    config: &Config)
  -> Result<BlockDiagSpinHamiltonian,CluEError>
{

  let Some(central_spin_mult) = config.detected_spin_multiplicity else{
    return Err(CluEError::NoDetectedSpinMultiplicity);
  };

  let s = (central_spin_mult as f64 - 1.0)/2.0;
  let spin_ms: Vec::<f64> = (0..central_spin_mult).map(|n| (n as f64 - s))
    .collect();

  let Some(transition) = config.detected_spin_transition else {
    return Err(CluEError::NoCentralSpinTransition);
  };

  let ms_beta = spin_ms[transition[0]];
  let ms_alpha = spin_ms[transition[1]];


  let spin_multiplicities: Vec::<usize> = 
    spin_indices.iter().map(|idx| tensors.spin_multiplicities[*idx]).collect();

  let mut dim: usize = 1;
  spin_multiplicities.iter().for_each(|spin_mul| dim *= spin_mul);

  let mut ham0 = CxMat::zeros((dim,dim));
  let mut ham_ms = CxMat::zeros((dim,dim));
  
  // electron Zeeman
  if let Some(vec) = tensors.spin1_tensors.get(0){
    ham_ms = CxMat::eye(dim)*vec.z();
  }

  let cluster_size = spin_indices.len();


  for (sop_idx0, &ten_idx0) in spin_indices.iter().enumerate(){

    let spin_mult0 = tensors.spin_multiplicities[ten_idx0];
    let sx0 = spin_ops.get(&SpinOp::Sx,spin_mult0,cluster_size,sop_idx0)?;
    let sy0 = spin_ops.get(&SpinOp::Sy,spin_mult0,cluster_size,sop_idx0)?;
    let sz0 = spin_ops.get(&SpinOp::Sz,spin_mult0,cluster_size,sop_idx0)?;

    // nuclear Zeeman
    if let Some(vec) = tensors.spin1_tensors.get(ten_idx0){
      ham0 = ham0 + sx0*vec.x();
      ham0 = ham0 + sy0*vec.y();
      ham0 = ham0 + sz0*vec.z();
    }

    // nuclear hyperfine
    if let Some(ten) = tensors.spin2_tensors.get(0,ten_idx0){
      ham_ms = ham_ms + sx0*ten.zx();
      ham_ms = ham_ms + sy0*ten.zy();
      ham_ms = ham_ms + sz0*ten.zz();
    }

    for (sop_idx1, &ten_idx1) in spin_indices.iter().enumerate().skip(sop_idx0){

      let spin_mult1 = tensors.spin_multiplicities[ten_idx1];
      let sx1 = spin_ops.get(&SpinOp::Sx,spin_mult1,cluster_size,sop_idx1)?;
      let sy1 = spin_ops.get(&SpinOp::Sy,spin_mult1,cluster_size,sop_idx1)?;
      let sz1 = spin_ops.get(&SpinOp::Sz,spin_mult1,cluster_size,sop_idx1)?;
      
      // dipole-dipole, and electric quadrupole
      if let Some(ten) = tensors.spin2_tensors.get(ten_idx0,ten_idx1){
        ham0 = ham0 + sx0.dot(sx1)*ten.xx();
        ham0 = ham0 + sx0.dot(sy1)*ten.xy();
        ham0 = ham0 + sx0.dot(sz1)*ten.xz();
        ham0 = ham0 + sy0.dot(sx1)*ten.yx();
        ham0 = ham0 + sy0.dot(sy1)*ten.yy();
        ham0 = ham0 + sy0.dot(sz1)*ten.yz();
        ham0 = ham0 + sz0.dot(sx1)*ten.zx();
        ham0 = ham0 + sz0.dot(sy1)*ten.zy();
        ham0 = ham0 + sz0.dot(sz1)*ten.zz();
      }
     
      
    }
  }

  let beta = ham0.clone() + ham_ms.clone()*ms_beta;
  let alpha = ham0 + ham_ms*ms_alpha;

  BlockDiagSpinHamiltonian::new(&beta,&alpha)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#[cfg(test)]
mod tests {
  use super::*;
  use crate::space_3d::{SymmetricTensor3D,Vector3D};
  use crate::quantum::tensors::*;
  use ndarray::array;
  use crate::cluster_methods::appa::{
    appa_hahn,
    appa_hahn_frequency,
  };
  use crate::quantum::spin_states::SpinStates;


  use rand_chacha::ChaCha20Rng;
  use rand::SeedableRng;


  //----------------------------------------------------------------------------
  #[test]
  fn test_propagate_pulse_sequence_block_diag(){

    let z0 = 33.0e9;
    let z1 = 80.0e6;
    let a1 = 10.0e6;
    let a2 = -10.0e6;
    let b = 10.0e3;
    let tensors = build_restricted_three_spin_tensors(z0, z1, a1, a2, b);
  
    let spin_indices = vec![1,2];
    
    let mut config = Config::new();
    config.set_defaults().unwrap();

    let spin_ops = ClusterSpinOperators::new(1,&vec![2],2,&config).unwrap();

    let nt = 21;
    config.number_timepoints = vec![nt];
    let delta_hf = a1 - a2;
    let freq = appa_hahn_frequency(delta_hf,b);
    config.tau_increments = vec![0.05/freq];
    config.pulse_sequence = Some(PulseSequence::CarrPurcell(1));

    config.set_defaults().unwrap();
    config.set_tau_axis().unwrap();
  
    let hamiltonian = build_block_diag_hamiltonian(&spin_indices,&spin_ops, &tensors,
        &config).unwrap();

    let mut rng = ChaCha20Rng::from_entropy();
    let states = SpinStates::generate(&mut rng, &tensors,&config).unwrap();
    let density_matrix = states.density_matrix_for(&spin_indices).unwrap()
        .unwrap();

    let signal = propagate_pulse_sequence_block_diag(
        &density_matrix, &hamiltonian, &config).unwrap();

    assert_eq!(signal.data.len(),nt);

    let ref_signal_opt = appa_hahn(
        &spin_indices,&tensors,&config).unwrap();
    let Some(ref_signal) = ref_signal_opt else{
      panic!("Could not calculate reference signal.");
    };

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
  fn test_build_block_diag_hamiltonian(){

    let z0 = 1000.0;
    let z1 = 100.0;
    let a1 = 2.0;
    let a2 = -1.0;
    let b = 0.1;

    let tensors = build_restricted_three_spin_tensors(z0, z1, a1, a2, b);
  

    let mut config = Config::new();
    config.set_defaults().unwrap();

    let spin_indices = vec![1,2];
    let spin_ops = ClusterSpinOperators::new(1,&vec![2],2,&config).unwrap();
    let mut config = Config::new();
    config.set_defaults().unwrap();

    let hamiltonian = build_block_diag_hamiltonian(&spin_indices,&spin_ops, &tensors,
        &config).unwrap();

    let ms = 0.5;
    let z0 = ms*z0*ONE;
    let z1 = z1*ONE;
    let a1 = a1*ONE;
    let a2 = a2*ONE;
    let b = b*ONE;
    
    let beta = array![
      [-z0 + z1 - (a1+a2)/4.0 + b/4.0,ZERO,ZERO,ZERO],
      [ZERO,-z0 -(a1-a2)/4.0 - b/4.0 ,-b/4.0,ZERO],
      [ZERO,-b/4.0, -z0 +(a1-a2)/4.0 - b/4.0,ZERO],
      [ZERO,ZERO,ZERO,-z0 -z1 + (a1+a2)/4.0 + b/4.0],
    ];
    
      let alpha = array![
      [z0 + z1 + (a1+a2)/4.0 + b/4.0,ZERO,ZERO,ZERO],
      [ZERO,z0 +(a1-a2)/4.0 - b/4.0 ,-b/4.0,ZERO],
      [ZERO,-b/4.0, z0 -(a1-a2)/4.0 - b/4.0,ZERO],
      [ZERO,ZERO,ZERO,z0 -z1 - (a1+a2)/4.0 + b/4.0],
    ];


    let halpha = hamiltonian.alpha();
    let hbeta = hamiltonian.beta();

    assert!(approx_eq(&hbeta, &beta, 1e-12));
    assert!(approx_eq(&halpha, &alpha, 1e-12));
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
  #[test]
  fn test_get_cluster_thermal_density_matrix(){
    let sz = spin_z(2);
    let delta_energy = 416732382466.5515; // kB*T/h at T = 20 K.
    let beta = (sz.clone() -  CxMat::eye(2))*delta_energy;
    let alpha = (sz.clone() +  CxMat::eye(2))*delta_energy;
    let spin_hamiltonian = BlockDiagSpinHamiltonian::new(&beta,&alpha).unwrap();

    let mut config = Config::new();
    config.temperature = Some(20.0);
    let density_matrix = get_cluster_thermal_density_matrix(
        &spin_hamiltonian,&config).unwrap();

    let e_inv = (-ONE).exp();
    let z = e_inv*e_inv + e_inv;
    let expected = Array2::from_diag(&array![e_inv*e_inv, e_inv])/z;
    assert!(approx_eq(&density_matrix, &expected, 1e-12));

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_propagators_complex_time() {
    let sy = spin_y(2);
    let ham = sy*2.0;

    let times = vec![ONE/8.0];

    let propagators = get_propagators_complex_time(&ham,&times).unwrap();
    let u = &propagators[0];

    for irow in  0..2 {
      for icol in 1..2 {
        let err: f64 = (u[[irow,icol]].conj()*u[[irow,icol]] - 0.5*ONE).norm();
        assert!(err < 1e-12);
      }
    }

  }  
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_propagators() {
    let sy = spin_y(2);
    let ham = sy*2.0;

    let times = vec![1.0/8.0];

    let propagators = get_propagators(&ham,&times).unwrap();
    let u = &propagators[0];

    for irow in  0..2 {
      for icol in 1..2 {
        let err: f64 = (u[[irow,icol]].conj()*u[[irow,icol]] - 0.5).norm();
        assert!(err < 1e-12);
      }
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
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

}



