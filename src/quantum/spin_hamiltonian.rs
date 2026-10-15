use crate::physical_constants::*;
use crate::clue_errors::*;
use crate::config::{Config,PulseSequence,DensityMatrixMethod};
use crate::HamiltonianTensors;
use crate::signal::Signal;

use std::fmt;
use ndarray::Array2;
use ndarray_linalg::{Eigh, UPLO, Trace};
use ndarray::linalg::kron;
use num_complex::Complex;

type CxMat = Array2::<Complex<f64>>;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This function evolves the quantum system according to the specified
/// pulse sequence and returns the time dependent expectation value of the
/// measurement.
pub fn propagate_pulse_sequence(
    density_matrix: &CxMat, hamiltonian: &SpinHamiltonian, config: &Config)
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
  let n_tot = number_timepoints.iter().sum::<usize>();

  let mut signal = Vec::<Complex<f64>>::with_capacity(n_tot);

  let du_betas = get_propagators(&hamiltonian.beta,tau_increments)?;
  let du_alphas = get_propagators(&hamiltonian.alpha,tau_increments)?;


  let mut u_beta = CxMat::eye(du_betas[0].dim().0);
  let mut u_alpha = CxMat::eye(du_alphas[0].dim().0);
  let mut u_beta_dag = CxMat::eye(du_betas[0].dim().0);
  let mut u_alpha_dag = CxMat::eye(du_alphas[0].dim().0);

  for (idt,_dt) in tau_increments.iter().enumerate(){
    let n_timepoints = number_timepoints[idt];
    for _inumt in 0..n_timepoints{
    
      let u: CxMat;
      match pulse_sequence{
        PulseSequence::CarrPurcell(0) => // FID
          u = u_alpha_dag.dot(&u_beta),

        PulseSequence::CarrPurcell(1) => // Hahn echo
          u = u_alpha_dag.dot(&u_beta_dag.dot(&u_alpha.dot(&u_beta))),  
        
        PulseSequence::CarrPurcell(2) => // CP-2
          u = u_beta_dag.dot(&u_alpha_dag.dot(&u_alpha_dag.dot(&u_beta_dag
                .dot(&u_alpha.dot(&u_beta.dot(&u_alpha.dot(&u_beta))))))),
        
        PulseSequence::CarrPurcell(n_pi) => { // CP-n
          let u_aa = u_alpha.dot(&u_alpha);
          let u_bb = u_beta.dot(&u_beta);
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
            u = ((u_alpha.dot(&bbaa.dot(&u_bb.dot(&u_alpha))))
                .t().map(|v| v.conj()))
              .dot( &(u_beta.dot(&aabb.dot(&u_aa.dot(&u_beta)))));
          }else{
            u = ((u_beta.dot(&aabb.dot(&u_alpha))).t().map(|v| v.conj()))
              .dot( &(u_alpha.dot(&bbaa.dot(&u_beta))) );
          }
        },
      }

      let it = std::iter::zip(density_matrix,&u);
      let v = it.map(|(rho_ij,u_ij)| rho_ij*u_ij).sum::<Complex<f64>>();
      signal.push(v);


      u_beta = du_betas[idt].dot(&u_beta);
      u_alpha = du_alphas[idt].dot(&u_alpha);

      u_beta_dag = u_beta.t().map(|u_ij| u_ij.conj() );
      u_alpha_dag = u_alpha.t().map(|u_ij| u_ij.conj() );
    }
  }

  Ok(Signal{data: signal})
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/*
/// This function is still under construction.
/// Please be patient.
pub fn propagate_pulse_sequence_gcce(
    density_matrix: &CxMat, hamiltonian: &CxMat, config: &Config)
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
  let n_tot = number_timepoints.iter().sum::<usize>();

  let mut signal = Vec::<Complex<f64>>::with_capacity(n_tot);

  let dus = get_propagators(&hamiltonian,tau_increments)?;
  let mut u_of_tau = CxMat::eye(dus[0].dim().0);

  let Some(spin_multiplicity) = config.detected_spin_multiplicity else{
    return Err(CluEError::NoDetectedSpinMultiplicity);
  };

  let Some(transition) = &config.detected_spin_transition else{
    return Err(CluEError::NoDetectedSpinTransition);
  };

  let u_pi = ideal_pulse(&SpinOp::Sx, PI, 
    spin_multiplicity, transition)?;

  let detection_op0 = CxMat::eye(dus[0].dim().0);

  for (idt,_dt) in tau_increments.iter().enumerate(){
    let n_timepoints = number_timepoints[idt];
    for _inumt in 0..n_timepoints{
    
      let mut u_sequence = CxMat::eye(dus[0].dim().0);

      match pulse_sequence{
        PulseSequence::CarrPurcell(n_pi) => { // CP-n
          let u_seg = u_of_tau.dot(&u_pi.dot(&u_of_tau));
          for _ii in 0..*n_pi{
            u_sequence = u_sequence.dot(&u_seg);
          }
        }
      }

      let u_sequence_dag = u_sequence.t().map(|u_ij| u_ij.conj() );
      
      let detection_op = u_sequence_dag.dot(&detection_op0.dot(&u_sequence));

      let it = std::iter::zip(density_matrix,&detection_op);
      let v = it.map(|(rho_ij,u_ij)| rho_ij*u_ij).sum::<Complex<f64>>();
      signal.push(v);


      u_of_tau = dus[idt].dot(&u_of_tau);
    }
  }

  panic!("CluE gCCE is still under construction.");
  //Ok(Signal{data: signal})
}
*/
//------------------------------------------------------------------------------
/// This function builds the density matrix.
pub fn get_density_matrix(hamiltonian: &SpinHamiltonian, config: &Config)
  -> Result<CxMat,CluEError>
{

  let Some(density_matrix_method) = &config.density_matrix else{
    return Err(CluEError::NoDensityMatrixMethod);
  };

  let mut density_matrix: CxMat;

  match density_matrix_method{
    // TODO: Deprecate ApproxThermal.
    DensityMatrixMethod::ApproxThermal(temperature) => {

      let beta = I/(temperature*BOLTZMANN/HBAR);
      
      let mean_hamiltonian = (&hamiltonian.beta + &hamiltonian.alpha)/2.0;
      
      let mut rhos = get_propagators_complex_time(&mean_hamiltonian,
          &[-beta])?;

      density_matrix = rhos.remove(0);
    },
    DensityMatrixMethod::Identity => {
      let dim = hamiltonian.beta.dim().0;
      density_matrix = CxMat::eye(dim);
    },
    DensityMatrixMethod::Thermal(temperature) => {

      let beta = I/(temperature*BOLTZMANN/HBAR);
      
      let rho_alpha = get_propagators_complex_time(&hamiltonian.alpha,
          &[-beta])?;

      let rho_beta = get_propagators_complex_time(&hamiltonian.beta,
          &[-beta])?;

      density_matrix = &rho_alpha[0] - &rho_beta[0];
    },
  }

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
/*
fn ideal_two_state_pulse(spin_op: &SpinOp, angle: f64)
  -> CxMat
{


  let e = CxMat::eye(2);
  let i_sigma = get_spin_operator(2,spin_op)*(I*2.0);

  e*( (angle/2.0).cos() ) + i_sigma*((angle/2.0).sin() )


}
*/
//------------------------------------------------------------------------------
/*
fn ideal_pulse(spin_op: &SpinOp, angle: f64, 
    spin_multipliciy: usize, transition: &[usize;2])
  -> Result<CxMat,CluEError>
{
  let u2 = ideal_two_state_pulse(spin_op, angle);

  let mut u = CxMat::eye(spin_multipliciy);

  for (ii,m) in transition.iter().enumerate(){
    let idx0 = spin_multipliciy - m - 1;
    for (jj,n) in transition.iter().enumerate(){
      let idx1 = spin_multipliciy - n - 1;

      u[[idx0,idx1]] = u2[[ii,jj]];
    }
  }

  Ok(u)
}
*/
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `SpinHamiltonian` defines a block diagonal spin Hamiltonian, where
/// the detected spin's eigenstates are beta and alpha.
pub struct SpinHamiltonian{
  beta: CxMat,
  alpha: CxMat,
}
//------------------------------------------------------------------------------
/// This function builds the cluster spin Hamiltonian assuming
/// < mS | H | mS' > = 0, for mS != mS',
pub fn build_hamiltonian(spin_indices: &[usize],
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors, 
    config: &Config)
  -> Result<SpinHamiltonian,CluEError>
{


  let Some(central_spin) = config.detected_spin_identity else{
    return Err(CluEError::NoCentralSpin);
  };

  let central_spin_mult = central_spin.spin_multiplicity();
  let s = (central_spin_mult as f64 - 1.0)/2.0;
  let spin_ms: Vec::<f64> = (0..central_spin_mult).map(|n| n as f64 - s)
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
  
Ok(SpinHamiltonian{beta,alpha})
}
//------------------------------------------------------------------------------
/// This function builds the cluster spin Hamiltonian without assuming
/// < mS | H | mS' > = 0, for mS != mS',
pub fn build_spin_hamiltonian(spin_indices: &[usize],
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors)
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

    let spin_mult0 = tensors.spin_multiplicities[ten_idx0];
    let sx0 = spin_ops.get(&SpinOp::Sx,spin_mult0,cluster_size,sop_idx0)?;
    let sy0 = spin_ops.get(&SpinOp::Sy,spin_mult0,cluster_size,sop_idx0)?;
    let sz0 = spin_ops.get(&SpinOp::Sz,spin_mult0,cluster_size,sop_idx0)?;

    // Zeeman
    if let Some(vec) = tensors.spin1_tensors.get(ten_idx0){
      ham = ham + sx0*vec.x();
      ham = ham + sy0*vec.y();
      ham = ham + sz0*vec.z();
    }

    for (sop_idx1, &ten_idx1) in spin_indices.iter().enumerate().skip(sop_idx0){

      let spin_mult1 = tensors.spin_multiplicities[ten_idx1];
      let sx1 = spin_ops.get(&SpinOp::Sx,spin_mult1,cluster_size,sop_idx1)?;
      let sy1 = spin_ops.get(&SpinOp::Sy,spin_mult1,cluster_size,sop_idx1)?;
      let sz1 = spin_ops.get(&SpinOp::Sz,spin_mult1,cluster_size,sop_idx1)?;
      
      // hyperfine, dipole-dipole, and electric quadrupole
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

  Ok(ham)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// 'ClusterSpinOperators' contains the spin operators used for clusters of 
/// spins.
/// `max_size` is the maximum cluster size specified.
/// 'spin_multiplicities' contains the allowed spin multiplicities of the
/// spins within each cluster (they must all be the ssame).
/// 'cluster_spin_ops' contains the matrices.
pub struct ClusterSpinOperators {
  max_size: usize, 
  spin_multiplicities: Vec<usize>, 
  cluster_spin_ops: Vec<KronSpinOpXYZ>,
}
//------------------------------------------------------------------------------
impl<'a> ClusterSpinOperators {
  /// This function builds 'ClusterSpinOperators' for clusters of 
  /// `spin_multiplicities` up to size `max_size`.
  pub fn new(spin_multiplicities: &[usize], max_size: usize) 
   -> Result<ClusterSpinOperators, CluEError> {
    
    let n_mults = spin_multiplicities.len();

    let mut cluster_spin_ops = Vec::<KronSpinOpXYZ>::with_capacity(n_mults);

    for spin_multiplicity in spin_multiplicities.iter() {
      let sops = KronSpinOpXYZ::new(*spin_multiplicity, max_size)?;
      cluster_spin_ops.push(sops);
    }

    Ok(ClusterSpinOperators{
        max_size,
        spin_multiplicities: spin_multiplicities.to_owned(),
        cluster_spin_ops,
        })
  }

  //----------------------------------------------------------------------------
  /// This function tries to find the matrix corresponding to the specified 
  /// spin operator for a spin of the specified multiplicity in a cluster
  /// of the indicated size, where the single spin operator has `op_pos`
  /// within the tensor product.
  pub fn get(&'a self, 
      sop: &SpinOp, spin_multiplicity: usize, cluster_size: usize,
      op_pos: usize) -> Result<&'a CxMat,CluEError> {
  
    if cluster_size > self.max_size {
      return Err(CluEError::NoSpinOpForClusterSize(cluster_size,self.max_size));
    }

    for (ii,&ispin_mult) in self.spin_multiplicities.iter().enumerate() {
      if ispin_mult == spin_multiplicity {
        
       let sop = self.cluster_spin_ops[ii]
         .get(sop, op_pos,cluster_size)?;
       return Ok(sop);
      }
    }
    Err(CluEError::NoSpinOpWithMultiplicity(spin_multiplicity))
  }


}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `KronSpinOpXYZ` contains the Sx, Sy, and Sz spin operators tensored
/// with various identity matrices on either side.
pub struct KronSpinOpXYZ {
  sx_list: KronSpinOpList,
  sy_list: KronSpinOpList,
  sz_list: KronSpinOpList,
}

impl<'a> KronSpinOpXYZ {

  /// This function builds `KronSpinOpXYZ` for spin with the input spin 
  /// multiplicity for clusters up to size `max_size`.
  pub fn new(spin_multiplicity: usize, max_size: usize) 
    -> Result<KronSpinOpXYZ,CluEError> 
  {

    let sx_list = KronSpinOpList::new(spin_multiplicity, SpinOp::Sx, max_size)?;
    let sy_list = KronSpinOpList::new(spin_multiplicity, SpinOp::Sy, max_size)?;
    let sz_list = KronSpinOpList::new(spin_multiplicity, SpinOp::Sz, max_size)?;

    Ok(KronSpinOpXYZ{
      sx_list,  
      sy_list,  
      sz_list,  
    })
  }
  //----------------------------------------------------------------------------
  /// This function tries to find the matrix corresponding to the specified 
  /// spin operator in a cluster of size `n_ops`,
  /// where the single spin operator has `op_pos` within the tensor product.
  pub fn get(&'a self, sop: &SpinOp, op_pos: usize, n_ops: usize) 
    -> Result<&'a CxMat,CluEError> 
  {

     match sop {
       SpinOp::Sx => self.sx_list.get(op_pos,n_ops),
       SpinOp::Sy => self.sy_list.get(op_pos,n_ops),
       SpinOp::Sz => self.sz_list.get(op_pos,n_ops),
       _ => Err(CluEError::CannotFindSpinOp(sop.to_string())),
     }
  } 
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `KronSpinOpList` contains a spin operator tensored with various identity 
/// matrices on either side.
pub struct KronSpinOpList {
  sop_list: Vec::<CxMat>,
}
//------------------------------------------------------------------------------
impl<'a> KronSpinOpList {

  /// This function builds `KronSpinOpList` with the specified 
  /// spin operator, for spins with the input spin 
  /// multiplicity for clusters up to size `max_size`.
  pub fn new(
      spin_multiplicity: usize,
      spin_operator: SpinOp,
      max_size: usize) -> Result<KronSpinOpList,CluEError> {

    let n_ops = (max_size*( max_size+1) )/2;
    let mut sop_list = Vec::<CxMat>::with_capacity(n_ops);

    for n_ops in 1..=max_size {
      let spin_mults = vec![spin_multiplicity; n_ops];

      for p_idx in 0..n_ops {

        let mut sops = vec![SpinOp::E; n_ops];
        sops[p_idx] = spin_operator;

        let sop = kron_spin_op(&spin_mults,&sops)?;
        sop_list.push(sop);
      }
    }

    Ok(KronSpinOpList {
      sop_list,
    })
  }
  //----------------------------------------------------------------------------
  /// This function tries to find the matrix of `n_ops` operators,
  /// with the `op_pos` operator being the non-identity within the 
  /// tensor product.
  pub fn get(&'a self, op_pos: usize, n_ops: usize) ->
    Result<&'a CxMat,CluEError> {


    if op_pos >= n_ops {
      return Err(CluEError::UnavailableSpinOp(op_pos,n_ops));
    }

    let idx = KronSpinOpList::get_index(op_pos, n_ops);

    if idx >= self.sop_list.len(){
      return Err(CluEError::UnavailableSpinOp(idx,self.sop_list.len()));
    }

    Ok(&self.sop_list[idx])
  }
  //----------------------------------------------------------------------------
  // This function translates the position within the tensor product to the
  // index within the list where that product is stored.
  fn get_index(op_pos: usize, n_ops: usize) -> usize {
    ( n_ops*(n_ops - 1) )/2 + op_pos
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//------------------------------------------------------------------------------
/// This function generates the matrix for spin operators, `sops`,
/// corresponding to spins with `spin_mults` spin multiplicities.
pub fn kron_spin_op(spin_mults: &[usize], sops: &[SpinOp]) ->
  Result<CxMat,CluEError> 
{

  if spin_mults.len() !=  sops.len() {
    return Err(CluEError::UnequalLengths(
          "spin_multiplicities".to_string(),
          spin_mults.len(),
          "spin operators".to_string(),
          sops.len(),
    ));
  }

  let mut sop = CxMat::eye(1);
  for ii in 0..spin_mults.len() {
    let s = get_spin_operator(spin_mults[ii],&sops[ii]);
    sop = kron(&sop,&s);
  }

  Ok(sop)
}
//------------------------------------------------------------------------------
/// This enum list possible spin operators. 
#[derive(Copy,Debug,Clone,PartialEq)]
pub enum SpinOp{
  E,
  Sx,
  Sy,
  Sz,
  Sp,
  Sm,
  S2,
}
impl fmt::Display for SpinOp {
    // This function translates `SpinOp` to strings.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
      match self{
        SpinOp::E => write!(f, "E"),
        SpinOp::Sx => write!(f, "Sx"),
        SpinOp::Sy => write!(f, "Sy"),
        SpinOp::Sz => write!(f, "Sz"),
        SpinOp::Sp => write!(f, "S+"),
        SpinOp::Sm => write!(f, "S-"),
        SpinOp::S2 => write!(f, "S^2"),
      }
    }
}
//------------------------------------------------------------------------------
// This function builds the matrix corresponding to the specified spin operator
// and multiplicity.
fn get_spin_operator(spin_multiplicity: usize, sop: &SpinOp) -> CxMat{
  match sop {
    SpinOp::E => CxMat::eye(spin_multiplicity),
    SpinOp::Sx => spin_x(spin_multiplicity),
    SpinOp::Sy => spin_y(spin_multiplicity),
    SpinOp::Sz => spin_z(spin_multiplicity),
    SpinOp::Sp => spin_plus(spin_multiplicity),
    SpinOp::Sm => spin_minus(spin_multiplicity),
    SpinOp::S2 => spin_squared(spin_multiplicity),
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//------------------------------------------------------------------------------
/// This function generates the Sx matrix:
/// <m'|Sx|m> = 1/2*(delta_{m',m+1} + delta_{m'+1,1})*sqrt(S*(S+1) - m'*m).
pub fn spin_x(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity - 1 {
    let n = ii as f64;
    let ms = spin - n;
    let value = 0.5*ONE*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt();
    op[[ii,ii+1]] = value;
    op[[ii+1,ii]] = value;
  }

  op
}
//------------------------------------------------------------------------------
/// This function generates the Sy matrix:
/// <m'|Sy|m> = -i/2*(delta_{m',m+1} - delta_{m'+1,1})*sqrt(S*(S+1) - m'*m).
pub fn spin_y(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity - 1 {
    let n = ii as f64;
    let ms = spin - n;
    let value = -0.5*I*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt();
    op[[ii,ii+1]] = value;
    op[[ii+1,ii]] = -value;
  }

  op
}
//------------------------------------------------------------------------------
/// This function generates the Sz matrix:
/// <m'|Sz|m> = delta_{m',m}*m.
pub fn spin_z(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity  {
    let n = ii as f64;
    let ms = spin - n;
    op[[ii,ii]] =  ms*ONE;
  }

  op
}
//------------------------------------------------------------------------------
/// This function generates the lowering ladder operator matrix:
/// <m'|S-|m> = delta_{m'+1,m} * sqrt(S*(S+1) - m'*m).
pub fn spin_minus(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity - 1 {
    let n = ii as f64;
    let ms = spin - n;
    op[[ii+1,ii]] = ONE*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt();
  }

  op
}
//------------------------------------------------------------------------------
/// This function generates the raising ladder operator matrix:
/// <m'|S+|m> = delta_{m',m+1} * sqrt(S*(S+1) - m'*m).
pub fn spin_plus(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity - 1 {
    let n = ii as f64;
    let ms = spin - n;
    op[[ii,ii+1]] = ONE*(spin*(spin+1.0) - (ms - 1.0)*ms).sqrt();
  }

  op
}
//------------------------------------------------------------------------------
/// This function generates the S^2 matrix:
/// <m'|S^2|m> = delta_{m',m}*S*(S+1).
pub fn spin_squared(spin_multiplicity: usize) -> CxMat {

  let spin: f64 = (spin_multiplicity as f64)/2.0 - 0.5;
  let mut op = CxMat::zeros((spin_multiplicity,spin_multiplicity));

  if spin_multiplicity == 0 {return op;}

  for ii in 0..spin_multiplicity  {
    let value = ONE*spin*(spin+1.0);
    op[[ii,ii]] = value;
  }

  op
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



#[cfg(test)]
mod tests {
  use super::*;
  use crate::space_3d::{SymmetricTensor3D,Vector3D};
  use crate::quantum::tensors::*;
  use ndarray::array;
  use crate::signal::calculate_analytic_restricted_2cluster_signals::{
    analytic_restricted_2cluster_signal,
    hahn_three_spin_modulation_frequency,
    //hahn_three_spin_modulation_depth
  };


  //----------------------------------------------------------------------------
  #[test]
  fn test_propagate_pulse_sequence(){

    let z0 = 33.0e9;
    let z1 = 80.0e6;
    let a1 = 10.0e6;
    let a2 = -10.0e6;
    let b = 10.0e3;
    let tensors = build_restricted_three_spin_tensors(z0, z1, a1, a2, b);
  
    let spin_indices = vec![1,2];
    
    let spin_ops = ClusterSpinOperators::new(&vec![2],2).unwrap();

    let mut config = Config::new();
    let nt = 21;
    config.number_timepoints = vec![nt];
    let delta_hf = a1 - a2;
    let freq = hahn_three_spin_modulation_frequency(delta_hf,b);
    config.tau_increments = vec![0.05/freq];
    config.pulse_sequence = Some(PulseSequence::CarrPurcell(1));

    config.set_defaults().unwrap();
    config.set_time_axis().unwrap();
  
    let hamiltonian = build_hamiltonian(&spin_indices,&spin_ops, &tensors,
        &config).unwrap();

    config.density_matrix = Some(DensityMatrixMethod::Identity);
    let density_matrix = get_density_matrix(&hamiltonian,&config).unwrap();

    let signal = propagate_pulse_sequence(&density_matrix, &hamiltonian, 
        &config).unwrap();

    assert_eq!(signal.data.len(),nt);

    let ref_signal_opt = analytic_restricted_2cluster_signal(
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
  /*
  #[test]
  fn test_propagate_pulse_sequence_gcce(){
    assert!(false);
  }
  */
  //----------------------------------------------------------------------------
  /*
  #[test]
  fn test_ideal_two_state_pulse(){
    assert!(false);
  }
  */
  //----------------------------------------------------------------------------
  /*
  #[test]
  fn test_ideal_pulse(){
    assert!(false);
  }
  */
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_hamiltonian(){

    let z0 = 1000.0;
    let z1 = 100.0;
    let a1 = 2.0;
    let a2 = -1.0;
    let b = 0.1;

    let tensors = build_restricted_three_spin_tensors(z0, z1, a1, a2, b);
  

    let spin_indices = vec![1,2];
    let spin_ops = ClusterSpinOperators::new(&vec![2],2).unwrap();
    let mut config = Config::new();
    config.set_defaults().unwrap();

    let hamiltonian = build_hamiltonian(&spin_indices,&spin_ops, &tensors,
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

    assert!(approx_eq(&hamiltonian.beta, &beta, 1e-12));
    assert!(approx_eq(&hamiltonian.alpha, &alpha, 1e-12));
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
  

    let spin_indices = vec![0,1,2];
    let spin_ops = ClusterSpinOperators::new(&vec![2],3).unwrap();

    let hamiltonian = build_spin_hamiltonian(&spin_indices,&spin_ops, &tensors,
        ).unwrap();

    let mut config = Config::new();
    config.set_defaults().unwrap();

    let spin_indices = vec![1,2];

    let block_hamiltonian = build_hamiltonian(&spin_indices,&spin_ops, &tensors,
        &config).unwrap();


    assert_eq!(hamiltonian.len(), 4*block_hamiltonian.beta.len());

    let beta = hamiltonian.slice(ndarray::s![4..,4..]).to_owned();
    assert!(approx_eq(&block_hamiltonian.beta, &beta, 1e-12));

    let alpha = hamiltonian.slice(ndarray::s![0..4,0..4]).to_owned();
    assert!(approx_eq(&block_hamiltonian.alpha, &alpha, 1e-12));


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
      } 
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_density_matrix(){
    let sz = spin_z(2);
    let delta_energy = 416732382466.5515; // kB*T/h at T = 20 K.
    let ham = sz*delta_energy;

    let spin_hamiltonian = SpinHamiltonian{beta: ham.clone(),alpha: ham};

    let mut config = Config::new();

    config.density_matrix = Some(DensityMatrixMethod::Identity);
    let density_matrix = get_density_matrix(&spin_hamiltonian,&config).unwrap();

    let mut expected = CxMat::eye(2)/2.0;

    assert!(approx_eq(&density_matrix, &expected, 1e-12));

    config.density_matrix = Some(DensityMatrixMethod::ApproxThermal(20.0));
    
    let density_matrix = get_density_matrix(&spin_hamiltonian,&config).unwrap();


    let z = (0.5*ONE).exp() + (-0.5*ONE).exp();
    expected[[0,0]] = (-0.5*ONE).exp()/z;
    expected[[1,1]] = (0.5*ONE).exp()/z;

    assert!(approx_eq(&density_matrix, &expected, 1e-9));

    let sx = spin_x(2);
    let ham = sx*delta_energy;
    let spin_hamiltonian = SpinHamiltonian{beta: ham.clone(),alpha: ham};

    let density_matrix = get_density_matrix(&spin_hamiltonian,&config).unwrap();

    let c = ( (0.5*ONE).exp() + (-0.5*ONE).exp())/2.0;
    let s = ( (0.5*ONE).exp() - (-0.5*ONE).exp())/2.0;
    
    expected[[0,0]] = 0.5*ONE;
    expected[[1,1]] = 0.5*ONE;
    expected[[0,1]] = -0.5*ONE*s/c;
    expected[[1,0]] = -0.5*ONE*s/c;

    assert!(approx_eq(&density_matrix, &expected, 1e-9));
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
  fn commutator(mat0: &CxMat, mat1: &CxMat) -> CxMat{
    mat0.dot(mat1) - mat1.dot(mat0)
  }
  //----------------------------------------------------------------------------
  #[test]
  #[allow(non_snake_case)]
  fn test_ClusterSpinOperators() {

    let spin_multiplicities = vec![2,3];
    let max_size = 3;
    let sops = ClusterSpinOperators::new(&spin_multiplicities,max_size).unwrap();


    for ispin_mult in &spin_multiplicities {
      let ispin_mult = *ispin_mult;
      for cluster_size in 1..=max_size{
        for op_pos in 0..cluster_size {
          let sx = sops.get(&SpinOp::Sx,ispin_mult,cluster_size,op_pos)
            .unwrap();

          let sy = sops.get(&SpinOp::Sy,ispin_mult,cluster_size,op_pos)
            .unwrap();

          let sz = sops.get(&SpinOp::Sz,ispin_mult,cluster_size,op_pos)
            .unwrap();

          let sp = sx + sy*I;
          let sm = sx - sy*I;
          let s2 = sx.dot(sx) + sy.dot(sy) + sz.dot(sz);

          assert_eq!(sx.ncols(), ispin_mult.pow(cluster_size as u32));
          assert!(check_spin_ops(sx,sy,sz,&sp,&sm,&s2));
        }
      }
    }

  }
  //----------------------------------------------------------------------------
  #[test]
  #[allow(non_snake_case)]
  fn test_KronSpinOpXYZ() {
    let spin_multiplicity = 2;
    let max_size = 2;
    let sops = KronSpinOpXYZ::new(spin_multiplicity, max_size).unwrap();

    for n_ops in 1..=max_size{
      for op_pos in 0..n_ops {
        let sx = sops.get(&SpinOp::Sx, op_pos,n_ops).unwrap();
        let sy = sops.get(&SpinOp::Sy, op_pos,n_ops).unwrap();
        let sz = sops.get(&SpinOp::Sz, op_pos,n_ops).unwrap();
        let sp = sx + sy*I;
        let sm = sx - sy*I;
        let s2 = sx.dot(sx) + sy.dot(sy) + sz.dot(sz);

        assert_eq!(sx.ncols(), spin_multiplicity.pow(n_ops as u32));
        assert!(check_spin_ops(sx,sy,sz,&sp,&sm,&s2));
      }
    }

  }
  //----------------------------------------------------------------------------
  #[test]
  #[allow(non_snake_case)]
  fn test_KronSpinOpList() {

    let spin_multiplicity = 2;
    let max_size = 2;

    let sx_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::Sx, max_size).unwrap();

    let sy_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::Sy, max_size).unwrap();

    let sz_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::Sz, max_size).unwrap();

    let sp_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::Sp, max_size).unwrap();

    let sm_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::Sm, max_size).unwrap();

    let s2_list = KronSpinOpList::new(
        spin_multiplicity, SpinOp::S2, max_size).unwrap();

    for n_ops in 1..=max_size{
      for op_pos in 0..n_ops {
        let sx = sx_list.get(op_pos,n_ops).unwrap();
        let sy = sy_list.get(op_pos,n_ops).unwrap();
        let sz = sz_list.get(op_pos,n_ops).unwrap();
        let sp = sp_list.get(op_pos,n_ops).unwrap();
        let sm = sm_list.get(op_pos,n_ops).unwrap();
        let s2 = s2_list.get(op_pos,n_ops).unwrap();

        assert_eq!(sx.ncols(), spin_multiplicity.pow(n_ops as u32));
        assert!(check_spin_ops(sx,sy,sz,sp,sm,s2));
      }
  }

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_kron_spin_op() {

    let spin_mults = Vec::<usize>::from([2,2]);

    let xx = kron_spin_op(&spin_mults,&vec![SpinOp::Sx,SpinOp::Sx]).unwrap();
    let yy = kron_spin_op(&spin_mults,&vec![SpinOp::Sy,SpinOp::Sy]).unwrap();
    let zz = kron_spin_op(&spin_mults,&vec![SpinOp::Sz,SpinOp::Sz]).unwrap();
    let pm = kron_spin_op(&spin_mults,&vec![SpinOp::Sp,SpinOp::Sm]).unwrap();
    let mp = kron_spin_op(&spin_mults,&vec![SpinOp::Sm,SpinOp::Sp]).unwrap();

    assert!(approx_eq(
          &( &(&xx + &yy) + &zz), 
          &(&zz + &(&(pm*(0.5*ONE))+&(mp*(0.5*ONE)))),1e-12)
        );
  }
//------------------------------------------------------------------------------
  #[test]
  fn test_spin_ops() {

    for spin_multiplicity in 0..14 {
      let sx = spin_x(spin_multiplicity);
      let sy = spin_y(spin_multiplicity);
      let sz = spin_z(spin_multiplicity);
      let sp = spin_plus(spin_multiplicity);
      let sm = spin_minus(spin_multiplicity);
      let s2 = spin_squared(spin_multiplicity);


      assert_eq!(sx.ncols(), sx.nrows());
      assert_eq!(sx.ncols(), spin_multiplicity);

      assert!(check_spin_ops(&sx,&sy,&sz,&sp,&sm,&s2));
      

    }
  }
  //----------------------------------------------------------------------------
  
  fn check_spin_ops(
      sx: &CxMat,
      sy: &CxMat,
      sz: &CxMat,
      sp: &CxMat,
      sm: &CxMat,
      s2: &CxMat,
      ) -> bool {

    let tol = 1e-12;

    let sx2 = sx.dot(sx);
    let sy2 = sy.dot(sy);
    let sz2 = sz.dot(sz);
    let spin2 = sx2 + sy2 + sz2;

    let mut pass: bool = true;  

    pass &= approx_eq( 
          &commutator(sx,sy), 
          &(I*sz), tol ) ; 
      
    pass &= approx_eq( 
          &commutator(sz,sx), 
          &(I*sy), tol ) ; 

    pass &= approx_eq( 
          &commutator(sy,sz), 
          &(I*sx), tol ) ; 
      
    pass &= approx_eq( 
          &sp, 
          &(sx + I*sy), tol ) ; 
      
    pass &= approx_eq( 
          &sm, 
          &(sx - I*sy), tol ) ; 
      
    pass &= approx_eq( 
          &s2, 
          &spin2, tol ) ; 
      
      
    pass 
  }
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

}



