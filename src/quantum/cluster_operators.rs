use crate::physical_constants::*;
use crate::clue_errors::*;
use crate::config::{Config, DetectedPopulation,DetectionOp,DetFrame};
use crate::quantum::pulse_sequences::{
  ideal_pulse,
  PI_OVER_2_PULSE_NAME,
  PI_PULSE_NAME
};
use crate::quantum::general_spin_hamiltonian::get_electron_cluster_thermal_density_matrix;
use crate::math::{
  anticommutator,
  commutator,
  cxmat_pow_n,
  hilbert_schmidt,
};
use crate::space_3d::{UnitSpherePoint,Vector3D};


use std::fmt;
use std::collections::HashMap;
use ndarray::{Array1,Array2,s};
use ndarray::linalg::kron;
use ndarray_linalg::UPLO;
use crate::ndarray_linalg_eigh_workaround::eigh_hermitian;
use num_complex::Complex64;


type CxMat = Array2::<Complex64>;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// 'ClusterSpinOperators' contains the spin operators used for clusters of 
/// spins of the form 'S⊗I⊗...I'.
/// Within a given 'S⊗I⊗...I', all the 'I's have the same dimension.
/// 'det_multiplicity' is the multiplicity of the detected spin:
/// 'dim(S) = det_multiplicity' (1 for CCE, which has no detected spin
/// operators).
/// 'bath_multiplicities' contains the allowed spin multiplicities of the
/// bath spins within each cluster:
/// 'dim(I) = 2s + 1', for 's' in 'bath_multiplicities'.  
/// `max_size` is the maximum number of spin operators in a product.
/// 'cluster_spin_ops' contains the matrices.
pub struct ClusterSpinOperators {
  det_multiplicity: usize,
  bath_multiplicities: Vec<usize>, 
  max_size: usize, 
  cluster_spin_ops: Vec<KronSpinOperators>,
  detection_ops: Vec<DetectionSpinOperators>,
}
impl<'a> ClusterSpinOperators {
  /// This function builds 'ClusterSpinOperators' for clusters of 
  /// `bath_multiplicities` up to size `max_size`.
  pub fn new(det_hamiltonian: (Array1::<f64>,CxMat),
      bath_multiplicities: &[usize], max_size: usize,
      config: &Config) 
   -> Result<Self, CluEError> 
  {

    let det_multiplicity = det_hamiltonian.0.len();
   
    let Some(ist_max_l) = config.max_spherical_tensor_rank else{
      return Err(CluEError::NoMaxISTRank);
    }; 

    let n_mults = bath_multiplicities.len();

    let mut cluster_spin_ops = Vec::<KronSpinOperators>::with_capacity(n_mults);

    let detection_ops = build_detection_operators(
        det_hamiltonian, 
        bath_multiplicities, max_size, config)?;

    for spin_multiplicity in bath_multiplicities.iter() {

      let sops = KronSpinOperators::new(det_multiplicity,
          *spin_multiplicity, max_size, ist_max_l)?;

      cluster_spin_ops.push(sops);
    }

    Ok(ClusterSpinOperators{
        det_multiplicity,
        bath_multiplicities: bath_multiplicities.to_owned(),
        max_size,
        cluster_spin_ops,
        detection_ops,
        })
  }

  //----------------------------------------------------------------------------
  /// This function rebuilds the detected-spin density matrix, detection
  /// operator, and pulses from `det_hamiltonian`, the eigenvalues and 
  /// eigenvectors of the detected-spin Hamiltonian at the current 
  /// orientation (`build_detected_spin_hamiltonian` of the rotated tensors).
  /// The bath spin operators do not depend on the orientation and are kept.
  /// The detected-spin multiplicity must equal the one used in `new`.
  pub fn set_detected_spin(&mut self, 
      det_hamiltonian: (Array1::<f64>,CxMat), config: &Config)
    -> Result<(),CluEError>
  {
    let det_multiplicity = det_hamiltonian.0.len();
    if det_multiplicity != self.det_multiplicity{
      return Err(CluEError::DetectedSpinMultiplicityChanged(
            self.det_multiplicity, det_multiplicity));
    }
    self.detection_ops = build_detection_operators(det_hamiltonian,
        &self.bath_multiplicities, self.max_size, config)?;
    Ok(())
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

    for (ii,&ispin_mult) in self.bath_multiplicities.iter().enumerate() {
      if ispin_mult == spin_multiplicity {
        
       let sop = self.cluster_spin_ops[ii]
         .get(sop, op_pos,cluster_size)?;
       return Ok(sop);
      }
    }
    Err(CluEError::NoSpinOpWithMultiplicity(spin_multiplicity))
  }
  //----------------------------------------------------------------------------
  pub fn get_detected_operator(&'a self,
      spin_multiplicity: usize, cluster_size: usize)
    -> Result<&'a CxMat,CluEError> 
  {
    for (ii,&ispin_mult) in self.bath_multiplicities.iter().enumerate() {
      if ispin_mult == spin_multiplicity {
        return self.detection_ops[ii].get_detection_operator(cluster_size);
      }
    }
    Err(CluEError::NoDetOpWithMultiplicity(spin_multiplicity))
  }
  //----------------------------------------------------------------------------
  pub fn get_density_matrix(&'a self,
      spin_multiplicity: usize, cluster_size: usize)
    -> Result<&'a CxMat,CluEError> 
  {
    for (ii,&ispin_mult) in self.bath_multiplicities.iter().enumerate() {
      if ispin_mult == spin_multiplicity {
        return self.detection_ops[ii].get_density_matrix(cluster_size);
      }
    }
    Err(CluEError::NoDensityMatrixWithMultiplicity(spin_multiplicity))
  }
  //----------------------------------------------------------------------------
  pub fn get_pulse_operator(&'a self, pulse_name: &str,
      spin_multiplicity: usize, cluster_size: usize)
    -> Result<&'a CxMat,CluEError> 
  {
    for (ii,&ispin_mult) in self.bath_multiplicities.iter().enumerate() {
      if ispin_mult == spin_multiplicity {
        return self.detection_ops[ii].get_pulse(pulse_name,cluster_size);
      }
    }
    Err(CluEError::NoPulseOpWithMultiplicity(pulse_name.to_string(),
          spin_multiplicity))
  }
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//------------------------------------------------------------------------------
// This function builds, for each bath multiplicity, the detected-spin
// density matrix, detection operator, and pulses, in the Zeeman basis, from
// `det_hamiltonian`, the eigenvalues (ascending) and eigenvectors V of the 
// detected-spin Hamiltonian at the current orientation.
// Level l (0 = ground state) is |v_l⟩ in both frames; `transition` and
// `detection_operator = [i, j]` name levels. Matrices and operator names
// are read in `detected_spin.frame` (see `DetFrame`). The thermal density
// matrix does not depend on the frame.
fn build_detection_operators(det_hamiltonian: (Array1::<f64>,CxMat),
    bath_multiplicities: &[usize],
    max_size: usize,
    config: &Config) 
    -> Result<Vec::<DetectionSpinOperators>,CluEError>
{

  let det_multiplicity = det_hamiltonian.0.len();

  match det_multiplicity {
    0 => return Err(CluEError::InvalidSpinMultiplicity(0)),
    1 => return Ok(Vec::<DetectionSpinOperators>::new()),
    _ => (),
  }

  // V: the eigenvectors (columns) in ascending energy order, so that level l
  // (0 = ground state) is column l.
  let eigvecs = &det_hamiltonian.1;

  { 
    // Check ordering of the eigenvalues.
    // It should not change unless ndarray changes it.
    let n = det_hamiltonian.0.len();
    let en_0 = det_hamiltonian.0[0];
    let en_n = det_hamiltonian.0[n-1];
    assert!( en_0 <= en_n);
  }

  // The eigen frame (plan_populations.md, S2-E (d)): W = V·J, the columns of
  // V in reverse order, so that matrix index i is level m-1-i and the ground
  // state is the bottom-right element. When the eigenvectors are the Zeeman
  // states of an electron (ν > 0), W is the identity, given the phases of
  // fix_eigenvector_phases().
  let eig_frame = eigvecs.slice(s![.., ..;-1]).to_owned();
  let eig_frame_dag = eig_frame.t().map(|v| v.conj());
  let to_zeeman_basis = |mat: &CxMat| -> CxMat {
    eig_frame.dot( &mat.dot( &eig_frame_dag) )
  };

  let Some(frame) = &config.detected_spin_frame else{
    return Err(CluEError::NoDetectionFrame);  
  }; 

  // Without detected_spin.density_matrix (allowed for ensemble gCCE, which
  // then uses the thermal state of the detected spin and cluster), no 
  // detected-spin density matrix is built.
  let density_matrix_opt = match &config.detected_population{
    None => None,
    Some(DetectedPopulation::Matrix(mat)) => {
      let rho = match frame{
        DetFrame::Eigen => to_zeeman_basis(mat),
        DetFrame::Zeeman => mat.clone(),  
      };
      Some(rho)
    },
    Some(DetectedPopulation::S(sop)) => {
      let mut rho = get_spin_operator(det_multiplicity,sop);
      match frame{
        DetFrame::Eigen => rho = to_zeeman_basis(&rho),
        DetFrame::Zeeman => (),  
      }
      Some(rho)
    },
    // exp(-βH_S)/Z is returned in the Zeeman basis and does not depend on
    // the basis or on the eigenvector phases, so `frame` does not apply.
    Some(DetectedPopulation::Thermal) => {
      Some(get_electron_cluster_thermal_density_matrix(
          &det_hamiltonian.0, &det_hamiltonian.1, config,
        )?)
    },   
  };

  let Some(det_op) =&config.detection_operator else {
    return Err(CluEError::NoDetectedSpinDetectionOperator);
  }; 

  let detection_operator = match det_op{
    DetectionOp::Matrix(mat) => { 
      
      match frame{
        DetFrame::Eigen => to_zeeman_basis(mat),
        DetFrame::Zeeman => mat.clone(),  
      }
    },
    DetectionOp::S(sop) => {
      let mut det_op = get_spin_operator(det_multiplicity,sop);
      match frame{
        DetFrame::Eigen => det_op = to_zeeman_basis(&det_op),
        DetFrame::Zeeman => (),  
      }
      det_op
    },
    // Levels name energy levels in both frames: |v_i⟩⟨v_j|.
    DetectionOp::Transition(level_0,level_1) => {
      let v_i = eigvecs.column(*level_0);
      let v_j = eigvecs.column(*level_1);
      let mut mat = CxMat::zeros([det_multiplicity,det_multiplicity]);
      for row in 0..det_multiplicity{
        for col in 0..det_multiplicity{
          mat[[row,col]] = v_i[row]*v_j[col].conj();
        }
      }
      mat
    }
  };

  let pulses = if config.pulses.is_empty(){
    match &config.detected_spin_transition{
      Some(transition) => {

        // Levels name energy levels in both frames. ideal_pulse() places
        // level l at matrix index m-1-l, which is column m-1-l of W, |v_l⟩.
        let u_half_pi = to_zeeman_basis(&ideal_pulse(&SpinOp::Sy, 0.5*PI,
            det_multiplicity, transition));
        let u_pi = to_zeeman_basis(&ideal_pulse(&SpinOp::Sy, PI,
            det_multiplicity, transition));

        HashMap::<String,CxMat>::from([
          (PI_OVER_2_PULSE_NAME.to_string(), u_half_pi),
          (PI_PULSE_NAME.to_string(), u_pi),
        ])},
      None => return Err(CluEError::NoDetectedSpinTransition),
    }  
  }else{
    let mut pulse_list = config.pulses.clone();
    match frame{
      DetFrame::Eigen => {
        for p in pulse_list.values_mut(){
          *p = to_zeeman_basis(p);
        }
      },
      DetFrame::Zeeman => (),  
    }
    pulse_list
  };


  let n_mults = bath_multiplicities.len();
  let mut detection_ops = Vec::<DetectionSpinOperators>::with_capacity(n_mults);

  for spin_multiplicity in bath_multiplicities.iter() {
    detection_ops.push(DetectionSpinOperators::new(
        &density_matrix_opt,
        &detection_operator,
        &pulses,    
        det_multiplicity, *spin_multiplicity, max_size,    
        )?);

  }

  Ok(detection_ops)
}

//------------------------------------------------------------------------------


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `KronSpinOperators` contains the Sx, Sy, and Sz spin operators tensored
/// with various identity matrices on either side.
pub struct KronSpinOperators {
  sx_list: KronSpinOpList,
  sy_list: KronSpinOpList,
  sz_list: KronSpinOpList,
  sp_list: KronSpinOpList,
  spherical_operators: HashMap::<(i32,i32),KronSpinOpList>,
}

impl<'a> KronSpinOperators {

  /// This function builds `KronSpinOperators` for spin with the input spin 
  /// multiplicity for clusters up to size `max_size`.
  pub fn new(det_multiplicity: usize,
      spin_multiplicity: usize, max_size: usize, ist_max_l: usize) 
    -> Result<KronSpinOperators,CluEError> 
  {

    let sx_list = KronSpinOpList::new(det_multiplicity,
        spin_multiplicity, SpinOp::Sx, max_size)?;
    let sy_list = KronSpinOpList::new(det_multiplicity,
        spin_multiplicity, SpinOp::Sy, max_size)?;
    let sz_list = KronSpinOpList::new(det_multiplicity,
        spin_multiplicity, SpinOp::Sz, max_size)?;
    let sp_list = KronSpinOpList::new(det_multiplicity,
        spin_multiplicity, SpinOp::Sp, max_size)?;
    
    let mut spherical_operators = HashMap::<(i32,i32), KronSpinOpList>::new();
    for il in 0..=ist_max_l{
      let l = il as i32;
      for im in 0..=2*l+1 {
        let m = im - l;
        let tlm_list = KronSpinOpList::new(det_multiplicity,
          spin_multiplicity, SpinOp::T(l,m), max_size)?;
        spherical_operators.insert( (l,m), tlm_list);
      }
    }
   

    Ok(KronSpinOperators{
      sx_list,  
      sy_list,  
      sz_list,  
      sp_list,  
      spherical_operators,
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
       SpinOp::Sp => self.sp_list.get(op_pos,n_ops),
       SpinOp::T(l,m) => {
         let Some(tlm_list) = self.spherical_operators.get(&(*l,*m)) else{
           return Err(CluEError::CannotFindSpinOp(sop.to_string()));
         };
           
         tlm_list.get(op_pos,n_ops)
       },
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
      det_multiplicity: usize,
      spin_multiplicity: usize,
      spin_operator: SpinOp,
      max_size: usize) -> Result<KronSpinOpList,CluEError> {

    let mut n_ops = (max_size*( max_size+1) )/2;
    if det_multiplicity > 1{ n_ops += 1;}


    let mut sop_list = Vec::<CxMat>::with_capacity(n_ops);

    for n_ops in 1..=max_size {
      let mut spin_mults = vec![spin_multiplicity; n_ops];
      if det_multiplicity > 1{
        spin_mults[0] = det_multiplicity;
      }

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
#[derive(Copy,Debug,Clone,PartialEq)]
pub enum DetOp{
  Detection,
  HalfPiPulse,  
  PiPulse,  
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `DetectionSpinOperators` contains the detection operator tensored
/// with various identity matrices on the right.
pub struct DetectionSpinOperators {
  density_matrix_list: DetectionSpinOpList,
  detection_list: DetectionSpinOpList,
  pulses_list: HashMap::<String,DetectionSpinOpList>,
}

impl<'a> DetectionSpinOperators {
  pub fn new(
    density_matrix_opt: &Option<CxMat>,
    detection_operator: &CxMat,
    pulses: &HashMap::<String,CxMat>,    
    det_multiplicity: usize, spin_multiplicity: usize, max_size: usize,    
      ) -> Result<Self,CluEError>
  {
    let density_matrix_list = match density_matrix_opt{
      Some(density_matrix) => DetectionSpinOpList::new(density_matrix,
          det_multiplicity,spin_multiplicity,max_size)?,
      None => DetectionSpinOpList::new(&CxMat::ones([1,1]),
          1,spin_multiplicity,max_size)?,
    };

    let detection_list = DetectionSpinOpList::new(detection_operator,
        det_multiplicity,spin_multiplicity,max_size)?;

    let mut pulses_list 
      = HashMap::<String,DetectionSpinOpList>::with_capacity(pulses.len());

    for (pulse_name,pulse_matrix) in pulses.iter(){
      let pl = DetectionSpinOpList::new(pulse_matrix,
        det_multiplicity,spin_multiplicity,max_size)?;
      pulses_list.insert(pulse_name.to_string(),pl);
    }

    Ok(Self{
      density_matrix_list,  
      detection_list,
      pulses_list,    
    })
  }
  //----------------------------------------------------------------------------
  pub fn get_detection_operator(&'a self, n_ops: usize) 
    -> Result<&'a CxMat,CluEError> {
    self.detection_list.get(n_ops) 
  }
  //----------------------------------------------------------------------------
  pub fn get_density_matrix(&'a self, n_ops: usize) 
    -> Result<&'a CxMat,CluEError> {
    self.density_matrix_list.get(n_ops) 
  }
  //----------------------------------------------------------------------------
  pub fn get_pulse(&'a self, pulse_name: &str, n_ops: usize)
    -> Result<&'a CxMat,CluEError> {
    let Some(pulse_list) = self.pulses_list.get(pulse_name) else{
      return Err(CluEError::UnknownPulse(pulse_name.to_string()));
    };   
   
    pulse_list.get(n_ops)
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub struct DetectionSpinOpList {
  sop_list: Vec::<CxMat>,
}
//------------------------------------------------------------------------------
impl<'a> DetectionSpinOpList {
  /// This function builds `KronSpinOpList` with the specified 
  /// spin operator, for spins with the input spin 
  /// multiplicity for clusters up to size `max_size`.
  pub fn new(
      det_operator: &CxMat,
      det_multiplicity: usize,
      spin_multiplicity: usize,
      max_size: usize) -> Result<Self,CluEError> {

    if det_multiplicity <= 1{
      return Ok(Self{sop_list: Vec::<CxMat>::new()});
    }
    let dims = det_operator.shape();
    if dims.len() != 2 || dims[0] != det_multiplicity || dims[1] != dims[0]{
      return Err(CluEError::InvalidDetectionOperator(format!(
              "expected a {}x{} matrix, but got shape {:?}",
              det_multiplicity, det_multiplicity, dims)));
    }
    if max_size == 0 {
      return Err(CluEError::ZeroMaxClusterSize)
    }
    
    let n_ops = (max_size*( max_size+1) )/2;

    let mut sop_list = Vec::<CxMat>::with_capacity(n_ops);

    for n_ops in 0..max_size {
      let spin_mults = vec![spin_multiplicity; n_ops];

      let mut sop = det_operator.clone();

      for &spin_mult in spin_mults.iter(){  
        let s = get_spin_operator(spin_mult,&SpinOp::E);
        sop = kron(&sop,&s);
      }
      sop_list.push(sop);
    }

    Ok(DetectionSpinOpList {
      sop_list,
    })
  }
  //----------------------------------------------------------------------------
  /// This function tries to find the matrix of O(1)⊗E(2)⊗...E(`n_ops`) 
  /// operators within the data structure.
  pub fn get(&'a self, n_ops: usize) ->
    Result<&'a CxMat,CluEError> {

    let idx = n_ops - 1;

    Ok(&self.sop_list[idx])
  }
  //----------------------------------------------------------------------------
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
  T(i32,i32),
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
        SpinOp::T(l,m) => write!(f, "T[{}][{}]",l,m),
      }
    }
}

impl std::str::FromStr for SpinOp{
  type Err = CluEError;
  /// This function parses a spin operator name ("Sx", "Sy", "Sz", "S+", "S-",
  /// or "S^2", with an optional lower-case "s").
  fn from_str(s: &str) -> Result<Self,CluEError>{
    // TODO: add ISTs 
    match s{
      "sx" | "Sx" => Ok(Self::Sx),  
      "sy" | "Sy" => Ok(Self::Sy),  
      "sz" | "Sz" => Ok(Self::Sz),  
      "s+" | "S+" => Ok(Self::Sp),  
      "s-" | "S-" => Ok(Self::Sm),  
      "s^2" | "S^2" => Ok(Self::S2),  
      _ => Err(CluEError::CannotParseSpinOp(s.to_string())),
    }
  }
}
//------------------------------------------------------------------------------
/// This function builds the matrix corresponding to the specified spin operator
/// and multiplicity.
pub fn get_spin_operator(spin_multiplicity: usize, sop: &SpinOp) -> CxMat{
  match sop {
    SpinOp::E => CxMat::eye(spin_multiplicity),
    SpinOp::Sx => spin_x(spin_multiplicity),
    SpinOp::Sy => spin_y(spin_multiplicity),
    SpinOp::Sz => spin_z(spin_multiplicity),
    SpinOp::Sp => spin_plus(spin_multiplicity),
    SpinOp::Sm => spin_minus(spin_multiplicity),
    SpinOp::S2 => spin_squared(spin_multiplicity),
    SpinOp::T(l,m) => spin_ist(spin_multiplicity,*l,*m),
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
//------------------------------------------------------------------------------
pub fn spin_ist(spin_multiplicity: usize, l: i32, m: i32) -> CxMat{

  assert!( l >= 0 );
  if l == 0{
    return CxMat::eye(spin_multiplicity);
  }
  if l >= spin_multiplicity as i32{
    return CxMat::zeros([spin_multiplicity,spin_multiplicity]);
  }
  let sp = spin_plus(spin_multiplicity);
  let sm = spin_minus(spin_multiplicity);
  let mut t = (-SQRT2_INV).powi(l) * ONE * cxmat_pow_n(&sp,l as usize);

  for n in 0..2*l{
    let mm = l-n;
    if mm  == m { break; }
    let denominator = ( (l*(l+1) - mm*(mm-1)) as f64 ).sqrt();
    t = commutator(&sm,&t)*(ONE/denominator);
  }

  t
}
//------------------------------------------------------------------------------
pub fn spherical_operator_decomposition(matrix: &CxMat, 
    zero_threshold: f64)
    -> Vec::<CxMat>
{

  let mut coefficients = Vec::<CxMat>::new();
  let dim  = matrix.dim().0;
  assert_eq!(matrix.dim().1, dim);

  
  for l in 0..=(dim as i32){
    
    let mut any_nonzero = false;
    let mut coefs = CxMat::zeros([2*l as usize + 1,1]);

    for (idx,m) in (-l..=l).enumerate(){
    
      let tlm = spin_ist(dim,l,m);
      let c = hilbert_schmidt(&tlm, matrix)*ONE/hilbert_schmidt(&tlm, &tlm);
      if c.norm() > zero_threshold{
        coefs[[idx,0]] = c;
        any_nonzero = true;
      }
    }
    if any_nonzero{
      coefficients.push(coefs);
    }
  }
  coefficients
}
//------------------------------------------------------------------------------
pub fn expmap_spin(multiplicity: usize, sop: &SpinOp, theta: f64) -> 
    Result<CxMat,CluEError>
{

  let s = get_spin_operator(multiplicity, sop);

  // TODO: This only works for Hermitian spinoperators.
  let Ok((eigvals, eigvecs)) = eigh_hermitian(&s, UPLO::Lower) else{
    return Err(CluEError::CannotDiagonalizeOperator(
          crate::math::cx_matrix_summary(&s)));
  };

  let inv_eigvecs = eigvecs.t().map(|v| v.conj());

  let u_eig = CxMat::from_diag(&eigvals.map(|nu|
    { 
      let i_phase: Complex64 = (I*nu)*theta;
          i_phase.exp()
    }
    )
  );

  let u = eigvecs.dot( &u_eig.dot( &inv_eigvecs) );

  Ok(u)
}
//------------------------------------------------------------------------------
pub fn expmap_spin_axis_angle(multiplicity: usize, axis: &Vector3D, angle: f64) 
    -> Result<CxMat,CluEError>
{

  let norm = axis.norm();
  if norm < 1e-12{
    return Err(CluEError::CannotNormalizeVector);
  }
  let n = axis.scale(1.0/norm);

  let sx = get_spin_operator(multiplicity, &SpinOp::Sx);
  let sy = get_spin_operator(multiplicity, &SpinOp::Sy);
  let sz = get_spin_operator(multiplicity, &SpinOp::Sz);
  let s = n.x()*ONE*sx + n.y()*ONE*sy + n.z()*ONE*sz;

  let Ok((eigvals, eigvecs)) = eigh_hermitian(&s, UPLO::Lower) else{
    return Err(CluEError::CannotDiagonalizeOperator(
          crate::math::cx_matrix_summary(&s)));
  };

  let inv_eigvecs = eigvecs.t().map(|v| v.conj());

  let u_eig = CxMat::from_diag(&eigvals.map(|nu|
    { 
      let i_phase: Complex64 = (I*nu)*angle;
          i_phase.exp()
    }
    )
  );

  let u = eigvecs.dot( &u_eig.dot( &inv_eigvecs) );

  Ok(u)
}
//------------------------------------------------------------------------------
pub fn wigner_dir(multiplicity: usize, dir: &UnitSpherePoint) 
    -> Result<CxMat,CluEError>
{
  let theta = dir.theta();
  let phi = dir.phi();

  let uz = expmap_spin(multiplicity, &SpinOp::Sz, theta)?;
  let uy = expmap_spin(multiplicity, &SpinOp::Sy, phi)?;

  let d = uz.dot(&uy);

  Ok(d)

}
//------------------------------------------------------------------------------
pub fn wigner_euler(multiplicity: usize, angles: &[f64; 3]) 
    -> Result<CxMat,CluEError>
{
  let ua = expmap_spin(multiplicity, &SpinOp::Sz, angles[0])?;
  let ub = expmap_spin(multiplicity, &SpinOp::Sy, angles[1])?;
  let uc = expmap_spin(multiplicity, &SpinOp::Sz, angles[2])?;

  let d = ua.dot(&ub.dot(&uc));

  Ok(d)

}
//------------------------------------------------------------------------------
pub fn spin_stevens(spin_multiplicity: usize, k: i32, q: i32) 
    -> Result<CxMat,CluEError>
{
  let spin = 0.5*(spin_multiplicity as f64 - 1.0);
  let s = spin*(spin + 1.0)*ONE;
  let cp = 0.5*ONE;
  let cm = -0.5*I;

  let (c,pm) = if q >= 0{
    (cp,ONE)
  }else{
    (cm,-ONE)
  };

  let e = CxMat::eye(spin_multiplicity);
  let sz = spin_z(spin_multiplicity);
  let sp = spin_plus(spin_multiplicity);
  let sm = spin_minus(spin_multiplicity);

  let a = anticommutator;
  let pow = cxmat_pow_n;
  let okq = match (k,q.abs()) {
    (2,0) => 3.0*ONE*sz.dot(&sz) - spin_squared(spin_multiplicity),
    (2,1) => 0.5*c*a(&sz, &(sp + pm*sm)),
    (2,2) => c*(sp.dot(&sp) + pm*sm.dot(&sm) ),
    (4,0) => 35.0*ONE*pow(&sz,4) 
        - (30.0*s - 25.0)*ONE*pow(&sz,2) 
        + (3.0*s*s - 6.0*s)*e,
    (4,1) => 0.5*c*a(
        &(7.0*ONE*pow(&sz,3) - (3.0*s + ONE)*sz),
        &(sp + pm*sm)),
    (4,2) => 0.5*c*a(
        &(7.0*ONE*pow(&sz,2) - (s + 5.0*ONE)*e),
        &(pow(&sp,2) + pm*pow(&sm,2))),
    (4,3) => 0.5*c*a(&sz,&(pow(&sp,3) + pm*pow(&sm,3) )),
    (4,4) => c*(pow(&sp,4) + pm*pow(&sm,4) ),
    (6,0) => 231.0*ONE*pow(&sz,6) 
        - (315.0*s - 735.0*ONE)*pow(&sz,4)
        + (105.0*s*s - 525.0*s + 294.0*ONE)*pow(&sz,2)
        - (5.0*s*s*s - 40.0*s*s + 60.0*s)*e,
    (6,1) => 0.5*c*a(
        &(
            33.0*ONE*pow(&sz,5) 
            - (30.0*s - 15.0*ONE)*pow(&sz,3)
            + (5.0*s*s - 10.0*s + 12.0*ONE)*sz
         ),
         &(sp + pm*sm)       
        ),
    (6,2) => 0.5*c*a(
        &(
            33.0*ONE*pow(&sz,4) 
            - (18.0*s + 123.0*ONE)*pow(&sz,2)
            + (s*s + 10.0*s + 102.0*ONE)*e
         ),    
        &(pow(&sp,2) + pm*pow(&sm,2))
        ),
    (6,3) => 0.5*c*a(
        &(
            11.0*ONE*pow(&sz,3) 
            - (3.0*s + 59.0*ONE)*sz
         ),    
        &(pow(&sp,3) + pm*pow(&sm,3))
        ),
    (6,4) => 0.5*c*a(
        &(
            11.0*ONE*pow(&sz,2) 
            - (s + 38.0*ONE)*e
         ),    
        &(pow(&sp,4) + pm*pow(&sm,4))
        ),
    (6,5) => 0.5*c*a(&sz,&(pow(&sp,5) + pm*pow(&sm,5)) ),
    (6,6) => c*(pow(&sp,6) + pm*pow(&sm,6)),
    _ => return Err(CluEError::NoStevensOp(k,q)),
  };

  Ok(okq)
}  
//------------------------------------------------------------------------------
pub fn operator_form_stevens_coefficients(spin_multiplicity: usize,
    coefficients: &[CxMat]) 
    -> Result<CxMat,CluEError>
{

  let mut t = CxMat::zeros([spin_multiplicity,spin_multiplicity]);

  for (index,coefs) in coefficients.iter().enumerate(){
    let k = (coefs.dim().0 as i32 - 1)/2;
    if coefs.dim().1 != 1{
      return Err(CluEError::InvalidStevensCoefficientShape(index,
            coefs.dim().0, coefs.dim().1));
    };

    for (ii,bkq) in coefs.iter().enumerate(){
      let q = -k + ii as i32;
      let okq = spin_stevens(spin_multiplicity, k, q)?;
      t = t + *bkq*okq;
    } 
  }

  Ok(t)
}  
//------------------------------------------------------------------------------
pub fn stevens_to_spherical_coefficients(spin_multiplicity: usize,
    coefficients: &[CxMat],tol: f64) 
    -> Result<Vec::<CxMat>,CluEError>
{
  let t = operator_form_stevens_coefficients(spin_multiplicity,coefficients)?;
  let ist_coefs = spherical_operator_decomposition(&t,tol);
  Ok(ist_coefs)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



#[cfg(test)]
mod tests {
  use super::*;
  use ndarray::array;
  use ndarray_linalg::Norm;

  
  //----------------------------------------------------------------------------
  fn assemble_spherical_tensor(mult: usize, coefficients: &CxMat)
      -> CxMat
  {
  
    let l = (coefficients.dim().0 as i32 - 1)/2;
    let mut out = CxMat::eye(mult);
    for (im,&c) in coefficients.iter().enumerate(){
      let m = im as i32 - l;
      let tlm = spin_ist(mult,l,m);
      out = out + tlm*c; 
    }
    out
  }
  //----------------------------------------------------------------------------
  fn assert_hermitian(h: &CxMat, tol: f64){
    let h_dag = h.t().map(|u_ij| u_ij.conj() );
    let err = (h-h_dag).norm();
    assert!(err < tol);
  }
  //----------------------------------------------------------------------------
  fn check_spherical_coefficients(mult: usize, coefs: &CxMat){
    let h = assemble_spherical_tensor(mult, coefs);
    assert_hermitian(&h, 5e-12);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_operator_form_stevens_coefficients_shape_error(){
    let column = array![[ONE],[ZERO],[ZERO],[ZERO],[ZERO]];
    let wide = CxMat::zeros([5,2]);
    let err = operator_form_stevens_coefficients(3,&[column,wide])
      .unwrap_err();
    assert_eq!(err, CluEError::InvalidStevensCoefficientShape(1,5,2));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_stevens_to_spherical_coefficients(){

    let tol = 1e-12;
    let mult = 3;

    let stevens = vec![
      array![[ONE],[ZERO],[ZERO],[ZERO],[ZERO]],
    ];   
    let coefs = stevens_to_spherical_coefficients(mult,&stevens,tol).unwrap();
      let expected = [array![[I],[ZERO],[ZERO],[ZERO],[-I]]];   
    for (ii,c) in coefs.iter().enumerate(){
      let c0 = &expected[ii];
      assert_eq!(c.shape(),c0.shape());
      assert!( (c-c0).norm() < tol );
      check_spherical_coefficients(mult,c);
    }

    let stevens = vec![
      array![[ZERO],[ONE],[ZERO],[ZERO],[ZERO]],
    ];   
    let coefs = stevens_to_spherical_coefficients(mult,&stevens,tol).unwrap();
      let expected = [array![[ZERO],[0.5*I],[ZERO],[0.5*I],[ZERO]]];   
    for (ii,c) in coefs.iter().enumerate(){
      let c0 = &expected[ii];
      assert_eq!(c.shape(),c0.shape());
      assert!( (c-c0).norm() < tol );
      check_spherical_coefficients(mult,c);
    }
  

    let stevens = vec![
      array![[ZERO],[ZERO],[ONE],[ZERO],[ZERO]],
    ];   
    let coefs = stevens_to_spherical_coefficients(mult,&stevens,tol).unwrap();
      let expected = [array![[ZERO],[ZERO],[SQRT2*SQRT3*ONE],[ZERO],[ZERO]]];   
    for (ii,c) in coefs.iter().enumerate(){
      let c0 = &expected[ii];
      assert_eq!(c.shape(),c0.shape());
      assert!( (c-c0).norm() < tol );
      check_spherical_coefficients(mult,c);
    }
  

    let stevens = vec![
      array![[ZERO],[ZERO],[ZERO],[ONE],[ZERO]],
    ];   
    let coefs = stevens_to_spherical_coefficients(mult,&stevens,tol).unwrap();
      let expected = [array![[ZERO],[0.5*ONE],[ZERO],[-0.5*ONE],[ZERO]]];   
    for (ii,c) in coefs.iter().enumerate(){
      let c0 = &expected[ii];
      assert_eq!(c.shape(),c0.shape());
      assert!( (c-c0).norm() < tol );
      check_spherical_coefficients(mult,c);
    }
  
    let stevens = vec![
      array![[ZERO],[ZERO],[ZERO],[ZERO],[ONE]],
    ];   
    let coefs = stevens_to_spherical_coefficients(mult,&stevens,tol).unwrap();
      let expected = [array![[ONE],[ZERO],[ZERO],[ZERO],[ONE]]];   
    for (ii,c) in coefs.iter().enumerate(){
      let c0 = &expected[ii];
      assert_eq!(c.shape(),c0.shape());
      assert!( (c-c0).norm() < tol );
      check_spherical_coefficients(mult,c);
    }
  
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_stevens_operators(){
  
    let tol = 1e-12;
    let mult = 3;
    let k = 2;

    let stevens = [
      spin_stevens(mult, k, -2).unwrap(),
      spin_stevens(mult, k, -1).unwrap(),
      spin_stevens(mult, k,  0).unwrap(),
      spin_stevens(mult, k,  1).unwrap(),
      spin_stevens(mult, k,  2).unwrap(),
    ];   

    let expected = [
      I*(spin_ist(mult,k, -2) - spin_ist(mult,k,  2)),
      0.5*I*(spin_ist(mult,k, -1) + spin_ist(mult,k,  1)),
      SQRT2*SQRT3*ONE*spin_ist(mult,k,  0),
      0.5*ONE*(spin_ist(mult,k, -1) - spin_ist(mult,k,  1)),
      spin_ist(mult,k, -2) + spin_ist(mult,k,  2),
    ];

    for (ii,o1) in stevens.iter().enumerate(){
      let o0 = &expected[ii];
      assert_eq!(o1.shape(),o0.shape());
      assert!( (o1-o0).norm() < tol );
      let o2 = o1.t().map(|u_ij| u_ij.conj() );
      assert!( (o1-o2).norm() < tol );
    }  

    for k in [4,6]{
      for q in-k..=k{
        let okq = spin_stevens(mult, k, q).unwrap();
        let okq_dag = okq.t().map(|u_ij| u_ij.conj() );
        assert!( (okq-okq_dag).norm() < tol );
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_spherical_operator_decomposition(){
  
    let tol = 1e-12;

    for mult in 2..=3{
      let e = CxMat::eye(mult);
      let coefs = spherical_operator_decomposition(&e, tol);
      let expected = [array![[ONE]]];
      assert_eq!(coefs.len(),expected.len());
      for (ii,c) in coefs.iter().enumerate(){
        let c0 = &expected[ii];
        assert_eq!(c.shape(),c0.shape());
        assert!( (c-c0).norm() < tol );
      }

      let sz = spin_z(mult);
      let coefs = spherical_operator_decomposition(&sz, tol);
      let expected = [array![[ZERO],[ONE],[ZERO]]];
      assert_eq!(coefs.len(),expected.len());
      for (ii,c) in coefs.iter().enumerate(){
        let c0 = &expected[ii];
        assert_eq!(c.shape(),c0.shape());
        assert!( (c-c0).norm() < tol );
      }

      let t = spin_z(mult) + CxMat::eye(mult);
      let coefs = spherical_operator_decomposition(&t, tol);
      let expected = [
        array![[ONE]], 
        array![[ZERO],[ONE],[ZERO]], 
      ];
      assert_eq!(coefs.len(),expected.len());
      for (ii,c) in coefs.iter().enumerate(){
        let c0 = &expected[ii];
        assert_eq!(c.shape(),c0.shape());
        assert!( (c-c0).norm() < tol );
      }

      let mut t = CxMat::zeros([mult,mult]);
      let mut x = ONE;
      for l in 0..=2_i32{
        for m in -l..=l{
          let tlm = spin_ist(mult,l,m);
          t = t + x*tlm;

          x += ONE;
        }
      }
      let coefs = spherical_operator_decomposition(&t, tol);
      let mut expected = vec![
        array![[ONE]], 
        array![[2.0*ONE],[3.0*ONE],[4.0*ONE]], 
      ];
      if mult >= 3{
        expected.push(
            array![[5.0*ONE],[6.0*ONE],[7.0*ONE],[8.0*ONE],[9.0*ONE]]
        );
      }
      assert_eq!(coefs.len(),expected.len());
      for (ii,c) in coefs.iter().enumerate(){
        let c0 = &expected[ii];
        assert_eq!(c.shape(),c0.shape());
        assert!( (c-c0).norm() < tol );
      }
    }

    for mult in 2..8{
      for k in [4,6]{
        for q in-k..=k{
          let okq = spin_stevens(mult, k, q).unwrap();
          let coefs = spherical_operator_decomposition(&okq, tol);
          for c in coefs.iter(){
            check_spherical_coefficients(mult,c);
          }
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_expmap_spin(){
    let d = expmap_spin(3, &SpinOp::Sz,0.0).unwrap();
    let u = CxMat::eye(3)*ONE;
    let err = (d - u).norm();
    assert!(err < 1e-12);

    let d = expmap_spin(3, &SpinOp::Sx,2.0*PI).unwrap();
    let u = CxMat::eye(3)*ONE;
    let err = (d - u).norm();
    assert!(err < 1e-12);

    let d = expmap_spin(2, &SpinOp::Sy,2.0*PI).unwrap();
    let u = -CxMat::eye(2)*ONE;
    let err = (d - u).norm();
    assert!(err < 1e-12);

    let d = expmap_spin(2, &SpinOp::Sx,4.0*PI).unwrap();
    let u = CxMat::eye(2)*ONE;
    let err = (d - u).norm();
    assert!(err < 1e-12);

    // expmap_spin(Sy, θ) = exp(iθSy). For spin 1/2 and θ = -π this is
    // cos(π/2) - i·sin(π/2)·σy = -iσy = [[0,-1],[1,0]].
    // (Before the eigh workaround, src/ndarray_linalg_eigh_workaround.rs,
    // eigh returned the eigenvectors of Sy* = -Sy, the code computed
    // exp(-iθSy), and this test expected [[0,1],[-1,0]].)
    let d = expmap_spin(2, &SpinOp::Sy,-PI).unwrap();
    let u = array![
      [ZERO,-ONE],
      [ONE,ZERO],
    ];
    let err = (d - u).norm();
    assert!(err < 1e-12);

    // Independent check: exp(iθS) from its Taylor series, for spin 1 and
    // each of Sx, Sy, Sz, and for an axis with all three components.
    let taylor = |s: &CxMat, theta: f64| -> CxMat {
      let n = s.dim().0;
      let a = s*(I*theta);
      let mut term = CxMat::eye(n);
      let mut sum = CxMat::eye(n);
      for k in 1..60{
        term = term.dot(&a)/(k as f64);
        sum += &term;
      }
      sum
    };
    for sop in [SpinOp::Sx, SpinOp::Sy, SpinOp::Sz]{
      let d = expmap_spin(3, &sop, 0.7).unwrap();
      let reference = taylor(&get_spin_operator(3, &sop), 0.7);
      assert!((d - reference).norm() < 1e-12);
    }
    let axis = Vector3D::from([0.3, -0.5, 0.8]);
    let n = axis.scale(1.0/axis.norm());
    let s = get_spin_operator(3, &SpinOp::Sx)*(ONE*n.x()) 
      + get_spin_operator(3, &SpinOp::Sy)*(ONE*n.y())
      + get_spin_operator(3, &SpinOp::Sz)*(ONE*n.z());
    let d = expmap_spin_axis_angle(3, &axis, 1.1).unwrap();
    assert!((d - taylor(&s, 1.1)).norm() < 1e-12);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_wigner_dir(){
    let d = wigner_dir(2, &UnitSpherePoint::new(0.0,0.0)).unwrap();
    let a = CxMat::eye(2);
    assert!( (d-a).norm() < 1e-12 );

    let d = wigner_dir(2, &UnitSpherePoint::new(2.0*PI,0.0)).unwrap();
    let a = -CxMat::eye(2);
    assert!( (d-a).norm() < 1e-12 );

    let d = wigner_dir(2, &UnitSpherePoint::new(4.0*PI,0.0)).unwrap();
    let a = CxMat::eye(2);
    assert!( (d-a).norm() < 1e-12 );

    let d = wigner_dir(2, &UnitSpherePoint::new(0.0,2.0*PI)).unwrap();
    let a = -CxMat::eye(2);
    assert!( (d-a).norm() < 1e-12 );

    let d = wigner_dir(2, &UnitSpherePoint::new(0.0,4.0*PI)).unwrap();
    let a = CxMat::eye(2);
    assert!( (d-a).norm() < 1e-12 );

    let d = wigner_dir(2, &UnitSpherePoint::new(2.0*PI,2.0*PI)).unwrap();
    let a = CxMat::eye(2);
    assert!( (d-a).norm() < 1e-12 );
  }
  //----------------------------------------------------------------------------
  #[test]
  #[allow(non_snake_case)]
  fn test_ClusterSpinOperators() {

    let mut config = Config::new();
    config.set_defaults().unwrap();

    let spin_multiplicities = vec![2,3];
    let max_size = 3;
    let det_h = (Array1::<f64>::zeros(1) ,CxMat::eye(1));
    let sops = ClusterSpinOperators::new(det_h,
        &spin_multiplicities,max_size,&config).unwrap();


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
  fn test_KronSpinOperators() {
    let spin_multiplicity = 2;
    let max_size = 2;
    let sops = KronSpinOperators::new(1,spin_multiplicity, max_size,3).unwrap();

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
  fn test1_KronSpinOpList() {

    let spin_multiplicity = 2;
    let max_size = 2;

    let sx_list = KronSpinOpList::new(1,
        spin_multiplicity, SpinOp::Sx, max_size).unwrap();

    let sy_list = KronSpinOpList::new(1,
        spin_multiplicity, SpinOp::Sy, max_size).unwrap();

    let sz_list = KronSpinOpList::new(1,
        spin_multiplicity, SpinOp::Sz, max_size).unwrap();

    let sp_list = KronSpinOpList::new(1,
        spin_multiplicity, SpinOp::Sp, max_size).unwrap();

    let sm_list = KronSpinOpList::new(1,
        spin_multiplicity, SpinOp::Sm, max_size).unwrap();

    let s2_list = KronSpinOpList::new(1,
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
  //----------------------------------------------------------------------------
  #[test]
  #[allow(non_snake_case)]
  fn test2_KronSpinOpList() {

    let det_multiplicity = 3;
    let spin_multiplicity = 2;
    let max_size = 3;

    let sx_list = KronSpinOpList::new(det_multiplicity,
        spin_multiplicity, SpinOp::Sx, max_size).unwrap();

    let sy_list = KronSpinOpList::new(det_multiplicity,
        spin_multiplicity, SpinOp::Sy, max_size).unwrap();

    let sz_list = KronSpinOpList::new(det_multiplicity,
        spin_multiplicity, SpinOp::Sz, max_size).unwrap();

    let sp_list = KronSpinOpList::new(det_multiplicity,
        spin_multiplicity, SpinOp::Sp, max_size).unwrap();

    let sm_list = KronSpinOpList::new(det_multiplicity,
        spin_multiplicity, SpinOp::Sm, max_size).unwrap();

    let s2_list = KronSpinOpList::new(det_multiplicity,
        spin_multiplicity, SpinOp::S2, max_size).unwrap();

    for n_ops in 1..=max_size{
      for op_pos in 0..n_ops {
        let sx = sx_list.get(op_pos,n_ops).unwrap();
        let sy = sy_list.get(op_pos,n_ops).unwrap();
        let sz = sz_list.get(op_pos,n_ops).unwrap();
        let sp = sp_list.get(op_pos,n_ops).unwrap();
        let sm = sm_list.get(op_pos,n_ops).unwrap();
        let s2 = s2_list.get(op_pos,n_ops).unwrap();

        assert!(check_spin_ops(sx,sy,sz,sp,sm,s2));

        assert_eq!(sx.ncols(), 
            det_multiplicity*spin_multiplicity.pow( (n_ops - 1) as u32));
      }
  }

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_kron_spin_op() {

    let spin_mults = Vec::<usize>::from([3,2]);

    let xx = kron_spin_op(&spin_mults,&[SpinOp::Sx,SpinOp::Sx]).unwrap();
    let yy = kron_spin_op(&spin_mults,&[SpinOp::Sy,SpinOp::Sy]).unwrap();
    let zz = kron_spin_op(&spin_mults,&[SpinOp::Sz,SpinOp::Sz]).unwrap();
    let pm = kron_spin_op(&spin_mults,&[SpinOp::Sp,SpinOp::Sm]).unwrap();
    let mp = kron_spin_op(&spin_mults,&[SpinOp::Sm,SpinOp::Sp]).unwrap();

    assert!(approx_eq(
          &( &(&xx + &yy) + &zz), 
          &(&zz + &(&(pm*(0.5*ONE))+&(mp*(0.5*ONE)))),1e-12)
        );
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_spin_ist(){
    let tol = 1e-12;

    for mult in 2..=3{
      let sz = spin_z(mult);
      let sp = spin_plus(mult);
      let sm = spin_minus(mult);

      let t = spin_ist(mult,1,0);
      let err = t - sz.clone();
      assert!( hilbert_schmidt(&err,&err).norm() < tol );

      let t = spin_ist(mult,1,-1);
      let err = t - SQRT2_INV*ONE*sm.clone();
      assert!( hilbert_schmidt(&err,&err).norm() < tol );

      let t = spin_ist(mult,1,1);
      let err = t - -SQRT2_INV*ONE*sp.clone();
      assert!( hilbert_schmidt(&err,&err).norm() < tol );

      let t = spin_ist(mult,2,0);
      let a = SQRT2*SQRT3_INV*ONE*(
          sz.dot(&sz) - 0.25*ONE*(sm.dot(&sp) + sp.dot(&sm)));
      let err = t - a;
      let eta = hilbert_schmidt(&err,&err).norm();
      assert!( eta < tol );

      let t = spin_ist(mult,2,1);
      let a = -0.5*ONE*(sz.dot(&sp) + sp.dot(&sz));
      let err = t - a;
      let eta = hilbert_schmidt(&err,&err).norm();
      assert!( eta < tol );

      let t = spin_ist(mult,2,-1);
      let a = 0.5*ONE*(sz.dot(&sm) + sm.dot(&sz));
      let err = t - a;
      let eta = hilbert_schmidt(&err,&err).norm();
      assert!( eta < tol );

      let t = spin_ist(mult,2,2);
      let a = 0.5*ONE*sp.dot(&sp);
      let err = t - a;
      let eta = hilbert_schmidt(&err,&err).norm();
      assert!( eta < tol );

      let t = spin_ist(mult,2,-2);
      let a = 0.5*ONE*sm.dot(&sm);
      let err = t - a;
      let eta = hilbert_schmidt(&err,&err).norm();
      assert!( eta < tol );

    }
    
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
          sp, 
          &(sx + I*sy), tol ) ; 
      
    pass &= approx_eq( 
          sm, 
          &(sx - I*sy), tol ) ; 
      
    pass &= approx_eq( 
          s2, 
          &spin2, tol ) ; 
      
      
    pass 
  }
  //----------------------------------------------------------------------------
  /*
  fn test_ideal_pulse(){
    let ref_pi_over_2 = ndarray::array![
      [SQRT2_INV,SQRT2_INV],
      [-SQRT2_INV,SQRT2_INV],
    ];
    let pi_over_2 = ideal_pulse(&SpinOp::Sy, 0.5*PI,2,&[0,1]).unwrap();
    let err = (pi_over_2 - ref_pi_over_2).norm();
    assert!(err< 1e-12);
  }
  */
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

  //----------------------------------------------------------------------------
  // Stage 2 of plan_populations.md: detected-spin operators per orientation,
  // the thermal state, the eigenvector phases, and the eigen-frame level
  // convention (S2-E (d)).
  //----------------------------------------------------------------------------
  use crate::quantum::tensors::{
    HamiltonianTensors, Spin1Tensors, Spin2Tensors, SpinSphericalTensors,
  };
  use crate::quantum::general_spin_hamiltonian::{
    build_bath_spin_hamiltonian,
    build_detected_spin_hamiltonian,
    get_gcce_product_density_matrix,
  };
  use crate::math::fix_eigenvector_phases;
  use crate::space_3d::SymmetricTensor3D;

  // Detected-spin gyromagnetic ratio of a free electron (rad/s/T).
  const GAMMA_E: f64 = -1.7609e11;
  //----------------------------------------------------------------------------
  // Tensors for a detected spin of multiplicity `m` (index 0) and one 
  // spin-1/2 bath nucleus (index 1) with hyperfine coupling `a_hf` (Hz, 
  // isotropic) to it. The detected spin has the Zeeman vector `nu` (Hz),
  // the Cartesian ZFS `zfs`, and the gamma matrix `gamma` (used by 
  // rotate_pasive).
  fn detected_spin_tensors(m: usize, nu: Vector3D, 
      zfs: Option<SymmetricTensor3D>, gamma: SymmetricTensor3D, a_hf: f64)
    -> HamiltonianTensors
  {
    let mut spin1_tensors = Spin1Tensors::new(2);
    spin1_tensors.set(0, nu);
    spin1_tensors.set(1, Vector3D::from([0.0, 0.0, -5.1e7]));
    let mut spin2_tensors = Spin2Tensors::new(2);
    if let Some(d) = zfs{
      spin2_tensors.set(0, 0, d);
    }
    spin2_tensors.set(0, 1, SymmetricTensor3D::from([a_hf, 0.0, 0.0,
                                                           a_hf, 0.0,
                                                                 a_hf]));
    HamiltonianTensors{
      spin_multiplicities: vec![m, 2],
      // No bath gamma matrices: the bath Zeeman vector is fixed.
      bath_gamma_matrices: vec![None; spin1_tensors.len()],
      spin1_tensors,
      spin2_tensors,
      spin_spherical_tensors: SpinSphericalTensors::new(),
      max_spherical_tensor_rank: 0,
      detected_gamma_matrix: gamma,
      magnetic_field: Vector3D::from([0.0, 0.0, 1.2]),
      mean_field_couplings: None,
    }
  }
  //----------------------------------------------------------------------------
  fn isotropic_gamma(g_factor: f64) -> SymmetricTensor3D{
    let gamma = GAMMA_E*g_factor/2.0023;
    SymmetricTensor3D::from([gamma, 0.0, 0.0,
                                   gamma, 0.0,
                                          gamma])
  }
  //----------------------------------------------------------------------------
  // A gCCE config with the given detected-spin settings.
  fn detected_spin_config(frame: DetFrame, 
      population: Option<DetectedPopulation>, det_op: DetectionOp,
      transition: Option<[usize;2]>, temperature: Option<f64>) -> Config
  {
    let mut config = Config::new();
    config.set_defaults().unwrap();
    config.detected_spin_frame = Some(frame);
    config.detected_population = population;
    config.detection_operator = Some(det_op);
    config.detected_spin_transition = transition;
    config.temperature = temperature;
    config
  }
  //----------------------------------------------------------------------------
  // The detected-spin objects in a one-spin cluster: (ρ_S, detection 
  // operator, π/2 pulse, π pulse). ρ_S is `None` without a density matrix.
  struct DetectedObjects{
    rho: Option<CxMat>,
    detection: CxMat,
    half_pi: CxMat,
    pi: CxMat,
  }
  fn detected_objects(spin_ops: &ClusterSpinOperators, config: &Config)
    -> DetectedObjects
  {
    let rho = config.detected_population.as_ref().map(|_| 
        spin_ops.get_density_matrix(2,1).unwrap().clone());
    DetectedObjects{
      rho,
      detection: spin_ops.get_detected_operator(2,1).unwrap().clone(),
      half_pi: spin_ops.get_pulse_operator(PI_OVER_2_PULSE_NAME,2,1)
        .unwrap().clone(),
      pi: spin_ops.get_pulse_operator(PI_PULSE_NAME,2,1).unwrap().clone(),
    }
  }
  fn objects_from(det_h: (Array1::<f64>,CxMat), config: &Config) 
    -> DetectedObjects
  {
    let spin_ops = ClusterSpinOperators::new(det_h, &[2], 1, config)
      .unwrap();
    detected_objects(&spin_ops, config)
  }
  //----------------------------------------------------------------------------
  fn dagger(m: &CxMat) -> CxMat{
    m.t().map(|z| z.conj())
  }
  //----------------------------------------------------------------------------
  fn assert_close(a: &CxMat, b: &CxMat, tol: f64, what: &str){
    let err = (a - b).norm();
    assert!(err <= tol, "{}: |a - b| = {:e}\na = {}\nb = {}", what, err, a, b);
  }
  //----------------------------------------------------------------------------
  // exp(a) by scaling and squaring with a Taylor series (no eigh).
  fn expm_taylor(a: &CxMat) -> CxMat{
    let k = 12;
    let scaled = a/Complex64::new(2f64.powi(k), 0.0);
    let n = a.dim().0;
    let mut term = CxMat::eye(n);
    let mut sum = CxMat::eye(n);
    for j in 1..30{
      term = term.dot(&scaled)/Complex64::new(j as f64, 0.0);
      sum += &term;
    }
    cxmat_pow_n(&sum, 1 << k)
  }
  //----------------------------------------------------------------------------
  // H = ν·S + S·D·S, built here from the spin matrices.
  fn single_spin_matrix(m: usize, nu: &Vector3D, 
      zfs: &Option<SymmetricTensor3D>) -> CxMat
  {
    let s = [spin_x(m), spin_y(m), spin_z(m)];
    let mut h = &s[0]*nu.x() + &s[1]*nu.y() + &s[2]*nu.z();
    if let Some(d) = zfs{
      let dm = [[d.xx(), d.xy(), d.xz()], [d.yx(), d.yy(), d.yz()],
        [d.zx(), d.zy(), d.zz()]];
      for (i, si) in s.iter().enumerate(){
        for (j, sj) in s.iter().enumerate(){
          h = h + si.dot(sj)*dm[i][j];
        }
      }
    }
    h
  }
  //----------------------------------------------------------------------------
  // S2-E (d) condition 3 and test (i): for an electron whose eigenstates
  // are its Zeeman states, W = 1 and every input gives the same operator 
  // in both frames. Test (iv): the transition pulses and the [i, j] 
  // detection operator equal those of the code before Stage 2 
  // (eigen frame), given the same V.
  #[test]
  fn test_detected_spin_frames_coincide(){
    for m in [2, 3]{
      let nu = Vector3D::from([0.0, 0.0, 3.363e10]);
      let tensors = detected_spin_tensors(m, nu, None, isotropic_gamma(2.0023),
          0.0);
      let det_h = build_detected_spin_hamiltonian(&tensors).unwrap();

      // V = J exactly, so W = V·J = 1.
      let v = &det_h.1;
      for i in 0..m{
        for j in 0..m{
          let expected = if i + j == m - 1 { ONE } else { ZERO };
          assert_eq!(v[[i,j]], expected);
        }
      }

      let mut off_diagonal = CxMat::zeros((m,m));
      off_diagonal[[0,1]] = Complex64::new(0.3, -0.2);
      off_diagonal[[1,0]] = Complex64::new(0.3, 0.2);
      off_diagonal[[m-1,m-1]] = ONE;
      let mut ground = CxMat::zeros((m,m));
      ground[[m-1,m-1]] = ONE;
      let mut user_pulse = CxMat::eye(m);
      user_pulse[[0,0]] = I;

      let populations = [
        DetectedPopulation::Matrix(ground.clone()),
        DetectedPopulation::Matrix(off_diagonal.clone()),
        DetectedPopulation::S(SpinOp::Sz),
        DetectedPopulation::S(SpinOp::Sp),
      ];
      let detection_ops = [
        DetectionOp::Matrix(off_diagonal.clone()),
        DetectionOp::S(SpinOp::Sp),
        DetectionOp::S(SpinOp::Sz),
        DetectionOp::Transition(0, 1),
      ];
      for pop in populations.iter(){
        for det_op in detection_ops.iter(){
          let objects = |frame: DetFrame, pulses: bool| {
            let mut config = detected_spin_config(frame, Some(pop.clone()),
                det_op.clone(), Some([0, 1]), None);
            if pulses{
              config.pulses = HashMap::from([
                (PI_OVER_2_PULSE_NAME.to_string(), user_pulse.clone()),
                (PI_PULSE_NAME.to_string(), user_pulse.clone()),
              ]);
            }
            let spin_ops = ClusterSpinOperators::new(det_h.clone(),
                &[2], 1, &config).unwrap();
            let rho = spin_ops.get_density_matrix(2,1).unwrap().clone();
            let det = spin_ops.get_detected_operator(2,1).unwrap().clone();
            let p1 = spin_ops.get_pulse_operator(PI_OVER_2_PULSE_NAME,2,1)
              .unwrap().clone();
            let p2 = spin_ops.get_pulse_operator(PI_PULSE_NAME,2,1)
              .unwrap().clone();
            (rho, det, p1, p2)
          };
          for pulses in [false, true]{
            let e = objects(DetFrame::Eigen, pulses);
            let z = objects(DetFrame::Zeeman, pulses);
            assert_close(&e.0, &z.0, 1e-12, "density matrix");
            assert_close(&e.1, &z.1, 1e-12, "detection operator");
            assert_close(&e.2, &z.2, 1e-12, "pi/2 pulse");
            assert_close(&e.3, &z.3, 1e-12, "pi pulse");
          }
        }
      }

      // In both frames, "Sz" is +S_z and the bottom-right matrix element is
      // the ground state m_S = -S.
      let config = detected_spin_config(DetFrame::Eigen, 
          Some(DetectedPopulation::S(SpinOp::Sz)), DetectionOp::S(SpinOp::Sp),
          Some([0, 1]), None);
      let spin_ops = ClusterSpinOperators::new(det_h.clone(), &[2], 1, 
          &config).unwrap();
      assert_close(spin_ops.get_density_matrix(2,1).unwrap(), &spin_z(m), 
          1e-15, "Sz");
      assert_close(spin_ops.get_detected_operator(2,1).unwrap(), 
          &spin_plus(m), 1e-15, "Sp");

      // (iv) The pulses and [i, j] of the code before Stage 2 
      // (eigen frame, same V).
      let config = detected_spin_config(DetFrame::Eigen, 
          Some(DetectedPopulation::Matrix(ground.clone())),
          DetectionOp::Transition(0, 1), Some([0, 1]), None);
      let new = objects_from(det_h.clone(), &config);
      let (row, col) = (m - 1, m - 2);
      let v_dag = dagger(v);
      let old_half_pi = v.dot(&ideal_pulse(&SpinOp::Sy, 0.5*PI, m, 
            &[row, col])).dot(&v_dag);
      let old_pi = v.dot(&ideal_pulse(&SpinOp::Sy, PI, m, &[row, col]))
        .dot(&v_dag);
      let mut old_det = CxMat::zeros((m,m));
      old_det[[0,1]] = ONE;
      let old_det = v.dot(&old_det).dot(&v_dag);
      assert_close(&new.half_pi, &old_half_pi, 1e-12, "old pi/2 pulse");
      assert_close(&new.pi, &old_pi, 1e-12, "old pi pulse");
      assert_close(&new.detection, &old_det, 1e-12, "old [0, 1]");
    }
  }
  //----------------------------------------------------------------------------
  // S2-B and S2-5 (b): the thermal ρ_S has ⟨S_z⟩ = -tanh(hν/2k_BT)/2 for
  // an isotropic spin-1/2 in both frames, and for a spin-1 with a tilted
  // ZFS it equals exp(-βH)/Z (Taylor series, no eigh) in both frames.
  #[test]
  fn test_detected_spin_thermal_state(){
    let temperature = 1.0;
    let beta = 1.0/(JOULES_TO_HERTZ*BOLTZMANN*temperature);

    // Spin-1/2.
    let nu_z = 3.363e10;
    let tensors = detected_spin_tensors(2, Vector3D::from([0.0, 0.0, nu_z]),
        None, isotropic_gamma(2.0023), 0.0);
    let det_h = build_detected_spin_hamiltonian(&tensors).unwrap();
    let sz_ref = -0.5*(0.5*beta*nu_z).tanh();
    assert!(sz_ref < -0.3);
    for frame in [DetFrame::Eigen, DetFrame::Zeeman]{
      let config = detected_spin_config(frame, 
          Some(DetectedPopulation::Thermal), DetectionOp::Transition(0,1),
          Some([0,1]), Some(temperature));
      let rho = objects_from(det_h.clone(), &config).rho.unwrap();
      let sz = rho.dot(&spin_z(2)).diag().sum();
      assert!((sz.re - sz_ref).abs() < 1e-12, "{} vs {}", sz, sz_ref);
      assert!(sz.im.abs() < 1e-15);
    }

    // Spin-1, ZFS comparable to the Zeeman energy and tilted from the field.
    let nu = Vector3D::from([1.0e9, -2.0e8, 9.0e9]);
    let (d, theta): (f64, f64) = (6.0e9, 0.6);
    let (st, ct) = (theta.sin(), theta.cos());
    let zfs = Some(SymmetricTensor3D::from([
        d*(st*st - 1.0/3.0), 0.0, d*st*ct,
                        -d/3.0, 0.0,
                                d*(ct*ct - 1.0/3.0)]));
    let h = single_spin_matrix(3, &nu, &zfs);
    let tensors = detected_spin_tensors(3, nu, zfs, isotropic_gamma(2.0023),
        0.0);
    let det_h = build_detected_spin_hamiltonian(&tensors).unwrap();
    let g = expm_taylor(&(&h*Complex64::new(-beta, 0.0)));
    let rho_ref = &g/g.diag().sum();
    let mut rhos = Vec::<CxMat>::new();
    for frame in [DetFrame::Eigen, DetFrame::Zeeman]{
      let config = detected_spin_config(frame, 
          Some(DetectedPopulation::Thermal), DetectionOp::Transition(0,1),
          Some([0,1]), Some(temperature));
      let spin_ops = ClusterSpinOperators::new(det_h.clone(), &[2], 1, 
          &config).unwrap();
      let rho = spin_ops.get_density_matrix(2,1).unwrap().clone();
      assert_close(&rho, &rho_ref, 1e-12, "thermal spin-1");
      assert_close(&commutator(&rho, &h), &CxMat::zeros((3,3)), 
          1e-12*h.norm(), "[ρ, H]");
      rhos.push(rho);
    }
    assert_close(&rhos[0], &rhos[1], 1e-14, "frame invariance");
  }
  //----------------------------------------------------------------------------
  // S2-5 (c): arbitrary phases on the eigenvectors, removed by
  // fix_eigenvector_phases, give the same ρ_S, detection operator and 
  // pulses. The thermal ρ_S does not depend on the phases at all.
  #[test]
  fn test_detected_spin_gauge(){
    let nu = Vector3D::from([2.0e9, 1.5e9, 9.0e9]);
    let (d, theta): (f64, f64) = (4.0e9, 0.9);
    let (st, ct) = (theta.sin(), theta.cos());
    let zfs = Some(SymmetricTensor3D::from([
        d*(st*st - 1.0/3.0), 0.0, d*st*ct,
                        -d/3.0, 0.0,
                                d*(ct*ct - 1.0/3.0)]));
    let tensors = detected_spin_tensors(3, nu, zfs, isotropic_gamma(2.0023),
        0.0);
    let det_h = build_detected_spin_hamiltonian(&tensors).unwrap();

    let mut phased = det_h.1.clone();
    for (k, phi) in [0.4, -1.9, 2.7].iter().enumerate(){
      let p = Complex64::from_polar(1.0, *phi);
      phased.column_mut(k).mapv_inplace(|z| z*p);
    }
    let mut regauged = phased.clone();
    fix_eigenvector_phases(&mut regauged);

    let mut mat = CxMat::zeros((3,3));
    mat[[0,1]] = Complex64::new(0.2, 0.1);
    mat[[1,0]] = Complex64::new(0.2, -0.1);
    mat[[2,2]] = ONE;
    let config = detected_spin_config(DetFrame::Eigen,
        Some(DetectedPopulation::Matrix(mat.clone())), 
        DetectionOp::Transition(0, 2), Some([0, 1]), Some(1.0));
    let a = objects_from(det_h.clone(), &config);
    let b = objects_from((det_h.0.clone(), regauged), &config);
    assert_close(a.rho.as_ref().unwrap(), b.rho.as_ref().unwrap(), 1e-15, 
        "density matrix");
    assert_close(&a.detection, &b.detection, 1e-15, "detection");
    assert_close(&a.half_pi, &b.half_pi, 1e-15, "pi/2");
    assert_close(&a.pi, &b.pi, 1e-15, "pi");

    // Without the gauge, the eigen-frame objects change.
    let c = objects_from((det_h.0.clone(), phased.clone()), &config);
    assert!((&a.detection - &c.detection).norm() > 0.1);
    assert!((a.rho.as_ref().unwrap() - c.rho.as_ref().unwrap()).norm() > 0.01);

    // The thermal state is the same for any phases.
    let config = detected_spin_config(DetFrame::Eigen,
        Some(DetectedPopulation::Thermal), DetectionOp::Transition(0, 2),
        Some([0, 1]), Some(1.0));
    let a = objects_from(det_h.clone(), &config);
    let c = objects_from((det_h.0.clone(), phased), &config);
    assert_close(a.rho.as_ref().unwrap(), c.rho.as_ref().unwrap(), 1e-15,
        "thermal");
  }
  //----------------------------------------------------------------------------
  // S2-5 (a) and S2-E test (ii), and S2-5 (f): with an anisotropic g, 
  // after set_detected_spin at a tilted orientation, the eigen-frame 
  // `[[0,0],[0,1]]` is the ground state (1 - n·σ)/2 of the rotated
  // H_S = ν·S, n = ν/|ν|, and the ensemble-gCCE product state built from the
  // stored ρ_S has Tr_C ρ = ρ_S.
  #[test]
  fn test_set_detected_spin_follows_orientation(){
    let gamma = SymmetricTensor3D::from([GAMMA_E*0.97, 0.01*GAMMA_E, 0.0,
                                                   GAMMA_E*1.0, 0.0,
                                                          GAMMA_E*1.06]);
    let mut tensors = detected_spin_tensors(2, Vector3D::zeros(), None, 
        gamma, 2.0e6);
    tensors.rotate_pasive(&UnitSpherePoint::from(
          Vector3D::from([0.0, 0.0, 1.0]))).unwrap();

    let mut ground = CxMat::zeros((2,2));
    ground[[1,1]] = ONE;
    let mut coherent = CxMat::zeros((2,2));
    coherent[[0,0]] = Complex64::new(0.3, 0.0);
    coherent[[1,1]] = Complex64::new(0.7, 0.0);
    coherent[[0,1]] = Complex64::new(0.1, 0.2);
    coherent[[1,0]] = Complex64::new(0.1, -0.2);

    let mut config_ground = detected_spin_config(DetFrame::Eigen,
        Some(DetectedPopulation::Matrix(ground.clone())), 
        DetectionOp::Transition(0, 1), Some([0, 1]), Some(10.0));
    config_ground.ensemble_cce = Some(true);
    let mut config_coherent = config_ground.clone();
    config_coherent.detected_population 
      = Some(DetectedPopulation::Matrix(coherent.clone()));

    // Built once from the unrotated tensors, as before the orientation loop.
    let det_h0 = build_detected_spin_hamiltonian(&tensors).unwrap();
    let mut spin_ops = ClusterSpinOperators::new(det_h0, &[2], 2, 
        &config_ground).unwrap();
    let rho_0 = spin_ops.get_density_matrix(2,1).unwrap().clone();

    let sigma = [spin_x(2)*2.0, spin_y(2)*2.0, spin_z(2)*2.0];
    let dirs = [Vector3D::from([0.5, -0.3, 0.8]), 
      Vector3D::from([-0.7, 0.6, 0.2])];
    let mut previous: Option<CxMat> = None;
    for dir in dirs.iter(){
      let mut rotated = tensors.clone();
      rotated.rotate_pasive(&UnitSpherePoint::from(dir.normalize()))
        .unwrap();
      let nu = rotated.spin1_tensors.get(0).unwrap().clone();
      let n = nu.normalize();
      let expected = (CxMat::eye(2) - &sigma[0]*n.x() - &sigma[1]*n.y() 
          - &sigma[2]*n.z())*Complex64::new(0.5, 0.0);

      let det_h = build_detected_spin_hamiltonian(&rotated).unwrap();
      spin_ops.set_detected_spin(det_h.clone(), &config_ground).unwrap();
      let rho = spin_ops.get_density_matrix(2,1).unwrap().clone();
      assert_close(&rho, &expected, 1e-12, "ground state");
      // The tilt is O(Δg/g), and it is resolved.
      assert!((&rho - &rho_0).norm() > 1e-3);

      // (f) Ensemble gCCE with a density matrix: the empty cluster state 
      // is ρ_S, and a one-nucleus cluster state has Tr_C ρ = ρ_S.
      spin_ops.set_detected_spin(det_h, &config_coherent).unwrap();
      let rho_s = spin_ops.get_density_matrix(2,1).unwrap().clone();
      if let Some(prev) = &previous{
        assert!((&rho_s - prev).norm() > 1e-3);
      }
      let empty_ham = build_bath_spin_hamiltonian(&[0], &spin_ops, &rotated)
        .unwrap();
      let rho_empty = get_gcce_product_density_matrix(&rho_s, &empty_ham,
          Some(10.0)).unwrap();
      assert_close(&rho_empty, &rho_s, 1e-12, "empty cluster");
      let bath_ham = build_bath_spin_hamiltonian(&[0, 1], &spin_ops, 
          &rotated).unwrap();
      let rho_c = get_gcce_product_density_matrix(&rho_s, &bath_ham,
          Some(10.0)).unwrap();
      let mut tr_c = CxMat::zeros((2,2));
      for i in 0..2{
        for j in 0..2{
          tr_c[[i,j]] = rho_c[[2*i,2*j]] + rho_c[[2*i+1,2*j+1]];
        }
      }
      assert_close(&tr_c, &rho_s, 1e-12, "Tr_C");
      previous = Some(rho_s);
    }
  }
  //----------------------------------------------------------------------------
  // S2-5 (d) and S2-E test (iii): a spin-1 with a ZFS that reorders the 
  // levels. `transition = [0, 1]` acts on the two lowest levels, in both 
  // frames and at a tilted orientation, and `detection_operator = [0, 1]` 
  // is |v_0⟩⟨v_1| in the Zeeman frame too.
  #[test]
  fn test_detected_spin_levels_with_zfs(){
    // D along z: E(m=0) = -2D/3 < E(m=-1) = D/3 - ν < E(m=+1) = D/3 + ν.
    let (nu_z, d) = (2.0e9, 5.0e9);
    let zfs_z = SymmetricTensor3D::from([-d/3.0, 0.0, 0.0,
                                                -d/3.0, 0.0,
                                                   2.0*d/3.0]);
    let gamma = isotropic_gamma(2.0023);
    let tensors = detected_spin_tensors(3, Vector3D::from([0.0, 0.0, nu_z]),
        Some(zfs_z.clone()), gamma.clone(), 0.0);
    let det_h = build_detected_spin_hamiltonian(&tensors).unwrap();
    // Zeeman indices: m = +1 → 0, m = 0 → 1, m = -1 → 2.
    let mut expected_det = CxMat::zeros((3,3));
    expected_det[[1,2]] = ONE;
    for frame in [DetFrame::Eigen, DetFrame::Zeeman]{
      let config = detected_spin_config(frame.clone(), None, 
          DetectionOp::Transition(0, 1), Some([0, 1]), None);
      let o = objects_from(det_h.clone(), &config);
      assert_close(&o.detection, &expected_det, 1e-15, "[0, 1]");
      // The π pulse swaps m = 0 and m = -1 and leaves m = +1 alone.
      assert!((o.pi[[0,0]] - ONE).norm() < 1e-15);
      assert!((o.pi[[2,1]].norm() - 1.0).abs() < 1e-15);
      assert!((o.pi[[1,2]].norm() - 1.0).abs() < 1e-15);
    }

    // Tilted orientations: the levels come from the rotated H_S, checked
    // against H_S built here from the rotated tensors.
    for dir in [Vector3D::from([0.4, 0.2, 0.9]), 
      Vector3D::from([0.9, -0.4, 0.1])]
    {
      let mut rotated = tensors.clone();
      rotated.rotate_pasive(&UnitSpherePoint::from(dir.normalize()))
        .unwrap();
      let h = single_spin_matrix(3, rotated.spin1_tensors.get(0).unwrap(),
          &rotated.spin2_tensors.get(0,0).cloned());
      let (e, v) = build_detected_spin_hamiltonian(&rotated).unwrap();
      for k in 0..3{
        let vk = v.column(k).to_owned();
        let hv = h.dot(&vk);
        assert!((&hv - &vk*Complex64::new(e[k], 0.0)).norm() 
            <= 1e-12*h.norm());
      }
      for frame in [DetFrame::Eigen, DetFrame::Zeeman]{
        let config = detected_spin_config(frame, None, 
            DetectionOp::Transition(0, 1), Some([0, 1]), None);
        let mut spin_ops = ClusterSpinOperators::new(det_h.clone(), &[2], 1,
            &config).unwrap();
        spin_ops.set_detected_spin((e.clone(), v.clone()), &config)
          .unwrap();
        let o = detected_objects(&spin_ops, &config);
        let in_level_basis = |u: &CxMat| dagger(&v).dot(u).dot(&v);
        for u in [&o.half_pi, &o.pi]{
          let u_l = in_level_basis(u);
          // Level 2 is untouched, and span{v_0, v_1} is invariant.
          assert!((u_l[[2,2]] - ONE).norm() < 1e-12);
          for k in 0..2{
            assert!(u_l[[2,k]].norm() < 1e-12);
            assert!(u_l[[k,2]].norm() < 1e-12);
          }
        }
        // The π pulse takes level 0 to level 1.
        assert!((in_level_basis(&o.pi)[[1,0]].norm() - 1.0).abs() < 1e-12);
        let mut v0v1 = CxMat::zeros((3,3));
        for r in 0..3{
          for c in 0..3{
            v0v1[[r,c]] = v[[r,0]]*v[[c,1]].conj();
          }
        }
        assert_close(&o.detection, &v0v1, 1e-15, "[0, 1] tilted");
      }
    }
  }
  //----------------------------------------------------------------------------
  // S2-5 (g): set_detected_spin rejects a Hamiltonian of another size.
  #[test]
  fn test_set_detected_spin_wrong_multiplicity(){
    let config = detected_spin_config(DetFrame::Eigen, None,
        DetectionOp::Transition(0, 1), Some([0, 1]), None);
    let t2 = detected_spin_tensors(2, Vector3D::from([0.0, 0.0, 3.0e10]), 
        None, isotropic_gamma(2.0023), 0.0);
    let t3 = detected_spin_tensors(3, Vector3D::from([0.0, 0.0, 3.0e10]), 
        None, isotropic_gamma(2.0023), 0.0);
    let mut spin_ops = ClusterSpinOperators::new(
        build_detected_spin_hamiltonian(&t2).unwrap(), &[2], 1, &config)
      .unwrap();
    let err = spin_ops.set_detected_spin(
        build_detected_spin_hamiltonian(&t3).unwrap(), &config).unwrap_err();
    assert_eq!(err, CluEError::DetectedSpinMultiplicityChanged(2, 3));
    assert!(spin_ops.set_detected_spin(
        build_detected_spin_hamiltonian(&t2).unwrap(), &config).is_ok());
  }
  //----------------------------------------------------------------------------
  // S2-5 (h) and S2-D: mean fields do not enter the detected-spin 
  // Hamiltonian, so ρ_S, the detection operator and the pulses are the 
  // same with and without them.
  #[test]
  fn test_detected_spin_ignores_mean_fields(){
    let gamma = SymmetricTensor3D::from([GAMMA_E*0.97, 0.0, 0.0,
                                                   GAMMA_E*1.0, 0.0,
                                                          GAMMA_E*1.06]);
    let mut tensors = detected_spin_tensors(2, Vector3D::zeros(), None, 
        gamma, 5.0e7);
    tensors.rotate_pasive(&UnitSpherePoint::from(
          Vector3D::from([0.3, 0.4, 0.8]).normalize())).unwrap();
    let det_h = build_detected_spin_hamiltonian(&tensors).unwrap();

    let mut with_mf = tensors.clone();
    let up = ndarray::array![[ONE],[ZERO]];
    with_mf.set_mean_field_couplings(&[up.clone(), up]).unwrap();
    // The Overhauser field on the detected spin is not zero.
    let overhauser = with_mf.get_mean_field_couplings(0, &[0]).unwrap();
    assert!(overhauser.norm() > 1e7);
    let det_h_mf = build_detected_spin_hamiltonian(&with_mf).unwrap();
    assert_eq!(det_h.0, det_h_mf.0);
    assert_eq!(det_h.1, det_h_mf.1);

    let mut mat = CxMat::zeros((2,2));
    mat[[0,1]] = Complex64::new(0.2, 0.1);
    mat[[1,0]] = Complex64::new(0.2, -0.1);
    mat[[1,1]] = ONE;
    let config = detected_spin_config(DetFrame::Eigen,
        Some(DetectedPopulation::Matrix(mat)), 
        DetectionOp::S(SpinOp::Sp), Some([0, 1]), None);
    let a = objects_from(det_h, &config);
    let b = objects_from(det_h_mf, &config);
    assert_eq!(a.rho, b.rho);
    assert_eq!(a.detection, b.detection);
    assert_eq!(a.half_pi, b.half_pi);
    assert_eq!(a.pi, b.pi);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_detection_spin_op_list_errors(){
    let det_op = CxMat::eye(2);
    assert_eq!(DetectionSpinOpList::new(&det_op, 2, 2, 0).err(),
        Some(CluEError::ZeroMaxClusterSize));
    assert_eq!(DetectionSpinOpList::new(&det_op, 3, 2, 2).err(),
        Some(CluEError::InvalidDetectionOperator(
            "expected a 3x3 matrix, but got shape [2, 2]".to_string())));
    assert!(DetectionSpinOpList::new(&det_op, 2, 2, 2).is_ok());
  }
}



