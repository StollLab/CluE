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
use ndarray::{Array1,Array2};
use ndarray::linalg::kron;
use ndarray_linalg::{Eigh, UPLO};
use num_complex::Complex64;


type CxMat = Array2::<Complex64>;
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// 'ClusterSpinOperators' contains the spin operators used for clusters of 
/// spins of the form 'S⊗I⊗...I'.
/// Within a given 'S⊗I⊗...I', all the 'I's have the same dimension.
/// 'det_multiplicity' is the multiplicity of the detected spin:
/// 'dim(S) = 2×det_multiplicity + 1'.  
/// 'bath_multiplicities' contains the allowed spin multiplicities of the
/// bath spins within each cluster:
/// 'dim(I) = 2s + 1', for 's' in 'bath_multiplicities'.  
/// `max_size` is the maximum number of spin operators in a product.
/// 'cluster_spin_ops' contains the matrices.
pub struct ClusterSpinOperators {
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
        bath_multiplicities: bath_multiplicities.to_owned(),
        max_size,
        cluster_spin_ops,
        detection_ops,
        })
  }

  //----------------------------------------------------------------------------
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

  let Some(detected_population) = &config.detected_population else{
    return Err(CluEError::NoDetectedSpinDensityMatrix);
  };

  let eigvecs = &det_hamiltonian.1;
  let inv_eigvecs = eigvecs.t().map(|v| v.conj());

  { 
    // Check ordering of the eigenvalues.
    // It should not change unless ndarray changes it.
    let n = det_hamiltonian.0.len();
    let en_0 = det_hamiltonian.0[0];
    let en_n = det_hamiltonian.0[n-1];
    assert!( en_0 <= en_n);
  }

  let Some(frame) = &config.detected_spin_frame else{
    return Err(CluEError::NoDetectionFrame);  
  }; 

  let density_matrix_opt = match detected_population{
    DetectedPopulation::Matrix(mat) => {
      let rho = match frame{
        DetFrame::Eigen => eigvecs.dot( &mat.dot( &inv_eigvecs) ),
        DetFrame::Zeeman => mat.clone(),  
      };
      Some(rho)
    },
    DetectedPopulation::S(sop) => {
      let mut rho = get_spin_operator(det_multiplicity,sop);
      match frame{
        DetFrame::Eigen => rho = eigvecs.dot( &rho.dot( &inv_eigvecs) ),
        DetFrame::Zeeman => (),  
      }
      Some(rho)
    },
    DetectedPopulation::Thermal => {
      let mut rho = get_electron_cluster_thermal_density_matrix(
          &det_hamiltonian.0, &det_hamiltonian.1, config,
        )?;
      match frame{
        DetFrame::Eigen => rho = eigvecs.dot( &rho.dot( &inv_eigvecs) ),
        DetFrame::Zeeman => (),  
      }
      Some(rho)
    },   
  };

  let Some(det_op) =&config.detection_operator else {
    return Err(CluEError::NoDetectedSpinDetectionOperator);
  }; 

  let detection_operator = match det_op{
    DetectionOp::Matrix(mat) => { 
      let det_op = match frame{
        DetFrame::Eigen => eigvecs.dot( &mat.dot( &inv_eigvecs) ),
        DetFrame::Zeeman => mat.clone(),  
      };
      det_op
    },
    DetectionOp::S(sop) => {
      let mut det_op = get_spin_operator(det_multiplicity,sop);
      match frame{
        DetFrame::Eigen => det_op = eigvecs.dot( &det_op.dot( &inv_eigvecs) ),
        DetFrame::Zeeman => (),  
      }
      det_op
    },
    DetectionOp::Transition(level_0,level_1) => {
      let mut mat = CxMat::zeros([det_multiplicity,det_multiplicity]);
      let (row,col) = match frame{
        DetFrame::Eigen => {
          (*level_0,*level_1)
        },
        DetFrame::Zeeman 
          => (det_multiplicity - *level_0 - 1,det_multiplicity - *level_1 - 1),
      };
      mat[[row,col]] = ONE;

      if *frame == DetFrame::Eigen{
        mat = eigvecs.dot( &mat.dot( &inv_eigvecs) );
      }

      mat
    }
  };

  let pulses = if config.pulses.is_empty(){
    match &config.detected_spin_transition{
      Some(transition) => {

      let level_0 = &transition[0];  
      let level_1 = &transition[1];  

      let (row,col) = match frame{
        DetFrame::Eigen 
          => (det_multiplicity - *level_0 - 1,det_multiplicity - *level_1 - 1),
        DetFrame::Zeeman 
          => (*level_0,*level_1)
      };

        let mut u_half_pi = ideal_pulse(&SpinOp::Sy, 0.5*PI,
            det_multiplicity, &[row,col]);
        let mut u_pi = ideal_pulse(&SpinOp::Sy, PI,det_multiplicity, &[row,col]);

        if *frame == DetFrame::Eigen{
          u_half_pi = eigvecs.dot( &u_half_pi.dot( &inv_eigvecs) );
          u_pi = eigvecs.dot( &u_pi.dot( &inv_eigvecs) );
        }

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
        for (_,p) in pulse_list.iter_mut(){
          *p = eigvecs.dot( &p.dot( &inv_eigvecs) );
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
        let m = im as i32 - l;
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
      return Err(CluEError::InvalidPulse(pulse_name.to_string()));
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
      return Err(CluEError::InvalidDetectionOperator)
    }
    if max_size == 0 {
      return Err(CluEError::InvalidDetectionOperator)
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

impl SpinOp{
  pub fn from_str(s: &str) -> Result<Self,CluEError>{
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
  let Ok((eigvals, eigvecs)) = s.eigh(UPLO::Lower) else{
    return Err(CluEError::CannotDiagonalizeOperator(s.to_string()));
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

  let Ok((eigvals, eigvecs)) = s.eigh(UPLO::Lower) else{
    return Err(CluEError::CannotDiagonalizeOperator(s.to_string()));
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

  for coefs in coefficients.iter(){
    let k = (coefs.dim().0 as i32 - 1)/2;
    if coefs.dim().1 != 1{
      return Err(CluEError::Generic(
            "Stevens coeeficients should be a list of k×1 matrices".to_string()));
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
  fn test_stevens_to_spherical_coefficients(){
  
    let tol = 1e-12;
    let mult = 3;

    let stevens = vec![
      array![[ONE],[ZERO],[ZERO],[ZERO],[ZERO]],
    ];   
    let coefs = stevens_to_spherical_coefficients(mult,&stevens,tol).unwrap();
      let expected = vec![
      array![[I],[ZERO],[ZERO],[ZERO],[-I]],
    ];   
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
      let expected = vec![
      array![[ZERO],[0.5*I],[ZERO],[0.5*I],[ZERO]],
    ];   
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
      let expected = vec![
      array![[ZERO],[ZERO],[SQRT2*SQRT3*ONE],[ZERO],[ZERO]],
    ];   
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
      let expected = vec![
      array![[ZERO],[0.5*ONE],[ZERO],[-0.5*ONE],[ZERO]],
    ];   
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
      let expected = vec![
      array![[ONE],[ZERO],[ZERO],[ZERO],[ONE]],
    ];   
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

    let stevens = vec![
      spin_stevens(mult, k, -2).unwrap(),
      spin_stevens(mult, k, -1).unwrap(),
      spin_stevens(mult, k,  0).unwrap(),
      spin_stevens(mult, k,  1).unwrap(),
      spin_stevens(mult, k,  2).unwrap(),
    ];   

    let expected = vec![
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
      let expected = vec![array![[ONE]] ];
      assert_eq!(coefs.len(),expected.len());
      for (ii,c) in coefs.iter().enumerate(){
        let c0 = &expected[ii];
        assert_eq!(c.shape(),c0.shape());
        assert!( (c-c0).norm() < tol );
      }

      let sz = spin_z(mult);
      let coefs = spherical_operator_decomposition(&sz, tol);
      let expected = vec![array![[ZERO],[ONE],[ZERO]] ];
      assert_eq!(coefs.len(),expected.len());
      for (ii,c) in coefs.iter().enumerate(){
        let c0 = &expected[ii];
        assert_eq!(c.shape(),c0.shape());
        assert!( (c-c0).norm() < tol );
      }

      let t = spin_z(mult) + CxMat::eye(mult);
      let coefs = spherical_operator_decomposition(&t, tol);
      let expected = vec![
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
      for l in 0..=(2 as i32){
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

    let d = expmap_spin(2, &SpinOp::Sy,-PI).unwrap();
    let u = array![
      [ZERO,ONE],
      [-ONE,ZERO],
    ];
    let err = (d - u).norm();
    assert!(err < 1e-12);
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

}



