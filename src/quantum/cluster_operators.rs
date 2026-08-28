use crate::physical_constants::*;
use crate::clue_errors::*;
use crate::config::{Config, DetectedPopulation};
use crate::quantum::pulse_sequences::{
  ideal_pulse,
  PI_OVER_2_PULSE_NAME,
  PI_PULSE_NAME
};
use crate::math::commutator;

use std::fmt;
use std::collections::HashMap;
use ndarray::Array2;
use ndarray::linalg::kron;
use num_complex::Complex;

type CxMat = Array2::<Complex<f64>>;
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
  pub fn new(det_multiplicity: usize, 
      bath_multiplicities: &[usize], max_size: usize,config: &Config) 
   -> Result<Self, CluEError> {
    

    let n_mults = bath_multiplicities.len();

    let mut cluster_spin_ops = Vec::<KronSpinOperators>::with_capacity(n_mults);

    let detection_ops = build_detection_operators(det_multiplicity, 
        bath_multiplicities, max_size, config)?;

    for spin_multiplicity in bath_multiplicities.iter() {

      let sops = KronSpinOperators::new(det_multiplicity,
          *spin_multiplicity, max_size)?;

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
fn build_detection_operators(det_multiplicity: usize, 
    bath_multiplicities: &[usize],
    max_size: usize,
    config: &Config) 
    -> Result<Vec::<DetectionSpinOperators>,CluEError>
{


  match det_multiplicity {
    0 => return Err(CluEError::InvalidSpinMultiplicity(0)),
    1 => return Ok(Vec::<DetectionSpinOperators>::new()),
    _ => (),
  }

  let Some(detected_population) = &config.detected_population else{
    return Err(CluEError::NoDetectedSpinDensityMatrix);
  };

  let density_matrix_opt = match detected_population{
    DetectedPopulation::Matrix(mat) => Some(mat.clone()),
    DetectedPopulation::S(sop) => Some(get_spin_operator(det_multiplicity,sop)),
    DetectedPopulation::Thermal => None,   
  };

  let Some(det_op) =&config.detection_operator else {
    return Err(CluEError::NoDetectedSpinDetectionOperator);
  }; 

  let pulses = if config.pulses.is_empty(){
    match &config.detected_spin_transition{
      Some(transition) => {
        HashMap::<String,CxMat>::from([
          (PI_OVER_2_PULSE_NAME.to_string(), ideal_pulse(&SpinOp::Sy, 0.5*PI,
                                 det_multiplicity, transition)),
          (PI_PULSE_NAME.to_string(), ideal_pulse(&SpinOp::Sy, PI,
                                 det_multiplicity, transition)),
        ])},
      None => return Err(CluEError::NoDetectedSpinTransition),
    }  
  }else{
    config.pulses.clone()
  };


  let n_mults = bath_multiplicities.len();
  let mut detection_ops = Vec::<DetectionSpinOperators>::with_capacity(n_mults);

  for spin_multiplicity in bath_multiplicities.iter() {
    detection_ops.push(DetectionSpinOperators::new(
        &density_matrix_opt,
        det_op,
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
}

impl<'a> KronSpinOperators {

  /// This function builds `KronSpinOperators` for spin with the input spin 
  /// multiplicity for clusters up to size `max_size`.
  pub fn new(det_multiplicity: usize,
      spin_multiplicity: usize, max_size: usize) 
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

    Ok(KronSpinOperators{
      sx_list,  
      sy_list,  
      sz_list,  
      sp_list,  
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

impl SpinOp{
  pub fn from_str(s: &str) -> Result<Self,CluEError>{
    match s{
      "sx" | "Sx" => Ok(Self::Sx),  
      "sy" | "Sy" => Ok(Self::Sy),  
      "sz" | "Sz" => Ok(Self::Sz),  
      "s+" | "S+" => Ok(Self::Sp),  
      "s-" | "S+" => Ok(Self::Sm),  
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

  let sp = spin_plus(spin_multiplicity);
  let sm = spin_minus(spin_multiplicity);
  let mut t = (-SQRT2_INV).powi(l as i32) * ONE * sp;

  for n in 0..2*l{
    if l - n  == m { break; }
    let a = ONE/( (l*(l+1) - m*(m-1)) as f64 ).sqrt(); 
    t = commutator(&sm,&t)*a;
  }

  t
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



#[cfg(test)]
mod tests {
  use super::*;
  use ndarray::array;

  //----------------------------------------------------------------------------
  #[test]
  #[allow(non_snake_case)]
  fn test_ClusterSpinOperators() {

    let mut config = Config::new();
    config.set_defaults().unwrap();

    let spin_multiplicities = vec![2,3];
    let max_size = 3;
    let sops = ClusterSpinOperators::new(1,
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
    let sops = KronSpinOperators::new(1,spin_multiplicity, max_size).unwrap();

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
    let t = spin_ist(2,1,1);
    assert!(approx_eq(&t, &(-SQRT2_INV*ONE*spin_plus(2)),tol));
    let t = spin_ist(2,1,0);
    assert!(approx_eq(&t, &spin_z(2),tol));
    let t = spin_ist(2,1,-1);
    assert!(approx_eq(&t, &(SQRT2_INV*ONE*spin_minus(2)),tol));

    let t = spin_ist(2,2,-2);
    assert!(approx_eq(&t, &(0.5*ONE*spin_minus(2)*spin_minus(2)),tol));
    let t = spin_ist(3,2,-2);
    assert!(approx_eq(&t, &(0.5*ONE*spin_minus(3)*spin_minus(3)),tol));
    let t = spin_ist(4,2,-2);
    assert!(approx_eq(&t, &(0.5*ONE*spin_minus(4)*spin_minus(4)),tol));
    
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



