use crate::physical_constants::PI;
use crate::CluEError;
use crate::misc::{
  are_all_same_type,
  mat_from_toml_array,
};

use rand::distributions::Uniform;
use rand_distr::Distribution;
use rand_chacha::ChaCha20Rng;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug, Clone,PartialEq)]
pub struct SymmetricTensor3D{
  pub elements: [f64; 6],
}

impl SymmetricTensor3D{
  //----------------------------------------------------------------------------
  pub fn from_toml_value(value: toml::Value) -> Result<Self,CluEError> 
  {
    match value{
      toml::Value::Array(array) => Self::from_toml_array(array),
      _ => Err(CluEError::TOMLArrayDoesNotSpecifyATensor)
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_array(array: Vec::<toml::Value>) -> Result<Self,CluEError>
  {

    if array.is_empty(){
      return Err(CluEError::TOMLArrayIsEmpty);
    }
    if !are_all_same_type(&array){
      return Err(CluEError::TOMLArrayContainsMultipleTypes);
    }
    let vec: Vec::<f64> = if array[0].is_array(){
        let m = mat_from_toml_array(array)?;
        if m.dim().0 != 3 || m.dim().1 != 3{
          return Err(CluEError::TOMLArrayDoesNotSpecifyATensor);
        }
        vec![ m[[0,0]], m[[0,1]], m[[0,2]],
              m[[1,0]], m[[1,1]], m[[1,2]],
              m[[2,0]], m[[2,1]], m[[2,2]],
         ]
    }else if array[0].is_float(){
      array.iter()
          .filter_map(|v| v.as_float()).collect()
    }else{      
      return Err(CluEError::TOMLArrayDoesNotSpecifyATensor);
    };


    match vec.len(){
      1 => Ok(Self::eye().scale(vec[0] )),
      3 => Ok(Self{ elements: 
            [vec[0],    0.0, 0.0,
                     vec[1], 0.0,
                             vec[2] ] 
      }),
      6 => Ok(Self{ elements: 
            [vec[0], vec[1], vec[2],
                     vec[3], vec[4],
                             vec[5] ] 
      }),
      9 => Ok(Self{ elements: 
            [vec[0], vec[1], vec[2],
                     vec[4], vec[5],
                             vec[8] ] 
      }),
      _ => Err(CluEError::TOMLArrayDoesNotSpecifyATensor)
    }
  }
  
  //----------------------------------------------------------------------------
  pub fn any_nan(&self) -> bool {
    for x in self.elements.iter(){
      if x.is_nan(){
        return true;
      }
    }
    false
  }
  //----------------------------------------------------------------------------
  pub fn zeros() -> Self{
    SymmetricTensor3D{elements: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}
  }
  //----------------------------------------------------------------------------
  pub fn from(elements: [f64;6]) -> Self{
    SymmetricTensor3D{elements}
  }
  //----------------------------------------------------------------------------
  pub fn eye() -> Self{
    SymmetricTensor3D{elements: [1.0, 0.0, 0.0, 
                                      1.0, 0.0, 
                                           1.0]}
  }
  //----------------------------------------------------------------------------
  pub fn rotate_pasive(&self, dir: &UnitSpherePoint) -> Self
  {
    let txx = self.xx();
    let txy = self.xy();
    let txz = self.xz();
    let tyy = self.yy();
    let tyz = self.yz();
    let tzz = self.zz();

    let cos_theta = dir.cos_theta();
    let cos_theta_squared = dir.cos_theta_squared();

    let sin_theta = dir.sin_theta();
    let sin_theta_squared = dir.sin_theta_squared();

    let cos_phi = dir.cos_phi();
    let cos_phi_squared = dir.cos_phi_squared();

    let sin_phi = dir.sin_phi();
    let sin_phi_squared = dir.sin_phi_squared();

    SymmetricTensor3D{ elements: [
      // xx
      txx*cos_theta_squared*cos_phi_squared 
      + 2.0*txy*sin_phi*cos_theta_squared*cos_phi 
      - 2.0*txz*sin_theta*cos_theta*cos_phi 
      + tyy*sin_phi_squared*cos_theta_squared 
      - 2.0*tyz*sin_theta*sin_phi*cos_theta 
      + tzz*sin_theta_squared
      ,
      // xy
      (-txx*sin_phi*cos_theta*cos_phi 
      - txy*sin_phi_squared*cos_theta 
      + txy*cos_theta*cos_phi_squared 
      + txz*sin_theta*sin_phi 
      + tyy*sin_phi*cos_theta*cos_phi 
      - tyz*sin_theta*cos_phi)
      ,
      // xz
      (txx*sin_theta*cos_theta*cos_phi_squared 
      + 2.0*txy*sin_theta*sin_phi*cos_theta*cos_phi 
      - txz*sin_theta_squared*cos_phi 
      + txz*cos_theta_squared*cos_phi 
      + tyy*sin_theta*sin_phi_squared*cos_theta 
      - tyz*sin_theta_squared*sin_phi 
      + tyz*sin_phi*cos_theta_squared 
      - tzz*sin_theta*cos_theta)
      ,
      // yy
      txx*sin_phi_squared 
      - 2.0*txy*sin_phi*cos_phi 
      + tyy*cos_phi_squared
      ,
      // yz
      -txx*sin_theta*sin_phi*cos_phi 
      - txy*sin_theta*sin_phi_squared 
      + txy*sin_theta*cos_phi_squared 
      - txz*sin_phi*cos_theta 
      + tyy*sin_theta*sin_phi*cos_phi 
      + tyz*cos_theta*cos_phi
      ,
      // zz
      txx*sin_theta_squared*cos_phi_squared 
      + 2.0*txy*sin_theta_squared*sin_phi*cos_phi 
      + 2.0*txz*sin_theta*cos_theta*cos_phi 
      + tyy*sin_theta_squared*sin_phi_squared 
      + 2.0*tyz*sin_theta*sin_phi*cos_theta 
      + tzz*cos_theta_squared
      ]}
  }
  //----------------------------------------------------------------------------
  /*
  pub fn rotate_active(&self, dir: &UnitSpherePoint) -> Self
  {

    let txx = self.xx();
    let txy = self.xy();
    let txz = self.xz();
    let tyy = self.yy();
    let tyz = self.yz();
    let tzz = self.zz();

    let cos_theta = dir.cos_theta();
    let cos_theta_squared = dir.cos_theta_squared();

    let sin_theta = dir.sin_theta();
    let sin_theta_squared = dir.sin_theta_squared();

    let cos_phi = dir.cos_phi();
    let cos_phi_squared = dir.cos_phi_squared();

    let sin_phi = dir.sin_phi();
    let sin_phi_squared = dir.sin_phi_squared();

    SymmetricTensor3D{ elements: [
      sin_phi_squared*tyy - 2.0*cos_phi*sin_phi*(cos_theta*txy + sin_theta*tyz) 
        +cos_phi_squared*(cos_theta_squared*txx + 2.0*cos_theta*sin_theta*txz 
            + sin_theta_squared*tzz)
      ,
      cos_phi_squared*(cos_theta*txy + sin_theta*tyz) 
        -sin_phi_squared*(cos_theta*txy + sin_theta*tyz) 
        +cos_phi*sin_phi*(cos_theta_squared*txx + 2.0*cos_theta*sin_theta*txz 
            - tyy + sin_theta_squared*tzz)
      ,
      -sin_theta*(cos_phi*cos_theta*txx - sin_phi*txy + cos_phi*sin_theta*txz) 
        +cos_theta*(cos_phi*cos_theta*txz - sin_phi*tyz + cos_phi*sin_theta*tzz)
      ,
      cos_theta_squared*sin_phi_squared*txx + 2.0*cos_theta*sin_phi*(cos_phi*txy 
          + sin_phi*sin_theta*txz) +cos_phi_squared*tyy 
        + 2.0*cos_phi*sin_phi*sin_theta*tyz 
        + sin_phi_squared*sin_theta_squared*tzz
      ,
      -sin_theta*(cos_theta*sin_phi*txx + cos_phi*txy + sin_phi*sin_theta*txz) 
        + cos_theta*(cos_theta*sin_phi*txz + cos_phi*tyz 
            + sin_phi*sin_theta*tzz)
      ,
      sin_theta_squared*txx - 2.0*cos_theta*sin_theta*txz 
        + cos_theta_squared*tzz
    ]}
  }
  */
  //----------------------------------------------------------------------------
  pub fn trace(&self) -> f64{
    self.xx() + self.yy() + self.zz()
  }
  //----------------------------------------------------------------------------
  pub fn scale(&self, a: f64) -> SymmetricTensor3D {
    let mut tensor = self.clone();
    tensor.scale_mut(a);
    tensor
  }
  //----------------------------------------------------------------------------
  pub fn scale_mut(&mut self, a: f64){
    for x in self.elements.iter_mut(){
      *x *= a;
    }
  }
  //----------------------------------------------------------------------------

  pub fn xx(&self) -> f64 { self.elements[0]}
  pub fn xy(&self) -> f64 { self.elements[1]}
  pub fn xz(&self) -> f64 { self.elements[2]}
  pub fn yx(&self) -> f64 { self.xy() }
  pub fn yy(&self) -> f64 { self.elements[3]}
  pub fn yz(&self) -> f64 { self.elements[4]}
  pub fn zx(&self) -> f64 { self.xz()}
  pub fn zy(&self) -> f64 { self.yz()}
  pub fn zz(&self) -> f64 { self.elements[5]}

  pub fn set_xx(&mut self, value: f64) { self.elements[0] = value;}
  pub fn set_xy(&mut self, value: f64) { self.elements[1] = value;}
  pub fn set_xz(&mut self, value: f64) { self.elements[2] = value;}
  pub fn set_yx(&mut self, value: f64) { self.elements[1] = value;}
  pub fn set_yy(&mut self, value: f64) { self.elements[3] = value;}
  pub fn set_yz(&mut self, value: f64) { self.elements[4] = value;}
  pub fn set_zx(&mut self, value: f64) { self.elements[2] = value;}
  pub fn set_zy(&mut self, value: f64) { self.elements[4] = value;}
  pub fn set_zz(&mut self, value: f64) { self.elements[5] = value;}

}


impl std::fmt::Display for SymmetricTensor3D {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    write!(f,"[{}, {}, {}, {}, {}, {}]", self.elements[0], self.elements[1],
        self.elements[2], self.elements[3],self.elements[4], self.elements[5])
  }
}
impl std::ops::Mul<f64> for &SymmetricTensor3D
{
  type Output = SymmetricTensor3D ;

  fn mul(self, rhs: f64) -> SymmetricTensor3D  
  {
   self.scale(rhs)
  }
}
impl std::ops::Mul<&SymmetricTensor3D> for f64
{
  type Output = SymmetricTensor3D;

  fn mul(self, rhs: &SymmetricTensor3D ) -> SymmetricTensor3D  
  {
    rhs.scale(self)
  }
}
impl std::ops::Add<&SymmetricTensor3D> for &SymmetricTensor3D{
  type Output = SymmetricTensor3D;
  fn add(self, rhs: &SymmetricTensor3D ) -> SymmetricTensor3D{
    let n = self.elements.len();
    assert_eq!(n, rhs.elements.len());
    let mut out = self.clone();

    for ii in 0..n{
      out.elements[ii] += rhs.elements[ii];  
    }
    out
  }
}
//------------------------------------------------------------------------------
impl std::ops::Sub<&SymmetricTensor3D> for &SymmetricTensor3D{
  type Output = SymmetricTensor3D;
  fn sub(self, rhs: &SymmetricTensor3D ) -> SymmetricTensor3D{
    let n = self.elements.len();
    assert_eq!(n, rhs.elements.len());
    let mut out = self.clone();

    for ii in 0..n{
      out.elements[ii] -= rhs.elements[ii];  
    }
    out
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug, Clone)]
pub struct Vector3D{
  pub elements: [f64; 3],
}

impl Vector3D{

  //----------------------------------------------------------------------------
  pub fn from_toml_array(array: Vec::<toml::Value>) -> Result<Self,CluEError>
  {
    if array.len() != 3 {
      return Err(CluEError::TOMLArrayDoesNotSpecifyAVector);
    }
    if !are_all_same_type(&array){
      return Err(CluEError::TOMLArrayContainsMultipleTypes);
    }
    if !array[0].is_float(){
      return Err(CluEError::TOMLArrayDoesNotSpecifyAVector);
    }

    let r = array.iter().filter_map(|x| x.as_float())
      .collect::<Vec::<f64>>();
    
    Ok(Self{elements: [r[0],r[1],r[2]]})

  }
  //----------------------------------------------------------------------------
  pub fn any_nan(&self) -> bool {
    for x in self.elements.iter(){
      if x.is_nan(){
        return true;
      }
    }
    false
  }
  //----------------------------------------------------------------------------
  pub fn from(elements: [f64;3]) -> Self{
    Vector3D{elements}
  }
  //----------------------------------------------------------------------------
  pub fn from_vec(v: Vec::<f64>) -> Result<Self,CluEError>{
    if v.len() != 3{
      return Err(CluEError::NotA3DVector(v.len()));
    }
    Ok(Self::from([ v[0],v[1],v[2] ]))
  }
  //----------------------------------------------------------------------------
  pub fn zeros() -> Self{
    Vector3D{elements: [0.0, 0.0, 0.0]}
  }
  //----------------------------------------------------------------------------
  pub fn random_direction(rng: &mut ChaCha20Rng) -> Self {
    let range = Uniform::new(0.0f64, 1.0);
    
    let random_number_for_phi = range.sample(rng);
    let phi: f64 = 2.0*PI*random_number_for_phi;

    let random_number_for_theta = range.sample(rng);
    let theta: f64 = (2.0*random_number_for_theta-1.0).acos();

    let x = theta.sin()*phi.cos();
    let y = theta.sin()*phi.sin();
    let z = theta.cos();

    Vector3D::from([x,y,z])
  }
  //----------------------------------------------------------------------------
  pub fn x(&self) -> f64 { self.elements[0]}
  //----------------------------------------------------------------------------
  pub fn y(&self) -> f64 { self.elements[1]}
  //----------------------------------------------------------------------------
  pub fn z(&self) -> f64 { self.elements[2]}
  //----------------------------------------------------------------------------

  pub fn set_x(&mut self, value: f64) { self.elements[0] = value;}
  //----------------------------------------------------------------------------
  pub fn set_y(&mut self, value: f64) { self.elements[1] = value;}
  //----------------------------------------------------------------------------
  pub fn set_z(&mut self, value: f64) { self.elements[2] = value;}

  //----------------------------------------------------------------------------
  pub fn r(&self) -> f64{
    self.elements.iter().map(|x| x*x).sum::<f64>().sqrt()
  }

  //----------------------------------------------------------------------------
  pub fn theta(&self) -> f64 {
    (self.z()/self.r()).acos()
  }

  //----------------------------------------------------------------------------
  pub fn phi(&self) -> f64 {
  
    let rho = (self.x()*self.x() + self.y()*self.y()).sqrt();
    if rho < 1e-12{ return 0.0};
    let mut phi = (self.x()/rho).acos();
    
    if self.y() < 0.0 {
      phi = -phi;
    }
    phi
  }
  //----------------------------------------------------------------------------
  pub fn norm_squared(&self) -> f64 {

    let mut l2 = 0.0;
    for el in self.elements.iter(){
      l2 += (*el)*(*el);
    }
    l2
  }
  //----------------------------------------------------------------------------
  pub fn norm(&self) -> f64 {

    self.norm_squared().sqrt()
  
  }
  //----------------------------------------------------------------------------
  pub fn scale(&self,scale: f64) -> Vector3D {                                   
    let mut out = self.clone();
    for el in out.elements.iter_mut(){
      *el *= scale;
    }
    out
  }                                               
  //----------------------------------------------------------------------------
  pub fn normalize(&self) -> Vector3D {                                        
    let r = self.norm();                                                        
    self.scale(1.0/r)                                                            
  }  
  //----------------------------------------------------------------------------
  pub fn dot(&self, other: &Vector3D) -> f64 {

    assert_eq!(self.elements.len(),other.elements.len());
    let mut out = 0.0;
    for ii in 0..self.elements.len(){
      out += self.elements[ii]*other.elements[ii];
    }

    out
  } 
  //----------------------------------------------------------------------------
  pub fn cross(&self, other: &Vector3D) -> Self{
    assert_eq!(self.elements.len(),3);
    assert_eq!(other.elements.len(),3);
    let mut out = self.clone();
    
    out.set_x(self.y()*other.z() - self.z()*other.y() );
    out.set_y(self.z()*other.x() - self.x()*other.z() );
    out.set_z(self.x()*other.y() - self.y()*other.x() );
    out
  }
  //----------------------------------------------------------------------------
  /// This function calculate the ___v___ trans(___v___).
  pub fn self_outer(&self) -> SymmetricTensor3D{
    let elements = [
      self.x()*self.x(), 
      self.x()*self.y(), 
      self.x()*self.z(), 
      self.y()*self.y(), 
      self.y()*self.z(), 
      self.z()*self.z(), 
    ];
    SymmetricTensor3D::from(elements)
  }
  //----------------------------------------------------------------------------
  /*
  pub fn rotate_active(&self, dir: &UnitSpherePoint) -> Self
  {
    let x = self.x();
    let y = self.y();
    let z = self.z();

    let cos_theta = dir.cos_theta();
    let sin_theta = dir.sin_theta();

    let cos_phi = dir.cos_phi();
    let sin_phi = dir.sin_phi();

    Vector3D{ elements: [
      cos_phi*cos_theta*x - sin_phi*y + cos_phi*sin_theta*z,
      cos_theta*sin_phi*x + cos_phi*y + sin_phi*sin_theta*z,
      -sin_theta*x + cos_theta*z
    ]}
  }
  */
  //----------------------------------------------------------------------------
  pub fn rotate_pasive(&self, dir: &UnitSpherePoint) -> Self
  {
    let x = self.x();
    let y = self.y();
    let z = self.z();

    let cos_theta = dir.cos_theta();
    let sin_theta = dir.sin_theta();

    let cos_phi = dir.cos_phi();
    let sin_phi = dir.sin_phi();

    Vector3D{ elements: [
      cos_phi*cos_theta*x + sin_phi*cos_theta*y - sin_theta*z,
      -sin_phi*x + cos_phi*y,
      sin_theta*cos_phi*x + sin_phi*sin_theta*y + cos_theta*z
    ]}
  }
  //----------------------------------------------------------------------------
}

impl std::ops::Add<&Vector3D> for &Vector3D{
  type Output = Vector3D;
  fn add(self, rhs: &Vector3D ) -> Vector3D{
    let n = self.elements.len();
    assert_eq!(n, rhs.elements.len());
    let mut out = self.clone();

    for ii in 0..n{
      out.elements[ii] += rhs.elements[ii];  
    }
    out
  }
}
//------------------------------------------------------------------------------
impl std::ops::Sub<&Vector3D> for &Vector3D{
  type Output = Vector3D;
  fn sub(self, rhs: &Vector3D ) -> Vector3D{
    let n = self.elements.len();
    assert_eq!(n, rhs.elements.len());
    let mut out = self.clone();

    for ii in 0..n{
      out.elements[ii] -= rhs.elements[ii];  
    }
    out
  }
}
//------------------------------------------------------------------------------
impl std::ops::Neg for Vector3D {

  type Output = Vector3D;
  fn neg(self) -> Vector3D{
    let mut out = self;
    let n = out.elements.len();
    for ii in 0..n{
      out.elements[ii] = -out.elements[ii];  
    }

    out
  }
}
//------------------------------------------------------------------------------
impl std::ops::Mul<f64> for &Vector3D
{
  type Output = Vector3D ;

  fn mul(self, rhs: f64) -> Vector3D  
  {
   self.scale(rhs)
  }
}
//------------------------------------------------------------------------------
impl std::ops::Mul<&Vector3D> for f64
{
  type Output = Vector3D;

  fn mul(self, rhs: &Vector3D ) -> Vector3D  
  {
    rhs.scale(self)
  }
}
//------------------------------------------------------------------------------
impl std::cmp::Ord for Vector3D{
  fn cmp(&self, other: &Self) -> std::cmp::Ordering {
    ( (self.norm()*1e12) as u64).cmp(&((other.norm()*1e12) as u64))
  }
}
//------------------------------------------------------------------------------
impl PartialOrd for Vector3D{
  fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
    Some(self.cmp(other))
  }
}
//------------------------------------------------------------------------------
impl PartialEq for Vector3D{
  fn eq(&self, other: &Self) -> bool {
    (self - other).norm() < 1e-12
  }
}
impl Eq for Vector3D {}
//------------------------------------------------------------------------------
impl std::fmt::Display for Vector3D {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    write!(f,"[{}, {}, {}]", self.elements[0], self.elements[1],
        self.elements[2])
  }
}
//------------------------------------------------------------------------------
//
// Let n be an integer and v,v0,v1, and step be real-valued vectors,
// This function returns the v1 := v + n*step founds by
//
//     v1 = argmin(|v0-v1|) = argmin( |v0 - (v + n*step)| ).
//
pub fn minimize_absolute_difference_for_vector3d_step(
    v: &Vector3D, v0: &Vector3D, step: &Vector3D,) -> Vector3D    
{


  let mut v1 = v.clone();

  let mut err0 = (&v1-v0).norm();
  let mut err = (&v1-&(v0 + step)).norm();

  let delta_v: Vector3D = if err < err0{
    step.clone()
  }else{
    -step.clone()
  };

  
  loop{
  
    let v2 = &v1 + &delta_v;
    err0 = err;
    err =  (&v2-v0).norm();
    if err <= err0{
      v1 = v2;
    }else{
      break;
    }
  
  }

  v1
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This function calulates ___T___***v***, where ___T___ is a 
/// `SymmetricTensor3D`. and ___v___ a `Vector3D`.
pub fn symten3d_times_vec3d(ten: &SymmetricTensor3D, vec: &Vector3D) 
  -> Vector3D
{
  Vector3D::from([
    vec.x()*ten.xx() + vec.y()*ten.xy() + vec.z()*ten.xz(),
    vec.x()*ten.yx() + vec.y()*ten.yy() + vec.z()*ten.yz(),
    vec.x()*ten.zx() + vec.y()*ten.zy() + vec.z()*ten.zz(),
  ])

} 
impl std::ops::Mul<&Vector3D> for &SymmetricTensor3D
{
  type Output = Vector3D;

  fn mul(self, rhs: &Vector3D) -> Vector3D {
    symten3d_times_vec3d(self,rhs)
  }
}
impl std::ops::Mul<&SymmetricTensor3D> for &Vector3D
{
  type Output = Vector3D;

  fn mul(self, rhs: &SymmetricTensor3D) -> Vector3D {
    symten3d_times_vec3d(rhs,self)
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug, Clone,PartialEq)]
pub struct UnitSpherePoint{ 
  theta: f64,
  phi: f64,
  cos_theta: f64, cos_theta_squared: f64,
  sin_theta: f64, sin_theta_squared: f64,
  cos_phi: f64, cos_phi_squared: f64, 
  sin_phi: f64, sin_phi_squared: f64
}

impl UnitSpherePoint{
  pub fn new(theta: f64, phi: f64) -> Self{

    let cos_theta = theta.cos();
    let sin_theta = theta.sin();
    
    let cos_theta_squared = cos_theta*cos_theta;
    let sin_theta_squared = sin_theta*sin_theta;

    let cos_phi = phi.cos();
    let sin_phi = phi.sin();

    let cos_phi_squared = cos_phi*cos_phi;
    let sin_phi_squared = sin_phi*sin_phi;

    UnitSpherePoint{ 
      theta,
      phi,
      cos_theta, cos_theta_squared,
      sin_theta, sin_theta_squared,
      cos_phi, cos_phi_squared, 
      sin_phi, sin_phi_squared
    }
  }
  //----------------------------------------------------------------------------
  pub fn from(r: Vector3D) -> Self{
    let theta = r.theta();
    let phi = r.phi();

    UnitSpherePoint::new(theta,phi)
  }
  //----------------------------------------------------------------------------
  pub fn theta(&self) -> f64{ self.theta }
  //----------------------------------------------------------------------------
  pub fn cos_theta(&self) -> f64{ self.cos_theta }
  //----------------------------------------------------------------------------
  pub fn cos_theta_squared(&self) -> f64{ self.cos_theta_squared }
  //----------------------------------------------------------------------------
  pub fn sin_theta(&self) -> f64{ self.sin_theta }
  //----------------------------------------------------------------------------
  pub fn sin_theta_squared(&self) -> f64{ self.sin_theta_squared }
  //----------------------------------------------------------------------------
  pub fn phi(&self) -> f64{ self.phi }
  //----------------------------------------------------------------------------
  pub fn cos_phi(&self) -> f64{ self.cos_phi }
  //----------------------------------------------------------------------------
  pub fn cos_phi_squared(&self) -> f64{ self.cos_phi_squared }
  //----------------------------------------------------------------------------
  pub fn sin_phi(&self) -> f64{ self.sin_phi }
  //----------------------------------------------------------------------------
  pub fn sin_phi_squared(&self) -> f64{ self.sin_phi_squared }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  use crate::physical_constants::PI;

  //----------------------------------------------------------------------------
  /*
  #[test]
  fn test_rotate_active(){

    let tol = 1e-12;
    let x_axis = Vector3D::from([1.0,0.0,0.0]);
    let y_axis = Vector3D::from([0.0,1.0,0.0]);
    let z_axis = Vector3D::from([0.0,0.0,1.0]);
    let xyz = Vector3D::from([1.0,1.0,1.0]);
    let ten = SymmetricTensor3D::from([2.0, 0.0, 0.0,
                                            73.0, 0.0,
                                                 547.0 ]);
    
    let vec4 = Vector3D::from([-1.0,4.1,16.3]);
    let ten4 = vec4.self_outer();


    let phi_list = (0..=10).map(|x| (x as f64)/10.0*PI*2.0)
      .collect::<Vec::<f64>>();

    let theta_list = (0..=10).map(|x| (2.0*(x as f64)/10.0 - 1.0f64).acos())
      .collect::<Vec::<f64>>();

    for &theta in theta_list.iter(){
      let ct = theta.cos();
      let st = theta.sin();

      for &phi in phi_list.iter(){
        let cp = phi.cos();
        let sp = phi.sin();

        let point = UnitSpherePoint::new(theta,phi);

        let rot_vec = x_axis.rotate_active(&point);
        let rot_x = Vector3D::from([ct*cp, ct*sp, -st]);
        assert!( (&rot_vec-&rot_x).norm() < tol );

        let rot_vec = y_axis.rotate_active(&point);
        let rot_y = Vector3D::from([-sp, cp, 0.0]);
        assert!( (&rot_vec-&rot_y).norm() < tol );

        let rot_vec = z_axis.rotate_active(&point);
        let rot_z = Vector3D::from([st*cp,st*sp,ct]);
        assert!( (&rot_vec-&rot_z).norm() < tol );

        let rot_vec = xyz.rotate_active(&point);
        let rot_xyz = Vector3D::from([
            ct*cp -sp + st*cp, 
            ct*sp + cp +st*sp, 
            -st + ct]);
        assert!( (&rot_vec-&rot_xyz).norm() < tol );

        let rot_ten = ten.rotate_active(&point);
        let ref_ten = &(
            &(&rot_x.self_outer()*2.0) + &(&rot_y.self_outer()*73.0))
          + &(&rot_z.self_outer()*547.0);

        assert!((rot_ten.xx()-ref_ten.xx()).abs() < tol );
        assert!((rot_ten.xy()-ref_ten.xy()).abs() < tol );
        assert!((rot_ten.xz()-ref_ten.xz()).abs() < tol );
        assert!((rot_ten.yy()-ref_ten.yy()).abs() < tol );
        assert!((rot_ten.yz()-ref_ten.yz()).abs() < tol );
        assert!((rot_ten.zz()-ref_ten.zz()).abs() < tol );

        let rot_vec4 = vec4.rotate_active(&point);
        let ref_ten4 = rot_vec4.self_outer();
        let rot_ten4 = ten4.rotate_active(&point);
        assert!((rot_ten4.xx()-ref_ten4.xx()).abs() < tol );
        assert!((rot_ten4.xy()-ref_ten4.xy()).abs() < tol );
        assert!((rot_ten4.xz()-ref_ten4.xz()).abs() < tol );
        assert!((rot_ten4.yy()-ref_ten4.yy()).abs() < tol );
        assert!((rot_ten4.yz()-ref_ten4.yz()).abs() < tol );
        assert!((rot_ten4.zz()-ref_ten4.zz()).abs() < tol );

      }
    }

  }
  */
  //----------------------------------------------------------------------------
  #[test]
  fn test_rotate_pasive(){

    let tol = 1e-12;
    let x_axis = Vector3D::from([1.0,0.0,0.0]);
    let y_axis = Vector3D::from([0.0,1.0,0.0]);
    let z_axis = Vector3D::from([0.0,0.0,1.0]);
    let xyz = Vector3D::from([1.0,1.0,1.0]);
    let ten = SymmetricTensor3D::from([2.0, 0.0, 0.0,
                                            73.0, 0.0,
                                                 547.0 ]);
    
    let vec4 = Vector3D::from([-1.0,4.1,16.3]);
    let ten4 = vec4.self_outer();


    let phi_list = (0..=10).map(|x| (x as f64)/10.0*PI*2.0)
      .collect::<Vec::<f64>>();

    let theta_list = (0..=10).map(|x| (2.0*(x as f64)/10.0 - 1.0f64).acos())
      .collect::<Vec::<f64>>();

    for &theta in theta_list.iter(){
      let ct = theta.cos();
      let st = theta.sin();

      for &phi in phi_list.iter(){
        let cp = phi.cos();
        let sp = phi.sin();

        let point = UnitSpherePoint::new(theta,phi);

        let rot_vec = x_axis.rotate_pasive(&point);
        assert!((rot_vec.norm() - x_axis.norm()) < tol);
        let rot_x = Vector3D::from([ct*cp, -sp, st*cp]);
        assert!( (&rot_vec-&rot_x).norm() < tol );

        let rot_vec = y_axis.rotate_pasive(&point);
        assert!((rot_vec.norm() - y_axis.norm()) < tol);
        let rot_y = Vector3D::from([sp*ct, cp, st*sp]);
        assert!( (&rot_vec-&rot_y).norm() < tol );

        let rot_vec = z_axis.rotate_pasive(&point);
        assert!((rot_vec.norm() - z_axis.norm()) < tol);
        let rot_z = Vector3D::from([-st, 0.0,ct]);
        assert!( (&rot_vec-&rot_z).norm() < tol );

        let rot_vec = xyz.rotate_pasive(&point);
        assert!((rot_vec.norm() - xyz.norm()) < tol);
        let rot_xyz = Vector3D::from([
            ct*cp  +sp*ct  - st, 
            -sp + cp, 
            st*cp + st*sp + ct]);
        assert!( (&rot_vec-&rot_xyz).norm() < tol );

        let rot_ten = ten.rotate_pasive(&point);
        let ref_ten = &(
            &(&rot_x.self_outer()*2.0) + &(&rot_y.self_outer()*73.0))
          + &(&rot_z.self_outer()*547.0);

        assert!((rot_ten.xx()-ref_ten.xx()).abs() < tol );
        assert!((rot_ten.xy()-ref_ten.xy()).abs() < tol );
        assert!((rot_ten.xz()-ref_ten.xz()).abs() < tol );
        assert!((rot_ten.yy()-ref_ten.yy()).abs() < tol );
        assert!((rot_ten.yz()-ref_ten.yz()).abs() < tol );
        assert!((rot_ten.zz()-ref_ten.zz()).abs() < tol );

        let rot_vec4 = vec4.rotate_pasive(&point);
        let ref_ten4 = rot_vec4.self_outer();
        let rot_ten4 = ten4.rotate_pasive(&point);
        assert!((rot_ten4.xx()-ref_ten4.xx()).abs() < tol );
        assert!((rot_ten4.xy()-ref_ten4.xy()).abs() < tol );
        assert!((rot_ten4.xz()-ref_ten4.xz()).abs() < tol );
        assert!((rot_ten4.yy()-ref_ten4.yy()).abs() < tol );
        assert!((rot_ten4.yz()-ref_ten4.yz()).abs() < tol );
        assert!((rot_ten4.zz()-ref_ten4.zz()).abs() < tol );

      }
    }

  }
  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_Vector3D_ord(){

    let x = Vector3D::from([1.0, 0.0, 0.0]);
    let y = Vector3D::from([0.0, 1.0, 0.0]);

    let mut vecs = Vec::<Vector3D>::with_capacity(9);
    for ix in -1..=1{
      for iy in -1..=1{
        let v = &x.scale(ix as f64) + &y.scale(iy as f64);
        vecs.push(v);
    }}

    vecs.sort();
    assert_eq!(vecs[0],Vector3D::from([0.0, 0.0, 0.0]));

    let sqrt2 = f64::sqrt(2.0);
    let radii = vec![
      0.0, 
      1.0, 1.0, 1.0, 1.0, 
      sqrt2, sqrt2, sqrt2, sqrt2];
    for ii in 0..9{
      assert_eq!(vecs[ii].norm(),radii[ii]);
    }

  }
  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_Vector3D_cross(){

    let x = Vector3D::from([1.0, 0.0, 0.0]);
    let y = Vector3D::from([0.0, 1.0, 0.0]);
    let z = Vector3D::from([0.0, 0.0, 1.0]);

    assert_eq!(x,y.cross(&z));
    assert_eq!(y,z.cross(&x));
    assert_eq!(z,x.cross(&y));
    assert_eq!(z.scale(-1.0),y.cross(&x));
    assert_eq!(x.scale(-1.0),z.cross(&y));
    assert_eq!(y.scale(-1.0),x.cross(&z));


    let v1 = Vector3D::from([1.0,1.0,1.0]).scale(1.0/3.0f64.sqrt());
    let v2 = Vector3D::from([-1.0,-1.0,2.0]).scale(1.0/6.0f64.sqrt());
    let v3 = v1.cross(&v2);

    assert!( (v1.norm() - 1.0) < 1e-12);
    assert!( (v2.norm() - 1.0) < 1e-12);
    assert!( (v3.norm() - 1.0) < 1e-12);
    assert!( v1.dot(&v2) < 1e-12);
    assert!( v1.dot(&v3) < 1e-12);
    assert!( v2.dot(&v3) < 1e-12);

  }
  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_Vector3D_norm(){
     let phi_list = (0..=10).map(|x| (x as f64)/10.0*PI*2.0)
      .collect::<Vec::<f64>>();

    let theta_list = (0..=10).map(|x| (2.0*(x as f64)/10.0 - 1.0f64).acos())
      .collect::<Vec::<f64>>();

    for &theta in theta_list.iter(){
      let ct = theta.cos();
      let st = theta.sin();

      for &phi in phi_list.iter(){
        let cp = phi.cos();
        let sp = phi.sin();


        let norm_vec = &Vector3D::from([st*cp,st*sp,ct]);

        assert!((norm_vec.norm()-1.0).abs()<1e-12);
      }
    }
  }
  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_Vector3D_normalize(){
     let phi_list = (0..=10).map(|x| (x as f64)/10.0*PI*2.0)
      .collect::<Vec::<f64>>();

    let theta_list = (0..=10).map(|x| (2.0*(x as f64)/10.0 - 1.0f64).acos())
      .collect::<Vec::<f64>>();

    let r_list = vec![0.01, 0.2 , 1.0, 2.1, 10.5];

    for &theta in theta_list.iter(){
      let ct = theta.cos();
      let st = theta.sin();

      for &phi in phi_list.iter(){
        let cp = phi.cos();
        let sp = phi.sin();


        for & r in r_list.iter(){
          let vec = &Vector3D::from([r*st*cp,r*st*sp,r*ct]);
          let norm_vec = vec.normalize();

          assert!((norm_vec.norm()-1.0).abs()<1e-12);
        }
      }
    }
  }
  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_SymmetricTensor3D(){
    let mut ten = SymmetricTensor3D{elements: [1.0, 2.0, 3.0, 
                                                4.0, 5.0, 
                                                     6.0]};
    assert_eq!(ten.xx(),1.0);
    assert_eq!(ten.xy(),2.0);
    assert_eq!(ten.xz(),3.0);

    assert_eq!(ten.yx(),2.0);
    assert_eq!(ten.yy(),4.0);
    assert_eq!(ten.yz(),5.0);
    
    assert_eq!(ten.zx(),3.0);
    assert_eq!(ten.zy(),5.0);
    assert_eq!(ten.zz(),6.0);

    let a = 2.0;
    ten = a*&ten;
    assert_eq!(ten.xx(),a*1.0);
    assert_eq!(ten.xy(),a*2.0);
    assert_eq!(ten.xz(),a*3.0);

    assert_eq!(ten.yx(),a*2.0);
    assert_eq!(ten.yy(),a*4.0);
    assert_eq!(ten.yz(),a*5.0);

    assert_eq!(ten.zx(),a*3.0);
    assert_eq!(ten.zy(),a*5.0);
    assert_eq!(ten.zz(),a*6.0);

    let a = 0.5;
    ten = &ten * a;
    assert_eq!(ten.xx(),1.0);
    assert_eq!(ten.xy(),2.0);
    assert_eq!(ten.xz(),3.0);

    assert_eq!(ten.yx(),2.0);
    assert_eq!(ten.yy(),4.0);
    assert_eq!(ten.yz(),5.0);
    
    assert_eq!(ten.zx(),3.0);
    assert_eq!(ten.zy(),5.0);
    assert_eq!(ten.zz(),6.0);

    ten.set_xx(-1.0);
    assert_eq!(ten.xx(),-1.0);

    ten.set_xy(-2.0);
    assert_eq!(ten.xy(),-2.0);

    ten.set_xz(-3.0);
    assert_eq!(ten.xz(),-3.0);

    ten.set_yx(-20.0);
    assert_eq!(ten.yx(),-20.0);

    ten.set_yy(-4.0);
    assert_eq!(ten.yy(),-4.0);

    ten.set_yz(-5.0);
    assert_eq!(ten.yz(),-5.0);
    
    ten.set_zx(-30.0);
    assert_eq!(ten.zx(),-30.0);

    ten.set_zy(-50.0);
    assert_eq!(ten.zy(),-50.0);

    ten.set_zz(-6.0);
    assert_eq!(ten.zz(),-6.0);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_symten3d_times_vec3d(){
    let vec = Vector3D::from([1.0, 2.0, 3.0]);
    let ten = SymmetricTensor3D{elements: [    1.0, 0.1, 0.01, 
                                           /*0.1*/  1.0, 0.1, 
                                          /*0.01,  0.1*/ 1.0]};
    let x = &ten * &vec;
    assert_eq!(x,Vector3D::from([1.23, 2.4, 3.21])); 
  }
  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_SymmetricTensor3D_self_outer(){
     let phi_list = (0..=10).map(|x| (x as f64)/10.0*PI*2.0)
      .collect::<Vec::<f64>>();

    let theta_list = (0..=10).map(|x| (2.0*(x as f64)/10.0 - 1.0f64).acos())
      .collect::<Vec::<f64>>();

    for &theta in theta_list.iter(){
      let ct = theta.cos();
      let st = theta.sin();

      for &phi in phi_list.iter(){
        let cp = phi.cos();
        let sp = phi.sin();


        let norm_vec = &Vector3D::from([st*cp,st*sp,ct]);

        let nnt = norm_vec.self_outer();

        let ref_nnt = SymmetricTensor3D::from(
            [cp*cp*st*st, cp*sp*st*st, cp*st*ct,
                          sp*sp*st*st, sp*st*ct,
                                          ct*ct]);

        let err_ten = &nnt - &ref_nnt;

        let tol = 1e-9;
        assert!(err_ten.xx().abs() < tol);
        assert!(err_ten.xy().abs() < tol);
        assert!(err_ten.xz().abs() < tol);
        assert!(err_ten.yy().abs() < tol);
        assert!(err_ten.yz().abs() < tol);
        assert!(err_ten.zz().abs() < tol);
      }
    }
  }
  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_UnitSpherePoint(){
  
    let theta = PI/6.0;
    let phi = PI/3.0;

    let point = UnitSpherePoint::new(theta,phi);

    let tol = 1e-12; 
    assert!( (point.theta() - theta).abs() < tol );
    assert!( (point.phi() - phi).abs() < tol );
    
    assert!( (point.cos_theta() - 0.8660254037844).abs() < tol );
    assert!( (point.cos_theta_squared() - 0.75).abs() < tol );
    assert!( (point.sin_theta() - 0.5).abs() < tol );
    assert!( (point.sin_theta_squared() - 0.25).abs() < tol );

    assert!( (point.sin_phi() - 0.8660254037844).abs() < tol );
    assert!( (point.sin_phi_squared() - 0.75).abs() < tol );
    assert!( (point.cos_phi() - 0.5).abs() < tol );
    assert!( (point.cos_phi_squared() - 0.25).abs() < tol );
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
