use crate::CluEError;

use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::{Hash,Hasher};
use std::collections::hash_map::DefaultHasher;

use num_complex::Complex64;
use ndarray::Array2;

type CxMat = Array2::<Complex64>;

/// This function compares two lists for equality.
pub fn are_vecs_equal<T: std::cmp::PartialEq>
(vec0: &[T], vec1: &[T])-> bool{
  if vec0.len() != vec1.len() {return false;}

  for ii in 0..vec0.len(){
    if vec0[ii] != vec1[ii]{ return false;}
  }
  true
}
//------------------------------------------------------------------------------
/// This function rounds a floating point number up to an interger. 
pub fn ceil(x: f64) -> f64{

  let mut a = x as i32;

  let err = (x-a as f64).abs();
  if err > 1e-12 && x >= 0.0{
    a += 1;
  }

  a as f64
}
//------------------------------------------------------------------------------
/// This function sorts a vector and removes duplicate entries.
pub fn unique<T>(vec: Vec::<T>) -> Vec::<T>
where T: std::cmp::Eq  + std::hash::Hash + std::cmp::Ord
{
  let mut out = vec.into_iter().collect::<HashSet<_>>()
      .into_iter()
      .collect::<Vec::<T>>();

  out.sort();

  out
}
//------------------------------------------------------------------------------
/// This function converts `input` to a `u64` via a `String` intermediary.
pub fn str_hash<T: Debug>(input: &T) -> u64 {
  let in_str = format!("{:?}",input);
  let mut out = DefaultHasher::new();
  in_str.hash(&mut out);
  out.finish()
}
//------------------------------------------------------------------------------
/*   
fn mat_pow_n(m: &Array2::<f64>, n: usize) -> Array2::<f64>
{
  let dim = m.dim();
  assert_eq!(dim.0,dim.1);

  let mut a = Array2::<f64>::eye(dim.0);
  for b in format!("{n:b}").chars(){
    if b == '1'{
      a = a.dot(m);
    }
  }
  a
}
*/
//------------------------------------------------------------------------------
/// This function calulates the nth power of a square complex matrix.
/// The function will panic if the matrix is not square.
pub fn cxmat_pow_n(m: &CxMat, n: usize) -> CxMat
{
  let dim = m.dim();
  assert_eq!(dim.0,dim.1);

  let mut a = CxMat::eye(dim.0);
  for b in format!("{n:b}").chars(){
    a = a.dot(&a);
    if b == '1'{
      a = a.dot(m);
    }
  }
  a
}
//------------------------------------------------------------------------------
pub fn commutator(mat0: &CxMat, mat1: &CxMat) -> CxMat{
  mat0.dot(mat1) - mat1.dot(mat0)
}
//------------------------------------------------------------------------------
pub fn hilbert_schmidt(a: &CxMat, b: &CxMat) -> Complex64
{
  let a_star = a.map(|a_ij| a_ij.conj() );
  let it = std::iter::zip(a_star,b);
  it.map(|(a_ij,b_ij)| a_ij*b_ij).sum::<Complex64>()
}
//------------------------------------------------------------------------------
// TODO: Rethink format and then add error checks.
pub fn expectation_value(op: &CxMat, v: &CxMat) -> Result<Complex64,CluEError> 
{
  let z = vectran_op_vec(v, op, v)?;
  Ok(z[[0,0]])
  
}
//------------------------------------------------------------------------------
// TODO: Rethink format and then add error checks.
pub fn vectran_op_vec(v0: &CxMat,op: &CxMat, v1: &CxMat) 
    -> Result<CxMat,CluEError>
{

  let op_ket = op.dot(v1);
  let bra = v0.t().map(|u| u.conj() );
  let z = bra.dot(&op_ket);
  Ok(z) 
}
//------------------------------------------------------------------------------
// TODO: pending
/*
pub fn vec_vec_transpose<T: Clone>(mat: &Vec::<Vec::<T>>) 
    -> Result<Vec::<Vec::<T>>,CluEError>
{
  if mat.is_empty(){
    return Ok(mat);
  }

  let n_col = mat[0].len();

  for v in mat.iter().skip(1){
    if mat.len() != n_col{
      return Err(CluEError)
    }
  }

  let n_row = mat.len();

  let mut mat_t: Vec::<Vec::<T>> = (0..n_col)
      .map(|_| Vec::<T>::with_capacity(n_row)).collect();

  for row in 0..n_row{
    for col in 0..n_col{
      let value = mat[row][col].clone();
      mat_t[col].push(value);
    }
  }

  Ok(mat_t)
}
*/
//------------------------------------------------------------------------------


#[cfg(test)]
mod tests{
  use super::*;
  use crate::physical_constants::{ONE,ZERO};
  use crate::quantum::cluster_operators::{
    spin_x,
    spin_y,
    spin_z,
  };
  use ndarray::array;

  #[test]
  fn test_are_vecs_equal(){
    let a = vec![1,2,3];
    let b = vec![1,2,3];
    let c = vec![1,2,3,4];

    assert!(are_vecs_equal(&a,&b));
    assert!(!are_vecs_equal(&a,&c));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_unique(){
    let a = vec![1,1,2,3,2];
    let a = unique(a);
    assert_eq!(a,vec![1,2,3]);
  }
  //----------------------------------------------------------------------------
  /*
  #[test]
  fn test_mat_pow_n(){
  
    let m = array![[1.0,1.0], [1.0,0.0]] ;
    let a = mat_pow_n(&m, 13);
    let r = array![[377.0,233.0], [233.0,144.0]] ;
    for row in 0..2{for col in 0..2{  
      assert!( (a[[row,col]]-r[[row,col]]).abs() < 1e-12  )
    }}
  }
  */
  //----------------------------------------------------------------------------
  #[test]
  fn test_cxmat_pow_n(){
  
    let m = array![[ONE,ONE], [ONE,ZERO]] ;
    let a = cxmat_pow_n(&m, 13);
    let r = array![[377.0*ONE,233.0*ONE], [233.0*ONE,144.0*ONE]] ;
    for row in 0..2{for col in 0..2{  
      assert!( (a[[row,col]]-r[[row,col]]).norm() < 1e-12  )
    }}
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_ceil(){
    assert_eq!(ceil(2.0),2.0);
    assert_eq!(ceil(2.1),3.0);
    assert_eq!(ceil(-2.1),-2.0);
    assert_eq!(ceil(-2.0),-2.0);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_hilbert_schmidt(){
    let e = CxMat::eye(2);
    let x = spin_x(2);
    let y = spin_y(2);
    let z = spin_z(2);
    let s2 = Complex64{re:0.5, im:0.0};
    assert!( (hilbert_schmidt(&e,&e) 
          - Complex64{re: 2.0,im: 0.0}).norm() < 1e-12);  
    assert!( (hilbert_schmidt(&x,&x) - s2).norm() < 1e-12);  
    assert!( (hilbert_schmidt(&y,&y) - s2).norm() < 1e-12);  
    assert!( (hilbert_schmidt(&z,&z) - s2).norm() < 1e-12);  
    assert!( (hilbert_schmidt(&x,&x) - s2).norm() < 1e-12);  
    assert!( hilbert_schmidt(&x,&e).norm() < 1e-12);  
    assert!( hilbert_schmidt(&y,&e).norm() < 1e-12);  
    assert!( hilbert_schmidt(&z,&e).norm() < 1e-12);  
    assert!( hilbert_schmidt(&x,&y).norm() < 1e-12);  
    assert!( hilbert_schmidt(&y,&x).norm() < 1e-12);  
    assert!( hilbert_schmidt(&x,&z).norm() < 1e-12);  
    assert!( hilbert_schmidt(&z,&x).norm() < 1e-12);  
    assert!( hilbert_schmidt(&y,&z).norm() < 1e-12);  
    assert!( hilbert_schmidt(&z,&y).norm() < 1e-12);  
  }
  //----------------------------------------------------------------------------

}
