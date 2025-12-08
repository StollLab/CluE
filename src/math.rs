use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::{Hash,Hasher};
use std::collections::hash_map::DefaultHasher;

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
  #[test]
  fn test_ceil(){
    assert_eq!(ceil(2.0),2.0);
    assert_eq!(ceil(2.1),3.0);
    assert_eq!(ceil(-2.1),-2.0);
    assert_eq!(ceil(-2.0),-2.0);
  }
}
