use crate::clue_errors::CluEError;

use toml::Value;
use num_complex::Complex;
use ndarray::Array2;
type CxMat = Array2::<Complex<f64>>;
type Mat = Array2::<f64>;

//------------------------------------------------------------------------------
pub fn eq_variant<T>(a: &T, b: &T) -> bool
{
  std::mem::discriminant(a) == std::mem::discriminant(b)
}
//------------------------------------------------------------------------------
pub fn are_all_same_type<T>(array: &Vec::<T>) 
  -> bool
{
  if array.len() <= 1 { return true; }

  for value in array.iter(){
    if !eq_variant(value,&array[0]){
      return false;
    } 
  }
  true
}
//------------------------------------------------------------------------------
pub fn cxmat_from_toml_array(array: Vec::<toml::Value>) 
    -> Result<CxMat,CluEError>
{
  let mat = vec_vec_f64_from_toml_array(array)?;

  let n_rows = mat.len();

  let n_cols = mat[0].len();

  for row in mat.iter().skip(1){
    if row.len() != n_cols{
      return Err(CluEError::TOMLArrayIsNotAMatrix);
    }
  }

  let mut out = CxMat::zeros((n_rows,n_cols));
  for (irow, row) in mat.iter().enumerate(){
    for (icol,el) in row.iter().enumerate(){
      out[[irow,icol]] += el;
    }
  } 

  Ok(out)
} 
//------------------------------------------------------------------------------
pub fn mat_from_toml_array(array: Vec::<toml::Value>) 
    -> Result<Mat,CluEError>
{
  let mat = vec_vec_f64_from_toml_array(array)?;

  let n_rows = mat.len();

  let n_cols = mat[0].len();

  for row in mat.iter().skip(1){
    if row.len() != n_cols{
      return Err(CluEError::TOMLArrayIsNotAMatrix);
    }
  }

  let mut out = Mat::zeros((n_rows,n_cols));
  for (irow, row) in mat.iter().enumerate(){
    for (icol,el) in row.iter().enumerate(){
      out[[irow,icol]] += el;
    }
  } 

  Ok(out)
} 
//------------------------------------------------------------------------------
pub fn vec_vec_f64_from_toml_array(array: Vec::<toml::Value>) 
    -> Result<Vec::<Vec::<f64>>,CluEError>
{

  let mut out = Vec::<Vec::<f64>>::with_capacity(array.len());

  for row in array{
    match row{
      toml::Value::Array(arr) => out.push(vec_f64_from_toml_array(arr)?),
      _ => return Err(CluEError::ExpectedTOMLArray(row.type_str().to_string())),
    }
  }
  Ok(out)
}    
//------------------------------------------------------------------------------
pub fn vec_f64_from_toml_array(array: Vec::<toml::Value>) 
    -> Result<Vec::<f64>,CluEError>
{
  let mut out = Vec::<f64>::with_capacity(array.len());
  for el in array.iter(){
    match el{
      toml::Value::Float(x) => out.push(*x as f64),
      _ => return Err(CluEError::ExpectedTOMLFloat(el.type_str().to_string())),
    }
  }
  Ok(out)
}    
//------------------------------------------------------------------------------

