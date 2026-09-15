use crate::clue_errors::CluEError;

use num_complex::Complex64;
use ndarray::Array2;
type CxMat = Array2::<Complex64>;
type Mat = Array2::<f64>;

//------------------------------------------------------------------------------
pub fn eq_variant<T>(a: &T, b: &T) -> bool
{
  std::mem::discriminant(a) == std::mem::discriminant(b)
}
//------------------------------------------------------------------------------
pub fn are_all_same_type<T>(array: &[T]) 
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
  let mat = vec_vec_complex_f64_from_toml_array(array)?;

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
// [x,y,z], [[x],[y],[z]] -> [[x],[y],[z]]
// [[x,y,z]]] -> [[x,y,z]]] 
// [[ [a,b], [c,d] ],[ [e,f], [g,h]]] -> [ [a+bi, c,di], [e+fi, g+hi] ]
pub fn vec_vec_complex_f64_from_toml_array(array: Vec::<toml::Value>) 
    -> Result<Vec::<Vec::<Complex64>>,CluEError>
{

  let mut out = Vec::<Vec::<Complex64>>::with_capacity(array.len());

  for row in array{
    match row{
      toml::Value::Array(arr) => out.push(vec_complex_f64_from_toml_array(arr)?),
      toml::Value::Float(x) => out.push(vec![Complex64{re: x, im:0.0}]),
      _ => return Err(CluEError::ExpectedTOMLArray(row.type_str().to_string())),
    }
  }
  Ok(out)
}    
//------------------------------------------------------------------------------
pub fn vec_cxmat_from_toml_array(array: Vec::<toml::Value>) 
    -> Result<Vec::<CxMat>,CluEError>
{

  let mut out = Vec::<CxMat>::with_capacity(array.len());

  for row in array{
    match row{
      toml::Value::Array(arr) => {
        out.push(cxmat_from_toml_array(arr)?)
      },
      _ => return Err(CluEError::ExpectedTOMLArray(row.type_str().to_string())),
    }
  }
  Ok(out)
}    
//------------------------------------------------------------------------------
pub fn vec_complex_f64_from_toml_array(array: Vec::<toml::Value>) 
    -> Result<Vec::<Complex64>,CluEError>
{
  let mut out = Vec::<Complex64>::with_capacity(array.len());
  for el in array{
    match el{
      toml::Value::Float(x) => out.push(Complex64{re: x, im:0.0}),
      toml::Value::Array(arr) => out.push(complex_f64_from_toml_array(arr)?),
      _ => return Err(CluEError::ExpectedTOMLFloat(el.type_str().to_string())),
    }
  }
  Ok(out)
}    
//------------------------------------------------------------------------------
pub fn complex_f64_from_toml_array(array: Vec::<toml::Value>) 
    -> Result<Complex64,CluEError>
{
  if array.len() > 2{
    return Err(CluEError::TOMLArrayIsNotAComplexNumber);
  }
  let mut z = Vec::<f64>::from([0.0, 0.0]);
  
  for (ii,el) in array.iter().enumerate(){
    match el{
      toml::Value::Float(x) => z[ii] +=*x,
      _ => return Err(CluEError::ExpectedTOMLArray(el.type_str().to_string())),
    }
  }
  Ok(Complex64{re: z[0], im: z[1]})
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
      toml::Value::Float(x) => out.push(*x),
      _ => return Err(CluEError::ExpectedTOMLFloat(el.type_str().to_string())),
    }
  }
  Ok(out)
}    
//------------------------------------------------------------------------------
pub fn vec_usize_from_toml_array(array: Vec::<toml::Value>) 
    -> Result<Vec::<usize>,CluEError>
{
  let mut out = Vec::<usize>::with_capacity(array.len());
  for el in array.iter(){
    match el{
      toml::Value::Integer(x) => {
        if *x < 0 {
          return Err(CluEError::ExpectedTOMLUInt(el.type_str().to_string()));
        }
        out.push(*x as usize)
      },
      _ => return Err(CluEError::ExpectedTOMLUInt(el.type_str().to_string())),
    }
  }
  Ok(out)
}    
//------------------------------------------------------------------------------

#[cfg(test)]
mod tests{
  use super::*;

  use ndarray::arr2;

  //----------------------------------------------------------------------------
  #[test]
  fn test_cxmat_from_toml_array(){
    let table = get_test_table();
    
    let toml::Value::Array(matrix) = &table["real_matrix"]else{
      panic!("failure");
    };

    let result = cxmat_from_toml_array(matrix.clone()).unwrap();
    assert_eq!(result, arr2(&[
          [Complex64{re: 1.0, im: 0.0} ,Complex64{re: 2.0,im: 0.0}],
          [Complex64{re: 3.0, im: 0.0}, Complex64{re: 4.0, im: 0.0}]
    ]));
    
    let toml::Value::Array(matrix) = &table["cx_matrix"]else{
      panic!("failure");
    };
    let result = cxmat_from_toml_array(matrix.clone()).unwrap();
    assert_eq!(result, arr2(&[
        [Complex64{re: 1.0, im: 0.0} ,Complex64{re: 2.0, im: 1.0} ],
        [Complex64{re: 3.0, im: 2.0} ,Complex64{re: 4.0, im: 3.0}]
    ]));
  }  
  //----------------------------------------------------------------------------
  #[test]
  fn test_vec_vec_complex_f64_from_toml_array(){
    let table = get_test_table();
    let toml::Value::Array(matrix) = &table["cx_matrix"]else{
      panic!("failure");
    };
    let result = vec_vec_complex_f64_from_toml_array(matrix.clone()).unwrap();
    assert_eq!(result, vec![
        vec![Complex64{re: 1.0, im: 0.0} ,Complex64{re: 2.0, im: 1.0} ],
        vec![Complex64{re: 3.0, im: 2.0} ,Complex64{re: 4.0, im: 3.0}]
    ]);
  }  
  //----------------------------------------------------------------------------
  #[test]
  fn test_vec_complex_f64_from_toml_array(){
    let table = get_test_table();
    let toml::Value::Array(z) = &table["cx_vector"]else{
      panic!("failure");
    };
    let result = vec_complex_f64_from_toml_array(z.clone()).unwrap();
    assert_eq!(result, vec![ Complex64{re: 1.0, im: 0.0}, Complex64{ re: 2.0, im: 1.0}] );

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_complex_f64_from_toml_array(){ 
    let table = get_test_table();
    let toml::Value::Array(z) = &table["cx_number"]else{
      panic!("failure");
    };
    let result = complex_f64_from_toml_array(z.clone()).unwrap(); 
    assert_eq!(result, Complex64{ re: 1.0, im: 2.0} );
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_mat_from_toml_array(){
    let table = get_test_table();
    let toml::Value::Array(real_matrix) = &table["real_matrix"]else{
      panic!("failure");
    };
    let result = mat_from_toml_array(real_matrix.clone()).unwrap();
    assert_eq!(result, arr2(&[[1.0,2.0],[3.0,4.0]]));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_vec_vec_f64_from_toml_array(){
    let table = get_test_table();
    let toml::Value::Array(real_matrix) = &table["real_matrix"]else{
      panic!("failure");
    };
    let result = vec_vec_f64_from_toml_array(real_matrix.clone()).unwrap();
    assert_eq!(result, vec![vec![1.0,2.0],vec![3.0,4.0]]);
  }  
  //----------------------------------------------------------------------------
  #[test]
  fn test_vec_f64_from_toml_array(){
    let table = get_test_table();
    let toml::Value::Array(vector) = &table["real_vector"]else{
      panic!("failure");
    };
    let result = vec_f64_from_toml_array(vector.clone()).unwrap();
    assert_eq!(result, vec![1.0,2.0]);
  }  
  //----------------------------------------------------------------------------
  fn get_test_table() -> toml::Table{
    r##"
      real_number = 1.0
      real_vector = [1.0,2.0]
      real_matrix = [
        [1.0,2.0],
        [3.0,4.0]
      ]
      cx_number = [1.0,2.0]  
      cx_vector = [1.0 , [2.0,1.0]]
      cx_matrix = [ 
        [ [1.0, 0.0], [2.0, 1.0] ],
        [ [3.0,2.0] , [4.0,3.0] ]
      ]
    "##.parse::<toml::Table>().unwrap()
  }
  //----------------------------------------------------------------------------
}
