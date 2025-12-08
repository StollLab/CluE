use pyo3::prelude::*;

use numpy::PyArray;

use clue_oxide::quantum::spin_hamiltonian as clue_spin;
use clue_oxide::{
  clue_errors::CluEError,
  physical_constants::ELEMENTARY_CHARGE,
};

use ndarray::{Array2,Ix2};
use num_complex::Complex;

use crate::py_clue_errors::PyCluEError;

type CxMat = Array2::<Complex<f64>>;
//------------------------------------------------------------------------------
#[pyfunction]
pub fn spin_x(spin_multiplicity: usize)
  -> Py<PyArray<Complex::<f64>,Ix2>>
{
  let s = clue_spin::spin_x(spin_multiplicity);
  Python::with_gil(|py|{
      PyArray::from_owned_array_bound(py, s).unbind() 
  })
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn spin_y(spin_multiplicity: usize)
  -> Py<PyArray<Complex::<f64>,Ix2>>
{
  let s = clue_spin::spin_y(spin_multiplicity);

  Python::with_gil(|py|{
      PyArray::from_owned_array_bound(py, s).unbind() 
  })
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn spin_z(spin_multiplicity: usize)
  -> Py<PyArray<Complex::<f64>,Ix2>>
{
  let s = clue_spin::spin_z(spin_multiplicity);

  Python::with_gil(|py|{
      PyArray::from_owned_array_bound(py, s).unbind() 
  })
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn spin_plus(spin_multiplicity: usize)
  -> Py<PyArray<Complex::<f64>,Ix2>>
{
  let s = clue_spin::spin_plus(spin_multiplicity);

  Python::with_gil(|py|{
      PyArray::from_owned_array_bound(py, s).unbind() 
  })
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn spin_minus(spin_multiplicity: usize)
  -> Py<PyArray<Complex::<f64>,Ix2>>
{
  let s = clue_spin::spin_minus(spin_multiplicity);

  Python::with_gil(|py|{
      PyArray::from_owned_array_bound(py, s).unbind() 
  })
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn spin_squared(spin_multiplicity: usize) 
  -> Py<PyArray<Complex::<f64>,Ix2>>
{
  let s = clue_spin::spin_squared(spin_multiplicity);

  Python::with_gil(|py|{
      PyArray::from_owned_array_bound(py, s).unbind() 
  })
}
//------------------------------------------------------------------------------
#[pyfunction]
pub fn nuclear_quadrupole_tensor_operator(
    spin_multiplicity: usize, 
    quadrupole_moment: f64
    ) -> Result<Vec::<Vec::<Vec::<Complex<f64>>>>,PyCluEError>
{
  let mut tensor = Vec::<Vec::<Vec::<Complex<f64>>>>::with_capacity(6);
  for ii in 0..3{
    for jj in ii..3{
      let q = get_nuclear_quadrupole_tensor_component(
          spin_multiplicity,quadrupole_moment,ii,jj)?;

      tensor.push(cx_matrix_to_vec_vec(q))

    }
  }
  Ok(tensor)
}
//------------------------------------------------------------------------------
fn get_nuclear_quadrupole_tensor_component(
    spin_multiplicity: usize,quadrupole_moment:f64,
    ax0: usize,ax1: usize)
  -> Result<CxMat,CluEError>
{
  let spin_op0 = axis_to_spin_operator(spin_multiplicity, ax0)?;
  let spin_op1 = axis_to_spin_operator(spin_multiplicity, ax1)?;

 let mut qij = (spin_op0.clone()*spin_op1.clone() + spin_op1*spin_op0)*1.5;
 if ax0 == ax1{
   qij = qij - clue_spin::spin_squared(spin_multiplicity);
 }

 let spin = (spin_multiplicity as f64 - 1.0)/2.0;
 qij = qij*quadrupole_moment*ELEMENTARY_CHARGE/ (spin*(2.0*spin - 1.0));

 Ok(qij)
}  
//------------------------------------------------------------------------------
fn axis_to_spin_operator(spin_multiplicity: usize,ax: usize) 
  -> Result<CxMat,CluEError>
{
  match ax{
    0 => Ok(clue_spin::spin_x(spin_multiplicity)),
    1 => Ok(clue_spin::spin_y(spin_multiplicity)),
    2 => Ok(clue_spin::spin_z(spin_multiplicity)),
    _ => Err(CluEError::InvalidAxes)
  }
}
//------------------------------------------------------------------------------
fn cx_matrix_to_vec_vec(cx_matrix: CxMat) -> Vec::<Vec<Complex<f64>>>{


  let n_rows = cx_matrix.len_of(ndarray::Axis(0));
  let n_cols = cx_matrix.len_of(ndarray::Axis(1));
  let mut out = Vec::<Vec<Complex<f64>>>::with_capacity(n_rows);
  {
    let zeros = (0..n_cols).map(|_| Complex::<f64>{re: 0.0,im:0.0})
      .collect::<Vec<Complex<f64>>>();
    for _ii in 0..n_rows{
      out.push(zeros.clone());
    }
  }

  for row in 0..n_rows{
    for col in 0..n_cols{
      out[row][col] = cx_matrix[[row,col]];
    }
  }

  out

}
