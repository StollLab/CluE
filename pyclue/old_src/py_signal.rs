use pyo3::prelude::*;

use clue_oxide::signal::Signal;
use clue_oxide::signal::load_csv_to_vec_signals;

use crate::py_clue_errors::PyCluEError;

use num_complex::Complex;
use std::ops::{Add,Sub,Mul,Div};

#[pyclass(name = "Signal")]
pub struct PySignal{
  signal: Signal,
}

#[pymethods]
impl PySignal{
  //#[new]
  //#[pyo3(signature = (data) )]

  //----------------------------------------------------------------------------
  #[staticmethod]
  pub fn read_from_csv(filename: String) -> Result<Self,PyCluEError>{

    let signal = Signal::read_from_csv(&filename)?;
  
    Ok(PySignal{signal})
  }
  //----------------------------------------------------------------------------
  pub fn get_data(&self) -> Vec::<Complex<f64>>{
    self.signal.data.clone()
  }
  //----------------------------------------------------------------------------
  #[staticmethod]
  pub fn ones(n: usize) -> Self{
    Self{signal: Signal::ones(n) }
  }
  //----------------------------------------------------------------------------
  #[staticmethod]
  pub fn zeros(n: usize) -> Self{
    Self{signal: Signal::zeros(n) }
  }
  //----------------------------------------------------------------------------
  pub fn write_to_csv(&self, filename: String) -> Result<(),PyCluEError>
  {
    Ok(self.signal.write_to_csv(&filename)?)
  }
  //----------------------------------------------------------------------------
  fn len(&self) -> usize{
    self.signal.len()
  }
  //----------------------------------------------------------------------------
  fn mut_scale(&mut self, value: Complex::<f64>){
    self.signal.mut_scale(value)
  }
  //----------------------------------------------------------------------------
  fn scale(&self, value: Complex::<f64>) -> Self{
    Self{signal: self.signal.scale(value)}
  }
  //----------------------------------------------------------------------------
  fn rmsd_with(&self, other: &PySignal) -> f64{
    
    let n = self.len(); 
    assert_eq!(n, other.len());

    let mut rmsd = 0.0;
    for (ii, z) in self.signal.data.iter().enumerate(){
      rmsd += (z - other.signal.data[ii]).norm_sqr();
    }

    rmsd /= n as f64;

    rmsd.sqrt()
  }
  //----------------------------------------------------------------------------
  fn __add__(&self, rhs: &PySignal) -> PySignal {
    self + rhs
  }
  //----------------------------------------------------------------------------
  fn __sub__(&self, rhs: &PySignal) -> PySignal {
    self - rhs
  }
  //----------------------------------------------------------------------------
  fn __mul__(&self, rhs: &PySignal) -> PySignal {
    self * rhs
  }
  //----------------------------------------------------------------------------
  fn __truediv__(&self, rhs: &PySignal) -> PySignal {
    self / rhs
  }
  //----------------------------------------------------------------------------
  fn __str__(&self) -> String{
    self.signal.to_string()
  }
  //----------------------------------------------------------------------------
  fn __repr__(&self) -> String{
    self.signal.to_string()
  }
  //----------------------------------------------------------------------------
  fn __len__(&self) -> usize{
    self.len()
  }
  //----------------------------------------------------------------------------
  fn __getitem__(&self, index: usize) -> Complex::<f64>{
    self.signal.data[index]
  }
  //----------------------------------------------------------------------------
  fn __setitem__(&mut self, index: usize, value: Complex::<f64>){
    self.signal.data[index] = value;
  }
  //----------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
impl Add for &PySignal{
  type Output = PySignal;

  fn add(self,rhs: &PySignal) -> PySignal{
    let signal = &self.signal + &rhs.signal;
    PySignal{signal}
  }
}
//------------------------------------------------------------------------------
impl Sub for &PySignal{
  type Output = PySignal;

  fn sub(self,rhs: &PySignal) -> PySignal{
    let signal = &self.signal - &rhs.signal;
    PySignal{signal}
  }
}
//------------------------------------------------------------------------------
impl Mul for &PySignal{
  type Output = PySignal;

  fn mul(self,rhs: &PySignal) -> PySignal{
    let signal = &self.signal * &rhs.signal;
    PySignal{signal}
  }
}
//------------------------------------------------------------------------------
impl Div for &PySignal{
  type Output = PySignal;

  fn div(self,rhs: &PySignal) -> PySignal{
    let signal = &self.signal / &rhs.signal;
    PySignal{signal}
  }
}
//------------------------------------------------------------------------------

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass(name = "Signals")]
pub struct PySignals{
  signals: Vec::<Signal>,
}

#[pymethods]
impl PySignals{
  #[staticmethod]
  pub fn read_from_csv(filename: String) -> Result<PySignals,PyCluEError>{
    let signals = load_csv_to_vec_signals(&filename)?;

    Ok(PySignals{signals})
  }
  //----------------------------------------------------------------------------
  pub fn get_data(&self) -> Vec< Vec::<Complex<f64>> >{

    let n_sigs = self.signals.len();
    let mut data = Vec::<Vec::<Complex<f64>>>::with_capacity(n_sigs);

    for sig in self.signals.iter(){
      data.push( sig.data.clone() )
    }
    data
  }
  //----------------------------------------------------------------------------
  fn len(&self) -> usize{
    self.signals.len()
  }
  //----------------------------------------------------------------------------
  fn __len__(&self) -> usize{
    self.len()
  }
  //----------------------------------------------------------------------------
  fn __getitem__(&self, index: i32) -> PySignal{
    let index = index.rem_euclid(self.len() as i32) as usize;
    PySignal{ signal: self.signals[index].clone() }
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

/*
#[pyfunction]
pub fn read_signals_file(filename: String) 
    -> Result< Vec< Vec::<Complex<f64>> > ,PyCluEError> 
{
  let signals = PySignals::read_from_csv(filename)?;

  Ok(signals.get_data())
}
*/
