use pyo3::prelude::*;

use serde::{Serialize,Deserialize};

use std::collections::HashMap;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass(name = "Value")]
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
pub enum PyValue{
  String(String),
  Integer(i64),
  Float(f64),
  Boolean(bool),
  Array(PyArray),
  Table(PyTable),
}

impl From<toml::Value> for PyValue {
    fn from(item: toml::Value) -> Self {
        match item{
          toml::Value::String(s) => Self::String(s),
          toml::Value::Integer(i) => Self::Integer(i),
          toml::Value::Float(f) => Self::Float(f),
          toml::Value::Boolean(b) => Self::Boolean(b),
          toml::Value::Datetime(_) => panic!("Cannot convert Datetime."),
          toml::Value::Array(array) => Self::Array(PyArray::from(array)),
          toml::Value::Table(table) => Self::Table(PyTable::from(table)),
        }
    }
}
impl From<i32> for PyValue {
    fn from(item: i32) -> Self {
      Self::Integer(item as i64)
    }
}
impl From<i64> for PyValue {
    fn from(item: i64) -> Self {
      Self::Integer(item)
    }
}
impl From<f64> for PyValue {
    fn from(item: f64) -> Self {
      Self::Float(item)
    }
}
impl From<bool> for PyValue {
    fn from(item: bool) -> Self {
      Self::Boolean(item)
    }
}

#[pymethods]
impl PyValue{
  //----------------------------------------------------------------------------
  pub fn db_print(&self){
    println!("{:?}",self);
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass(name = "Array")]
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
pub struct PyArray{
  array: Vec::<PyValue>
}
#[pymethods]
impl PyArray{
  //----------------------------------------------------------------------------
  pub fn __getitem__(&self,index: usize) -> PyValue{
    self.array[index].clone()
  }
  //----------------------------------------------------------------------------
  pub fn __setitem__(&mut self, index: usize, value: PyValue){
    self.array[index] = value;
  }
  //----------------------------------------------------------------------------
  pub fn db_print(&self){
    println!("{:?}",self);
  }
  //----------------------------------------------------------------------------
  #[staticmethod]
  pub fn with_capacity(capacity: usize) -> Self{
    let array = Vec::<PyValue>::with_capacity(capacity);
    Self{array}
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn append(&mut self, value: PyValue){
    self.array.push(value);
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn push(&mut self, value: PyValue){
    self.array.push(value);
  }
  //----------------------------------------------------------------------------
}

impl From<Vec::<toml::Value>> for PyArray {
    fn from(item: Vec::<toml::Value>) -> Self {
      let mut out = PyArray::with_capacity(item.len());
      for value in item.iter(){
        let py_value = PyValue::from(value.clone());
        out.push(py_value);
      }
      out
    }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[pyclass(name = "Table")]
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
pub struct PyTable{
  table: HashMap<String, PyValue>,
}
#[pymethods]
impl PyTable{
  //----------------------------------------------------------------------------
  pub fn __getitem__(&self,key: String) -> PyValue{
    self.table[&key].clone()
  }
  //----------------------------------------------------------------------------
  pub fn __setitem__(&mut self,key: String, value: PyValue){
    self.table.insert(key,value);
  }
  //----------------------------------------------------------------------------
  pub fn db_print(&self){
    println!("{:?}",self);
  }
}
impl PyTable{
  //----------------------------------------------------------------------------
  pub fn with_capacity(capacity: usize) -> Self{
    let table = HashMap::<String, PyValue>::with_capacity(capacity);
    Self{table}
  }
  //----------------------------------------------------------------------------
  pub fn insert<T>(&mut self, key: String,value: T)
    where PyValue: From<T>
  {
    let py_value: PyValue = value.into();
    self.table.insert(key,py_value);
  }
}

impl From<toml::Table> for PyTable {
    fn from(item: toml::Table) -> Self {
      let mut out = PyTable::with_capacity(item.len());
      for (key,value) in item.iter(){
        let py_value = PyValue::from(value.clone());
        out.insert(key.to_string(),py_value);
      }
      out
    }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/*
impl PyValue{
  fn from_toml_string(toml_str: &str) -> Result<Self,PyCluEError>{
    let decoded: Result<ConfigTOML,_> = toml::from_str(toml_str);
    match decoded {
      Ok(mut config) => {
        config.set_default_units();
        Ok(config)
      },
      //Err(err) => Err(CluEError::CannotReadConfigTOML( format!("{}",err) )),
      Err(err) => panic!("TODO: implement error: {}.",err)
    }
  }
}
*/

