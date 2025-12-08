use pyo3::prelude::*;

use clue_oxide::structure::particle::Particle;
use clue_oxide::elements::Element;
use clue_oxide::isotopes::Isotope;

use crate::py_clue_errors::PyCluEError;

#[pyclass(name = "Particle")]
#[derive(Debug,Clone)]
pub struct PyParticle{
  pub particle: Particle,
}


#[pymethods]
impl PyParticle{
  //----------------------------------------------------------------------------
  pub fn get_element(&self) -> PyElement {
    PyElement{element: self.particle.element.clone() }
  }
  //----------------------------------------------------------------------------
  pub fn get_coordinates(&self) -> [f64;3] {
    self.particle.coordinates.elements.clone()
  }
  //----------------------------------------------------------------------------
  pub fn set_coordinates(&mut self, value: [f64;3]){
    self.particle.coordinates.elements = value;
  }
  //----------------------------------------------------------------------------
  pub fn g_value(&self) -> f64{
    self.particle.isotope.g_value()
  }
  //----------------------------------------------------------------------------
  pub fn spin_multiplicity(&self) -> usize{
    self.particle.isotope.spin_multiplicity()
  }
  //----------------------------------------------------------------------------
}


#[pyclass(name = "Element")]
#[derive(Debug,Clone)]
pub struct PyElement{
  pub element: Element,
}

#[pymethods]
impl PyElement{
  //----------------------------------------------------------------------------
  #[new]
  fn new(element_symbol: String) -> Result<Self,PyCluEError>{
    let element = Element::from(&element_symbol)?;
    Ok(Self{element})
  }
  //----------------------------------------------------------------------------
  fn __str__(&self) -> String{
    self.element.to_string()
  }
  //----------------------------------------------------------------------------
  fn __repr__(&self) -> String{
    self.element.to_string()
  }
  //----------------------------------------------------------------------------
  fn __eq__(&self, other: PyElement) -> bool{
    self.element == other.element
  }
  //----------------------------------------------------------------------------
  
}

#[pyclass(name = "Isotope")]
#[derive(Debug,Clone)]
pub struct PyIsotope{
  pub isotope: Isotope,
}

#[pymethods]
impl PyIsotope{
  //----------------------------------------------------------------------------
  #[new]
  fn new(isotope_symbol: String) -> Result<Self,PyCluEError>{
    let isotope = Isotope::from(&isotope_symbol)?;
    Ok(Self{isotope})
  }
  //----------------------------------------------------------------------------
  pub fn g_value(&self) -> f64{
    self.isotope.g_value()
  }
  //----------------------------------------------------------------------------
  pub fn spin_multiplicity(&self) -> usize{
    self.isotope.spin_multiplicity()
  }
  //----------------------------------------------------------------------------
  fn __str__(&self) -> String{
    self.isotope.to_string()
  }
  //----------------------------------------------------------------------------
  fn __repr__(&self) -> String{
    self.isotope.to_string()
  }
  //----------------------------------------------------------------------------
  
}
