use crate::Config;
use crate::config::pulse_sequence::PulseSequence;
use crate::signal::Signal;
use crate::HamiltonianTensors;
use crate::CluEError;
use crate::physical_constants::{ONE,PI};

use num_complex::Complex;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub fn appa(
    vertices: &[usize],
    tensors: &HamiltonianTensors, config: &Config) 
  -> Result<Option<Signal>,CluEError>
{
  let Some(pulse_sequence_name) = &config.pulse_sequence else{
    return Err(CluEError::NoPulseSequence);
  };

  match pulse_sequence_name{
      PulseSequence::CarrPurcell(1) => appa_hahn(vertices,tensors,config),
      _ => Err(
          CluEError::PulseSequenceNotSupported(
            "pulse sequence not supproted".to_string(),
            "".to_string(),
          )),
  }  
}
/// This function applies Witzel Das Sarma's analytic solution to calculate
/// the 2-CCE auxiliary signal the specified cluster.
pub fn appa_hahn(
    vertices: &[usize],
    tensors: &HamiltonianTensors, config: &Config) 
  -> Result<Option<Signal>,CluEError>
{
  if vertices.len() != 2 {
    return Err(CluEError::WrongClusterSizeForAnalyticCCE(vertices.len()));
  }
  let idx0 = vertices[0];
  let idx1 = vertices[1];

  let Some(dipdip) = &tensors.spin2_tensors.get(idx0,idx1) else {
    return Ok(Some(Signal::new()));
  };
  let Some(hf0) =  &tensors.spin2_tensors.get(0,idx0) else {
    return Ok(Some(Signal::new()));
  }; 
  let Some(hf1) =  &tensors.spin2_tensors.get(0,idx1) else {
    return Ok(Some(Signal::new()));
  }; 

  // For a point dipole b_zz = -(b_xx + b_yy), but
  // -(b_xx + b_yy) account for isotropic coupling as well.
  let b = -(dipdip.xx() + dipdip.yy() );
  let delta_hf = (hf0.zz() - hf1.zz()).abs();
  
  let omega = 2.0*PI*appa_hahn_frequency(delta_hf,b);
  let k = appa_hahn_modulation_depth(delta_hf,b);

  let tau_axis = config.get_tau_axis_as_ref()?;
  let nt = tau_axis.len();
  let mut data = Vec::<Complex<f64>>::with_capacity(nt);

  for t in tau_axis.iter(){
    let twotau = 2.0*(*t);
    let s4 = (omega*twotau).sin().powi(4);
    data.push(ONE - k*s4 );
  }

  Ok(Some(Signal{data}))

}
//------------------------------------------------------------------------------
/// This function calculate the modulation depth for a single 2-cluster.
pub fn appa_hahn_modulation_depth(delta_hf: f64,b:f64) -> f64{
  (2.0*b*delta_hf/( delta_hf*delta_hf + b*b) ).powi(2)
}
//------------------------------------------------------------------------------
/// This function return the modulation frequency to the three spin Hahn echo
/// in unit matching the inputs, which are assumed to have the same units as
/// each other.
pub fn appa_hahn_frequency(delta_hf: f64,b:f64) -> f64{
  (delta_hf*delta_hf + b*b).sqrt()/8.0
}
//------------------------------------------------------------------------------
/// This function calculates the 4th order Taylor series coefficient 
/// for a single 2-cluster.
pub fn appa_hahn_fourth_order_coefficient(delta_hf: f64,b:f64) -> f64{
  let k = appa_hahn_modulation_depth(delta_hf,b);
  let omega = 2.0*PI*appa_hahn_frequency(delta_hf,b);
  k*(omega.powi(4))
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#[cfg(test)]
mod tests{
  use super::*;
  use crate::physical_constants::SQRT2;

  //----------------------------------------------------------------------------
  #[test]
  fn test_appa_hahn_frequency(){
    let delta_hfs = [1e-2, 1e-1, 1.0, 1e1, 1e2];

    for &delta_hf in  delta_hfs.iter(){
      assert!( 
         (appa_hahn_frequency(delta_hf,delta_hf)- 0.125*SQRT2*delta_hf).abs()
          < 1e-12);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_appa_hahn_modulation_depth(){

    let delta_hfs = [1e-2, 1e-1, 1.0, 1e1, 1e2];
    let bs = [1e-2, 1e-1, 1.0, 1e1, 1e2];

    for &delta_hf in  delta_hfs.iter(){
      assert_eq!(appa_hahn_modulation_depth(delta_hf,delta_hf), 1.0);

      for &b in bs.iter(){
        assert!(appa_hahn_modulation_depth(delta_hf,b) <= 1.0);
      }
    }
  }
  //----------------------------------------------------------------------------
}







