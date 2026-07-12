use crate::Config;
use crate::config::pulse_sequence::PulseSequence;
use crate::signal::Signal;
use crate::HamiltonianTensors;
use crate::CluEError;
use crate::physical_constants::PI;

use num_complex::Complex;


pub fn lce(
    vertices: &[usize],
    tensors: &HamiltonianTensors, config: &Config)
  -> Result<Option<Signal>,CluEError>
{
  let Some(pulse_sequence_name) = &config.pulse_sequence else{
    return Err(CluEError::NoPulseSequence);
  };

  match pulse_sequence_name{
      PulseSequence::CarrPurcell(1) 
          => lce_hahn_2spin_2order(vertices,tensors,config),
      _ => Err(
          CluEError::PulseSequenceNotSupported(
            "pulse sequence not supproted".to_string(),
            "".to_string(),
          )),
  }
}

/// This function calculates the 2-spin 2nd order LCE signal
/// v(2τ) = exp(D2mn),
/// where `vertices = [m,n]`, and in energy units
/// D_2mn(2\tau) = -(2b_mn/ΔA_mn)^2 sin^4( ΔA_mn/8ħ 2τ).
/// This function assumes that `tensors` is in frequency units
pub fn lce_hahn_2spin_2order(
    vertices: &[usize],
    tensors: &HamiltonianTensors, 
    config: &Config
    ) ->Result<Option<Signal>,CluEError>
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
  
  let omega = 2.0*PI*delta_hf/8.0;
  let k = (2.0*b/delta_hf).powi(2);

  let tau_axis = config.get_tau_axis_as_ref()?;
  let nt = tau_axis.len();
  let mut data = Vec::<Complex<f64>>::with_capacity(nt);

  for t in tau_axis.iter(){
    let twotau = 2.0*(*t);
    let s4 = (omega*twotau).sin().powi(4);
    data.push( Complex::<f64>{ re: (-k*s4).exp(), im: 0.0} );
  }

  Ok(Some(Signal{data}))
}
