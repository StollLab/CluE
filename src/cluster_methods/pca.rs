use crate::Config;
use crate::signal::Signal;
use crate::HamiltonianTensors;
use crate::CluEError;

use num_complex::Complex;

/// TODO: under construction
pub fn pca(
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

  let h2 = 0.25*(3.0*b*b + 0.25*delta_hf*delta_hf + 2.0*delta_hf*b);
  
  let omega = h2.sqrt()*0.25;
  let k = 0.5*(b*delta_hf/h2).powi(2);

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
