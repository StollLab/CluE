use crate::config::pulse_sequence::PulseStepSpecifier;
use crate::clue_errors::CluEError;
use crate::quantum::cluster_operators::{
  ClusterSpinOperators,
  get_spin_operator,
  SpinOp,
};
use crate::physical_constants::I;

use ndarray::Array2;
use num_complex::Complex;

type CxMat = Array2::<Complex<f64>>;

pub const PI_OVER_2_PULSE_NAME: &str = "pi/2";
pub const PI_PULSE_NAME: &str = "pi";

#[derive(Debug,Clone,PartialEq)]
pub enum PulseStep<'a>{
  Pulse(&'a CxMat),
  FixedDelay(usize, Option<usize>),
  FixedDelay2(usize, Option<usize>),
  InvFixedDelay(usize, Option<usize>),
  InvFixedDelay2(usize, Option<usize>),
  TauDelay,
  Tau2Delay,
  InvTauDelay,
  InvTau2Delay,
  Detect(&'a CxMat),  
}


pub fn get_standard_pulses<'a>(
    spin_operators: &'a ClusterSpinOperators, 
    spin_multiplicity: usize, cluster_size: usize,
    ) 
    -> Result<Vec::<PulseStep<'a>>,CluEError>
{
  let pi_over_2 = spin_operators.get_pulse_operator(PI_OVER_2_PULSE_NAME,
      spin_multiplicity, cluster_size)?;

  let pi = spin_operators.get_pulse_operator(PI_PULSE_NAME,
      spin_multiplicity, cluster_size)?;
  
  let detect = spin_operators.get_detected_operator(
      spin_multiplicity,cluster_size)?;

  Ok(vec![
      PulseStep::Pulse(pi_over_2),
      PulseStep::Pulse(pi),
      PulseStep::Detect(detect),
  ])
}

//------------------------------------------------------------------------------
pub fn generate_pulse_sequence<'a>(
    pulse_sequence_specifier: &[PulseStepSpecifier],
    spin_operators: &'a ClusterSpinOperators, 
    spin_multiplicity: usize, cluster_size: usize,
    )
    -> Result<Vec::<PulseStep<'a>>,CluEError>
{
  let mut pulse_sequence = Vec::<PulseStep<'a>>::with_capacity(
      pulse_sequence_specifier.len());

  for step in pulse_sequence_specifier.iter(){
    let seq_step = match step{
      PulseStepSpecifier::Pulse(pulse_name) => PulseStep::Pulse(
        spin_operators.get_pulse_operator(pulse_name,
        spin_multiplicity, cluster_size)?
      ),
      PulseStepSpecifier::TauDelay => PulseStep::TauDelay,
      PulseStepSpecifier::Tau2Delay => PulseStep::Tau2Delay,
      PulseStepSpecifier::InvTauDelay => PulseStep::InvTauDelay,
      PulseStepSpecifier::InvTau2Delay => PulseStep::InvTau2Delay,
      PulseStepSpecifier::FixedDelay(number,index_opt) =>{
        if *number >= 0{
          PulseStep::FixedDelay(*number as usize, *index_opt )
        }else{
          PulseStep::InvFixedDelay((*number).unsigned_abs() as usize, *index_opt)
        }  
      },
      PulseStepSpecifier::FixedDelay2(number,index_opt) =>{
        if *number >= 0{
          PulseStep::FixedDelay2(*number as usize, *index_opt )
        }else{
          PulseStep::InvFixedDelay2((*number).unsigned_abs() as usize, *index_opt)
        }  
      },
      PulseStepSpecifier::Detect =>PulseStep::Detect(
        spin_operators.get_detected_operator(
            spin_multiplicity,cluster_size)?
      ),  
    };
    pulse_sequence.push(seq_step);
  }
  Ok(pulse_sequence)
}
//------------------------------------------------------------------------------
fn ideal_two_state_pulse(spin_op: &SpinOp, angle: f64)
  -> CxMat
{


  let e = CxMat::eye(2);
  let i_sigma = get_spin_operator(2,spin_op)*(I*2.0);

  e*( (angle/2.0).cos() ) + i_sigma*((angle/2.0).sin() )


}
//------------------------------------------------------------------------------
pub fn ideal_pulse(spin_op: &SpinOp, angle: f64,
    spin_multipliciy: usize, transition: &[usize;2])
  -> CxMat
{
  let u2 = ideal_two_state_pulse(spin_op, angle);

  let mut u = CxMat::eye(spin_multipliciy);
  if spin_multipliciy <= 1{
    return u
  }

  for (ii,m) in transition.iter().enumerate(){
    let idx0 = spin_multipliciy - m - 1;
    for (jj,n) in transition.iter().enumerate(){
      let idx1 = spin_multipliciy - n - 1;

      u[[idx0,idx1]] = u2[[ii,jj]];
    }
  }

  u
}
//------------------------------------------------------------------------------

#[cfg(test)]
mod tests{
  use super::*;

  use crate::physical_constants::{PI,SQRT2_INV};

  use ndarray_linalg::Norm;
  //----------------------------------------------------------------------------
  #[test]
  fn test_ideal_pulse(){
    let ref_pi_over_2 = ndarray::array![
      [SQRT2_INV,-SQRT2_INV],
      [SQRT2_INV,SQRT2_INV],
    ];
    let pi_over_2 = ideal_pulse(&SpinOp::Sy, 0.5*PI,2,&[0,1]);
    let err = (pi_over_2 - ref_pi_over_2).norm();
    assert!(err< 1e-12);
  }
  //----------------------------------------------------------------------------
}
