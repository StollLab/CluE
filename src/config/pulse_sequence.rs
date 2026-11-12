use crate::CluEError;

use std::str::FromStr;
use substring::Substring;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq)]
pub enum PulseStepSpecifier{
  Pulse(String),
  FixedDelay(i32,Option<usize>),  
  FixedDelay2(i32,Option<usize>),  
  TauDelay,
  Tau2Delay,
  InvTauDelay,
  InvTau2Delay,
  Detect,
}

impl PulseStepSpecifier{
  //----------------------------------------------------------------------------
  pub fn from_toml_value(pulse_step: &toml::Value) -> Result<Self,CluEError>{
    
    let toml::Value::Array(pulse_step) = &pulse_step else{
      return Err(CluEError::TomlPulseStepSpecifier(format!(
            "each step must be an array, for example [\"pulse\", \"pi\"], \
but got {}", pulse_step)));
    }; 
    if pulse_step.is_empty(){
      return Err(CluEError::TomlPulseStepSpecifier(
            "a step cannot be empty".to_string()));
    }
    if pulse_step[0].as_str() == Some("detect"){
      return Ok(PulseStepSpecifier::Detect);
    }
    if pulse_step.len() < 2{
      return Err(CluEError::TomlPulseStepSpecifier(format!(
        "{} needs a second element", toml::Value::Array(pulse_step.clone()))));
    }

    let out = match &pulse_step[0]{
      toml::Value::String(ps) => {
        match ps.as_str(){
          "pulse" => {
            let toml::Value::String(pulse_name) = &pulse_step[1] else{
              return Err(CluEError::TomlPulseStepSpecifier(format!(
                "pulse names must be strings, but got {}", pulse_step[1])));
            };
            PulseStepSpecifier::Pulse(pulse_name.to_owned())
          },
          "delay" => parse_delay(pulse_step,1)?,
          "delay2" => parse_delay(pulse_step,2)?,
          _ => return Err(CluEError::TomlPulseStepSpecifier(format!(
                "cannot parse \"{}\"; the first element of a step must be \
\"pulse\", \"delay\", \"delay2\", or \"detect\"",ps) )),
        }
      } 
      _ => return Err(CluEError::TomlPulseStepSpecifier(format!(
            "cannot parse {}; the first element of a step must be \
\"pulse\", \"delay\", \"delay2\", or \"detect\"",pulse_step[0]))),
    };

    Ok(out)  
  }
  //----------------------------------------------------------------------------
}

fn parse_delay(pulse_step: &[toml::Value],delay_dim: usize) 
    -> Result<PulseStepSpecifier,CluEError>
{
  if pulse_step.len() < 2{
    return Err(CluEError::TomlPulseStepSpecifier(
          "delay step arrays must have at least 2 elements".to_string()));
  }
  if pulse_step[0] != toml::Value::String("delay".to_string())
      && pulse_step[0] != toml::Value::String("delay2".to_string()){
    return Err(CluEError::TomlPulseStepSpecifier(
      "delay steps must be specified with \"delay\" or \"delay2\" as the first element".to_string()));
  }

  if pulse_step[1] == toml::Value::String("tau".to_string()){
    if delay_dim != 1{ 
      return Err(CluEError::WrongTauDimension(delay_dim,1));
    }
    return Ok(PulseStepSpecifier::TauDelay);
  }
  if pulse_step[1] == toml::Value::String("-tau".to_string()){
    if delay_dim != 1{ 
      return Err(CluEError::WrongTauDimension(delay_dim,1));
    }
    return Ok(PulseStepSpecifier::InvTauDelay);
  }

  if pulse_step[1] == toml::Value::String("tau2".to_string()){
    if delay_dim != 2{ 
      return Err(CluEError::WrongTauDimension(delay_dim,2));
    }
    return Ok(PulseStepSpecifier::Tau2Delay);
  }

  if pulse_step[1] == toml::Value::String("-tau2".to_string()){
    if delay_dim != 2{ 
      return Err(CluEError::WrongTauDimension(delay_dim,2));
    }
    return Ok(PulseStepSpecifier::InvTau2Delay);
  }

  let toml::Value::Integer(number) = pulse_step[1] else{
    return Err(CluEError::TomlPulseStepSpecifier(
      format!("a delay needs \"tau\", \"-tau\", \"tau2\", \"-tau2\", \
or an integer number of time steps, but got {}", pulse_step[1])));
  };
  

  let delay_index: Option::<usize> = if pulse_step.len() >=3 {
    match pulse_step[2]{
    toml::Value::Integer(idx) => {
      if idx < 0{ return Err(CluEError::TomlPulseStepSpecifier(
          "dt indices must be integers >= 0".to_string()));
      }
      Some(idx as usize)
    },
    _ => return Err(CluEError::TomlPulseStepSpecifier(
          "dt indices must be integers >= 0".to_string())),
    }
  }else{ None };

  match delay_dim{
    1 => Ok(PulseStepSpecifier::FixedDelay(number as i32,delay_index)),
    2 => Ok(PulseStepSpecifier::FixedDelay2(number as i32,delay_index)),
    _ => Err(CluEError::WrongDelayDimension(delay_dim)),
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `PulseSequence` lists the options for pulse sequences to simulate.
#[derive(Debug,Clone,PartialEq)]
pub enum PulseSequence{
  FreeEvolution,
  CarrPurcell(usize),
  RefocusedHahnEcho,
  Custom(Vec::<PulseStepSpecifier>),
}

impl PulseSequence{
  pub fn from_toml_value(pulse_seq: &toml::Value) -> Result<Self,CluEError>
  {
    match pulse_seq{
      toml::Value::String(ps) => Self::from_str(ps),
      toml::Value::Array(ps) => Self::from_toml_array(ps),  
      other => Err(CluEError::InvalidTOMLValue(
            "pulse_sequence.sequence".to_string(),
            "a string or an array of steps".to_string(),
            other.type_str().to_string())),
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_array(pulse_seq: &[toml::Value]) -> Result<Self,CluEError>{
    let mut pulse_sequence 
      = Vec::<PulseStepSpecifier>::with_capacity(pulse_seq.len());

    for (ii,step) in pulse_seq.iter().enumerate(){
      let pulse_step = PulseStepSpecifier::from_toml_value(step)
        .map_err(|err| CluEError::InvalidPulseStep(ii, Box::new(err)))?;
      pulse_sequence.push(pulse_step)
    }
  
    Ok(Self::Custom(pulse_sequence))
  }
  //----------------------------------------------------------------------------
}

impl FromStr for PulseSequence {
  type Err = CluEError;
  fn from_str(pulse_seq: &str) -> Result<Self, Self::Err> {
    if pulse_seq.substring(0,3) == "cp-"{
      let Ok(n_pi) = pulse_seq.substring(3,pulse_seq.len()).parse::<usize>()else{
        return Err(CluEError::CannotParsePulseSequence(pulse_seq.to_string()));
      };
      return Ok(Self::CarrPurcell(n_pi)); 
    }
    match pulse_seq{
      "free_evolution" => Ok(PulseSequence::FreeEvolution),
      "fid" => Ok(PulseSequence::CarrPurcell(0)),
      "hahn" => Ok(PulseSequence::CarrPurcell(1)),
      "refocused_hahn_echo" => Ok(PulseSequence::RefocusedHahnEcho),
      _ => Err(CluEError::CannotParsePulseSequence(pulse_seq.to_string())),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



#[cfg(test)]
mod tests{
  use super::*;
  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_PulseStepSpecifier_from_toml_value(){
    let expected = vec![
      PulseStepSpecifier::Pulse("pi/2".to_string()),
      PulseStepSpecifier::FixedDelay(4,None),  
      PulseStepSpecifier::FixedDelay(6,Some(0) ),  
      PulseStepSpecifier::Pulse("pi".to_string()),
      PulseStepSpecifier::TauDelay,
      PulseStepSpecifier::Detect,
    ];
    let table = r##"
      pulse_sequence = [
        ["pulse","pi/2"],
        ["delay", 4 ],
        ["delay", 6, 0],
        ["pulse","pi"],
        ["delay","tau"],
        ["detect"]
      ]
    "##.parse::<toml::Table>().unwrap();

    let toml::Value::Array(ps) = &table["pulse_sequence"] else{
      todo!();
    };
    for (ii,p) in ps.iter().enumerate(){
      let step = PulseStepSpecifier::from_toml_value(p).unwrap();
      assert_eq!(step,expected[ii]); 
    }
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  #[test]
  fn test_pulse_step_errors(){
    let table = r##"
      steps = ["pi", [], ["wait", 1], [1, 2], ["delay", "t"]]
    "##.parse::<toml::Table>().unwrap();
    let steps = table["steps"].as_array().unwrap();
    let messages = steps.iter()
      .map(|step| PulseStepSpecifier::from_toml_value(step).unwrap_err()
          .to_string())
      .collect::<Vec::<String>>();
    assert_eq!(messages, vec![
        "each step must be an array, for example [\"pulse\", \"pi\"], \
but got \"pi\"".to_string(),
        "a step cannot be empty".to_string(),
        "cannot parse \"wait\"; the first element of a step must be \"pulse\", \
\"delay\", \"delay2\", or \"detect\"".to_string(),
        "cannot parse 1; the first element of a step must be \"pulse\", \
\"delay\", \"delay2\", or \"detect\"".to_string(),
        "a delay needs \"tau\", \"-tau\", \"tau2\", \"-tau2\", or an integer \
number of time steps, but got \"t\"".to_string(),
    ]);

    // The step index is added by `PulseSequence::from_toml_array`.
    let err = PulseSequence::from_toml_array(steps).unwrap_err();
    assert_eq!(err, CluEError::InvalidPulseStep(0, Box::new(
            PulseStepSpecifier::from_toml_value(&steps[0]).unwrap_err())));

    let err = PulseSequence::from_toml_value(&toml::Value::Integer(3))
      .unwrap_err();
    assert_eq!(err.to_string(),
        "pulse_sequence.sequence must be a string or an array of steps, but \
got integer");
  }
}
