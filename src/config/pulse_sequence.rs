use crate::CluEError;

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
      return Err(CluEError::TomlPulseStepSpecifier(
            "pulse_step must be of type toml::Value::Array".to_string()));
    }; 
    if pulse_step.is_empty(){
      return Err(CluEError::TomlPulseStepSpecifier(
            "pulse_step cannot be empty".to_string()));
    }
    if pulse_step[0].as_str() == Some("detect"){
      return Ok(PulseStepSpecifier::Detect);
    }
    if pulse_step.len() < 2{
      return Err(CluEError::TomlPulseStepSpecifier(
        "missing elements from pulse_step".to_string()));
    }

    let out = match &pulse_step[0]{
      toml::Value::String(ps) => {
        match ps.as_str(){
          "pulse" => {
            let toml::Value::String(pulse_name) = &pulse_step[1] else{
              return Err(CluEError::TomlPulseStepSpecifier(
                "pulse step names must be strings".to_string()));
            };
            PulseStepSpecifier::Pulse(pulse_name.to_owned())
          },
          "delay" => parse_delay(pulse_step,1)?,
          "delay2" => parse_delay(pulse_step,2)?,
          _ => return Err(CluEError::TomlPulseStepSpecifier(
                format!("cannot parse \"{}\"",ps) )),
        }
      } 
      _ => return Err(CluEError::TomlPulseStepSpecifier(
            format!("cannot parse \"{:?}\"",pulse_step[0]))),
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
      "fixed delays require the number of dt steps".to_string()));
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
      _ => return Err(CluEError::ErrorPulseSequence(
            "cannot parse pule sequence".to_string())), 
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_array(pulse_seq: &[toml::Value]) -> Result<Self,CluEError>{
    let mut pulse_sequence 
      = Vec::<PulseStepSpecifier>::with_capacity(pulse_seq.len());

    for step in pulse_seq.iter(){
      let pulse_step = PulseStepSpecifier::from_toml_value(step)?;
      pulse_sequence.push(pulse_step)
    }
  
    Ok(Self::Custom(pulse_sequence))
  }
  //----------------------------------------------------------------------------
  pub fn from_str(pulse_seq: &str) -> Result<Self,CluEError>
  {
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
}
