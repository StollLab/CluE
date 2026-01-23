pub mod calculate_analytic_restricted_2cluster_signals;
pub mod calculate_signal;
pub mod cluster_correlation_expansion;

use crate::CluEError;
use crate::cluster::Cluster;
use crate::structure::Structure;

use std::ops::{Add,Sub,Mul,Div};
use num_complex::Complex;
use std::error::Error;

/// `Signal` is used to hold a calculated signal. 
#[derive(PartialEq,Debug,Clone,Default)]
pub struct Signal{
  pub data: Vec::<Complex<f64>>,
}
impl Signal{
  //----------------------------------------------------------------------------
  /// This function return the number of data points.
  pub fn len(&self) -> usize{
    self.data.len()
  }
  //----------------------------------------------------------------------------
  /// This function return `true` iff there are no data.
  pub fn is_empty(&self) -> bool{
    self.data.is_empty()
  }
  //----------------------------------------------------------------------------
  /// This function creates a new instance of `Signal` with no data.
  pub fn new() -> Self {Signal::default()}
  //----------------------------------------------------------------------------
  /// This function creates a new instance of `Signal` with `n` data points,
  /// all equal to one.
  pub fn ones(n: usize) -> Self {
    let mut data = Vec::<Complex<f64>>::with_capacity(n);

    for _ii in 0..n{
      data.push(Complex::<f64>{re:1.0,im: 0.0});
    }
    Signal{data}
  }
  //----------------------------------------------------------------------------
  /// This function creates a new instance of `Signal` with `n` data points,
  /// all equal to zero.
  pub fn zeros(n: usize) -> Self {
    let mut data = Vec::<Complex<f64>>::with_capacity(n);

    for _ii in 0..n{
      data.push(Complex::<f64>{re:0.0,im: 0.0});
    }
    Signal{data}
  }
  //----------------------------------------------------------------------------
  pub fn norm(&self) -> f64{
    self.data.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt()
  }
  //----------------------------------------------------------------------------
  /// This function scales every element by `scale_factor`.
  pub fn mut_scale(&mut self, scale_factor: Complex<f64>){

    for z in self.data.iter_mut(){
      *z *= scale_factor;
    }

  }
  //----------------------------------------------------------------------------
  /// This function returns a scaled `Signal`.
  pub fn scale(&self, scale_factor: Complex<f64>) -> Self{

    let mut out_signal = self.clone();
    
    out_signal.mut_scale(scale_factor);

    out_signal

  }
  //----------------------------------------------------------------------------
  /*
  // TODO: delete or fix.
  /// This function tries to read a csv file into an `Ok(Signal)`.  
  /// It will return an `Err(CluEError)` if it fails. 
  pub fn read_from_csv(filename: &str) -> Result<Self,CluEError>
  {
    match Self::read_single_signal_from_csv(filename){
      Ok(data) => Ok(Signal{data}),
      Err(_) => Err(CluEError::CannotOpenFile(filename.to_string())),
    }

  }
  //----------------------------------------------------------------------------
  // This function tries to read a csv file into an `Ok(Vec::<Complex<f64>>)`.
  // This is the back end to `read_from_csv()`. 
  fn read_single_signal_from_csv(filename: &str)
    -> Result<Vec::<Complex<f64>>,Box<dyn Error>>
  {
    // Count data points.
    let mut rdr = csv::Reader::from_path(filename)?;
    let mut num_data = 0;
    for _result in rdr.records() {
        num_data += 1;
    }

    // Load signal.
    let mut signal = Vec::<Complex<f64>>::with_capacity(num_data);
    let mut rdr = csv::Reader::from_path(filename)?;
    for result in rdr.records() {
        let record = result?;
        if let Some(v) = record.get(0){
          let v = v.parse::<Complex<f64>>()?;
          signal.push(v);
        }
    }
    Ok(signal)
  }
  */
  //----------------------------------------------------------------------------
  /// This function tries to write a `Signal` ro csv file.  
  /// It will return an `Err(CluEError)` if it fails. 
  pub fn write_to_csv(&self, filename: &str) -> Result<(),CluEError>{
    match self.write_single_signal_to_csv(filename){
      Ok(()) => Ok(()),
      Err(_) => Err(CluEError::CannotWriteFile(filename.to_string())),
    }
  }
  //----------------------------------------------------------------------------
  // This function tries to write a `Signal` ro csv file.  
  // This is the back end to `write_to_csv()`.
  fn write_single_signal_to_csv(&self, filename: &str) 
    -> Result<(),Box<dyn Error>>
  {
    let mut wtr = csv::Writer::from_path(filename)?;
    wtr.write_record(["signal"])?;
    for v in self.data.iter(){
      wtr.write_record(&[v.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
  }
  //----------------------------------------------------------------------------
}


impl std::fmt::Display for Signal {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    let mut string = format!("[{}", self.data[0]);

    for v in self.data.iter().skip(1){
      string = format!("{},{}",string,v);
    }

    write!(f,"{}]",string)
  }
}


//------------------------------------------------------------------------------
/// This function saves a `Vec::<Signal>` to a csv file.  
/// Each `Signal` is saved as a coloumn with the supplied header.
/// The function will return an error if the number of `Signal`s is different
/// form the number of headers, if the `Signal`s are not all the same length,
/// or if the csv file cannot be written.
pub fn write_vec_signals(signals: &[Signal], 
    headers: Vec::<String>, filename: &str) 
  -> Result<(),CluEError>
{
  
   if headers.len() != signals.len(){
     return Err(CluEError::MissingHeader(headers.len(),signals.len(),
           filename.to_string()));
   }

   for signal in signals.iter(){
     if signal.data.len() != signals[0].len(){
       return Err(CluEError::AllSignalsNotSameLength(filename.to_string()));
     }
   }

    match write_vec_signals_to_csv(signals,filename,headers){
      Ok(()) => Ok(()),
      Err(_) => Err(CluEError::CannotWriteFile(filename.to_string())),
    }
}
//------------------------------------------------------------------------------
// This function is the back end to `write_vec_signals()`, and saves a 
// list of `Signal's to a csv file.  
// Each `Signal` is saved as a coloumn with the supplied header.
// The function will return an error if the number of `Signal`s is different
// form the number of headers, if the `Signal`s are not all the same length,
// or if the csv file cannot be written.
fn write_vec_signals_to_csv(signals: &[Signal],filename: &str,
    headers: Vec::<String>) -> Result<(),Box<dyn Error>>
{
  let n_data = signals[0].data.len();

  let mut wtr = csv::Writer::from_path(filename)?;
  
  wtr.write_record(&headers)?;

  for ii in 0..n_data{
    
    let mut rec = Vec::<String>::with_capacity(n_data);
   
    for signal in signals.iter(){
      rec.push(signal.data[ii].to_string());
    }

    wtr.write_record(&rec)?;
  }

  wtr.flush()?;
  
  Ok(())
}
//----------------------------------------------------------------------------
/// This function attempts to load a csv file into a `Vec::<Signal>`.
/// This function will return an error if it cannot find/parse the file.
pub fn load_csv_to_vec_signals(filename: &str)
  -> Result<Vec::<Signal>, CluEError>
{
  // Count data points.
  let Ok(mut rdr) = csv::Reader::from_path(filename) else{
    return Err(CluEError::CannotOpenFile(filename.to_string()));
  };

  let mut num_data = 0;
  let mut num_sigs = 0;
  for result in rdr.records() {
      if let Ok(str_rec) = result{
        num_sigs = str_rec.len();
      }else{
        return Err(CluEError::CannotOpenFile(filename.to_string()));
      };
      num_data += 1;
  }

  // Load signal.
  let mut signals = (0..num_sigs).map(|_| Signal::zeros(num_data))
    .collect::<Vec::<Signal>>();

  let Ok(mut rdr) = csv::Reader::from_path(filename) else{
    return Err(CluEError::CannotOpenFile(filename.to_string()));
  };

  for (idata, result) in rdr.records().enumerate() {
      let Ok(record) = result else{
        return Err(CluEError::CannotOpenFile(filename.to_string()));
      };

      for (isig,entry) in record.iter().enumerate(){
        let Ok(v) = entry.parse::<Complex<f64>>() else{
          return Err(CluEError::CannotOpenFile(filename.to_string()));
        };
        signals[isig].data[idata] = v;
      }
  }
  Ok(signals)
}
//------------------------------------------------------------------------------
/// This function attempts to load a csv file, and attach each `Signal`
/// to the corresponding cluster.
/// The file is assummed to be saved from CluE Oxide as part of the
/// checkpointing process.
/// This function will return an error if it cannot find/parse the file.
pub fn load_batch_signals(clusters: &mut [Cluster], 
    idx: usize, batch_size: usize, 
    filename: &str) -> Result<(),CluEError>
{

  let signals = load_csv_to_vec_signals(filename)?;

  let mut ii = 0;
  for signal in signals{
    clusters[idx + ii].signal = Ok(Some(signal)); 
    ii += 1;
    if ii == batch_size{
      break;
    }
  } 

  Ok(())
}
//------------------------------------------------------------------------------
/// This function attempts to write a csv file containing the `Signal` data
/// from a list of `Cluster`s.
/// The file is assummed to be saved from CluE Oxide as part of the
/// checkpointing process.
/// This function will return an error if it cannot write the file.
pub fn write_batch_signals(clusters: &[Cluster], 
    n_data: usize, idx: usize, batch_size: usize, 
    filename: &str, structure: &Structure) -> Result<(),CluEError>
{

  let mut headers = Vec::<String>::with_capacity(batch_size);

  for cluster in clusters.iter().skip(idx).take(batch_size){
    let header = cluster.to_header(structure)?;
    headers.push(header);
  }

  let Ok(mut wtr) = csv::Writer::from_path(filename) else{
    return Err(CluEError::CannotWriteFile(filename.to_string()));
  };

  if wtr.write_record(&headers).is_err(){
    return Err(CluEError::CannotWriteFile(filename.to_string()));
  }

  for ii in 0..n_data{
    
    let mut rec = Vec::<String>::with_capacity(n_data);

    for cluster in clusters.iter()
      .skip(idx).take(batch_size){

        let signal  = match &cluster.signal{
          Err(err)  => return Err(err.clone()),
          Ok(None) 
            => return Err(CluEError::ClusterHasNoSignal(cluster.to_string())),
          Ok(Some(sig)) => sig,
        };

        if signal.data.len() != n_data{
          return Err(CluEError::AllSignalsNotSameLength(filename.to_string()));
        }

        rec.push(signal.data[ii].to_string());
    }

    if wtr.write_record(&rec).is_err() {
      return Err(CluEError::CannotWriteFile(filename.to_string()));
    }
  }

  if wtr.flush().is_err(){
    return Err(CluEError::CannotWriteFile(filename.to_string()));
  }

  Ok(())
}
//------------------------------------------------------------------------------
impl Add for &Signal{
  type Output = Signal;
 
  // The implementation of "+" is elementi-wise.
  fn add(self,rhs: &Signal) -> Signal{
    assert_eq!(self.len(), rhs.len());
    let mut data = Vec::<Complex<f64>>::with_capacity(self.len());
    for (ii,v) in self.data.iter().enumerate(){
      data.push(v + rhs.data[ii]);
    }
    Signal{data}
  }
}

impl Sub for &Signal{
  type Output = Signal;
 
  // The implementation of "-" is elementi-wise.
  fn sub(self,rhs: &Signal) -> Signal{
    assert_eq!(self.len(), rhs.len());
    let mut data = Vec::<Complex<f64>>::with_capacity(self.len());
    for (ii,v) in self.data.iter().enumerate(){
      data.push(v - rhs.data[ii]);
    }
    Signal{data}
  }
}


impl Mul for &Signal{
  // The implementation of "*" is elementi-wise.
  type Output = Signal;
 
  fn mul(self, rhs: &Signal) -> Signal{
    assert_eq!(self.len(), rhs.len());
    let mut data = Vec::<Complex<f64>>::with_capacity(self.len());
    for (ii,v) in self.data.iter().enumerate(){
      data.push(v*rhs.data[ii]);
    }
    Signal{data}
  }
}
impl Div for &Signal{
  type Output = Signal;
 
  // The implementation of "/" is elementi-wise.
  fn div(self,rhs: &Signal) -> Signal{
    assert_eq!(self.len(), rhs.len());
    let mut data = Vec::<Complex<f64>>::with_capacity(self.len());
    for (ii,v) in self.data.iter().enumerate(){
      data.push(v/rhs.data[ii]);
    }
    Signal{data}
  }
}

#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_signal_ops(){
    const ONE: Complex<f64> = Complex::<f64>{re:1.0,im:0.0};
    let signal0 = Signal{ data: vec![ONE,ONE],};
    let signal1 = Signal{ data: vec![ONE,0.5*ONE],};

    let signal2 = &signal0 + &signal1;
    assert_eq!(signal2.data,vec![2.0*ONE,1.5*ONE]);

    let signal3 = &signal0 - &signal1;
    assert_eq!(signal3.data,vec![0.0*ONE,0.5*ONE]);

    let signal4 = &signal0 * &signal1;
    assert_eq!(signal4.data,vec![ONE,0.5*ONE]);

    let signal5 = &signal0 / &signal1;
    assert_eq!(signal5.data,vec![ONE,2.0*ONE]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_load_csv_to_vec_signals(){
    let filename = "assets/auxiliary_signals.csv";

    let signals = load_csv_to_vec_signals(&filename).unwrap();

    let ref_signals = vec![
      Signal{data: vec![
        Complex::<f64>{re: 1.0, im: 0.0},  
        Complex::<f64>{re: 0.99, im: 0.01},
        Complex::<f64>{re: 0.8, im: 0.2},
        Complex::<f64>{re: 0.7, im: 0.3}
      ]},
      Signal{data: vec![
        Complex::<f64>{re: 1.0, im: 0.0},  
        Complex::<f64>{re: 0.9, im: 0.1},
        Complex::<f64>{re: 0.5, im: 0.5},
        Complex::<f64>{re: 0.6, im: 0.4}
      ]},
      Signal{data: vec![
        Complex::<f64>{re: 1.0, im: 0.0},  
        Complex::<f64>{re: 0.999, im: 0.001},
        Complex::<f64>{re: 0.98, im: 0.02},
        Complex::<f64>{re: 0.97, im: 0.03}
      ]},
    ];

    for (ii,ref_sig) in ref_signals.iter().enumerate(){
      assert_eq!(signals[ii], *ref_sig);
    }

  }
}


