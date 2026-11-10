pub mod cluster_signals;
pub mod calculate_signal;

use crate::CluEError;
use crate::cluster::Cluster;
use crate::structure::Structure;

use std::ops::{Add,Sub,Mul,Div};
use num_complex::Complex64;
use std::error::Error;

/// `Signal` is used to hold a calculated signal. 
#[derive(PartialEq,Debug,Clone,Default)]
pub struct Signal{
  pub data: Vec::<Complex64>,
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
    let mut data = Vec::<Complex64>::with_capacity(n);

    for _ii in 0..n{
      data.push(Complex64{re:1.0,im: 0.0});
    }
    Signal{data}
  }
  //----------------------------------------------------------------------------
  /// This function creates a new instance of `Signal` with `n` data points,
  /// all equal to zero.
  pub fn zeros(n: usize) -> Self {
    let mut data = Vec::<Complex64>::with_capacity(n);

    for _ii in 0..n{
      data.push(Complex64{re:0.0,im: 0.0});
    }
    Signal{data}
  }
  //----------------------------------------------------------------------------
  pub fn norm(&self) -> f64{
    self.data.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt()
  }
  //----------------------------------------------------------------------------
  /// This function scales every element by `scale_factor`.
  pub fn scale_mut(&mut self, scale_factor: Complex64){

    for z in self.data.iter_mut(){
      *z *= scale_factor;
    }

  }
  //----------------------------------------------------------------------------
  /// This function returns a scaled `Signal`.
  pub fn scale(&self, scale_factor: Complex64) -> Self{

    let mut out_signal = self.clone();
    
    out_signal.scale_mut(scale_factor);

    out_signal

  }
  //----------------------------------------------------------------------------
  pub fn read_from_csv(filename: &str) -> Result<Self,CluEError>
  {
    let signals = load_csv_to_vec_signals(filename)?;

    let n = signals.len();

    let signal = Self{data: signals[n-1].data.clone()};

    Ok(signal)
  }
  //----------------------------------------------------------------------------
  /// This function tries to write a `Signal` ro csv file.  
  /// It will return an `Err(CluEError)` if it fails. 
  pub fn write_to_csv(&self, filename: &str) -> Result<(),CluEError>{
    match self.write_single_signal_to_csv(filename){
      Ok(()) => Ok(()),
      Err(err) => Err(CluEError::CannotWriteCSV(filename.to_string(),
            err.to_string())),
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
      Err(err) => Err(CluEError::CannotWriteCSV(filename.to_string(),
            err.to_string())),
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
  let mut rdr = match csv::Reader::from_path(filename){
    Ok(rdr) => rdr,
    Err(err) => return Err(CluEError::CannotOpenFile(
        filename.to_string(), err.to_string())),
  };

  let mut num_data = 0;
  let mut num_sigs = 0;
  for result in rdr.records() {
      match result{
        Ok(str_rec) => num_sigs = str_rec.len(),
        Err(err) => return Err(CluEError::CannotReadFile(
            filename.to_string(), err.to_string())),
      };
      num_data += 1;
  }

  // Load signal.
  let mut signals = (0..num_sigs).map(|_| Signal::zeros(num_data))
    .collect::<Vec::<Signal>>();

  let mut rdr = match csv::Reader::from_path(filename){
    Ok(rdr) => rdr,
    Err(err) => return Err(CluEError::CannotOpenFile(
        filename.to_string(), err.to_string())),
  };

  for (idata, result) in rdr.records().enumerate() {
      let record = match result{
        Ok(record) => record,
        Err(err) => return Err(CluEError::CannotReadFile(
            filename.to_string(), err.to_string())),
      };

      for (isig,entry) in record.iter().enumerate(){
        let Ok(v) = entry.parse::<Complex64>() else{
          return Err(CluEError::CannotParseSignalEntry(filename.to_string(),
                idata + 1, entry.to_string()));
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
/// The batch starts at `clusters[idx]` and holds up to `batch_size` clusters.
/// This function will return an error if it cannot find/parse the file,
/// if `batch_size` is zero, if the file does not have exactly one signal
/// per cluster in the batch, or if the signals do not have `n_data` points.
pub fn load_batch_signals(clusters: &mut [Cluster],
    idx: usize, batch_size: usize, n_data: usize,
    filename: &str) -> Result<(),CluEError>
{
  if batch_size == 0{
    return Err(CluEError::ZeroClusterBatchSize);
  }

  let signals = load_csv_to_vec_signals(filename)?;

  // Number of clusters in this batch; the last batch can be shorter.
  let n_expected = std::cmp::min(batch_size,
      clusters.len().saturating_sub(idx));

  if signals.len() != n_expected{
    return Err(CluEError::WrongNumberOfBatchSignals(filename.to_string(),
          signals.len(), n_expected));
  }

  // All columns of a csv file have the same length.
  if let Some(signal) = signals.first(){
    if signal.len() != n_data{
      return Err(CluEError::WrongBatchSignalLength(filename.to_string(),
            signal.len(), n_data));
    }
  }

  for (ii, signal) in signals.into_iter().enumerate(){
    clusters[idx + ii].signal = Ok(Some(signal));
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

  let mut wtr = match csv::Writer::from_path(filename){
    Ok(wtr) => wtr,
    Err(err) => return Err(CluEError::CannotCreateFile(filename.to_string(),
          err.to_string())),
  };

  if let Err(err) = wtr.write_record(&headers){
    return Err(CluEError::CannotWriteCSV(filename.to_string(),
          err.to_string()));
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

    if let Err(err) = wtr.write_record(&rec){
      return Err(CluEError::CannotWriteCSV(filename.to_string(),
            err.to_string()));
    }
  }

  if let Err(err) = wtr.flush(){
    return Err(CluEError::CannotWriteCSV(filename.to_string(),
          err.to_string()));
  }

  Ok(())
}
//------------------------------------------------------------------------------
impl Add for &Signal{
  type Output = Signal;
 
  // The implementation of "+" is elementi-wise.
  fn add(self,rhs: &Signal) -> Signal{
    assert_eq!(self.len(), rhs.len());
    let mut data = Vec::<Complex64>::with_capacity(self.len());
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
    let mut data = Vec::<Complex64>::with_capacity(self.len());
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
    let mut data = Vec::<Complex64>::with_capacity(self.len());
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
    let mut data = Vec::<Complex64>::with_capacity(self.len());
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
    const ONE: Complex64 = Complex64{re:1.0,im:0.0};
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
        Complex64{re: 1.0, im: 0.0},  
        Complex64{re: 0.99, im: 0.01},
        Complex64{re: 0.8, im: 0.2},
        Complex64{re: 0.7, im: 0.3}
      ]},
      Signal{data: vec![
        Complex64{re: 1.0, im: 0.0},  
        Complex64{re: 0.9, im: 0.1},
        Complex64{re: 0.5, im: 0.5},
        Complex64{re: 0.6, im: 0.4}
      ]},
      Signal{data: vec![
        Complex64{re: 1.0, im: 0.0},  
        Complex64{re: 0.999, im: 0.001},
        Complex64{re: 0.98, im: 0.02},
        Complex64{re: 0.97, im: 0.03}
      ]},
    ];

    for (ii,ref_sig) in ref_signals.iter().enumerate(){
      assert_eq!(signals[ii], *ref_sig);
    }

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_load_csv_signal_entry_error(){
    let filename = std::env::temp_dir().join(format!(
          "clue_oxide_test_signal_{}.csv", std::process::id()));
    std::fs::write(&filename, "signal\n1+0i\nabc\n").unwrap();
    let filename = filename.to_str().unwrap().to_string();
    let result = load_csv_to_vec_signals(&filename);
    std::fs::remove_file(&filename).unwrap();
    assert_eq!(result.err(), Some(CluEError::CannotParseSignalEntry(
            filename.clone(), 2, "abc".to_string())));
  }
  //----------------------------------------------------------------------------
  // load_batch_signals attaches one signal per cluster of the batch, and
  // rejects files that do not match the batch or the time axis, and a
  // zero batch size.
  #[test]
  fn test_load_batch_signals_checks(){
    let filename = std::env::temp_dir().join(format!(
          "clue_oxide_test_batch_{}.csv", std::process::id()));
    let filename = filename.to_str().unwrap().to_string();
    std::fs::write(&filename, 
        "clu_1,clu_2\n1+0i,1+0i\n0.9+0i,0.8+0i\n0.7+0i,0.6+0i\n").unwrap();
    let new_clusters = || (1..=5).map(|v| Cluster::from(vec![v]))
      .collect::<Vec::<Cluster>>();

    // With batches of 2 over 5 clusters, the batch at index 3 holds
    // clusters 3 and 4.
    let mut clusters = new_clusters();
    load_batch_signals(&mut clusters, 3, 2, 3, &filename).unwrap();
    assert_eq!(clusters[3].signal, Ok(Some(Signal{data: vec![
        Complex64::new(1.0,0.0),Complex64::new(0.9,0.0),
        Complex64::new(0.7,0.0)]})));
    assert_eq!(clusters[4].signal, Ok(Some(Signal{data: vec![
        Complex64::new(1.0,0.0),Complex64::new(0.8,0.0),
        Complex64::new(0.6,0.0)]})));
    assert_eq!(clusters[2].signal, Ok(None));

    // The batch at index 4 has one cluster, but the file has two signals.
    let mut clusters = new_clusters();
    assert_eq!(load_batch_signals(&mut clusters, 4, 2, 3, &filename),
        Err(CluEError::WrongNumberOfBatchSignals(filename.clone(),2,1)));
    // A batch of three clusters.
    assert_eq!(load_batch_signals(&mut clusters, 0, 3, 3, &filename),
        Err(CluEError::WrongNumberOfBatchSignals(filename.clone(),2,3)));
    // The run has a different number of time points.
    assert_eq!(load_batch_signals(&mut clusters, 0, 2, 4, &filename),
        Err(CluEError::WrongBatchSignalLength(filename.clone(),3,4)));
    assert_eq!(load_batch_signals(&mut clusters, 0, 0, 3, &filename),
        Err(CluEError::ZeroClusterBatchSize));
    assert!(clusters.iter().all(|c| c.signal == Ok(None)));

    std::fs::remove_file(&filename).unwrap();
  }
  //----------------------------------------------------------------------------
  // Writing into a directory that does not exist gives CannotWriteCSV with
  // the file name and the reason from the operating system.
  #[test]
  fn test_write_csv_errors_give_the_reason(){
    let dir = std::env::temp_dir().join(format!(
          "clue_oxide_test_no_such_dir_{}", std::process::id()));
    assert!(!dir.exists());
    let filename = dir.join("signal.csv").to_str().unwrap().to_string();
    let signal = Signal{data: vec![Complex64::new(1.0,0.0)]};

    for result in [signal.write_to_csv(&filename),
        write_vec_signals(&[signal.clone()], vec!["a".to_string()], &filename)]{
      match result{
        Err(CluEError::CannotWriteCSV(file,reason)) => {
          assert_eq!(file, filename);
          assert!(!reason.is_empty());
        },
        other => panic!("expected CannotWriteCSV, got {:?}", other),
      }
    }
    assert!(!dir.exists());
  }
}


