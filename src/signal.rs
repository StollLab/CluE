pub mod cluster_signals;
pub mod calculate_signal;

use crate::CluEError;
use crate::cluster::Cluster;
use crate::config::Config;
use crate::io::write_atomically;
use crate::structure::Structure;

use std::ops::{Add,Sub,Mul,Div};
use num_complex::Complex64;
use std::error::Error;

// The headers of the time columns of a signal file.
const TAU_HEADER: &str = "tau";
const TAU2_HEADER: &str = "tau2";

/// `TauAxes` holds the time columns of a signal file, in the output unit of
/// time (units.time): one entry per time point, with tau in the outer loop
/// and tau2 in the inner loop, so that `tau[i]` and `tau2[i]` belong to 
/// point `i` of a signal.
#[derive(PartialEq,Debug,Clone)]
pub struct TauAxes{
  pub tau: Vec::<f64>,
  pub tau2: Vec::<f64>,
}
impl TauAxes{
  //----------------------------------------------------------------------------
  /// This function builds the time columns of a run from the tau and tau2
  /// axes of `config` (`Config::set_tau_axis()` must have been called),
  /// converted from seconds to units.time.
  pub fn from_config(config: &Config) -> Result<Self,CluEError>{
    let seconds_to_unit_of_time = match config.unit_of_time_to_seconds{
      Some(unit_of_time_to_seconds) => 1.0/unit_of_time_to_seconds,
      None => return Err(CluEError::NoUnitOfTime),
    };

    let tau_axis = config.get_tau_axis_as_ref()?;
    let tau2_axis = config.get_tau2_axis_as_ref()?;

    let n = tau_axis.len()*tau2_axis.len();
    let mut tau = Vec::<f64>::with_capacity(n);
    let mut tau2 = Vec::<f64>::with_capacity(n);

    // This is the order in which the propagators fill a signal.
    for t1 in tau_axis.iter(){
      for t2 in tau2_axis.iter(){
        tau.push(*t1*seconds_to_unit_of_time);
        tau2.push(*t2*seconds_to_unit_of_time);
      }
    }

    Ok(TauAxes{tau, tau2})
  }
  //----------------------------------------------------------------------------
}

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
  /// This function reads the last signal column of a csv file; the time 
  /// columns, if the file has them, are ignored.
  /// It returns an error if the file cannot be read or has no signal column.
  pub fn read_from_csv(filename: &str) -> Result<Self,CluEError>
  {
    let (_tau_axes, mut signals) = load_csv_to_vec_signals(filename)?;

    match signals.pop(){
      Some(signal) => Ok(signal),
      None => Err(CluEError::CannotReadFile(filename.to_string(),
            "the file has no signal column".to_string())),
    }
  }
  //----------------------------------------------------------------------------
  /// This function tries to write a `Signal` to a csv file, with the columns
  /// `tau`, `tau2`, and `signal`.
  /// It will return an `Err(CluEError)` if it fails, or if `tau_axes` does
  /// not have one entry per point of the signal.
  pub fn write_to_csv(&self, tau_axes: &TauAxes, filename: &str) 
    -> Result<(),CluEError>
  {
    write_vec_signals(tau_axes, std::slice::from_ref(self), 
        vec!["signal".to_string()], filename)
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
/// The file starts with the time columns `tau` and `tau2` from `tau_axes`;
/// each `Signal` is then saved as a column with the supplied header.
/// The function will return an error if the number of `Signal`s is different
/// form the number of headers, if the `Signal`s are not all the same length,
/// if the time columns do not have one entry per point of the `Signal`s,
/// or if the csv file cannot be written.
pub fn write_vec_signals(tau_axes: &TauAxes, signals: &[Signal], 
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

   // Without signals, the file has only the time columns.
   let n_data = signals.first().map_or(tau_axes.tau.len(), |s| s.len());
   if tau_axes.tau.len() != n_data || tau_axes.tau2.len() != n_data{
     return Err(CluEError::SignalTauAxesLength(filename.to_string(),
           tau_axes.tau.len(), tau_axes.tau2.len(), n_data));
   }

    write_atomically(filename, |tmp_file|
      match write_vec_signals_to_csv(tau_axes,signals,tmp_file,headers){
        Ok(()) => Ok(()),
        Err(err) => Err(CluEError::CannotWriteCSV(filename.to_string(),
              err.to_string())),
      }
    )
}
//------------------------------------------------------------------------------
// This function is the back end to `write_vec_signals()`, and saves the
// time columns and a list of `Signal's to a csv file.  
// Each `Signal` is saved as a column with the supplied header.
// The lengths are checked by `write_vec_signals()`.
fn write_vec_signals_to_csv(tau_axes: &TauAxes, signals: &[Signal],
    filename: &str, headers: Vec::<String>) -> Result<(),Box<dyn Error>>
{
  let n_data = tau_axes.tau.len();

  let mut wtr = csv::Writer::from_path(filename)?;
  
  wtr.write_record([TAU_HEADER, TAU2_HEADER].into_iter()
      .chain(headers.iter().map(|h| h.as_str())))?;

  for ii in 0..n_data{
    
    let mut rec = Vec::<String>::with_capacity(signals.len() + 2);

    rec.push(tau_axes.tau[ii].to_string());
    rec.push(tau_axes.tau2[ii].to_string());
   
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
/// If the first two columns are the time columns `tau` and `tau2`, they are
/// returned as `Some(TauAxes)`, and the other columns are the signals;
/// otherwise, every column is a signal and the time columns are `None`.
/// This function will return an error if it cannot find/parse the file,
/// or if a column `tau` or `tau2` is not in its place.
pub fn load_csv_to_vec_signals(filename: &str)
  -> Result<(Option<TauAxes>, Vec::<Signal>), CluEError>
{
  let mut rdr = match csv::Reader::from_path(filename){
    Ok(rdr) => rdr,
    Err(err) => return Err(CluEError::CannotOpenFile(
        filename.to_string(), err.to_string())),
  };

  let headers = match rdr.headers(){
    Ok(headers) => headers.clone(),
    Err(err) => return Err(CluEError::CannotReadFile(
        filename.to_string(), err.to_string())),
  };

  let n_axes = match count_time_columns(&headers){
    Ok(n_axes) => n_axes,
    Err(reason) => return Err(CluEError::CannotReadFile(
        filename.to_string(), reason)),
  };

  let mut tau = Vec::<f64>::new();
  let mut tau2 = Vec::<f64>::new();
  let mut signals = (n_axes..headers.len()).map(|_| Signal::new())
    .collect::<Vec::<Signal>>();

  for (idata, result) in rdr.records().enumerate() {
      let record = match result{
        Ok(record) => record,
        Err(err) => return Err(CluEError::CannotReadFile(
            filename.to_string(), err.to_string())),
      };

      if record.len() != headers.len(){
        return Err(CluEError::CannotReadFile(filename.to_string(), 
              format!("row {} has {} entries, but there are {} headers",
                idata + 1, record.len(), headers.len())));
      }

      for (icol,entry) in record.iter().enumerate(){
        if icol < n_axes{
          let Ok(t) = entry.parse::<f64>() else{
            return Err(CluEError::CannotParseSignalEntry(filename.to_string(),
                  idata + 1, entry.to_string()));
          };
          if icol == 0{
            tau.push(t);
          }else{
            tau2.push(t);
          }
        }else{
          let Ok(v) = entry.parse::<Complex64>() else{
            return Err(CluEError::CannotParseSignalEntry(filename.to_string(),
                  idata + 1, entry.to_string()));
          };
          signals[icol - n_axes].data.push(v);
        }
      }
  }

  let tau_axes = if n_axes == 2{
    Some(TauAxes{tau, tau2})
  }else{
    None
  };

  Ok((tau_axes, signals))
}
//------------------------------------------------------------------------------
// This function returns the number of time columns of a signal file with
// `headers`: 2 if the file starts with `tau` and `tau2`, and 0 otherwise.
// It returns the reason as an error if `tau` or `tau2` is anywhere else.
fn count_time_columns(headers: &csv::StringRecord) -> Result<usize,String>{
  let n_axes = if headers.get(0) == Some(TAU_HEADER) 
    && headers.get(1) == Some(TAU2_HEADER)
  {
    2
  }else{
    0
  };

  for (icol, header) in headers.iter().enumerate().skip(n_axes){
    if header == TAU_HEADER || header == TAU2_HEADER{
      return Err(format!("column {} has the header \"{}\", but the time \
columns must be the first two columns, \"{}\" and \"{}\"", 
            icol + 1, header, TAU_HEADER, TAU2_HEADER));
    }
  }

  Ok(n_axes)
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

  // Batch files have no time columns.
  let (_tau_axes, signals) = load_csv_to_vec_signals(filename)?;

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
/// The file is written through a temporary file (`io::write_atomically()`),
/// so an interrupted run never leaves a partly written batch file.
pub fn write_batch_signals(clusters: &[Cluster], 
    n_data: usize, idx: usize, batch_size: usize, 
    filename: &str, structure: &Structure) -> Result<(),CluEError>
{
  write_atomically(filename, |tmp_file| 
    write_batch_signals_to_file(clusters, n_data, idx, batch_size, filename,
        tmp_file, structure)
  )
}
//------------------------------------------------------------------------------
// This function is the back end of `write_batch_signals()`: it writes the
// batch to `tmp_file`; errors name `filename`.
fn write_batch_signals_to_file(clusters: &[Cluster], 
    n_data: usize, idx: usize, batch_size: usize, 
    filename: &str, tmp_file: &str, structure: &Structure) 
  -> Result<(),CluEError>
{

  let mut headers = Vec::<String>::with_capacity(batch_size);

  for cluster in clusters.iter().skip(idx).take(batch_size){
    let header = cluster.to_header(structure)?;
    headers.push(header);
  }

  let mut wtr = match csv::Writer::from_path(tmp_file){
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

    // The file has no time columns.
    let (tau_axes, signals) = load_csv_to_vec_signals(&filename).unwrap();
    assert_eq!(tau_axes, None);
    assert_eq!(signals.len(), 3);

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
  fn test_signal_csv_round_trip_is_exact(){
    // Resuming an output loads saved signals instead of computing them, 
    // so a written and read signal must have the same bits.
    let values = [0.0, -0.0, 1.0/3.0, 0.1 + 0.2, -1e300, 5e-324, 
      -2.2250738585072014e-308, 1.6562997521730228e-32, f64::MAX];
    let mut data = Vec::<Complex64>::new();
    for &re in values.iter(){
      for &im in values.iter(){
        data.push(Complex64::new(re, im));
      }
    }
    // The time columns must also read back to the same bits.
    let tau = data.iter().map(|z| z.re).collect::<Vec::<f64>>();
    let tau2 = data.iter().map(|z| z.im).collect::<Vec::<f64>>();
    let tau_axes = TauAxes{tau, tau2};
    let signal = Signal{data};
    let filename = std::env::temp_dir().join(format!(
          "clue_oxide_test_round_trip_{}.csv", std::process::id()));
    let filename = filename.to_str().unwrap().to_string();
    write_vec_signals(&tau_axes, &[signal.clone()], 
        vec!["signal_0".to_string()], &filename).unwrap();
    let (read_axes, read) = load_csv_to_vec_signals(&filename).unwrap();
    std::fs::remove_file(&filename).unwrap();
    assert_eq!(read.len(), 1);
    for (a, b) in signal.data.iter().zip(read[0].data.iter()){
      assert_eq!((a.re.to_bits(), a.im.to_bits()), 
          (b.re.to_bits(), b.im.to_bits()), "{} read as {}", a, b);
    }
    let read_axes = read_axes.unwrap();
    for (axis, read_axis) in [(&tau_axes.tau, &read_axes.tau), 
      (&tau_axes.tau2, &read_axes.tau2)]
    {
      assert_eq!(axis.len(), read_axis.len());
      for (a, b) in axis.iter().zip(read_axis.iter()){
        assert_eq!(a.to_bits(), b.to_bits(), "{} read as {}", a, b);
      }
    }
  }
  //----------------------------------------------------------------------------
  // The layout of a signal file: tau, tau2, then the signals, one row per
  // time point; Signal::write_to_csv uses the header "signal", and 
  // Signal::read_from_csv returns the last signal column.
  #[test]
  fn test_signal_file_columns(){
    let tau_axes = TauAxes{tau: vec![0.0, 0.0, 2.5], tau2: vec![0.0, 8.0, 0.0]};
    let a = Signal{data: vec![Complex64::new(1.0,0.0), Complex64::new(0.5,0.25),
      Complex64::new(0.0,-1.0)]};
    let b = a.scale(Complex64::new(2.0,0.0));
    let filename = std::env::temp_dir().join(format!(
          "clue_oxide_test_columns_{}.csv", std::process::id()));
    let filename = filename.to_str().unwrap().to_string();

    write_vec_signals(&tau_axes, &[a.clone(), b.clone()], 
        vec!["signal_0".to_string(), "signal_1".to_string()], &filename)
      .unwrap();
    assert_eq!(std::fs::read_to_string(&filename).unwrap(),
        "tau,tau2,signal_0,signal_1\n\
0,0,1+0i,2+0i\n\
0,8,0.5+0.25i,1+0.5i\n\
2.5,0,0-1i,0-2i\n");
    assert_eq!(load_csv_to_vec_signals(&filename).unwrap(),
        (Some(tau_axes.clone()), vec![a.clone(), b.clone()]));
    assert_eq!(Signal::read_from_csv(&filename).unwrap(), b);

    a.write_to_csv(&tau_axes, &filename).unwrap();
    assert_eq!(std::fs::read_to_string(&filename).unwrap(),
        "tau,tau2,signal\n0,0,1+0i\n0,8,0.5+0.25i\n2.5,0,0-1i\n");
    assert_eq!(Signal::read_from_csv(&filename).unwrap(), a);

    // Without signals, only the time columns are written.
    write_vec_signals(&tau_axes, &[], Vec::<String>::new(), &filename)
      .unwrap();
    assert_eq!(std::fs::read_to_string(&filename).unwrap(),
        "tau,tau2\n0,0\n0,8\n2.5,0\n");
    assert_eq!(load_csv_to_vec_signals(&filename).unwrap(),
        (Some(tau_axes.clone()), Vec::<Signal>::new()));
    assert_eq!(Signal::read_from_csv(&filename), 
        Err(CluEError::CannotReadFile(filename.clone(),
            "the file has no signal column".to_string())));

    std::fs::remove_file(&filename).unwrap();
  }
  //----------------------------------------------------------------------------
  // Time columns that do not match the signals are not written, and time
  // columns out of place, or a non-numeric time, are not read.
  #[test]
  fn test_time_column_errors(){
    let filename = std::env::temp_dir().join(format!(
          "clue_oxide_test_time_columns_{}.csv", std::process::id()));
    let filename = filename.to_str().unwrap().to_string();
    let signal = Signal::ones(3);

    let short = TauAxes{tau: vec![0.0, 1.0], tau2: vec![0.0, 0.0, 0.0]};
    assert_eq!(write_vec_signals(&short, &[signal.clone()], 
          vec!["signal_0".to_string()], &filename),
        Err(CluEError::SignalTauAxesLength(filename.clone(), 2, 3, 3)));
    let short2 = TauAxes{tau: vec![0.0, 1.0, 2.0], tau2: vec![0.0]};
    assert_eq!(signal.write_to_csv(&short2, &filename),
        Err(CluEError::SignalTauAxesLength(filename.clone(), 3, 1, 3)));
    assert!(!std::path::Path::new(&filename).exists());

    for (contents, reason) in [
      ("tau,signal_0\n0,1+0i\n", "column 1 has the header \"tau\", but the \
time columns must be the first two columns, \"tau\" and \"tau2\""),
      ("tau2,tau,signal_0\n0,0,1+0i\n", "column 1 has the header \"tau2\", \
but the time columns must be the first two columns, \"tau\" and \"tau2\""),
      ("tau,tau2,signal_0,tau\n0,0,1+0i,0\n", "column 4 has the header \
\"tau\", but the time columns must be the first two columns, \"tau\" and \
\"tau2\""),
      ("signal_0,tau2\n1+0i,0\n", "column 2 has the header \"tau2\", but \
the time columns must be the first two columns, \"tau\" and \"tau2\""),
    ]{
      std::fs::write(&filename, contents).unwrap();
      assert_eq!(load_csv_to_vec_signals(&filename).err(), 
          Some(CluEError::CannotReadFile(filename.clone(), reason.to_string())),
          "{}", contents);
    }

    std::fs::write(&filename, "tau,tau2,signal_0\n0,0,1+0i\n1,x,1+0i\n")
      .unwrap();
    assert_eq!(load_csv_to_vec_signals(&filename).err(), 
        Some(CluEError::CannotParseSignalEntry(filename.clone(), 2, 
            "x".to_string())));

    std::fs::remove_file(&filename).unwrap();
  }
  //----------------------------------------------------------------------------
  // TauAxes::from_config: tau in the outer loop, tau2 in the inner loop, in
  // units.time, with the arithmetic of the former tau_axis.csv writer.
  #[test]
  fn test_tau_axes_from_config(){
    let mut config = Config::new();
    config.number_timepoints = vec![3];
    config.tau_increments = vec![2e-9];
    config.number_timepoints2 = vec![2];
    config.tau2_increments = vec![8e-9];
    config.set_tau_axis().unwrap();
    assert_eq!(TauAxes::from_config(&config), Err(CluEError::NoUnitOfTime));

    config.unit_of_time_to_seconds = Some(1e-9);
    let tau_axes = TauAxes::from_config(&config).unwrap();
    assert_eq!(tau_axes.tau.len(), config.get_total_number_timesteps());
    assert_eq!(tau_axes.tau2.len(), config.get_total_number_timesteps());

    let tau_axis = config.get_tau_axis().unwrap();
    let tau2_axis = config.get_tau2_axis().unwrap();
    let to_ns = 1.0/1e-9;
    for ii in 0..6{
      assert_eq!(tau_axes.tau[ii].to_bits(), (tau_axis[ii/2]*to_ns).to_bits());
      assert_eq!(tau_axes.tau2[ii].to_bits(), 
          (tau2_axis[ii%2]*to_ns).to_bits());
    }
    let expected_tau = [0.0, 0.0, 2.0, 2.0, 4.0, 4.0];
    let expected_tau2 = [0.0, 8.0, 0.0, 8.0, 0.0, 8.0];
    for ii in 0..6{
      assert!((tau_axes.tau[ii] - expected_tau[ii]).abs() < 1e-12);
      assert!((tau_axes.tau2[ii] - expected_tau2[ii]).abs() < 1e-12);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_write_csv_errors_give_the_reason(){
    let dir = std::env::temp_dir().join(format!(
          "clue_oxide_test_no_such_dir_{}", std::process::id()));
    assert!(!dir.exists());
    let filename = dir.join("signal.csv").to_str().unwrap().to_string();
    let signal = Signal{data: vec![Complex64::new(1.0,0.0)]};
    let tau_axes = TauAxes{tau: vec![0.0], tau2: vec![0.0]};

    for result in [signal.write_to_csv(&tau_axes, &filename),
        write_vec_signals(&tau_axes, &[signal.clone()], vec!["a".to_string()],
          &filename)]
    {
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


