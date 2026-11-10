use crate::CluEError;
use crate::cluster::cluster_set::ClusterSet;
use crate::config::Config;
use crate::math;
use crate::signal::load_batch_signals;

use std::{
  collections::HashMap,
  error::Error,
  path::Path,
  str::FromStr,
};
//------------------------------------------------------------------------------
pub trait FromTOMLString{
  fn from_toml_string(s: &str) -> Result<Self,CluEError> where Self: Sized;
} 
//------------------------------------------------------------------------------
/// This function reads in the time_axis.csv file.
pub fn read_time_axis(filename: &str) -> Result<Vec::<f64>, CluEError >
{
  match read_single_vec_of_floats_from_csv(filename){
    Ok(data) => Ok(data),
    Err(err) => Err(CluEError::CannotReadFile(
        filename.to_string(), err.to_string())),
  }
}
//----------------------------------------------------------------------------
/// This function loads auxiliary signals from 'load_dir' into a `ClusterSet`.
/// The files are those written by `calculate_auxiliary_signals`:
/// `cluster_size_{n}_batch_{b}.csv` holds batch `b` of
/// `cluster_set.clusters[n]`, where `n` is the index into
/// `cluster_set.clusters` (index 0 is the empty cluster).
/// `config` must have the `cluster_batch_size` and time axis of the run
/// that wrote the files.
/// Batches without a file are left unchanged.
pub fn load_auxiliary_signals(cluster_set: &mut ClusterSet, 
    config: &Config, load_dir: &str) 
    -> Result<(),CluEError>
{
  let clusters = &mut cluster_set.clusters;

  let Some(batch_size) = config.cluster_batch_size else{
    return Err(CluEError::NoClusterBatchSize);
  };
  if batch_size == 0{
    return Err(CluEError::ZeroClusterBatchSize);
  }

  // The number of time points, as in Config::set_tau_axis();
  // `config` may not have had set_tau_axis() called.
  let n_data = {
    let (result, n1) = Config::construct_tau_axis(
        &config.number_timepoints, &config.tau_increments, "");
    result?;
    let (result, n2) = Config::construct_tau_axis(
        &config.number_timepoints2, &config.tau2_increments, "2");
    result?;
    n1*n2
  };

  // Loop over cluster sizes.
  for (cluster_size, clusters_of_size) in clusters.iter_mut().enumerate(){

    // Split the cluster into batches.
    let n_clusters = clusters_of_size.len();
    let n_batches 
      = math::ceil( (n_clusters as f64)/(batch_size as f64)) as usize;

    // Loop over batches
    for ibatch in 0..n_batches{
      let idx = ibatch*batch_size;

      // Check for saved data.
      let aux_filename = format!("{}/cluster_size_{}_batch_{}.csv",
          load_dir, cluster_size,ibatch);

      if Path::new(&aux_filename).exists(){
         load_batch_signals(clusters_of_size,
             idx,batch_size,n_data,&aux_filename)?;
      }
    }
  }
  Ok(())
}
//----------------------------------------------------------------------------
/// This function tries to read a csv file into an `Ok(Signal)`.
/// It will return an `Err(CluEError)` if it fails.
pub fn read_csv<T>(filename: &str)
  -> Result<HashMap::<String,Vec::<T>>,CluEError>
  where T: std::str::FromStr, <T as FromStr>::Err: std::error::Error, 
  <T as FromStr>::Err: 'static
{
  match read_csv_to_hashmap(filename){
    Ok(data) => Ok(data),
    Err(err) => Err(CluEError::CannotReadFile(
        filename.to_string(), err.to_string())),
  }

}
//----------------------------------------------------------------------------
// This function tries to read a csv file into an `Ok(Vec::<Complex<f64>>)`.
// This is the back end to `read_from_csv()`.
fn read_csv_to_hashmap<T>(filename: &str)
  -> Result<HashMap::<String,Vec::<T>>,Box<dyn Error>>
  where T: std::str::FromStr, <T as FromStr>::Err: std::error::Error, 
  <T as FromStr>::Err: 'static
{
  // Count data points.
  let mut rdr = csv::Reader::from_path(filename)?;
  let mut num_data = 0;
  for _result in rdr.records() {
      num_data += 1;
  }

  // Load signal.
  let headers = rdr.headers()?;
  let mut data = HashMap::<String, Vec::<T>>::with_capacity(headers.len());
  let mut rdr = csv::Reader::from_path(filename)?;
  for key in headers.iter(){
    data.insert(key.to_string(),Vec::<T>::with_capacity(num_data));
  }

  for result in rdr.records() {
    let record = result?;
    for (ii,s) in record.iter().enumerate(){
      let v = s.parse::<T>()?;
      let col = data.get_mut(&headers[ii]).ok_or("missing header")?;
      col.push(v);
    }
  }
  Ok(data)
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// This function tries to read a csv file into an `Ok(Vec::<Complex<f64>>)`.
// This is the back end to `read_from_csv()`. 
fn read_single_vec_of_floats_from_csv(filename: &str)
  -> Result<Vec::<f64>,Box<dyn Error>>
{
  // Count data points.
  let mut rdr = csv::Reader::from_path(filename)?;
  let mut num_data = 0;
  for _result in rdr.records() {
      num_data += 1;
  }

  // Load signal.
  let mut time_axis = Vec::<f64>::with_capacity(num_data);
  let mut rdr = csv::Reader::from_path(filename)?;
  for result in rdr.records() {
      let record = result?;
      if let Some(t) = record.get(0){
        let t = t.parse::<f64>()?;
        time_axis.push(t);
      }
  }
  Ok(time_axis)
}
//------------------------------------------------------------------------------
/// This function writes signal data to a csv.
pub fn write_data<T>(data: &[Vec::<T>], filename: &str,
    headers: Vec::<String>) -> Result<(),CluEError> where T: std::fmt::Display
{
   if headers.len() != data.len(){
     return Err(CluEError::MissingHeader(headers.len(), data.len(),
           filename.to_string()));
   }

   for datum in data.iter(){
     if datum.len() != data[0].len(){
       return Err(CluEError::AllVectorsNotSameLength(filename.to_string()));
     }
   }

    match write_vecs_to_csv(data,filename,headers){
      Ok(()) => Ok(()),
      Err(err) => Err(CluEError::CannotWriteCSV(filename.to_string(),
            err.to_string())),
    }
}
//------------------------------------------------------------------------------
// This function writes signal data to a csv.
// The function write_data calls this function, but has an extra error handling
// layer.
fn write_vecs_to_csv<T>(data: &[Vec::<T>],filename: &str,
    headers: Vec::<String>) 
  -> Result<(),Box<dyn Error>> where T: std::fmt::Display
{
  let n_data = data[0].len();

  let mut wtr = csv::Writer::from_path(filename)?;
  
  wtr.write_record(&headers)?;

  for ii in 0..n_data{
    
    let mut rec = Vec::<String>::with_capacity(n_data);
   
    for datum in data.iter(){
      rec.push(datum[ii].to_string());
    }

    wtr.write_record(&rec)?;
  }

  wtr.flush()?;
  
  Ok(())
}
//------------------------------------------------------------------------------


#[cfg(test)]
mod tests{
  use super::*;
  //----------------------------------------------------------------------------
  // Writing into a directory that does not exist gives CannotWriteCSV with
  // the file name and the reason from the operating system.
  #[test]
  fn test_write_data_error_gives_the_reason(){
    let dir = std::env::temp_dir().join(format!(
          "clue_oxide_test_io_no_such_dir_{}", std::process::id()));
    assert!(!dir.exists());
    let filename = dir.join("data.csv").to_str().unwrap().to_string();
    let result = write_data(&[vec![1.0, 2.0]], &filename,
        vec!["a".to_string()]);
    match result{
      Err(CluEError::CannotWriteCSV(file,reason)) => {
        assert_eq!(file, filename);
        assert!(!reason.is_empty());
      },
      other => panic!("expected CannotWriteCSV, got {:?}", other),
    }
    assert!(!dir.exists());
  }
  //----------------------------------------------------------------------------
  // A gCCE run that saves its auxiliary signals (in batches of 2) is
  // reloaded with load_auxiliary_signals. Every cluster, including the empty
  // cluster, must get a signal, and their product must be the signal of the
  // run. Before the fix, files were loaded one cluster size too low, which
  // panicked.
  #[test]
  fn test_load_auxiliary_signals_round_trip(){
    use crate::cluster::read_clusters::read_cluster_file;
    use crate::config::{SAVE_DIR_AUXILIARY_SIGNALS,SAVE_DIR_ORIENTATION_SIGNALS,
      SAVE_FILE_CLUSTERS};
    use crate::signal::{Signal,load_csv_to_vec_signals};
    use crate::structure::Structure;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    let out_dir = std::env::temp_dir().join(format!(
          "clue_oxide_test_load_aux_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out_dir);
    let out = out_dir.to_str().unwrap().to_string();

    let toml = std::fs::read_to_string(
        "assets/1omp_K26R1_0.003228966616048703_gcce.toml").unwrap()
      .replace("radius = 300","radius = 8\ncluster_batch_size = 2")
      .replace("number_timepoints = [301]","number_timepoints = [11]")
      .replace("auxiliary_signals = false","auxiliary_signals = true");
    let toml = format!("output_directory = \"{}\"\n{}",out,toml);
    let config = Config::from_toml_string(&toml).unwrap();
    assert_eq!(config.output_directory, Some(out.clone()));
    let (_, signal_first) = crate::run(config.clone()).unwrap();
    // A second run into the same directory reloads every batch, including
    // the empty cluster's, through the checkpoint path.
    let (_, signal_second) = crate::run(config.clone()).unwrap();
    assert_eq!(signal_first, signal_second);

    // The run_name of the asset is "system"; there is one orientation.
    let ori_root = format!("{}/system/{}",out,SAVE_DIR_ORIENTATION_SIGNALS);
    let entries = std::fs::read_dir(&ori_root).unwrap()
      .map(|e| e.unwrap().path()).collect::<Vec::<_>>();
    assert_eq!(entries.len(),1);
    let ori_dir = entries[0].to_str().unwrap().to_string();

    let mut rng = ChaCha20Rng::seed_from_u64(0);
    let structure = Structure::build_structure(&mut rng,&config).unwrap();
    let mut cluster_set = read_cluster_file(
        &format!("{}/{}.toml",ori_dir,SAVE_FILE_CLUSTERS), &structure).unwrap();
    assert!(cluster_set.clusters[1].len() > 2, "need several batches");

    load_auxiliary_signals(&mut cluster_set, &config, 
        &format!("{}/{}",ori_dir,SAVE_DIR_AUXILIARY_SIGNALS)).unwrap();

    let n_tot = 11;
    let mut product = Signal::ones(n_tot);
    for clusters_of_size in cluster_set.clusters.iter(){
      for cluster in clusters_of_size.iter(){
        let Ok(Some(signal)) = &cluster.signal else{
          panic!("cluster {} has no signal", cluster);
        };
        product = &product * signal;
      }
    }
    // The empty cluster's auxiliary signal (the bare detected-spin signal
    // in gCCE) is now saved and loaded.
    assert!(std::path::Path::new(&format!("{}/{}/cluster_size_0_batch_0.csv",
        ori_dir,SAVE_DIR_AUXILIARY_SIGNALS)).exists());
    assert_eq!(cluster_set.clusters[0][0].signal.as_ref().unwrap()
        .as_ref().unwrap().len(), n_tot);

    let run_signals = load_csv_to_vec_signals(
        &format!("{}/signal.csv",ori_dir)).unwrap();
    let reference = &run_signals[run_signals.len()-1];
    for (v,v0) in product.data.iter().zip(reference.data.iter()){
      assert!((v-v0).norm() < 1e-12, "{} vs {}", v, v0);
    }

    std::fs::remove_dir_all(&out_dir).unwrap();
  }
}
