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
    Err(_) => Err(CluEError::CannotOpenFile(filename.to_string())),
  }
}
//----------------------------------------------------------------------------
// TODO: Test
/// This function loads auxiliary signals from 'load_dir' into a `ClusterSet`.
pub fn load_auxiliary_signals(cluster_set: &mut ClusterSet, 
    config: &Config, load_dir: &str) 
    -> Result<(),CluEError>
{
  let clusters = &mut cluster_set.clusters;

  let Some(batch_size) = config.cluster_batch_size else{
    return Err(CluEError::NoClusterBatchSize);
  };

  let max_size = clusters.len();

  // Loop over cluster sizes.
  for cluster_size in 1..=max_size{

    // Split the cluster into batches.
    let n_clusters = clusters[cluster_size-1].len();
    let n_batches 
      = math::ceil( (n_clusters as f64)/(batch_size as f64)) as usize;

    // Loop over batches
    for ibatch in 0..n_batches{
      let idx = ibatch*batch_size;

      // Check for saved data.
      let aux_filename = format!("{}/cluster_size_{}_batch_{}.csv",
          load_dir, cluster_size,ibatch);

      if Path::new(&aux_filename).exists(){
         load_batch_signals(&mut clusters[cluster_size-1],
             idx,batch_size,&aux_filename)?;

         continue;
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
    Err(_) => Err(CluEError::CannotOpenFile(filename.to_string())),
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
      Err(err) => Err(CluEError::CannotWriteFile(
            format!("{}, {}", filename.to_string(),err.to_string() ))),
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
