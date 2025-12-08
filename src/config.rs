use crate::clue_errors::*;

use crate::cluster::{
  partition::PartitioningMethod,
  unit_of_clustering::UnitOfClustering,
};
use crate::config::config_toml::*;
use crate::config::particle_config::set_particle_configs_from_toml_table;
use crate::config::command_line_input::CommandLineInput;
use crate::config::particle_config::{
  CellType, ParticleConfig,EigSpecifier,TensorSpecifier
};
use crate::isotopes::Isotope;
use crate::physical_constants::*;
use crate::space_3d::Vector3D;
use crate::structure::particle_filter::VectorSpecifier;
use crate::integration_grid::IntegrationGrid;

use crate::misc::are_all_same_type;


use crate::io::{FromTOMLString,read_csv};


use substring::Substring;

use toml::Value;
use serde::{Serialize,Deserialize};

use std::collections::HashMap;
use std::fmt;
use std::fs;
use std::mem;

pub mod config_toml;
pub mod particle_config;
pub mod command_line_input;


// Define the directory and file names used to save outputs.
pub const SAVE_DIR_AUXILIARY_SIGNALS: &str = "auxiliary_signals";
pub const SAVE_FILE_BATH: &str = "bath";
pub const SAVE_FILE_CLUSTERS: &str = "clusters";
pub const SAVE_DIR_INFO: &str = "info";
pub const SAVE_FILE_DET_SPIN_COOR: &str = "detected_spin";
pub const SAVE_FILE_EXCHANGE_GROUPS: &str = "exchange_groups";
pub const SAVE_FILE_METHYL_PARTITIONS: &str = "methyl_partitions";
pub const SAVE_DIR_ORIENTATION_SIGNALS: &str = "orientations";
pub const SAVE_FILE_ORIENTATIONS: &str = "orientation_grid";
pub const SAVE_FILE_SANS_SPIN_SIGNALS: &str = "sans_spin_signals";
pub const SAVE_FILE_STRUCTURE_PDB: &str = "structure_pdb";
pub const SAVE_FILE_TENSORS: &str = "tensors";

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// Config contains all the setting for CluE.
#[derive(Debug,Clone,Default)]
pub struct Config{

  // Detected Spin
  pub density_matrix: Option<DensityMatrixMethod>,
  pub detected_spin_g_matrix: Option<TensorSpecifier>,
  pub detected_spin_identity: Option<Isotope>,
  pub detected_spin_multiplicity: Option<usize>,
  pub detected_spin_position: Option<DetectedSpinCoordinates>,
  pub detected_spin_transition: Option<[usize;2]>,

  // Structure
  pub replicate_unit_cell: Option<ReplicateUnitCell>,
  pub clash_distance: Option<f64>, 
  pub clash_distance_pbc: Option<f64>,
  pub extracell_particles: Vec::<ParticleConfig>,
  pub input_structure_file: Option<String>,
  pub particles: Vec::<ParticleConfig>,
  pub pdb_model_index: Option<usize>,
  pub radius: Option<f64>,
  
  // Pair Cutoffs
  pub neighbor_cutoff_coupling_xx_yy: Option<f64>,
  pub neighbor_cutoff_delta_hyperfine_zz: Option<f64>,
  pub neighbor_cutoff_point_dipole_perpendicular: Option<f64>,
  pub neighbor_cutoff_distance: Option<f64>,
  pub neighbor_cutoff_3_spin_hahn_mod_depth: Option<f64>,
  pub neighbor_cutoff_3_spin_hahn_taylor_4: Option<f64>,

  // CCE
  pub cluster_batch_size: Option<usize>, 
  pub cluster_source: Option<ClusterSource>,
  pub cluster_method: Option<ClusterMethod>,
  pub min_cell_size: Option<usize>,
  pub max_cell_size: Option<usize>,
  pub max_cluster_size: Option<usize>,
  pub max_spins: Option<usize>,
  pub partitioning: Option<PartitioningMethod>,
  pub partition_table: Option<PartitionTableConfig>,
  pub rng_seed: Option<u64>,
  pub unit_of_clustering: Option<UnitOfClustering>,

  // Experimental Details
  pub magnetic_field: Option<Vector3D>,
  pub orientation_grid: Option<OrientationAveraging>,
  pub pulse_sequence: Option<PulseSequence>,
  time_axis: Vec::<f64>,
  pub tau_increments: Vec::<f64>,
  pub number_runs: Option<usize>, 
  pub number_timepoints: Vec::<usize>,

  // Output
  pub output_directory: Option<String>,
  pub run_name: Option<String>,
  pub write_auxiliary_signals: Option<bool>, 
  pub write_bath: Option<bool>,
  pub write_clusters: Option<bool>,
  pub write_config: Option<bool>,
  pub write_detected_spin: Option<bool>,
  pub write_info: Option<bool>,
  pub write_exchange_groups: Option<bool>,
  pub write_methyl_partitions: Option<bool>,
  pub write_orientation_grid: Option<bool>,
  pub write_orientation_signals: Option<bool>,
  pub write_sans_spin_signals: Option<bool>,
  pub write_structure_pdb: Option<bool>,
  pub write_tensors: Option<bool>,

  // Units
  pub unit_of_time_to_seconds: Option<f64>,
}


impl FromTOMLString for Config{
  fn from_toml_string(toml_str: &str) -> Result<Self,CluEError>{

    let mut config = Config::new();

    let config_toml = ConfigTOML::from_toml_string(toml_str)?;
    config.set_from_config_toml(config_toml)?;

    config.set_extracell_particles();

    config.set_defaults()?;


    Ok(config)
  }
}

impl Config{
  //----------------------------------------------------------------------------
  /// This function generates a empty `Config`.
  pub fn new() -> Self{
    Default::default()
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_file(filename: &str) -> Result<Self,CluEError>
  {
  
    let toml_str: String = match fs::read_to_string(filename){
      Ok(lines) => lines,
      Err(_) => return Err(CluEError::CannotOpenFile(filename.to_string())),
    };

    let mut config = Self::from_toml_string(&toml_str)?;

    match config.set_time_axis(){
      Ok(_) => (),
      Err(err) => return Err(err),
    }

    Ok(config)
  }
  //----------------------------------------------------------------------------
  /// This function sets config setting that are necessary to be `Some`,
  /// but will likely be the same for most simulations.  
  /// Any field that is already `Some` will be left alone.
  pub fn set_defaults(&mut self) -> Result<(),CluEError> {

    if self.replicate_unit_cell.is_none(){
      self.replicate_unit_cell = Some(ReplicateUnitCell::Auto);
    }
    if self.clash_distance.is_none(){
      self.clash_distance = Some(1e-12);
    }


    if self.cluster_batch_size.is_none(){
      self.cluster_batch_size = Some(1000);
    }

    if self.cluster_source.is_none(){
      self.cluster_source = Some(ClusterSource::Structure);
    }

    if self.partitioning.is_none(){
      self.partitioning = Some(PartitioningMethod::Particles);
    }

    if self.pdb_model_index.is_none(){
      self.pdb_model_index = Some(0);
    }
    if self.density_matrix.is_none(){
      self.density_matrix = Some(DensityMatrixMethod::Identity);
    }

    if self.number_runs.is_none(){
      self.number_runs = Some(1);
    }

    if self.unit_of_clustering.is_none(){
      self.unit_of_clustering = Some(UnitOfClustering::Set);
    }
    if self.write_bath.is_none(){
      self.write_bath = Some(true);
    }

    if self.write_info.is_none(){
      self.write_info = Some(true);
    }
    if self.write_detected_spin.is_none(){
      self.write_detected_spin = Some(true);
    }

    if self.write_exchange_groups.is_none(){
      self.write_exchange_groups = Some(true);
    }

    if self.write_orientation_grid.is_none(){
      self.write_orientation_grid = Some(true);
    }

    if self.write_orientation_signals.is_none(){
      self.write_orientation_signals = Some(true);
    }


    if self.write_structure_pdb.is_none(){
      self.write_structure_pdb = Some(true);
    }
    

    self.set_detected_spin()
  }
  //----------------------------------------------------------------------------
  // This function sets the detected spin.
  fn set_detected_spin(&mut self) -> Result<(),CluEError>{

    if self.detected_spin_identity.is_none(){
      self.detected_spin_identity = Some(Isotope::Electron);
    }

    // Set g-matrix
    if self.detected_spin_g_matrix.is_some(){
     
    }else{
      self.detected_spin_g_matrix = Some(TensorSpecifier::Eig(EigSpecifier{
        values: Some([ELECTRON_G,ELECTRON_G,ELECTRON_G]),
        x_axis: Some(VectorSpecifier::Vector(Vector3D::from([1.0, 0.0, 0.0]))),
        y_axis: Some(VectorSpecifier::Vector(Vector3D::from([0.0, 1.0, 0.0]))),
        z_axis: None,
          }));
    }
    if self.detected_spin_multiplicity.is_none(){
      self.detected_spin_multiplicity = match self.detected_spin_identity{
        Some(isotope) => Some(isotope.spin_multiplicity()),
        None => return Err(CluEError::NoDetectedSpinIdentity),
      };
    }
    if self.detected_spin_transition.is_none(){
      self.detected_spin_transition = Some([0,1]);
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function the highest spin multiplicity for a "#[group(...)]". 
  pub fn max_spin_multiplicity_for_particle_config(&self, id: usize) 
    -> Option<usize>
  {


   let mult0 = if id < self.particles.len() {
     self.particles[id].max_possible_spin_multiplicity()
   } else{
     None
   };

   let mult1 = if id < self.extracell_particles.len() {
     self.extracell_particles[id].max_possible_spin_multiplicity()
   } else{
     None
   };

   match (mult0,mult1){
     (Some(m0),Some(m1)) => { if m0 >= m1 {Some(m0)}else{Some(m1)}  },
     (Some(m0),None) => Some(m0),
     (None,Some(m1)) => Some(m1),
     (None,None) => None,
   }

  }
  //----------------------------------------------------------------------------
  /// This function gets the simulated experiment's time axis.
  pub fn get_time_axis_as_ref(&self) -> Result<&Vec::<f64>,CluEError>
  {
    if self.time_axis.is_empty() {
      return Err(CluEError::NoTimeAxis);
    };
    Ok(&self.time_axis)
  }
  //----------------------------------------------------------------------------
  /// This function gets the simulated experiment's time axis.
  pub fn get_time_axis(&self) -> Result<Vec::<f64>,CluEError>
  {
    if self.time_axis.is_empty() {
      return self.construct_time_axis();
    };
    Ok(self.time_axis.clone())
  }
  //----------------------------------------------------------------------------
  pub fn set_time_axis(&mut self) -> Result<(),CluEError>
  {
    self.time_axis  = self.construct_time_axis()?;

    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function constructs the simulated experiment's time axis.
  pub fn construct_time_axis(&self) -> Result<Vec::<f64>,CluEError>
  {
    let dts = &self.tau_increments;
    if dts.is_empty(){
      return Err(CluEError::NoTimeIncrements);
    }

    let n_dts = &self.number_timepoints;
    if n_dts.is_empty(){
      return Err(CluEError::NoTimepoints);
    }

    if dts.len() != n_dts.len(){
      return Err(CluEError::LenghMismatchTimepointsIncrements(
            n_dts.len(),dts.len()));
    }

    let Some(pulse_sequence) = &self.pulse_sequence else{
      return Err(CluEError::NoPulseSequence);
    };

    let n_tot = n_dts.iter().sum::<usize>();
    let mut time_axis = Vec::<f64>::with_capacity(n_tot);
    match pulse_sequence{
      PulseSequence::CarrPurcell(n_pi_pulses) => {
        let mut t = 0.0;
        for (idx, &n_dt) in n_dts.iter().enumerate(){
          let dt = (*n_pi_pulses as f64 + 1.0)*dts[idx];
          for _ii in 0..n_dt{
            time_axis.push(t);
            t += dt;
          }
        }
      }
    }

    Ok(time_axis)
  }
  //----------------------------------------------------------------------------
  pub fn get_partition_table(&self) 
      -> Result<Vec::<(usize,usize)>,CluEError>
  {
    Ok(match &self.partition_table{
      Some(PartitionTableConfig::CSV(file_name)) => {
        let data: HashMap::<String, Vec::<usize>> = read_csv(file_name)?;
        
        let Some(indices) = data.get("index") else{
          return Err(CluEError::MissingFieldInCSVFile(
                "index".to_string(),file_name.to_string()));
        };

        let Some(blocks) = data.get("block") else{
          return Err(CluEError::MissingFieldInCSVFile(
                "block".to_string(),file_name.to_string()));
        };

        let mut p_table = Vec::<(usize,usize)>::with_capacity(indices.len());

        for (ii,idx) in indices.iter().enumerate(){
          p_table.push((*idx,blocks[ii]));
        }

        p_table
      },
      Some(PartitionTableConfig::IndexBlockList(p_table)) => p_table.clone(),
      None => Vec::<(usize,usize)>::new(),
    })

  }
  //----------------------------------------------------------------------------
  /// This function finds the `ParticleConfig` that has `label` if it exists.
  pub fn find_particle_config(&self, label: &str) 
    -> Option<(usize,&ParticleConfig)>
  {
    for (idx, p_cfg) in self.particles.iter().enumerate(){
      if *label == p_cfg.label{
        return Some((idx,p_cfg));
      }
    }
    None
  }
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `DensityMatrixMethod` specifies different methods for determining the
#[derive(Debug,Clone,PartialEq)]
pub enum DensityMatrixMethod{
  ApproxThermal(f64),
  Identity,
  Thermal(f64),
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ClusterSource` specifies where to search for clusters
#[derive(Debug,Clone,PartialEq)]
pub enum ClusterSource{
  Structure,
  File(String),
  OrientationsDir(String)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `DetectedSpinCoordinates` lists options for indicating the coordinates
/// of the detected spin.
#[derive(Debug,Clone,PartialEq)]
pub enum DetectedSpinCoordinates{
  CentroidOverGroup(Vec::<String>),
  CentroidOverSerials(Vec::<u32>),
  XYZ(Vector3D),
  ProbabilityDistribution(IntegrationGrid),
}
impl DetectedSpinCoordinates{
  pub fn from_toml_value(value: toml::Value) -> Result<Self,CluEError> 
  {
    match value{
      Value::Array(array) => Self::from_toml_array(array),
      Value::String(filename) => { 
        let int_grid = IntegrationGrid::read_from_csv(&filename)?;
        Ok(Self::ProbabilityDistribution(int_grid))
      },
      _ => Err(CluEError::TOMLArrayDoesNotSpecifyAVector),
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_array(array: Vec::<toml::Value>) -> Result<Self,CluEError>
  {

    if array.is_empty(){
      return Err(CluEError::TOMLArrayIsEmpty);
    }
    if !are_all_same_type(&array){
      return Err(CluEError::TOMLArrayContainsMultipleTypes);
    }
    match array[0]{
      Value::Array(_) => {
        let n = array.len();
        let mut points = Vec::<f64>::with_capacity(3*n);
        let mut weights = Vec::<f64>::with_capacity(n);
        for value in array.iter(){
          let toml::Value::Array(a) = value else{
            return Err(CluEError::TOMLArrayDoesNotSpecifyAVector);
          };
          if a.len() != 4{
            return Err(CluEError::TOMLArrayDoesNotSpecifyAVector);
          }
          if !are_all_same_type(a){
            return Err(CluEError::TOMLArrayContainsMultipleTypes);
          }
          let (Some(x),Some(y),Some(z),Some(w)) 
              = (a[0].as_float(), a[1].as_float(),a[2].as_float()
              ,a[3].as_float()) else{
            return Err(CluEError::TOMLArrayDoesNotSpecifyAVector);
          };
          points.push(x);
          points.push(y);
          points.push(z);
          weights.push(w);
        }
        Ok(Self::ProbabilityDistribution(
              IntegrationGrid{dim: 3, points,weights}))
      }
      Value::Float(_) => {
        let vec: Vec::<f64> = array.iter()
          .filter_map(|v| v.as_float()).collect();

        let vec3d = Vector3D::from_vec(vec)?;

        Ok(Self::XYZ(vec3d))

      },
      Value::Integer(_) => {
        let vec: Vec::<u32> = array.iter()
          .filter_map(|v| v.as_integer()).map(|n| n as u32).collect();

        Ok(Self::CentroidOverSerials(vec))
      },
      Value::String(_) => {
        let vec: Vec::<String> = array.iter()
          .filter_map(|v| v.as_str()).map(|s| s.to_string()).collect();

        Ok(Self::CentroidOverGroup(vec))
      },
      _ => Err(CluEError::TOMLArrayDoesNotSpecifyAVector),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq)]
pub enum ReplicateUnitCell{
  No,
  Auto,
  Specified([usize;3]),
}

impl ReplicateUnitCell{
  pub fn from_toml_value(value: toml::Value) -> Result<Self,CluEError> 
  {
    match value{
      Value::Array(array) => Self::from_toml_array(array),
      Value::Boolean(b) => {
        if b { Ok(Self::Auto)}
        else { Ok(Self::No)  }
      },
      _ => Err(CluEError::TOMLArrayDoesNotSpecifyAVector),
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_array(array: Vec::<toml::Value>) -> Result<Self,CluEError>
  {

    if array.is_empty(){
      return Err(CluEError::TOMLArrayIsEmpty);
    }
    if !are_all_same_type(&array){
      return Err(CluEError::TOMLArrayContainsMultipleTypes);
    }
    match array[0]{
      Value::Integer(_) => {
        let vec: Vec::<usize> = array.iter()
          .filter_map(|v| v.as_integer()).map(|n| n as usize).collect();

        Ok(Self::Specified([vec[0],vec[1],vec[2]]))
      },
      _ => Err(CluEError::TOMLArrayDoesNotSpecifyAVector),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ClusterMethod` list different cluster simulation methods.
#[derive(Debug,Clone,PartialEq)]
pub enum ClusterMethod{
  AnalyticRestricted2CCE,
  CCE,
  GCCE,
}
impl ClusterMethod{
  pub fn from(method_str: &str) -> Result<Self,CluEError>{
    match method_str.to_lowercase().as_str(){
      "r2cce" => Ok(Self::AnalyticRestricted2CCE),
      "cce" => Ok(Self::CCE),
      "gcce" => Ok(Self::GCCE),
      _ => Err(CluEError::CannotParseClusterMethod(method_str.to_string())),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq)]
pub enum PartitionTableConfig{
  CSV(String),
  IndexBlockList(Vec::<(usize,usize)>),
}
impl PartitionTableConfig{
  pub fn from_toml_value(value: toml::Value) -> Result<Self,CluEError> 
  {
    match value{
      Value::Array(array) => Self::from_toml_array(array),
      Value::String(s) => Ok(Self::CSV(s.to_string())),
      _ => Err(CluEError::TOMLValueDoesNotSpecifyAPartitionTable),
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_array(array: Vec::<toml::Value>) -> Result<Self,CluEError>
  {

    if array.is_empty(){
      return Err(CluEError::TOMLArrayIsEmpty);
    }
    if !are_all_same_type(&array){
      return Err(CluEError::TOMLArrayContainsMultipleTypes);
    }
    match array[0]{
      Value::Integer(_) => {
        if array.len() % 2 != 0{
          return Err(CluEError::TOMLValueDoesNotSpecifyAPartitionTable);
        }
        let mut p_table = Vec::<(usize,usize)>::with_capacity(array.len()/2);
        for ii in 0..array.len()/2{
          let Some(idx) = array[2*ii].as_integer() else{
            return Err(CluEError::TOMLValueDoesNotSpecifyAPartitionTable);
          };
          let Some(blk) = array[2*ii + 1].as_integer() else{
            return Err(CluEError::TOMLValueDoesNotSpecifyAPartitionTable);
          };
         p_table.push((idx as usize,blk as usize)) 
        }

        Ok(Self::IndexBlockList(p_table))
      },
      _ => Err(CluEError::TOMLValueDoesNotSpecifyAPartitionTable),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `PulseSequence` lists the options for pulse sequences to simulate.
#[derive(Debug,Clone,PartialEq)]
pub enum PulseSequence{
  CarrPurcell(usize),
  //RefocusedEcho,
}
impl PulseSequence{
  pub fn from(pulse_seq: &str) -> Result<Self,CluEError>
{
  if pulse_seq.substring(0,3) == "cp-"{
    let Ok(n_pi) = pulse_seq.substring(3,pulse_seq.len()).parse::<usize>()else{
      return Err(CluEError::CannotParsePulseSequence(pulse_seq.to_string()));
    };
    return Ok(Self::CarrPurcell(n_pi)); 
  }
  match pulse_seq{
    "hahn" => Ok(PulseSequence::CarrPurcell(1)),
    _ => Err(CluEError::CannotParsePulseSequence(pulse_seq.to_string())),
  }
}

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `OrientationAveraging` lists the different orientation averaging options.
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
pub enum OrientationAveraging{
  Grid(IntegrationGrid),
  Lebedev(usize),
  Random(usize),
}
impl FromTOMLString for OrientationAveraging{
  fn from_toml_string(toml_str: &str) -> Result<Self,CluEError>{
    let decoded: Result<OrientationAveraging,_> = toml::from_str(toml_str);
    match decoded {
      Ok(out) => Ok(out),
      Err(err) => Err(CluEError::CannotReadTOML( format!("{}",err) )),
    }   
  }
}
impl fmt::Display for OrientationAveraging{
   fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
     let toml = toml::to_string(self).unwrap();
     write!(f,"{}",toml)
   }
}
impl OrientationAveraging{
  pub fn from_orientations_toml(ori_toml: OrientationsTOML) 
    -> Result<Self,CluEError>
  {
    let Some(grid) = ori_toml.grid else{
      return Err(CluEError::NoOrientationGrid);
    };

    if grid == KEY_ORI_LEBEDEV{

      let Some(n_ori) = ori_toml.number else{
        return Err(CluEError::CannotParseOrientations(grid.to_string()));
      };
      return Ok(OrientationAveraging::Lebedev(n_ori));

    } else if grid == KEY_ORI_RANDOM{

      let Some(n_ori) = ori_toml.number else{
        return Err(CluEError::CannotParseOrientations(grid.to_string()));
      };
      return Ok(OrientationAveraging::Random(n_ori));

    } else if grid == KEY_ORI_FILE{
      
      let Some(filename) = &ori_toml.file else{
        return Err(CluEError::CannotParseOrientations(grid.to_string()));
      };
      let int_grid = IntegrationGrid::read_from_csv(filename)?;
      return Ok(OrientationAveraging::Grid(int_grid));

    } else if grid == KEY_ORI_VECTOR{
      let Some(v) = ori_toml.vector else{
        return Err(CluEError::CannotParseOrientations(grid.to_string()));
      };
      if v.len() != 3{
        return Err(CluEError::CannotParseOrientations(grid.to_string()));
      }
      return Ok(OrientationAveraging::Grid(IntegrationGrid{
          dim: 3,
          points: v,
          weights: vec![1.0],
      }));


    } else if grid == KEY_ORI_VECTORGRID{
      let Some(v_grid) = ori_toml.vector_grid else{
        return Err(CluEError::CannotParseOrientations(grid.to_string()));
      };
      let mut points = Vec::<f64>::with_capacity(3*v_grid.len());
      let mut weights = Vec::<f64>::with_capacity(v_grid.len());
      for v in v_grid.iter(){
        if v.len() != 4{
          return Err(CluEError::CannotParseOrientations(grid.to_string()));
        }
        points.push(v[0]);
        points.push(v[1]);
        points.push(v[2]);
        weights.push(v[3]);
      }

      return Ok(OrientationAveraging::Grid(IntegrationGrid{
          dim: 3,
          points,
          weights,
      }));

    } else{
      return Err(CluEError::CannotParseOrientations(grid.to_string()));
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl Config{
  //----------------------------------------------------------------------------
  /// This method modifies a `Config` from a `ConfigTOML`.
  /// For field that is `Some` in the `ConfigTOML`, 
  /// the corresponding field in the `Config` will be overwriten,
  /// but the `None`-valued fields in the `ConfigTOML` will not alter the
  /// `Config`.
  pub fn set_from_config_toml(&mut self, mut config_toml: ConfigTOML) 
    -> Result<(),CluEError>
  {
    // Get units. 
    let unit_of_distance = match &config_toml.unit_of_distance{
      Some(u) => distance_unit_to_meters(u)?,
      None => return Err(CluEError::NoUnitOfDistance),
    }; 
    let unit_of_energy = match &config_toml.unit_of_energy {
      Some(u) => energy_unit_to_hertz(u)?,
      None => return Err(CluEError::NoUnitOfEnergy),
    }; 
    let unit_of_magnetic_field = match &config_toml.unit_of_magnetic_field{
      Some(u) => magnetic_field_unit_to_tesla(u)?,
      None => return Err(CluEError::NoUnitOfMagneticField),
    }; 
    let unit_of_time = match &config_toml.unit_of_time {
      Some(u) => time_unit_to_seconds(u)?,
      None => return Err(CluEError::NoUnitOfTime),
    };
    self.unit_of_time_to_seconds = Some(unit_of_time);
 

    // A--C
    if let Some(clash_distance) = config_toml.clash_distance{
      self.clash_distance = Some(clash_distance*unit_of_distance);
    }

    if let Some(clash_distance_pbc) = config_toml.clash_distance_pbc{
      self.clash_distance_pbc = Some(clash_distance_pbc*unit_of_distance);
    }

    if config_toml.cluster_batch_size.is_some(){
      self.cluster_batch_size = config_toml.cluster_batch_size;
    }


    if let Some(clu_src) = config_toml.cluster_source{
      if clu_src.ends_with(".toml") || clu_src.ends_with(".txt"){
        self.cluster_source = Some(ClusterSource::File(clu_src.to_string()));
      }else {
        self.cluster_source 
          = Some(ClusterSource::OrientationsDir(clu_src.to_string()));
      }
    }

    if let Some(method_str) = &config_toml.cluster_method{
      let cluster_method = ClusterMethod::from(method_str)?;
      self.cluster_method = Some(cluster_method);
    }

    // D
    if let Some(mut detected_spin) = config_toml.detected_spin{
      if let Some(toml_value) = detected_spin.g_matrix{
        self.detected_spin_g_matrix 
          = Some(TensorSpecifier::from_toml_value(toml_value,1.0)?);
      }
      
      if let Some(spin_id) = &detected_spin.identity{                          
        self.detected_spin_identity =  Some(Isotope::from(spin_id)?);
      }
      if let Some(spin_mult) = detected_spin.multiplicity{
        self.detected_spin_multiplicity = Some(spin_mult);
      }

      if let Some(toml_value) = detected_spin.position{
        self.detected_spin_position 
          = Some(DetectedSpinCoordinates::from_toml_value(toml_value)? );
      }

      if detected_spin.transition.is_some(){
        mem::swap(
            &mut self.detected_spin_transition, 
            &mut detected_spin.transition);
      }
    
    }
    //E--G
    if let Some(groups) = config_toml.groups{
      for value in groups.iter(){
        let Some(group) = value.as_table() else{
          return Err(CluEError::ExpectedTOMLTable(
                value.type_str().to_string()));
        };
        set_particle_configs_from_toml_table(&mut self.particles,
            group.clone(),unit_of_distance,unit_of_energy)?;
      }
    }
    //I
    if config_toml.input_structure_file.is_some(){
      self.input_structure_file = config_toml.input_structure_file;
    }


    // J--M
    if let Some(bz) = config_toml.magnetic_field{
     let mag_field = Vector3D::from([0.0,0.0,bz*unit_of_magnetic_field]);
     self.magnetic_field = Some(mag_field);
    }

    if config_toml.max_cell_size.is_some(){
      self.max_cell_size = config_toml.max_cell_size;
    }

    if config_toml.max_cluster_size.is_some(){
      self.max_cluster_size = config_toml.max_cluster_size;
    }

    if config_toml.max_spins.is_some(){
      self.max_spins = config_toml.max_spins;
    }

    if config_toml.min_cell_size.is_some(){
      self.min_cell_size = config_toml.min_cell_size;   
    }


    // N
    if let Some(pair_cutoffs) = &config_toml.pair_cutoffs{
      if let Some(&cutoff) = pair_cutoffs.get(KEY_CUTOFF_COUPLING){ 
        self.neighbor_cutoff_coupling_xx_yy = Some(cutoff*unit_of_energy);
      }
      if let Some(&cutoff) = pair_cutoffs.get(KEY_CUTOFF_DELTA_HF){ 
        self.neighbor_cutoff_delta_hyperfine_zz = Some(cutoff*unit_of_energy);
      }
      if let Some(&cutoff) = pair_cutoffs.get(KEY_CUTOFF_DIPOLE_PERP){ 
        self.neighbor_cutoff_point_dipole_perpendicular 
            = Some(cutoff*unit_of_energy);
      }
      if let Some(&cutoff) = pair_cutoffs.get(KEY_CUTOFF_DISTANCE){ 
        self.neighbor_cutoff_distance = Some(cutoff*unit_of_distance);
      }
      if let Some(&cutoff) = pair_cutoffs.get(KEY_CUTOFF_HAHN_MOD_DEPTH){ 
        self.neighbor_cutoff_3_spin_hahn_mod_depth = Some(cutoff)
      }
      if let Some(&cutoff) = pair_cutoffs.get(KEY_CUTOFF_HAHN_TAYLOR_4){ 
        // The expected unit of neighbor_cutoff_3_spin_hahn_taylor_4 is
        // (rad/s)^4, while the input unit is in unit_of_energy^4.

        // Define the converion unit_of_energy -> (rad/s).
        let u = unit_of_energy*HZ_TO_RAD_PER_S;

        // Define the converion unit_of_energy^4 -> (rad/s)^4.
        let u4 = u*u*u*u;

        self.neighbor_cutoff_3_spin_hahn_taylor_4 = Some(cutoff*u4)
      }

    }

    if config_toml.number_runs.is_some(){
      self.number_runs = config_toml.number_runs;
    }

    if let Some(number_timepoints) = &mut config_toml.number_timepoints{
      mem::swap(&mut self.number_timepoints, number_timepoints);
    }

    // O
    if let Some(ori_toml) = config_toml.orientations{
      self.orientation_grid 
          = Some(OrientationAveraging::from_orientations_toml(ori_toml)?); 
    }

    // P--Q
    if config_toml.partitioning 
        == Some(KEY_PARTITION_PARTICLE.to_string())
    {
      self.partitioning = Some(PartitioningMethod::Particles);

    }else if config_toml.partitioning 
        == Some(KEY_PARTITION_EX_GROUPS.to_string())
    {
      self.partitioning 
          = Some(PartitioningMethod::ExchangeGroupsAndParticles);

    }else if let Some(s) = config_toml.partitioning
    {
      return Err(CluEError::CannotParsePartitioningMethod(0,s.to_string()));  
    }

    if let Some(p_table) = config_toml.partition_table{
      self.partition_table 
          = Some(PartitionTableConfig::from_toml_value(p_table)?); 
    }

    if config_toml.pdb_model_index.is_some(){
      self.pdb_model_index = config_toml.pdb_model_index;
    }
    if let Some(pulse_seq) = &config_toml.pulse_sequence{
      self.pulse_sequence = Some(PulseSequence::from(pulse_seq)?);
    }

    // R--S
    if let Some(radius) = config_toml.radius{
      self.radius = Some(radius*unit_of_distance);
    }

    if let Some(pbc) = config_toml.replicate_unit_cell{
      self.replicate_unit_cell = Some(ReplicateUnitCell::from_toml_value(pbc)?);
    }

    if config_toml.rng_seed.is_some(){
      self.rng_seed = config_toml.rng_seed;
    }

    if config_toml.output_directory.is_some(){
      mem::swap(&mut self.output_directory, &mut config_toml.output_directory);
    }

    if config_toml.run_name.is_some(){
      mem::swap(&mut self.run_name, &mut config_toml.run_name);
    }

    // T--V
    if config_toml.populations 
        == Some(KEY_DENSITY_MATRIX_THERMAL.to_string())
    {
      let Some(t) = config_toml.temperature else {
        return Err(CluEError::NoTemperature);
      };
      self.density_matrix = Some(DensityMatrixMethod::Thermal(t));  
    }else if config_toml.populations 
        == Some(KEY_DENSITY_MATRIX_ID.to_string())
    {
      self.density_matrix = Some(DensityMatrixMethod::Identity); 
    }


    if let Some(tau_increments) = &mut config_toml.tau_increments{ 
      for t in tau_increments.iter_mut(){
        *t *= unit_of_time;
      }
      mem::swap(&mut self.tau_increments, tau_increments);
    }

    // W--Z
    if let Some(output) = config_toml.output{
      if let Some(&b) = output.get(KEY_OUT_AUX_SIGS){ 
        self.write_auxiliary_signals = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_BATH){ 
        self.write_bath = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_DET_SPIN){ 
        self.write_detected_spin = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_CLUSTERS){ 
        self.write_clusters = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_CONFIG){ 
        self.write_config = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_INFO){ 
        self.write_info = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_EXCHANGE_GROUPS){ 
        self.write_exchange_groups = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_METHYL_PARTITIONS){ 
        self.write_methyl_partitions = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_ORI_GRID){ 
        self.write_orientation_grid = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_ORI_SIGS){ 
        self.write_orientation_signals = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_SANS_SPIN_SIGS){ 
        self.write_sans_spin_signals = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_STRUC_PDB){ 
        self.write_structure_pdb = Some(b)
      }
      if let Some(&b) = output.get(KEY_OUT_TENSORS){ 
        self.write_tensors = Some(b)
      }
    }

    Ok(())
  }
  //----------------------------------------------------------------------------

  /// This function takes the parsed command line input and develops the main
  /// configuration data structure.
  pub fn read_input(input: CommandLineInput) -> Result<Self,CluEError>{
    let Some(filename) = &input.config_file else{
      return Err(CluEError::NoInputFile);
    };

    return Self::from_toml_file(filename); 
  }
  //----------------------------------------------------------------------------
  // This functions sets properties for particles outside the primary cell.
  fn set_extracell_particles(&mut self){
    
    let mut n_extra = 0; 
    let mut n_primary = 0; 
    for particle_config in self.particles.iter(){
      match particle_config.cell_type{
      CellType::AllCells => {
        n_extra += 1;
        n_primary += 1;
      },
      CellType::PrimaryCell => n_primary += 1,
      CellType::Extracells => n_extra += 1,
      }
    }

    let mut particles = Vec::<ParticleConfig>::with_capacity(n_primary);
    self.extracell_particles = Vec::<ParticleConfig>::with_capacity(n_extra);

    for particle_config in self.particles.iter(){
      match particle_config.cell_type{
      CellType::AllCells => {
        self.extracell_particles.push(particle_config.clone());
        particles.push(particle_config.clone());
      },
      CellType::PrimaryCell => particles.push(particle_config.clone()),
      CellType::Extracells 
        => self.extracell_particles.push(particle_config.clone()),
      }
    }

    self.particles = particles;

  }
  //----------------------------------------------------------------------------
  /*
  /// This function updates config from a `Vec<TokenExpression>`.
  pub fn parse_token_stream(&mut self,token_stream: Vec<TokenExpression>) 
    -> Result<(),CluEError>
  {
    let mut mode = ModeAttribute::new();

    for expression in token_stream.iter(){
      if expression.lhs.is_empty(){ continue; }

      if let Token::Mode(new_mode) = &expression.lhs[0]{
        mode = new_mode.clone();
        continue;
      }

      match mode.mode{
        ConfigMode::Config =>  self.parse_config_line(expression)?,
        ConfigMode::Filter => {
          let Some(_label) = &mode.label else{
            return Err(CluEError::FilterNeedsALabel);
          };
          self.parse_filter_line(expression,&mode.label)?;
        },
        ConfigMode::SpinProperties => {
          let Some(label) = &mode.label else{
            return Err(CluEError::SpinPropertiesNeedsALabel);
          };
          if mode.isotope.is_none(){
            return Err(CluEError::SpinPropertiesNeedsAnIsotope(label.clone()));
          }
          self.parse_properties_line(expression,&mode.label,mode.isotope)?;
        },
        ConfigMode::StructureProperties => {
          let Some(label) = &mode.label else{
            return Err(CluEError::StructurePropertiesNeedsALabel);
          };
          if mode.isotope.is_some(){
            return Err(
               CluEError::StructurePropertiesDoesNotNeedAnIsotope(
                 label.clone()));
          }
          self.parse_properties_line(expression,&mode.label,mode.isotope)?;
        },
      }

    }


    Ok(())
  }
  */
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[cfg(test)]
mod tests{
  use super::*;
  use crate::elements::Element;
  use crate::structure::particle_filter::SecondaryParticleFilter;

  //----------------------------------------------------------------------------
  #[test]
  fn test_config_from_toml(){
  
    let config = Config::from_toml_string(r##"
      clash_distance = 0.1
      clash_distance_pbc = 0.1
      cluster_batch_size = 20000
      populations = "thermal"
      cluster_method = "cce"
      cluster_source = "clusters_file.toml"
      input_structure_file = "../../assets/TEMPO_wat_gly_70A.pdb"
      magnetic_field = 1.2 # T
      max_cell_size = 2
      max_cluster_size = 4
      max_spins = 8
      min_cell_size = 1
      number_runs = 2
      replicate_unit_cell = false    
      partitioning = "exchange_groups"
      partition_table = [
        1,1,
        2,1,
        3,2,
        4,2,
      ]
      radius = 80 # Å
      rng_seed = 42
      output_directory = "save_directory"
      # TODO: remove run_name and automate hash remove
      # TODO: rework hash vs run number system
      run_name = "system"
      temperature = 20

      pulse_sequence = "cp-1"
      tau_increments = [1e-3,5e-1] # μs
      number_timepoints = [40,60]


      [pair_cutoffs]
        coupling_xx_yy = 1e-3 # MHz
        delta_hyperfine_zz = 1e-2 # MHz
        point_dipole_perpendicular = 1e-4 # MHz
        distance = 10 # Å
        hahn_mod_depth = 1e-10 # unitless
        hahn_taylor_4 = 6.416238909177711e-37 # MHz^4

      [orientations]
        grid = "lebedev"
        number = 170

      [detected_spin]
      position = [28,29] # centroid over serials
      g_matrix.values = [2.0097, 2.0064, 2.0025]
      g_matrix.axes.x = {from = "tempo_c1", to = "tempo_c19"}
      g_matrix.axes.y = [-1.1500, -0.4700, 0.7100]

      [output]
      auxiliary_signals = true
      bath = true
      clusters = true
      exchange_groups = true
      info = true
      methyl_partitions = true
      orientation_signals = true
      sans_spin_signals = true
      structure_pdb = true
      tensors = false

      [[groups]]
      name = "tempo_c1"
      selection = {elements = ["C"], serials = [1]}

      [[groups]]
      name = "tempo_c19"
      selection = {elements = ["C"], serials = [19]}

      [[groups]]
      name = "tempo_n"
      selection = {elements = ["N"]}

      [groups.14N]
      hyperfine.values = [20.0,20.0,100.0] # MHz
      hyperfine.axes.x = {to_bonded_to = "tempo_o" }
      hyperfine.axes.y = { from_bonded_to="tempo_c1", to_bonded_to="tempo_c19"}

      electric_quadrupole.values = [-0.28, -1.47, 1.75] # MHz
      electric_quadrupole.axes.x = {to_bonded_to = "tempo_o" }
      electric_quadrupole.axes.y = {from_bonded_to="tempo_c1", to_bonded_to="tempo_c19"}

      [[groups]]
        name = "tempo_o"
        selection = {elements = ["O"]}

      [[groups]]
        name = "tempo_h"
        selection = {elements = ["H"]}
        1H.c3_tunnel_splitting = 0.120 # MHz

      [[groups]]
        name = "glycerol_oh"
        drop_probability = 0.5
        [groups.selection]
          elements = ["H"]
          residues = ["GLY"]
          bonded_elements = ["O"]

      [[groups]]
        name = "glycerol_ch"
        
        cosubstitute = "same_molecule"

        1H.abundace = 0.5
        2H.abundace = 0.5

        [groups.selection]
          elements = ["H"]
          residues = ["GLY"]
          bonded_elements = ["C"]

    "##).unwrap();

    assert_eq!(config.replicate_unit_cell, Some(ReplicateUnitCell::No));
    let clash_distance = config.clash_distance.unwrap();
    assert!( (clash_distance - 0.1e-10).abs()/(0.1e-10) < 1e-12 );
    let clash_distance_pbc = config.clash_distance_pbc.unwrap();
    assert!( (clash_distance_pbc - 0.1e-10).abs()/(0.1e-10) < 1e-12 );
    assert_eq!(config.cluster_batch_size,Some(20000));
    assert_eq!(config.density_matrix, Some(DensityMatrixMethod::Thermal(20.0)));
    assert_eq!(config.cluster_method,Some(ClusterMethod::CCE));

    let g_matrix = match config.detected_spin_g_matrix.unwrap(){
      TensorSpecifier::Eig(eig_specifier) => eig_specifier.clone(),
      _ => panic!("Expected TensorSpecifier::Eig(eig_specifier)."),
    };

    assert_eq!(g_matrix.values,
        Some([2.0097, 2.0064, 2.0025] ) );
    assert_eq!(g_matrix.z_axis, None);
    assert_eq!(g_matrix.y_axis,
        Some(VectorSpecifier::Vector(
            Vector3D::from([-1.1500, -0.4700, 0.7100]) ) ) );
    assert_eq!(g_matrix.x_axis,
        Some(VectorSpecifier::Diff(
            SecondaryParticleFilter::Filter, "tempo_c1".to_string(),
            SecondaryParticleFilter::Filter, "tempo_c19".to_string(),
            )));

    assert_eq!(config.detected_spin_identity, Some(Isotope::Electron));
    assert_eq!(config.detected_spin_multiplicity, Some(2));
    assert_eq!(config.detected_spin_position,
        Some(DetectedSpinCoordinates::CentroidOverSerials(vec![28,29])));
    assert_eq!(config.input_structure_file,
         Some("../../assets/TEMPO_wat_gly_70A.pdb".to_string()));
    assert_eq!(config.cluster_source,
         Some(ClusterSource::File("clusters_file.toml".to_string())));

    assert_eq!(config.max_cell_size, Some(2));
    assert_eq!(config.magnetic_field, Some(Vector3D::from([0.0,0.0,1.2])));
    assert_eq!(config.max_cluster_size, Some(4));
    assert_eq!(config.max_spins, Some(8));
    assert_eq!(config.min_cell_size, Some(1));
    assert_eq!(config.neighbor_cutoff_delta_hyperfine_zz, Some(1e4));
    assert_eq!(config.neighbor_cutoff_coupling_xx_yy, Some(1e3));
    assert_eq!(config.neighbor_cutoff_point_dipole_perpendicular, Some(1e2));
    assert_eq!(config.neighbor_cutoff_3_spin_hahn_mod_depth, Some(1e-10));
    assert_eq!(config.neighbor_cutoff_3_spin_hahn_taylor_4, Some(1e-9));
    assert_eq!(config.number_runs, Some(2));
    assert_eq!(config.number_timepoints, vec![40,60]);
    assert_eq!(config.partitioning, 
        Some(PartitioningMethod::ExchangeGroupsAndParticles));
    assert_eq!(config.partition_table, 
        Some(PartitionTableConfig::IndexBlockList(vec![
        (1,1),(2,1),(3,2),(4,2)])));
    assert_eq!(config.pulse_sequence, Some(PulseSequence::CarrPurcell(1)));
    assert_eq!(config.radius, Some(80.0e-10));
    assert_eq!(config.rng_seed, Some(42));
    assert_eq!(config.output_directory, Some(String::from("save_directory")));
    assert_eq!(config.run_name, Some(String::from("system")));
    assert_eq!(config.tau_increments, vec![1e-9,5e-7]);
    assert_eq!(config.write_auxiliary_signals, Some(true));
    assert_eq!(config.write_bath, Some(true));
    assert_eq!(config.write_clusters, Some(true));
    assert_eq!(config.write_exchange_groups, Some(true));
    assert_eq!(config.write_info, Some(true));
    assert_eq!(config.write_methyl_partitions, Some(true));
    assert_eq!(config.write_orientation_signals, Some(true));
    assert_eq!(config.write_structure_pdb, Some(true));
    assert_eq!(config.write_sans_spin_signals, Some(true));
    assert_eq!(config.write_tensors, Some(false));


    let groups = config.particles;
    assert_eq!(groups.len(),7);
    assert_eq!(groups[0].label, "tempo_c1".to_string());
    let filter = groups[0].filter.clone().unwrap();
    assert_eq!(filter.elements, vec![Element::Carbon]);
    assert_eq!(filter.serials, vec![1]);

    assert_eq!(groups[1].label, "tempo_c19".to_string());
    let filter = groups[1].filter.clone().unwrap();
    assert_eq!(filter.elements, vec![Element::Carbon]);
    assert_eq!(filter.serials, vec![19]);

    assert_eq!(groups[2].label, "tempo_n".to_string());
    let filter = groups[2].filter.clone().unwrap();
    assert_eq!(filter.elements, vec![Element::Nitrogen]);
    assert!(filter.serials.is_empty());
    let properties = groups[2].properties.clone().unwrap();
    let nitrogen = &properties
        .isotope_properties[&Isotope::Nitrogen14.to_string()];

    let x = VectorSpecifier::Diff(
        SecondaryParticleFilter::Particle, "self".to_string(),
        SecondaryParticleFilter::Bonded, "tempo_o".to_string());
    let y = VectorSpecifier::Diff(
        SecondaryParticleFilter::Bonded, "tempo_c1".to_string(),
        SecondaryParticleFilter::Bonded, "tempo_c19".to_string());

    assert_eq!(nitrogen.hyperfine_coupling,
        Some(TensorSpecifier::Eig(EigSpecifier{
          values: Some([20.0e6, 20.0e6, 100.0e6]),
          x_axis: Some(x.clone()),    
          y_axis: Some(y.clone()),    
          z_axis: None,    
        })));

    assert_eq!(nitrogen.electric_quadrupole_coupling,
        Some(TensorSpecifier::Eig(EigSpecifier{
          values: Some([-0.28e6, -1.47e6, 1.75e6]),
          x_axis: Some(x.clone()),    
          y_axis: Some(y.clone()),    
          z_axis: None,    
        })));

    assert_eq!(groups[3].label, "tempo_o".to_string());
    let filter = groups[3].filter.clone().unwrap();
    assert_eq!(filter.elements, vec![Element::Oxygen]);
    assert!(filter.serials.is_empty());

    assert_eq!(groups[4].label, "tempo_h".to_string());
    let filter = groups[4].filter.clone().unwrap();
    assert_eq!(filter.elements, vec![Element::Hydrogen]);
    assert!(filter.serials.is_empty());
    let properties = groups[4].properties.clone().unwrap();
    let hydrogen = &properties
        .isotope_properties[&Isotope::Hydrogen1.to_string()];
    let ex1 = hydrogen.exchange_coupling.clone().unwrap();
    let ex0 = -80.0e3;
    assert!(2.0*(ex0-ex1).abs()/(ex0+ex1).abs() < 1e-12);

    assert_eq!(groups[5].label, "glycerol_oh".to_string());
    let filter = groups[5].filter.clone().unwrap();
    assert_eq!(filter.elements, vec![Element::Hydrogen]);
    assert_eq!(filter.residues, vec!["GLY".to_string()]);
    assert_eq!(filter.bonded_elements, vec![Element::Oxygen]);
    let properties = groups[5].properties.clone().unwrap();
    assert_eq!(properties.isotopic_distribution.void_probability,
        Some(0.5));

    assert_eq!(groups[6].label, "glycerol_ch".to_string());
    let filter = groups[6].filter.clone().unwrap();
    assert_eq!(filter.elements, vec![Element::Hydrogen]);
    assert_eq!(filter.residues, vec!["GLY".to_string()]);
    assert_eq!(filter.bonded_elements, vec![Element::Carbon]);
    let properties = groups[6].properties.clone().unwrap();
    assert_eq!(properties.cosubstitute,
        Some(SecondaryParticleFilter::SameMolecule));


  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_set_time_axis(){
    let mut config = Config::new();
    config.pulse_sequence = Some(PulseSequence::CarrPurcell(1));
    config.number_timepoints = vec![100,91];
    config.tau_increments = vec![0.5e-8, 0.5e-7];
    config.set_time_axis().unwrap();
    let time_axis = config.get_time_axis().unwrap();
    assert_eq!(time_axis.len(),191);
    assert_eq!(time_axis[0],0.0);
    let t = 1e-8;
    assert!((time_axis[1]-t).abs()/t<1e-12);
    let t = 1e-6;
    assert!((time_axis[100]-t).abs()/t<1e-12);
    let t = 1.1e-6;
    assert!((time_axis[101]-t).abs()/t<1e-12);
    let t = 10e-6;
    assert!((time_axis[190]-t).abs()/t<1e-12);
    
  }
}

