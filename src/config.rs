use crate::clue_errors::*;

use crate::cluster::{
  partition::PartitioningMethod,
};
use crate::config::config_toml::*;
use crate::config::particle_config::set_particle_configs_from_toml_table;
use crate::config::command_line_input::CommandLineInput;
use crate::config::particle_config::{
  CellType, ParticleConfig,EigSpecifier,TensorSpecifier,
  IsotopeAbundance, DEFAULT_ABUNDANCE_TOLERANCE
};
use crate::config::pulse_sequence::*;
use crate::isotopes::Isotope;
use crate::physical_constants::*;
use crate::space_3d::Vector3D;
use crate::structure::particle_filter::VectorSpecifier;
use crate::integration_grid::IntegrationGrid;

use crate::quantum::cluster_operators::SpinOp;

use crate::misc::{
  are_all_same_type,
  cxmat_from_toml_array,  
  vec_usize_from_toml_array,
};


use crate::io::{FromTOMLString,read_csv};



use toml::Value;
use serde::{Serialize,Deserialize};

use std::collections::HashMap;
use std::fmt;
use std::fs;
use std::mem;

use num_complex::Complex64;
use ndarray::Array2;

type CxMat = Array2::<Complex64>;

pub mod command_line_input;
pub mod config_toml;
pub mod particle_config;
pub mod toml_keys;
pub mod pulse_sequence;
pub mod validate;


// Define the directory and file names used to save outputs.
pub const SAVE_DIR_AUXILIARY_SIGNALS: &str = "auxiliary_signals";
pub const SAVE_DIR_CLUSTER_SIGNALS: &str = "cluster_signals";
pub const SAVE_FILE_BATH: &str = "bath";
pub const SAVE_FILE_CLUSTERS: &str = "clusters";
pub const SAVE_DIR_INFO: &str = "info";
pub const SAVE_FILE_DET_SPIN_COOR: &str = "detected_spin";
pub const SAVE_FILE_EXCHANGE_GROUPS: &str = "exchange_groups";
pub const SAVE_FILE_METHYL_PARTITIONS: &str = "methyl_partitions";
pub const SAVE_DIR_ORIENTATION_SIGNALS: &str = "orientations";
pub const SAVE_FILE_ORIENTATIONS: &str = "orientation_grid";
pub const SAVE_FILE_PARTITION_TABLE: &str = "partition_table";
pub const SAVE_FILE_SANS_SPIN_SIGNALS: &str = "sans_spin_signals";
pub const SAVE_FILE_STRUCTURE_PDB: &str = "structure_pdb";
pub const SAVE_FILE_TENSORS: &str = "tensors";

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// Config contains all the setting for CluE.
#[derive(Debug,Clone,Default)]
pub struct Config{

  // Detected Spin
  pub cluster_populations: Option<ClusterPopulations>,
  pub detected_population: Option<DetectedPopulation>,
  pub thermalize_cluster_populations: Option<bool>,
  pub detected_spin_frame: Option<DetFrame>,
  pub detected_spin_g_matrix: Option<TensorSpecifier>,
  pub detected_spin_identity: Option<Isotope>,
  pub detected_spin_multiplicity: Option<usize>,
  pub detected_spin_position: Option<DetectedSpinCoordinates>,
  pub detected_spin_weight_tolerance: Option<f64>,
  pub detected_spin_transition: Option<[usize;2]>,
  pub detection_operator: Option<DetectionOp>,
  pub detected_spin_zerofield_coupling: Option<TensorSpecifier>,
  pub detected_spin_quadrupole_coupling: Option<TensorSpecifier>,
  pub mean_fields: Option<bool>,
  pub pulses: HashMap::<String,CxMat>,

  // Structure
  pub clash_distance: Option<f64>, 
  pub clash_distance_pbc: Option<f64>,
  pub extracell_particles: Vec::<ParticleConfig>,
  pub input_structure_file: Option<String>,
  pub max_spherical_tensor_rank: Option<usize>,
  pub particles: Vec::<ParticleConfig>,
  pub pdb_model_index: Option<usize>,
  pub replicate_unit_cell: Option<ReplicateUnitCell>,
  pub radius: Option<f64>,
  pub temperature: Option<f64>,
  
  // Pair Cutoffs
  pub neighbor_cutoff_coupling_xx_yy: Option<f64>,
  pub neighbor_cutoff_delta_hyperfine_zz: Option<f64>,
  pub neighbor_cutoff_point_dipole_perpendicular: Option<f64>,
  pub neighbor_cutoff_distance: Option<f64>,
  pub neighbor_cutoff_3_spin_hahn_mod_depth: Option<f64>,
  pub neighbor_cutoff_3_spin_hahn_taylor_4: Option<f64>,
  pub neighbor_cutoff_delta_zeeman: Option<f64>,
  pub neighbor_cutoff_average_zeeman: Option<f64>,

  // CCE
  pub cluster_batch_size: Option<usize>, 
  pub cluster_source: Option<ClusterSource>,
  pub cluster_method: Option<ClusterMethod>,
  pub ensemble_cce: Option<bool>,
  pub connect_exchange_groups: Option<bool>,
  pub min_cell_size: Option<usize>,
  pub max_cell_size: Option<usize>,
  pub max_cluster_size: Option<usize>,
  pub max_spins: Option<usize>,
  pub kmeans_iter_cap: Option<usize>,
  pub partitioning: Option<PartitioningMethod>,
  pub partition_table: Option<PartitionTableConfig>,
  pub rng_seed: Option<u64>,
  pub run_in_parallel: Option<bool>,

  // Experimental Details
  pub magnetic_field: Option<Vector3D>,
  pub orientation_grid: Option<OrientationAveraging>,
  pub orientation_weight_tolerance: Option<f64>,
  pub pulse_sequence: Option<PulseSequence>,
  tau_axis: Vec::<f64>,
  tau2_axis: Vec::<f64>,
  total_number_timepoints: usize,
  pub tau_increments: Vec::<f64>,
  pub tau2_increments: Vec::<f64>,
  pub number_runs: Option<usize>, 
  pub number_timepoints: Vec::<usize>,
  pub number_timepoints2: Vec::<usize>,

  // Output
  pub output_directory: Option<String>,
  pub run_name: Option<String>,
  pub write_auxiliary_signals: Option<bool>, 
  pub write_bath: Option<bool>,
  pub write_clusters: Option<bool>,
  pub write_cluster_signals: Option<bool>, 
  pub write_config: Option<bool>,
  pub write_detected_spin: Option<bool>,
  pub write_info: Option<bool>,
  pub write_exchange_groups: Option<bool>,
  pub write_methyl_partitions: Option<bool>,
  pub write_orientation_grid: Option<bool>,
  pub write_orientation_signals: Option<bool>,
  pub write_partition_table: Option::<bool>,
  pub write_sans_spin_signals: Option<bool>,
  pub write_structure_pdb: Option<bool>,
  pub write_tensors: Option<bool>,

  // Units
  pub unit_of_time_to_seconds: Option<f64>,
}


/// The default allowed difference between one and the sum of the weights of 
/// a detected spin probability grid.
pub const DEFAULT_SPIN_GRID_WEIGHT_TOLERANCE: f64 = 1e-9;

/// The default allowed difference between one and the sum of the weights of 
/// an orientation grid.
pub const DEFAULT_ORIENTATION_GRID_WEIGHT_TOLERANCE: f64 = 1e-9;

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
      Err(err) => return Err(CluEError::CannotOpenFile(
        filename.to_string(), err.to_string())),
    };

    let mut config = Self::from_toml_string(&toml_str)?;

    match config.set_tau_axis(){
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
    /*
    // Slow for large systems.   
    if self.clash_distance_pbc.is_none(){
      self.clash_distance_pbc = Some(1e-11);
    }
    */


    if self.cluster_batch_size.is_none(){
      self.cluster_batch_size = Some(1000);
    }

    if self.cluster_source.is_none(){
      self.cluster_source = Some(ClusterSource::Structure);
    }

    if self.ensemble_cce.is_none(){
      self.ensemble_cce = Some(false);
    }

    if self.max_spherical_tensor_rank.is_none(){
      self.max_spherical_tensor_rank = Some(0);
    }

    if self.mean_fields.is_none(){
      self.mean_fields = Some(false);
    }

    if self.partitioning.is_none(){
      self.partitioning = Some(PartitioningMethod::Particles);
    }

    if self.pdb_model_index.is_none(){
      self.pdb_model_index = Some(0);
    }
    if self.run_in_parallel.is_none(){
      self.run_in_parallel = Some(true);
    }

    if self.cluster_populations. is_none(){
      self.cluster_populations = Some(ClusterPopulations::Uniform);
    }

    if self.detected_spin_frame.is_none(){
      self.detected_spin_frame = Some(DetFrame::Eigen);
    }

    if self.neighbor_cutoff_delta_zeeman.is_none(){
      self.neighbor_cutoff_delta_zeeman = Some(1e-12);
    }

    if self.neighbor_cutoff_average_zeeman.is_none(){
      self.neighbor_cutoff_average_zeeman = Some(0.0);
    }

    if self.number_runs.is_none(){
      self.number_runs = Some(1);
    }

    if self.kmeans_iter_cap.is_none(){
      self.kmeans_iter_cap = Some(1000);
    }

    if self.number_timepoints2.is_empty(){
      self.number_timepoints2 = vec![1];
    }
    if self.tau2_increments.is_empty(){
      self.tau2_increments = vec![0.0];
    }

    if self.thermalize_cluster_populations.is_none(){
      self.thermalize_cluster_populations = Some(false);
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
    if self.detected_spin_multiplicity.is_none(){
      let Some(spin) = self.detected_spin_identity else{
        return Err(CluEError::NoDetectedSpinIdentity);
      };
      self.detected_spin_multiplicity = Some(spin.spin_multiplicity());
    }

    // Set g-matrix
    if self.detected_spin_g_matrix.is_none(){
      self.detected_spin_g_matrix = Some(TensorSpecifier::Eig(EigSpecifier{
        values: Some([ELECTRON_G,ELECTRON_G,ELECTRON_G]),
        x_axis: Some(VectorSpecifier::Vector(Vector3D::from([1.0, 0.0, 0.0]))),
        y_axis: Some(VectorSpecifier::Vector(Vector3D::from([0.0, 1.0, 0.0]))),
        z_axis: None,
        spherical_tensor_coefficients: Vec::<CxMat>::new(),
        stevens: false,
      }));
    }
    if self.detected_spin_multiplicity.is_none() {
      self.detected_spin_multiplicity = match self.detected_spin_identity{
        Some(isotope) => Some(isotope.spin_multiplicity()),
        None => return Err(CluEError::NoDetectedSpinIdentity),
      };
    };

    if self.detected_spin_transition.is_none(){
      let Some(s) = self.detected_spin_multiplicity else{
        return Err(CluEError::NoDetectedSpinMultiplicity);
      };

      if s == 2{
        self.detected_spin_transition = Some([s/2 - 1,s/2]);
      }
    }

    self.check_particle_configs()?;

    self.check_detected_spin_grid()?;

    self.check_orientation_grid()?;

    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function checks that the weights of a detected spin probability
  // grid sum to one within detected_spin.weight_tolerance, and normalizes
  // them.
  // A message is printed when normalizing changes the weights by more 
  // than floating point rounding.
  fn check_detected_spin_grid(&mut self) -> Result<(),CluEError>{
    let Some(DetectedSpinCoordinates::ProbabilityDistribution(grid)) 
      = &mut self.detected_spin_position else{
      return Ok(());
    };

    let tolerance = self.detected_spin_weight_tolerance
      .unwrap_or(DEFAULT_SPIN_GRID_WEIGHT_TOLERANCE);
    if !(0.0..1.0).contains(&tolerance){
      return Err(CluEError::InvalidSpinGridWeightTolerance(tolerance));
    }

    let sum = grid.weights.iter().fold(0.0, |sum, w| sum + w);
    if !((sum - 1.0).abs() <= tolerance){
      return Err(CluEError::SpinGridWeightsDoNotSumToOne(sum, tolerance));
    }

    if sum == 1.0{ return Ok(()); }
    for weight in grid.weights.iter_mut(){
      *weight /= sum;
    }

    // Differences up to this size are floating point rounding.
    const ROUNDING: f64 = 1e-12;
    if (sum - 1.0).abs() > ROUNDING{
      println!("The weights of the detected spin grid sum to {}; \
they were normalized to sum to one.", sum);
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function checks that the weights of an orientation grid sum to one
  // within orientations.weight_tolerance, and normalizes them.
  // The signal is the weighted sum of the orientation signals and is not 
  // normalized afterwards, so weights that do not sum to one scale it.
  // A message is printed when normalizing changes the weights by more 
  // than floating point rounding.
  fn check_orientation_grid(&mut self) -> Result<(),CluEError>{
    let Some(OrientationAveraging::Grid(grid)) 
      = &mut self.orientation_grid else{
      return Ok(());
    };

    let tolerance = self.orientation_weight_tolerance
      .unwrap_or(DEFAULT_ORIENTATION_GRID_WEIGHT_TOLERANCE);
    if !(0.0..1.0).contains(&tolerance){
      return Err(CluEError::InvalidOrientationGridWeightTolerance(tolerance));
    }

    let sum = grid.weights.iter().fold(0.0, |sum, w| sum + w);
    if !((sum - 1.0).abs() <= tolerance){
      return Err(CluEError::OrientationGridWeightsDoNotSumToOne(sum, 
            tolerance));
    }

    if sum == 1.0{ return Ok(()); }
    for weight in grid.weights.iter_mut(){
      *weight /= sum;
    }

    // Differences up to this size are floating point rounding.
    const ROUNDING: f64 = 1e-12;
    if (sum - 1.0).abs() > ROUNDING{
      println!("The weights of the orientation grid sum to {}; \
they were normalized to sum to one.", sum);
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function checks the isotope abundances and drop probabilities of 
  // the groups, and normalizes isotope abundances that sum to one within
  // the group's abundance_tolerance.
  // A message is printed when normalizing changes the abundances by more 
  // than floating point rounding.
  fn check_particle_configs(&mut self) -> Result<(),CluEError>{
    // Groups that apply to all unit cells are in both lists, 
    // so messages are printed once per group.
    let mut reported = std::collections::HashSet::<String>::new();

    for particle_config in self.particles.iter_mut()
      .chain(self.extracell_particles.iter_mut())
    {
      let label = particle_config.label.clone();
      let Some(properties) = &mut particle_config.properties else{ continue; };
      let distribution = &mut properties.isotopic_distribution;

      if let Some(p) = distribution.void_probability{
        if !(0.0..=1.0).contains(&p){
          return Err(CluEError::InvalidDropProbability(label, p));
        }
      }

      let abundances = &mut distribution.isotope_abundances;
      if abundances.is_empty(){ continue; }

      let tolerance = distribution.abundance_tolerance
        .unwrap_or(DEFAULT_ABUNDANCE_TOLERANCE);
      if !(0.0..1.0).contains(&tolerance){
        return Err(CluEError::InvalidAbundanceTolerance(label, tolerance));
      }

      for isotope_abundance in abundances.iter(){
        let abundance = isotope_abundance.abundance;
        if !abundance.is_finite() || abundance < 0.0{
          return Err(CluEError::InvalidIsotopeAbundance(label,
                isotope_abundance.isotope.to_string(), abundance));
        }
      }

      let list = |abundances: &[IsotopeAbundance]| abundances.iter()
        .map(|a| format!("{} {}", a.isotope, a.abundance))
        .collect::<Vec::<String>>().join(", ");

      let sum = abundances.iter().map(|a| a.abundance).sum::<f64>();
      if (sum - 1.0).abs() > tolerance{
        return Err(CluEError::IsotopeAbundancesDoNotSumToOne(label,
              list(abundances), sum, tolerance));
      }

      if sum == 1.0{ continue; }
      let original = list(abundances);
      for isotope_abundance in abundances.iter_mut(){
        isotope_abundance.abundance /= sum;
      }

      // Differences up to this size are floating point rounding.
      const ROUNDING: f64 = 1e-12;
      if (sum - 1.0).abs() > ROUNDING && reported.insert(label.clone()){
        println!("Group \"{}\": the isotope abundances ({}) sum to {}; \
they were normalized to ({}).", label, original, sum, list(abundances));
      }
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function gets the simulated experiment's time axis.
  pub fn get_tau_axis_as_ref(&self) -> Result<&Vec::<f64>,CluEError>
  {
    if self.tau_axis.is_empty() {
      return Err(CluEError::NoTimeAxis);
    };
    Ok(&self.tau_axis)
  }
  //----------------------------------------------------------------------------
  /// This function gets the simulated experiment's time axis.
  pub fn get_tau2_axis_as_ref(&self) -> Result<&Vec::<f64>,CluEError>
  {
    if self.tau2_axis.is_empty() {
      return Err(CluEError::NoTimeAxis);
    };
    Ok(&self.tau2_axis)
  }
  //----------------------------------------------------------------------------
  /// This function gets the simulated experiment's time axis.
  pub fn get_tau_axis(&self) -> Result<Vec::<f64>,CluEError>
  {
    if self.tau_axis.is_empty() {
      return Err(CluEError::NoTimeAxis);
    };
    Ok(self.tau_axis.clone())
  }
  //----------------------------------------------------------------------------
  /// This function gets the simulated experiment's time axis.
  pub fn get_tau2_axis(&self) -> Result<Vec::<f64>,CluEError>
  {
    if self.tau2_axis.is_empty() {
      return Err(CluEError::NoTimeAxis);
    };
    Ok(self.tau2_axis.clone())
  }
  //----------------------------------------------------------------------------
  pub fn set_tau_axis(&mut self) -> Result<(),CluEError>
  {
    let (result, n1) = Self::construct_tau_axis(
        &self.number_timepoints, &self.tau_increments, "");

    self.tau_axis = result?;

    let (result, n2) = Self::construct_tau_axis(
        &self.number_timepoints2, &self.tau2_increments, "2");

    self.tau2_axis = result?;

    self.total_number_timepoints = n1*n2;

    Ok(())
  }
  //----------------------------------------------------------------------------
  /// `suffix` is "" for number_timepoints and tau_increments, and "2" for
  /// number_timepoints2 and tau2_increments; it is used in error messages.
  pub fn construct_tau_axis(
    number_timepoints: &[usize], 
    tau_increments: &[f64],
    suffix: &str,
      ) -> (Result<Vec::<f64>,CluEError>,usize)
  {
    if tau_increments.len() != number_timepoints.len(){
      return (Err(CluEError::LengthMismatchTimepointsIncrements(
            suffix.to_string(),
            number_timepoints.len(),tau_increments.len())), 0);
    }
  
    let n_tot = number_timepoints.iter().sum::<usize>();
    
    let mut tau_axis = Vec::<f64>::with_capacity(n_tot);
    
    let mut t = 0.0;

    for (idx, &n_dt) in number_timepoints.iter().enumerate(){
      let dt = tau_increments[idx];
      for _ii in 0..n_dt{
        tau_axis.push(t);
        t += dt;
      }
    }

    (Ok(tau_axis), n_tot)
  }
  //----------------------------------------------------------------------------
  pub fn get_total_number_timesteps(&self) -> usize{
    self.total_number_timepoints
  }
  //----------------------------------------------------------------------------
  /*
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
      return Err(CluEError::LengthMismatchTimepointsIncrements(
            n_dts.len(),dts.len()));
    }

    let Some(pulse_sequence) = &self.pulse_sequence else{
      return Err(CluEError::NoPulseSequence);
    };

    let n_tot = n_dts.iter().sum::<usize>();
    let mut time_axis = Vec::<f64>::with_capacity(n_tot);
    match pulse_sequence{
      PulseSequence::FreeEvolution =>{
        let mut t = 0.0;
        for (idx, &n_dt) in n_dts.iter().enumerate(){
          let dt = dts[idx];
          for _ii in 0..n_dt{
            time_axis.push(t);
            t += dt;
          }
        }
      }
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
      PulseSequence::RefocusedHahnEcho => panic!("TODO"),
      PulseSequence::Custom(_) => panic!("TODO"),
    }

    Ok(time_axis)
  }
  */
  //----------------------------------------------------------------------------
  /// This function tries to build a `Vec::<(usize,usize)>`, where
  /// the first `usize` is the PDB serial number and the second `usize` 
  /// is the partinioning block id.  
  /// Block 0 is reserved for the detected spin.
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
#[derive(Debug,Clone,PartialEq)]
pub enum ClusterPopulations{
  Uniform,
  Random,
  Thermal,
  Zeeman,
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq)]
pub enum DetectionOp{
  Matrix(CxMat),
  S(SpinOp),
  Transition(usize,usize),  
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq)]
pub enum DetectedPopulation{
  Matrix(CxMat),
  S(SpinOp),
  Thermal,  
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
  pub fn from_toml_value(value: toml::Value, unit_of_distance: f64) 
      -> Result<Self,CluEError> 
  {
    match value{
      Value::Array(array) => Self::from_toml_array(array,unit_of_distance),
      Value::String(filename) => { 
        let mut int_grid = IntegrationGrid::read_from_csv(&filename)?;
        int_grid.scale_mut(unit_of_distance);
        Ok(Self::ProbabilityDistribution(int_grid))
      },
      other => Err(CluEError::InvalidTOMLValue(
            "detected_spin.position".to_string(),
            "an array or a csv file name".to_string(),
            other.type_str().to_string())),
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_array(array: Vec::<toml::Value>, unit_of_distance: f64) 
      -> Result<Self,CluEError>
  {

    if array.is_empty(){
      return Err(CluEError::InvalidTOMLValue("detected_spin.position".to_string(),
            "a non-empty array".to_string(), "an empty array".to_string()));
    }
    if !are_all_same_type(&array){
      return Err(CluEError::InvalidTOMLValue("detected_spin.position".to_string(),
            "an array whose elements have the same type".to_string(),
            "an array of mixed types".to_string()));
    }
    match array[0]{
      Value::Array(_) => {
        let n = array.len();
        let mut points = Vec::<f64>::with_capacity(3*n);
        let mut weights = Vec::<f64>::with_capacity(n);
        for value in array.iter(){
          let toml::Value::Array(a) = value else{
            return Err(CluEError::InvalidTOMLValue(
                  "detected_spin.position".to_string(),
                  "a row [x, y, z, weight] of four floats".to_string(),
                  value.to_string()));
          };
          if a.len() != 4{
            return Err(CluEError::InvalidTOMLValue(
                  "detected_spin.position".to_string(),
                  "a row [x, y, z, weight] of four floats".to_string(),
                  value.to_string()));
          }
          if !are_all_same_type(a){
            return Err(CluEError::InvalidTOMLValue(
                  "detected_spin.position".to_string(),
                  "a row [x, y, z, weight] of four floats".to_string(),
                  value.to_string()));
          }
          let (Some(x),Some(y),Some(z),Some(w)) 
              = (a[0].as_float(), a[1].as_float(),a[2].as_float()
              ,a[3].as_float()) else{
            return Err(CluEError::InvalidTOMLValue(
                  "detected_spin.position".to_string(),
                  "a row [x, y, z, weight] of four floats".to_string(),
                  value.to_string()));
          };
          points.push(x*unit_of_distance);
          points.push(y*unit_of_distance);
          points.push(z*unit_of_distance);
          weights.push(w);
        }
        Ok(Self::ProbabilityDistribution(
              IntegrationGrid{dim: 3, points,weights}))
      }
      Value::Float(_) => {
        let vec: Vec::<f64> = array.iter()
          .filter_map(|v| v.as_float()).collect();

        let vec3d = Vector3D::from_vec(vec)?.scale(unit_of_distance);

        Ok(Self::XYZ(vec3d))

      },
      Value::Integer(_) => {
        // Serial numbers are u32 (`Structure::map_serial_id_to_bath_index`).
        let mut vec = Vec::<u32>::with_capacity(array.len());
        for (ii,v) in array.iter().enumerate(){
          let Some(n) = v.as_integer().and_then(|n| u32::try_from(n).ok())
          else{
            return Err(CluEError::InvalidTOMLValue(
                  format!("detected_spin.position[{}]",ii),
                  "a serial number 0 <= n <= 4294967295".to_string(),
                  v.to_string()));
          };
          vec.push(n);
        }

        // Three integers may have been meant as coordinates.
        if let Some(message) = serial_position_message(&vec){
          println!("{}", message);
        }

        Ok(Self::CentroidOverSerials(vec))
      },
      Value::String(_) => {
        let vec: Vec::<String> = array.iter()
          .filter_map(|v| v.as_str()).map(|s| s.to_string()).collect();

        Ok(Self::CentroidOverGroup(vec))
      },
      _ => Err(CluEError::InvalidTOMLValue(
            "detected_spin.position".to_string(),
            "an array of floats [x, y, z], of serial numbers, of group names, \
or of rows [x, y, z, weight]".to_string(),
            format!("an array of {}", array[0].type_str()))),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//------------------------------------------------------------------------------
// This function returns a message for a detected spin position given as 
// three integers, which are read as serial numbers, but may have been meant
// as coordinates.
fn serial_position_message(serials: &[u32]) -> Option<String>{
  let [a, b, c] = serials else{
    return None;
  };
  Some(format!("detected_spin.position = [{a}, {b}, {c}] contains integers, \
so it is read as the centroid of the atoms with serial numbers {a}, {b}, and \
{c}; to give coordinates, write them with decimal points, for example \
[{a}.0, {b}.0, {c}.0]."))
}
//------------------------------------------------------------------------------

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,PartialEq)]
pub enum DetFrame{
  Eigen,
  Zeeman,  
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
      other => Err(CluEError::InvalidTOMLValue(
            "replicate_unit_cell".to_string(),
            "true, false, or [n1, n2, n3]".to_string(),
            other.to_string())),
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_array(array: Vec::<toml::Value>) -> Result<Self,CluEError>
  {

    if array.is_empty(){
      return Err(CluEError::InvalidTOMLValue("replicate_unit_cell".to_string(),
            "a non-empty array".to_string(), "an empty array".to_string()));
    }
    if !are_all_same_type(&array){
      return Err(CluEError::InvalidTOMLValue("replicate_unit_cell".to_string(),
            "an array whose elements have the same type".to_string(),
            "an array of mixed types".to_string()));
    }
    match array[0]{
      Value::Integer(_) => {
        // The number of extra cells per dimension is used as an i32 when
        // the cell offsets are built (`build_cell_shifts()`), so it must be
        // in 0..=i32::MAX.
        let expected = "true, false, or [n1, n2, n3] with integers \
0 <= n <= 2147483647";
        if array.len() != 3{
          return Err(CluEError::InvalidTOMLValue(
                "replicate_unit_cell".to_string(), expected.to_string(),
                format!("an array of length {}", array.len())));
        }
        let mut nc = [0_usize;3];
        for (ii,v) in array.iter().enumerate(){
          let Some(n) = v.as_integer()
            .and_then(|n| i32::try_from(n).ok())
            .and_then(|n| usize::try_from(n).ok())
          else{
            return Err(CluEError::InvalidTOMLValue(
                  format!("replicate_unit_cell[{}]",ii),
                  "an integer 0 <= n <= 2147483647".to_string(),
                  v.to_string()));
          };
          nc[ii] = n;
        }

        Ok(Self::Specified(nc))
      },
      _ => Err(CluEError::InvalidTOMLValue(
            "replicate_unit_cell".to_string(),
            "true, false, or [n1, n2, n3]".to_string(),
            format!("an array of {}", array[0].type_str()))),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ClusterMethod` list different cluster simulation methods.
#[derive(Debug,Clone,PartialEq)]
pub enum ClusterMethod{
  //AnalyticRestricted2CCE,
  CCE,
  GCCE,
  APPA,
  LCE,
}
impl ClusterMethod{
  pub fn from(method_str: &str) -> Result<Self,CluEError>{
    match method_str.to_lowercase().as_str(){
      //"r2cce" => Ok(Self::AnalyticRestricted2CCE),
      "cce" => Ok(Self::CCE),
      "gcce" => Ok(Self::GCCE),
      "appa" => Ok(Self::APPA),
      "lce" => Ok(Self::LCE),
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
      other => Err(CluEError::InvalidTOMLValue("partition_table".to_string(),
            "a csv file name or an array [index0, block0, index1, block1, ...] \
of integers".to_string(), other.type_str().to_string())),
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_array(array: Vec::<toml::Value>) -> Result<Self,CluEError>
  {

    if array.is_empty(){
      return Err(CluEError::InvalidTOMLValue("partition_table".to_string(),
            "a non-empty array".to_string(), "an empty array".to_string()));
    }
    if !are_all_same_type(&array){
      return Err(CluEError::InvalidTOMLValue("partition_table".to_string(),
            "an array whose elements have the same type".to_string(),
            "an array of mixed types".to_string()));
    }
    match array[0]{
      Value::Integer(_) => {
        if array.len() % 2 != 0{
          return Err(CluEError::InvalidTOMLValue("partition_table".to_string(),
                "a csv file name or an array [index0, block0, index1, block1, ...] \
of integers".to_string(),
                format!("an array of length {}", array.len())));
        }
        let mut p_table = Vec::<(usize,usize)>::with_capacity(array.len()/2);
        for ii in 0..array.len()/2{
          let (Some(idx),Some(blk)) = (array[2*ii].as_integer(),
              array[2*ii + 1].as_integer()) else{
            return Err(CluEError::InvalidTOMLValue("partition_table".to_string(),
                  "a csv file name or an array [index0, block0, index1, block1, ...] \
of integers".to_string(),
                  format!("the pair [{}, {}]", array[2*ii], array[2*ii + 1])));
          };
          // Serial numbers and blocks are checked further in
          // `get_partition_pretable()`, which also reads the csv form.
          let (Ok(idx),Ok(blk)) = (usize::try_from(idx),usize::try_from(blk))
          else{
            return Err(CluEError::InvalidTOMLValue(
                  format!("partition_table[{}..={}]", 2*ii, 2*ii + 1),
                  "a pair [index, block] of non-negative integers".to_string(),
                  format!("[{}, {}]", idx, blk)));
          };
          p_table.push((idx,blk));
        }

        Ok(Self::IndexBlockList(p_table))
      },
      _ => Err(CluEError::InvalidTOMLValue("partition_table".to_string(),
            "a csv file name or an array [index0, block0, index1, block1, ...] \
of integers".to_string(),
            format!("an array of {}", array[0].type_str()))),
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
  // TODO: Direction(VectorSpecifier),
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
        return Err(CluEError::MissingTOMLKey(
              "orientations.number".to_string(),
              format!("orientations.grid = \"{}\"", grid)));
      };
      Ok(OrientationAveraging::Lebedev(n_ori))

    } else if grid == KEY_ORI_RANDOM{

      let Some(n_ori) = ori_toml.number else{
        return Err(CluEError::MissingTOMLKey(
              "orientations.number".to_string(),
              format!("orientations.grid = \"{}\"", grid)));
      };
      return Ok(OrientationAveraging::Random(n_ori));

    } else if grid == KEY_ORI_FILE{
      
      let Some(filename) = &ori_toml.file else{
        return Err(CluEError::MissingTOMLKey(
              "orientations.file".to_string(),
              format!("orientations.grid = \"{}\"", grid)));
      };
      let int_grid = IntegrationGrid::read_from_csv(filename)?;
      return Ok(OrientationAveraging::Grid(int_grid));

    } else if grid == KEY_ORI_VECTOR{
      let Some(v) = ori_toml.vector else{
        return Err(CluEError::MissingTOMLKey(
              "orientations.vector".to_string(),
              format!("orientations.grid = \"{}\"", grid)));
      };
      if v.len() != 3{
        return Err(CluEError::InvalidTOMLValue(
              "orientations.vector".to_string(),
              "an array [x, y, z] of three floats".to_string(),
              format!("an array of length {}", v.len())));
      }
      return Ok(OrientationAveraging::Grid(IntegrationGrid{
          dim: 3,
          points: v,
          weights: vec![1.0],
      }));


    } else if grid == KEY_ORI_CIRCLE{
      let Some(v) = ori_toml.vector else{
        return Err(CluEError::MissingTOMLKey(
              "orientations.vector".to_string(),
              format!("orientations.grid = \"{}\"", grid)));
      };
      if v.len() != 3{
        return Err(CluEError::InvalidTOMLValue(
              "orientations.vector".to_string(),
              "an array [x, y, z] of three floats".to_string(),
              format!("an array of length {}", v.len())));
      }
      let Some(n_ori) = ori_toml.number else{
        return Err(CluEError::MissingTOMLKey(
              "orientations.number".to_string(),
              format!("orientations.grid = \"{}\"", grid)));
      };
      return Ok(OrientationAveraging::Grid(
          IntegrationGrid::circle_3d(&Vector3D::from([v[0],v[1],v[2]]), n_ori)
      ));
    } else if grid == KEY_ORI_VECTORGRID{
      let Some(v_grid) = ori_toml.vector_grid else{
        return Err(CluEError::MissingTOMLKey(
              "orientations.vector_grid".to_string(),
              format!("orientations.grid = \"{}\"", grid)));
      };
      let mut points = Vec::<f64>::with_capacity(3*v_grid.len());
      let mut weights = Vec::<f64>::with_capacity(v_grid.len());
      for v in v_grid.iter(){
        if v.len() != 4{
          return Err(CluEError::InvalidTOMLValue(
                "orientations.vector_grid".to_string(),
                "rows [x, y, z, weight] of four floats".to_string(),
                format!("a row of length {}", v.len())));
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

    if config_toml.connect_exchange_groups.is_some(){
      self.connect_exchange_groups = config_toml.connect_exchange_groups;
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

      if let Some(frame) = detected_spin.frame{
        match frame.as_str(){
          KEY_DET_FRAME_EIGEN 
              => self.detected_spin_frame = Some(DetFrame::Eigen),
          KEY_DET_FRAME_ZEEMAN 
              => self.detected_spin_frame = Some(DetFrame::Zeeman),
          _ => return Err(CluEError::InvalidDetectionFrame(frame.to_string())),
        }
      }
      if let Some(toml_value) = detected_spin.g_matrix{
        self.detected_spin_g_matrix 
          = Some(TensorSpecifier::from_toml_value(toml_value,1.0,
                "detected_spin.g_matrix")?);
      }
      
      if let Some(toml_value) = detected_spin.electric_quadrupole{
        self.detected_spin_quadrupole_coupling 
          = Some(TensorSpecifier::from_toml_value(toml_value,unit_of_energy,
                "detected_spin.electric_quadrupole")?);
      }
      
      if let Some(toml_value) = detected_spin.zerofield{
        self.detected_spin_zerofield_coupling 
          = Some(TensorSpecifier::from_toml_value(toml_value,unit_of_energy,
                "detected_spin.zerofield")?);
      }
      
      if let Some(spin_id) = &detected_spin.identity{                          
        self.detected_spin_identity =  Some(Isotope::from(spin_id)?);
      }
      if let Some(spin_mult) = detected_spin.multiplicity{
        self.detected_spin_multiplicity = Some(spin_mult);
      }

      if let Some(toml_value) = detected_spin.position{
        self.detected_spin_position 
          = Some(DetectedSpinCoordinates::from_toml_value(
                toml_value,unit_of_distance)? );
      }

      if detected_spin.weight_tolerance.is_some(){
        self.detected_spin_weight_tolerance = detected_spin.weight_tolerance;
      }

      if detected_spin.transition.is_some(){
        mem::swap(
            &mut self.detected_spin_transition, 
            &mut detected_spin.transition);
      }
      if let Some(det_op) = detected_spin.detection_operator{
        match det_op{
          toml::Value::Array(array) => {

            if let toml::Value::Integer(_) = array[0]{
              if array.len() != 2{
                return Err(CluEError::InvalidTOMLValue(
                      "detected_spin.detection_operator".to_string(),
                      "a transition [i, j] of two integers".to_string(),
                      format!("an array of length {}", array.len())));
              }
              let transition = vec_usize_from_toml_array(array)?;
              self.detection_operator = Some(
                  DetectionOp::Transition(transition[0],transition[1]));
            }
            else{    
              let op = cxmat_from_toml_array(array.clone())?;
              self.detection_operator = Some(DetectionOp::Matrix(op));
            }
          },
          toml::Value::String(s) => {
            let sop = SpinOp::from_str(&s)?;
            self.detection_operator = Some(DetectionOp::S(sop));
          },
          other => return Err(CluEError::InvalidTOMLValue(
                "detected_spin.detection_operator".to_string(),
                "a transition [i, j], a matrix, or a spin operator name"
                .to_string(), other.type_str().to_string())),
        }
      }
      if let Some(rho) = detected_spin.density_matrix{
        match rho{
          toml::Value::Array(array) => {
            let op = cxmat_from_toml_array(array.clone())?;
            self.detected_population = Some(DetectedPopulation::Matrix(op));
          },
          toml::Value::String(s) => {
            if s == KEY_DENSITY_MATRIX_THERMAL{
               self.detected_population = Some(DetectedPopulation::Thermal);
            }else{
              let sop = SpinOp::from_str(&s)?;
              self.detected_population = Some(DetectedPopulation::S(sop));
            }
          },
          other => return Err(CluEError::InvalidTOMLValue(
                "detected_spin.density_matrix".to_string(),
                format!("a matrix, \"{}\", or a spin operator name",
                  KEY_DENSITY_MATRIX_THERMAL), other.type_str().to_string())),
        }
      }
    
    }
    //E--G
    if config_toml.ensemble_cce.is_some(){
      self.ensemble_cce = config_toml.ensemble_cce;
    }
    if let Some(groups) = config_toml.groups{
      for value in groups.iter(){
        let Some(group) = value.as_table() else{
          return Err(CluEError::InvalidTOMLValue("groups".to_string(),
                "an array of tables ([[groups]])".to_string(),
                format!("an array containing a {}", value.type_str())));
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

    if config_toml.mean_fields.is_some(){
      self.mean_fields = config_toml.mean_fields;
    }

    if config_toml.max_cluster_size.is_some(){
      self.max_cluster_size = config_toml.max_cluster_size;
    }

    if config_toml.max_spins.is_some(){
      self.max_spins = config_toml.max_spins;
    }

    if config_toml.max_spherical_tensor_rank.is_some(){
      self.max_spherical_tensor_rank = config_toml.max_spherical_tensor_rank;
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
      if let Some(&cutoff) = pair_cutoffs.get(KEY_CUTOFF_DELTA_ZEEMAN){ 
        self.neighbor_cutoff_delta_zeeman = Some(cutoff)
      }

    }

    if config_toml.number_runs.is_some(){
      self.number_runs = config_toml.number_runs;
    }

    if let Some(number_timepoints) = &mut config_toml.number_timepoints{
      mem::swap(&mut self.number_timepoints, number_timepoints);
    }
    if let Some(number_timepoints2) = &mut config_toml.number_timepoints2{
      mem::swap(&mut self.number_timepoints2, number_timepoints2);
    }

    // O
    if let Some(ori_toml) = config_toml.orientations{
      if ori_toml.weight_tolerance.is_some(){
        self.orientation_weight_tolerance = ori_toml.weight_tolerance;
      }
      self.orientation_grid 
          = Some(OrientationAveraging::from_orientations_toml(ori_toml)?); 
    }

    // P--Q
    if config_toml.partitioning 
        == Some(KEY_PARTITION_PARTICLE.to_string())
    {
      self.partitioning = Some(PartitioningMethod::Particles);

    }
    else if config_toml.partitioning 
        == Some(KEY_PARTITION_EX_GROUPS.to_string())
    {
      self.partitioning 
          = Some(PartitioningMethod::ExchangeGroupsAndParticles);

    }
    else if config_toml.partitioning
        == Some(KEY_PARTITION_KMEANS.to_string())
    { 
      let Some(kmeans_size) = config_toml.kmeans_size else{
        return Err(CluEError::NoKMeansSize);
      };
      if kmeans_size < 1{
        return Err(CluEError::InvalidKMeansSize);
      }    
      self.partitioning 
          = Some(PartitioningMethod::KMeans(kmeans_size));
    
    }
    else if config_toml.partitioning 
        == Some(KEY_PARTITION_RESTRICTED_KMEANS.to_string())
    {
      let Some(kmeans_size) = config_toml.kmeans_size else{
        return Err(CluEError::NoKMeansSize);
      };
      if kmeans_size < 1{
        return Err(CluEError::InvalidKMeansSize);
      }    
      self.partitioning 
          = Some(PartitioningMethod::RestrictedKMeans(kmeans_size));
    }
    else if let Some(s) = config_toml.partitioning
    {
      return Err(CluEError::CannotParsePartitioningMethod(s.to_string()));
    }

    if let Some(kmeans_iter_cap) = config_toml.kmeans_iter_cap{
      if kmeans_iter_cap < 1{
        return Err(CluEError::InvalidKMeansIterCap);
      }
      self.kmeans_iter_cap = Some(kmeans_iter_cap);
    }

    if let Some(p_table) = config_toml.partition_table{
      self.partition_table 
          = Some(PartitionTableConfig::from_toml_value(p_table)?); 
    }

    if config_toml.pdb_model_index.is_some(){
      self.pdb_model_index = config_toml.pdb_model_index;
    }
    if let Some(pulse_seq) = &config_toml.pulse_sequence{
      self.pulse_sequence = Some(PulseSequence::from_toml_value(pulse_seq)?);
    }

    // R--S
    if let Some(radius) = config_toml.radius{
      self.radius = Some(radius*unit_of_distance);
    }
    if let Some(b) = config_toml.run_in_parallel{
      self.run_in_parallel = Some(b);
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
    if config_toml.temperature.is_some(){
      self.temperature = config_toml.temperature;
    }

    if config_toml.populations 
        == Some(KEY_DENSITY_MATRIX_THERMAL.to_string())
    {
      self.cluster_populations = Some(ClusterPopulations::Thermal);
    }else if config_toml.populations 
        == Some(KEY_DENSITY_MATRIX_ID.to_string())
    {
      self.cluster_populations = Some(ClusterPopulations::Uniform);
    }else if config_toml.populations 
        == Some(KEY_DENSITY_MATRIX_RANDOM.to_string())
    {
      self.cluster_populations = Some(ClusterPopulations::Random);
    }else if config_toml.populations 
        == Some(KEY_DENSITY_MATRIX_ZEEMAN.to_string())
    {
      self.cluster_populations = Some(ClusterPopulations::Zeeman);
    }


    if let Some(tau_increments) = &mut config_toml.tau_increments{ 
      for t in tau_increments.iter_mut(){
        *t *= unit_of_time;
      }
      mem::swap(&mut self.tau_increments, tau_increments);
    }
    if let Some(tau2_increments) = &mut config_toml.tau2_increments{ 
      for t in tau2_increments.iter_mut(){
        *t *= unit_of_time;
      }
      mem::swap(&mut self.tau2_increments, tau2_increments);
    }
    
    if config_toml.thermalize_cluster_populations.is_some(){
      self.thermalize_cluster_populations 
          = config_toml.thermalize_cluster_populations;
    } 

    // W--Z
    if let Some(output) = config_toml.output{
      if let Some(&b) = output.get(KEY_OUT_AUX_SIGS){ 
        self.write_auxiliary_signals = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_BATH){ 
        self.write_bath = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_CLU_SIGS){ 
        self.write_cluster_signals = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_DET_SPIN){ 
        self.write_detected_spin = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_CLUSTERS){ 
        self.write_clusters = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_CONFIG){ 
        self.write_config = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_INFO){ 
        self.write_info = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_EXCHANGE_GROUPS){ 
        self.write_exchange_groups = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_METHYL_PARTITIONS){ 
        self.write_methyl_partitions = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_ORI_GRID){ 
        self.write_orientation_grid = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_ORI_SIGS){ 
        self.write_orientation_signals = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_PART_TAB){ 
        self.write_partition_table = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_SANS_SPIN_SIGS){ 
        self.write_sans_spin_signals = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_STRUC_PDB){ 
        self.write_structure_pdb = Some(b);
      }
      if let Some(&b) = output.get(KEY_OUT_TENSORS){ 
        self.write_tensors = Some(b);
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

    Self::from_toml_file(filename)
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
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[cfg(test)]
mod tests{
  use super::*;
  use crate::elements::Element;
  use crate::structure::particle_filter::SecondaryParticleFilter;
  use crate::config::particle_config::ParticleProperties;

  //----------------------------------------------------------------------------
  #[test]
  fn test_config_from_toml(){
  
    let config = Config::from_toml_string(r##"
      clash_distance = 0.1
      clash_distance_pbc = 0.1
      cluster_batch_size = 20000
      connect_exchange_groups = true
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
      #position = [28,29] # centroid over serials
      position = [1.0, 2.0, 3.0]
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

        1H.abundance = 0.5
        2H.abundance = 0.5

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
    //assert_eq!(config.detected_spin_position,
    //    Some(DetectedSpinCoordinates::CentroidOverSerials(vec![28,29])));
    assert_eq!(config.detected_spin_position,
        Some(DetectedSpinCoordinates::XYZ(Vector3D { 
            elements: [1.0*ANGSTROM, 2.0*ANGSTROM, 3.0*ANGSTROM] }))
    );
    assert_eq!(config.input_structure_file,
         Some("../../assets/TEMPO_wat_gly_70A.pdb".to_string()));
    assert_eq!(config.cluster_source,
         Some(ClusterSource::File("clusters_file.toml".to_string())));
    assert_eq!(config.connect_exchange_groups,Some(true));

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
          spherical_tensor_coefficients: Vec::<CxMat>::new(),
          stevens: false,
        })));

    assert_eq!(nitrogen.electric_quadrupole_coupling,
        Some(TensorSpecifier::Eig(EigSpecifier{
          values: Some([-0.28e6, -1.47e6, 1.75e6]),
          x_axis: Some(x.clone()),    
          y_axis: Some(y.clone()),    
          z_axis: None,    
          spherical_tensor_coefficients: Vec::<CxMat>::new(),
          stevens: false,
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
        Some(crate::config::particle_config::Cosubstitute::SameMolecule));


  }
  //----------------------------------------------------------------------------
  // This function reads a config with one group "g" with the settings 
  // `settings`.
  fn config_with_group(settings: &str) -> Result<Config,CluEError>{
    Config::from_toml_string(&format!(r##"
        [[groups]]
        name = "g"
        selection.elements = ["C"]
        {}
        "##, settings))
  }
  // This function returns the isotope abundances of the group "g".
  fn abundances_of_group(config: &Config) -> Vec::<(Isotope,f64)>{
    config.particles[0].properties.as_ref().unwrap().isotopic_distribution
      .isotope_abundances.iter().map(|a| (a.isotope, a.abundance)).collect()
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_isotope_abundances_sum_to_one(){
    // Abundances within the tolerance of one are normalized.
    let config = config_with_group(
        "1H.abundance = 0.1\n2H.abundance = 0.2\n3H.abundance = 0.7")
      .unwrap();
    let abundances = abundances_of_group(&config);
    let sum = 0.1 + 0.2 + 0.7;
    for ((_, a), expected) in abundances.iter().zip([0.1, 0.2, 0.7]){
      assert!((a - expected/sum).abs() < 1e-15);
    }

    // Rounded table values need a larger tolerance.
    let settings = "12C.abundance = 0.9893\n13C.abundance = 0.0106";
    assert!(matches!(config_with_group(settings),
          Err(CluEError::IsotopeAbundancesDoNotSumToOne(group,_,_,tol))
          if group == "g" && tol == 1e-9));
    let config = config_with_group(&format!(
          "{}\nabundance_tolerance = 1e-4", settings)).unwrap();
    let abundances = abundances_of_group(&config);
    assert_eq!(abundances[0].0, Isotope::Carbon12);
    assert!((abundances[0].1 - 0.9893/0.9999).abs() < 1e-12);
    assert!((abundances[1].1 - 0.0106/0.9999).abs() < 1e-12);
    assert!((abundances[0].1 + abundances[1].1 - 1.0).abs() < 1e-15);

    // Abundances that do not sum to one are errors.
    for settings in ["13C.abundance = 0.3", 
        "13C.abundance = 0.8\n12C.abundance = 0.8"]
    {
      assert!(matches!(config_with_group(settings),
            Err(CluEError::IsotopeAbundancesDoNotSumToOne(group,_,_,_))
            if group == "g"), "{}", settings);
    }

    // Configs built in code are checked by set_defaults().
    let mut config = Config::new();
    let mut particle_config = ParticleConfig::new("g".to_string());
    let mut properties = ParticleProperties::new();
    properties.isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen2, abundance: 0.5});
    particle_config.properties = Some(properties);
    config.particles.push(particle_config);
    assert!(matches!(config.set_defaults(),
          Err(CluEError::IsotopeAbundancesDoNotSumToOne(..))));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_invalid_isotope_abundances(){
    assert_eq!(config_with_group(
          "12C.abundance = 1.2\n13C.abundance = -0.2").err(),
        Some(CluEError::InvalidIsotopeAbundance("g".to_string(),
            "13C".to_string(), -0.2)));
    assert!(matches!(config_with_group(
          "12C.abundance = nan\n13C.abundance = 0.5"),
          Err(CluEError::InvalidIsotopeAbundance(..))));

    for tolerance in ["-1.0", "1.0", "1.5", "nan"]{
      assert!(matches!(config_with_group(&format!(
            "12C.abundance = 1.0\nabundance_tolerance = {}", tolerance)),
            Err(CluEError::InvalidAbundanceTolerance(..))), "{}", tolerance);
    }

    // A tolerance without abundances is ignored, and the group keeps the 
    // isotopes of the structure file.
    let config = config_with_group("abundance_tolerance = 1.5").unwrap();
    assert!(abundances_of_group(&config).is_empty());
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_detected_spin_grid_weights(){
    let config_with_grid = |weights: [&str;2], tolerance: &str|{
      Config::from_toml_string(&format!(r##"
          [detected_spin]
          position = [[0.0, 0.0, 0.0, {}], [1.0, 0.0, 0.0, {}]]
          {}
          "##, weights[0], weights[1], tolerance))
    };
    let weights_of = |config: &Config| -> Vec::<f64>{
      let Some(DetectedSpinCoordinates::ProbabilityDistribution(grid))
        = &config.detected_spin_position else{ panic!("no grid"); };
      grid.weights.clone()
    };

    // Weights that sum to one, including negative weights, are accepted.
    let config = config_with_grid(["0.5","0.5"],"").unwrap();
    assert_eq!(weights_of(&config), vec![0.5, 0.5]);
    let config = config_with_grid(["1.2","-0.2"],"").unwrap();
    assert_eq!(weights_of(&config), vec![1.2, -0.2]);

    // Weights that do not sum to one are errors, 
    // unless they are within the tolerance, which normalizes them.
    assert_eq!(config_with_grid(["0.45","0.45"],"").err(),
        Some(CluEError::SpinGridWeightsDoNotSumToOne(0.9, 1e-9)));
    let config = config_with_grid(["0.45","0.45"],"weight_tolerance = 0.2")
      .unwrap();
    for w in weights_of(&config){
      assert!((w - 0.5).abs() < 1e-15);
    }
    assert_eq!(config.detected_spin_weight_tolerance, Some(0.2));

    for tolerance in ["-1.0", "1.0", "nan"]{
      assert!(matches!(config_with_grid(["0.5","0.5"],
            &format!("weight_tolerance = {}",tolerance)),
            Err(CluEError::InvalidSpinGridWeightTolerance(_))), "{}",
          tolerance);
    }

    // Positions that are not grids ignore the tolerance.
    let config = Config::from_toml_string(r##"
        [detected_spin]
        position = [1.0, 2.0, 3.0]
        weight_tolerance = 1.5
        "##).unwrap();
    assert!(matches!(config.detected_spin_position,
          Some(DetectedSpinCoordinates::XYZ(_))));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_orientation_grid_weights(){
    let config_with_grid = |weights: [&str;2], tolerance: &str|{
      Config::from_toml_string(&format!(r##"
          [orientations]
          grid = "vector_grid"
          vector_grid = [[0.0, 0.0, 1.0, {}], [1.0, 0.0, 0.0, {}]]
          {}
          "##, weights[0], weights[1], tolerance))
    };
    let weights_of = |config: &Config| -> Vec::<f64>{
      let Some(OrientationAveraging::Grid(grid))
        = &config.orientation_grid else{ panic!("no grid"); };
      grid.weights.clone()
    };

    // Weights that sum to one are accepted unchanged.
    let config = config_with_grid(["0.25","0.75"],"").unwrap();
    assert_eq!(weights_of(&config), vec![0.25, 0.75]);

    // Weights that do not sum to one are errors, 
    // unless they are within the tolerance, which normalizes them.
    assert_eq!(config_with_grid(["0.45","0.45"],"").err(),
        Some(CluEError::OrientationGridWeightsDoNotSumToOne(0.9, 1e-9)));
    let config = config_with_grid(["0.45","0.45"],"weight_tolerance = 0.2")
      .unwrap();
    for w in weights_of(&config){
      assert!((w - 0.5).abs() < 1e-15);
    }
    assert_eq!(config.orientation_weight_tolerance, Some(0.2));

    for tolerance in ["-1.0", "1.0", "nan"]{
      assert!(matches!(config_with_grid(["0.5","0.5"],
            &format!("weight_tolerance = {}",tolerance)),
            Err(CluEError::InvalidOrientationGridWeightTolerance(_))), "{}",
          tolerance);
    }

    // Built-in grids pass the check.
    for grid in [r#"grid = "vector"
        vector = [0.0, 0.0, 1.0]"#, r#"grid = "circle"
        vector = [0.0, 0.0, 1.0]
        number = 7"#]
    {
      let config = Config::from_toml_string(&format!(
            "[orientations]\n{}", grid)).unwrap();
      let sum: f64 = weights_of(&config).iter().sum();
      assert!((sum - 1.0).abs() < 1e-15, "{}", grid);
    }

    // Grids that are made later ignore the tolerance.
    let config = Config::from_toml_string(r##"
        [orientations]
        grid = "lebedev"
        number = 6
        weight_tolerance = 1.5
        "##).unwrap();
    assert_eq!(config.orientation_grid, Some(OrientationAveraging::Lebedev(6)));

    // A grid file without a header keeps all of its points, so its weights
    // sum to one.
    let filename = std::env::temp_dir().join(format!(
          "clue_oxide_test_orientations_{}.csv", std::process::id()));
    std::fs::write(&filename, "0,0,1,0.25\n1,0,0,0.75\n").unwrap();
    let config = Config::from_toml_string(&format!(r##"
        [orientations]
        grid = "file"
        file = "{}"
        "##, filename.to_str().unwrap()));
    std::fs::remove_file(&filename).unwrap();
    assert_eq!(weights_of(&config.unwrap()), vec![0.25, 0.75]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_detected_spin_position_unit_of_distance(){
    // Coordinates are in unit_of_distance, as grid points are.
    let config = Config::from_toml_string(r##"
        unit_of_distance = "nm"
        detected_spin.position = [1.0, 2.0, 3.0]
        "##).unwrap();
    assert_eq!(config.detected_spin_position,
        Some(DetectedSpinCoordinates::XYZ(
            Vector3D::from([10.0, 20.0, 30.0]).scale(ANGSTROM))));

    let config = Config::from_toml_string(r##"
        unit_of_distance = "nm"
        detected_spin.position = [[1.0, 2.0, 3.0, 1.0]]
        "##).unwrap();
    let Some(DetectedSpinCoordinates::ProbabilityDistribution(grid))
      = config.detected_spin_position else{ panic!("no grid"); };
    assert_eq!(grid.xyz(0).unwrap(),
        Vector3D::from([10.0, 20.0, 30.0]).scale(ANGSTROM));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_serial_position_message(){
    // Three integers are read as serial numbers, with a message.
    let message = serial_position_message(&[25, 25, 25]).unwrap();
    assert!(message.contains("serial numbers 25, 25, and 25"));
    assert!(message.contains("[25.0, 25.0, 25.0]"));
    let config = Config::from_toml_string(
        "detected_spin.position = [25, 25, 25]").unwrap();
    assert_eq!(config.detected_spin_position,
        Some(DetectedSpinCoordinates::CentroidOverSerials(vec![25,25,25])));

    // Other numbers of serials get no message.
    assert_eq!(serial_position_message(&[28, 29]), None);
    assert_eq!(serial_position_message(&[1, 2, 3, 4]), None);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_drop_probability_range(){
    for p in ["1.5", "-0.5", "nan"]{
      assert!(matches!(config_with_group(&format!("drop_probability = {}",p)),
            Err(CluEError::InvalidDropProbability(group,_)) if group == "g"),
          "{}", p);
    }
    for p in ["0.0", "1.0", "0.5"]{
      assert!(config_with_group(&format!("drop_probability = {}",p)).is_ok());
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_kmeans_iter_cap(){
    // Default.
    let config = Config::from_toml_string("").unwrap();
    assert_eq!(config.kmeans_iter_cap, Some(1000));

    // Set from TOML.
    let config = Config::from_toml_string(r##"
        partitioning = "restricted_kmeans"
        kmeans_size = 3
        kmeans_iter_cap = 25
        "##).unwrap();
    assert_eq!(config.kmeans_iter_cap, Some(25));
    assert_eq!(config.partitioning, 
        Some(PartitioningMethod::RestrictedKMeans(3)));

    // Zero is rejected.
    assert_eq!(Config::from_toml_string("kmeans_iter_cap = 0").err(),
        Some(CluEError::InvalidKMeansIterCap));
    assert_eq!(CluEError::InvalidKMeansIterCap.to_string(),
        "kmeans_iter_cap must be at least 1");
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_cosubstitute(){
    use crate::config::particle_config::Cosubstitute;
    for (value, expected) in [("group", Cosubstitute::Group),
        ("filter", Cosubstitute::Group), ("particle", Cosubstitute::Particle),
        ("same_molecule", Cosubstitute::SameMolecule)]
    {
      assert_eq!(Cosubstitute::from(value,"methyl"), Ok(expected.clone()));
      let config = Config::from_toml_string(&format!(r##"
          [[groups]]
          name = "methyl"
          selection.elements = ["H"]
          cosubstitute = "{}"
          "##, value)).unwrap();
      assert_eq!(config.particles[0].properties.as_ref().unwrap().cosubstitute,
          Some(expected));
    }

    // "bonded" is not supported for cosubstitution, 
    // and the error names the group.
    for value in ["bonded", "nonsense"]{
      let result = Config::from_toml_string(&format!(r##"
          [[groups]]
          name = "methyl"
          selection.elements = ["H"]
          cosubstitute = "{}"
          "##, value));
      assert_eq!(result.err(), Some(CluEError::CannotParseCosubstitute(
              "methyl".to_string(), value.to_string())));
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_isotope_selections_are_rejected(){
    // Groups are selected before isotopes are drawn, so selections by
    // isotope were removed; the keys are invalid.
    for (key, value) in [("isotopes", "2H"), ("not_isotopes", "1H")]{
      let result = Config::from_toml_string(&format!(r##"
          [[groups]]
          name = "deuterons"
          selection.{} = ["{}"]
          "##, key, value));
      assert_eq!(result.err(), Some(CluEError::InvalidTOMLKeys(vec![
              format!("groups[0].selection.{} (group \"deuterons\"): \
unknown key", key)])));
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_set_tau_axis(){
    let mut config = Config::new();
    config.pulse_sequence = Some(PulseSequence::CarrPurcell(1));
    config.number_timepoints = vec![100,91];
    config.tau_increments = vec![1.0e-8, 1e-7];
    config.set_tau_axis().unwrap();
    let tau_axis = config.get_tau_axis().unwrap();
    assert_eq!(tau_axis.len(),191);
    assert_eq!(tau_axis[0],0.0);
    let t = 1e-8;
    assert!((tau_axis[1]-t).abs()/t<1e-12);
    let t = 1e-6;
    assert!((tau_axis[100]-t).abs()/t<1e-12);
    let t = 1.1e-6;
    assert!((tau_axis[101]-t).abs()/t<1e-12);
    let t = 10e-6;
    assert!((tau_axis[190]-t).abs()/t<1e-12);
    
  }
  //----------------------------------------------------------------------------
  // Errors for invalid config values name the key (and group), and parse 
  // errors list the valid options.
  #[test]
  fn test_config_errors_name_the_key(){
    let s = |x: &str| x.to_string();
    let cases = vec![
      ("replicate_unit_cell = \"auto\"",
       CluEError::InvalidTOMLValue(s("replicate_unit_cell"),
         s("true, false, or [n1, n2, n3]"), s("\"auto\""))),
      ("orientations.grid = \"lebedev\"",
       CluEError::MissingTOMLKey(s("orientations.number"),
         s("orientations.grid = \"lebedev\""))),
      ("orientations.grid = \"vector\"\norientations.vector = [1.0, 0.0]",
       CluEError::InvalidTOMLValue(s("orientations.vector"),
         s("an array [x, y, z] of three floats"), s("an array of length 2"))),
      ("orientations.grid = \"spiral\"",
       CluEError::CannotParseOrientations(s("spiral"))),
      ("detected_spin.frame = \"lab\"",
       CluEError::InvalidDetectionFrame(s("lab"))),
      ("detected_spin.detection_operator = true",
       CluEError::InvalidTOMLValue(s("detected_spin.detection_operator"),
         s("a transition [i, j], a matrix, or a spin operator name"),
         s("boolean"))),
      ("partitioning = \"pairs\"",
       CluEError::CannotParsePartitioningMethod(s("pairs"))),
      ("partition_table = [1, 2, 3]",
       CluEError::InvalidTOMLValue(s("partition_table"),
         s("a csv file name or an array [index0, block0, index1, block1, ...] \
of integers"), s("an array of length 3"))),
      ("pulse_sequence = [[\"pulse\", \"pi/2\"], [\"delay2\", \"tau\"]]",
       CluEError::InvalidPulseStep(1,
         Box::new(CluEError::WrongTauDimension(2,1)))),
      ("[[groups]]\nname = \"g\"\ndrop_probability = \"high\"",
       CluEError::InvalidTOMLValue(s("group \"g\": drop_probability"),
         s("a float"), s("string"))),
      ("[[groups]]\nname = \"g\"\n1H.active = 1",
       CluEError::InvalidTOMLValue(s("group \"g\": 1H.active"),
         s("a bool"), s("integer"))),
      ("[[groups]]\nname = \"g\"\n1H.hyperfine.axes.x = [1.0, 0.0, 0.0]",
       CluEError::MissingTOMLKey(s("group \"g\": 1H.hyperfine.values"),
         s("group \"g\": 1H.hyperfine given as a table"))),
    ];
    for (toml_str, expected) in cases{
      let result = Config::from_toml_string(toml_str);
      assert_eq!(result.err(), Some(expected), "{}", toml_str);
    }

    // Messages.
    let err = Config::from_toml_string("[[groups]]\nname = \"g\"\n\
        drop_probability = \"high\"").unwrap_err();
    assert_eq!(err.to_string(),
        "group \"g\": drop_probability must be a float, but got string");
    let err = Config::from_toml_string("pulse_sequence = [[\"delay2\", \"tau\"]]")
      .unwrap_err();
    assert_eq!(err.to_string(), "pulse_sequence[0]: \"delay2\" cannot be used \
with \"tau\"; use \"delay\" with \"tau\" or \"-tau\", and \"delay2\" with \
\"tau2\" or \"-tau2\"");
    let err = Config::from_toml_string("cluster_method = \"dmrg\"").unwrap_err();
    assert_eq!(err.to_string(), "cannot parse \"dmrg\" as a cluster method; \
the options are \"cce\", \"gcce\", \"appa\", and \"lce\"");

    // The second time axis names its own keys.
    let (result, _) = Config::construct_tau_axis(&[10, 10], &[1e-8], "2");
    assert_eq!(result, Err(CluEError::LengthMismatchTimepointsIncrements(
            s("2"), 2, 1)));
    assert_eq!(result.unwrap_err().to_string(), "number_timepoints2 has 2 \
entries, but tau2_increments has 1; they must have the same length");
  }
  //----------------------------------------------------------------------------
  // Integer arrays in replicate_unit_cell, detected_spin.position, and
  // partition_table are range checked instead of being cast with `as`.
  #[test]
  fn test_config_integer_arrays_are_checked(){
    let s = |x: &str| x.to_string();
    let rep = "true, false, or [n1, n2, n3] with integers 0 <= n <= 2147483647";
    let rep_el = "an integer 0 <= n <= 2147483647";
    let cases = vec![
      // Wrong lengths (index out of bounds, or entries ignored).
      ("replicate_unit_cell = [2, 2]",
       CluEError::InvalidTOMLValue(s("replicate_unit_cell"), s(rep),
         s("an array of length 2"))),
      ("replicate_unit_cell = [1, 1, 1, 7]",
       CluEError::InvalidTOMLValue(s("replicate_unit_cell"), s(rep),
         s("an array of length 4"))),
      // Negative and too large (wrapped around).
      ("replicate_unit_cell = [1, -1, 1]",
       CluEError::InvalidTOMLValue(s("replicate_unit_cell[1]"), s(rep_el),
         s("-1"))),
      ("replicate_unit_cell = [0, 0, 2147483648]",
       CluEError::InvalidTOMLValue(s("replicate_unit_cell[2]"), s(rep_el),
         s("2147483648"))),
      ("detected_spin.position = [12, -5]",
       CluEError::InvalidTOMLValue(s("detected_spin.position[1]"),
         s("a serial number 0 <= n <= 4294967295"), s("-5"))),
      ("detected_spin.position = [4294967296]",
       CluEError::InvalidTOMLValue(s("detected_spin.position[0]"),
         s("a serial number 0 <= n <= 4294967295"), s("4294967296"))),
      ("partition_table = [10, 1, -1, 5]",
       CluEError::InvalidTOMLValue(s("partition_table[2..=3]"),
         s("a pair [index, block] of non-negative integers"), s("[-1, 5]"))),
      ("partition_table = [10, -1]",
       CluEError::InvalidTOMLValue(s("partition_table[0..=1]"),
         s("a pair [index, block] of non-negative integers"), s("[10, -1]"))),
    ];
    for (toml_str, expected) in cases{
      let result = Config::from_toml_string(toml_str);
      assert_eq!(result.err(), Some(expected), "{}", toml_str);
    }

    // Valid input, including the bounds, is read as before.
    let config = Config::from_toml_string(
        "replicate_unit_cell = [0, 3, 2147483647]").unwrap();
    assert_eq!(config.replicate_unit_cell,
        Some(ReplicateUnitCell::Specified([0, 3, 2147483647])));
    let config = Config::from_toml_string(
        "detected_spin.position = [0, 28, 4294967295]").unwrap();
    assert_eq!(config.detected_spin_position,
        Some(DetectedSpinCoordinates::CentroidOverSerials(
            vec![0, 28, 4294967295])));
    let config = Config::from_toml_string(
        "partition_table = [10, 1, 0, 0, 4294967296, 2]").unwrap();
    assert_eq!(config.partition_table,
        Some(PartitionTableConfig::IndexBlockList(
            vec![(10, 1), (0, 0), (4294967296, 2)])));

    let err = Config::from_toml_string("replicate_unit_cell = [2, 2]")
      .unwrap_err();
    assert_eq!(err.to_string(), "replicate_unit_cell must be true, false, \
or [n1, n2, n3] with integers 0 <= n <= 2147483647, but got an array of \
length 2");
  }

}
