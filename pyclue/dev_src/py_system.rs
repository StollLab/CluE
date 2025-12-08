use pyo3::prelude::*;

use num_complex::Complex;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

use clue_oxide::{
  config::{
    Config, 
      DensityMatrixMethod,
      OrientationAveraging,
      command_line_input::CommandLineInput,
    },
  info,
  space_3d::Vector3D,
};
use clue_oxide::io::FromTOMLString;
use crate::py_clue_errors::PyCluEError;

#[pyclass(name = "System")]
#[derive(Debug)]
pub struct PySystem{
  pub config: Config,
  pub rng: ChaCha20Rng,
}

#[pymethods]
impl PySystem{
 
  fn db_print(&self){
    println!("{:?}",self);
  }
  //----------------------------------------------------------------------------
  fn run(&self) -> Result<(Vec::<f64>, Vec::<Complex::<f64>>),PyCluEError>{
    Ok(clue_oxide::run(self.config.clone())?)
  }
  //----------------------------------------------------------------------------
  #[staticmethod]
  fn from_input(args: String) -> Result<Self,PyCluEError>{

    let mut input = Vec::<String>::new();
    input.push("clue_oxixe".to_string());

    for arg in args.split(" "){
      input.push(arg.to_string());
    }

    let input = CommandLineInput::new(input)?;

    // Decide what information to display.
    if input.show_help{
      info::help::print_help();
    }

    if input.show_license{
      info::license::print_license();
    }

    if input.show_warrenty{
      info::warrenty::print_warrenty();
    }

    if input.show_version{
      info::version::print_version();
    }

    if input.show_title{
      info::title::print_title();
    }

    let mut config = Config::read_input(input)?;

    config.set_defaults()?;

    config.construct_time_axis()?;

    let rng = match config.rng_seed{
      Some(seed) => ChaCha20Rng::seed_from_u64(seed),
      None => ChaCha20Rng::from_entropy(),
    };

    Ok(Self{
      config,
      rng,
    })
  }
  //----------------------------------------------------------------------------
  pub fn get_apply_pbc(&self) -> Option<bool>{
    self.config.apply_pbc.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_apply_pbc(&mut self, value:  bool){
    self.config.apply_pbc = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_clash_distance(&self) -> Option<f64>{
    self.config.clash_distance.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_clash_distance(&mut self, value: f64){
    self.config.clash_distance = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_clash_distance_pbc(&self) -> Option<f64>{
    self.config.clash_distance_pbc.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_clash_distance_pbc(&mut self, value: f64 ){
    self.config.clash_distance_pbc = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_cluster_batch_size(&self) -> Option<usize>{
    self.config.cluster_batch_size.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_cluster_batch_size(&mut self, value: usize){
    self.config.cluster_batch_size = Some(value);
  }
  //----------------------------------------------------------------------------
  // TODO: get/set cluster_method
  //----------------------------------------------------------------------------
  pub fn get_clusters_file(&self) -> Option<String>{
    self.config.clusters_file.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_clusters_file(&mut self, value: String ){
    self.config.clusters_file = Some(value);
  }
  //----------------------------------------------------------------------------
  // TODO: get/set density_matrix
  //----------------------------------------------------------------------------
  // TODO: get/set detected_spin_g_matrix
  //----------------------------------------------------------------------------
  // TODO: get/set detected_spin_identity
  //----------------------------------------------------------------------------
  pub fn get_detected_spin_multiplicity(&self) -> Option<usize>{
    self.config.detected_spin_multiplicity.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_detected_spin_multiplicity(&mut self, value:  usize){
    self.config.detected_spin_multiplicity = Some(value);
  }
  //----------------------------------------------------------------------------
  // TODO: get/set detected_spin_position
  //----------------------------------------------------------------------------
  pub fn get_detected_spin_transition(&self) -> Option<[usize;2]>{
    self.config.detected_spin_transition.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_detected_spin_transition(&mut self, value: [usize;2] ){
    self.config.detected_spin_transition = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_input_structure_file(&self) -> Option<String>{
    self.config.input_structure_file.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_input_structure_file(&mut self, value:  Option<String>){
    self.config.input_structure_file = value;
  }
  //----------------------------------------------------------------------------
  // TODO: get/set load_geometry
  //----------------------------------------------------------------------------
  pub fn get_magnetic_field(&self) -> Option<f64>{
    match &self.config.magnetic_field{
      Some(vector) => Some(vector.z()),
      None => None,
    }
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_magnetic_field(&mut self, value: f64){
    let value = Vector3D::from([ 0.0, 0.0, value]);
    self.config.magnetic_field = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_max_cell_size(&self) -> Option<usize>{
    self.config.max_cell_size.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_max_cell_size(&mut self, value:  usize){
    self.config.max_cell_size = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_max_cluster_size(&self) -> Option<usize>{
    self.config.max_cluster_size.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_max_cluster_size(&mut self, value:  usize){
    self.config.max_cluster_size = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_max_spin_order(&self) -> Option<usize>{
    self.config.max_spin_order.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_max_spin_order(&mut self, value: usize){
    self.config.max_spin_order = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_min_cell_size(&self) -> Option<usize>{
    self.config.min_cell_size.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_min_cell_size(&mut self, value: usize){
    self.config.min_cell_size = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_neighbor_cutoff_coupling(&self) -> Option<f64>{
    self.config.neighbor_cutoff_coupling.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_neighbor_cutoff_coupling(&mut self, value:  f64){
    self.config.neighbor_cutoff_coupling = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_neighbor_cutoff_delta_hyperfine(&self) -> Option<f64>{
    self.config.neighbor_cutoff_delta_hyperfine.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_neighbor_cutoff_delta_hyperfine(&mut self, value: f64){
    self.config.neighbor_cutoff_delta_hyperfine = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_neighbor_cutoff_dipole_perpendicular(&self) -> Option<f64>{
    self.config.neighbor_cutoff_dipole_perpendicular.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_neighbor_cutoff_dipole_perpendicular(&mut self, value: f64){
    self.config.neighbor_cutoff_dipole_perpendicular = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_neighbor_cutoff_distance(&self) -> Option<f64>{
    self.config.neighbor_cutoff_distance.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_neighbor_cutoff_distance(&mut self, value: f64){
    self.config.neighbor_cutoff_distance = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_neighbor_cutoff_3_spin_hahn_mod_depth(&self) -> Option<f64>{
    self.config.neighbor_cutoff_3_spin_hahn_mod_depth.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_neighbor_cutoff_3_spin_hahn_mod_depth(&mut self, value: f64){
    self.config.neighbor_cutoff_3_spin_hahn_mod_depth = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_neighbor_cutoff_3_spin_hahn_taylor_4(&self) -> Option<f64>{
    self.config.neighbor_cutoff_3_spin_hahn_taylor_4.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_neighbor_cutoff_3_spin_hahn_taylor_4(&mut self, value: f64){
    self.config.neighbor_cutoff_3_spin_hahn_taylor_4 = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_number_timepoints(&self) -> Vec::<usize>{
    self.config.number_timepoints.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_number_timepoints(&mut self, value:  Vec::<usize>){
    self.config.number_timepoints = value;
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (distribution,n_points))]
  pub fn set_orientations(&mut self, 
      distribution: String, n_points: Option<usize>)
    -> Result<(),PyCluEError>
  {

    let toml_str = match n_points{
      Some(n) => format!("{} = {}",distribution, n), 
      None => format!("{}",distribution),
    };

    let ori = OrientationAveraging::from_toml_string(&toml_str)?;
    self.config.orientation_grid = Some(ori);
    Ok(())
  }
  //----------------------------------------------------------------------------
  // TODO get/set partitioning_method
  //----------------------------------------------------------------------------
  // TODO: get/set particles, and consider better name
  //----------------------------------------------------------------------------
  pub fn get_pdb_model_index(&self) -> Option<usize>{
    self.config.pdb_model_index.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_pdb_model_index(&mut self, value: usize ){
    self.config.pdb_model_index = Some(value);
  }
  //----------------------------------------------------------------------------
  // TODO: get/set pulse_sequence.
  //----------------------------------------------------------------------------
  pub fn get_root_dir(&self) -> Option<String>{
    self.config.root_dir.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_root_dir(&mut self, value: String){
    self.config.root_dir = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_radius(&self) -> Option<f64>{
    self.config.radius.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_radius(&mut self, value: f64){
    self.config.radius = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_rng_seed(&self) -> Option<u64>{
    self.config.rng_seed.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_(&mut self, value: u64){
    self.config.rng_seed = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_save_name(&self) -> Option<String>{
    self.config.save_name.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_save_name(&mut self, value:  Option<String>){
    self.config.save_name = value;
  }
  //----------------------------------------------------------------------------
  pub fn get_system_name(&self) -> Option<String>{
    self.config.system_name.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_system_name(&mut self, value:  Option<String>){
    self.config.system_name = value;
  }
  //----------------------------------------------------------------------------
  pub fn get_temperature(&self) -> Option<f64>{
    match self.config.density_matrix{
      Some(DensityMatrixMethod::Thermal(value)) => Some(value),
      _ => None,  
    }
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_temperature(&mut self, value:  f64){
    self.config.density_matrix = Some(DensityMatrixMethod::Thermal(value));
  }
  //----------------------------------------------------------------------------
  pub fn get_time_increments(&self) -> Vec::<f64>{
    self.config.time_increments.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_time_increments(&mut self, value: Vec::<f64>){
    self.config.time_increments = value;
  }
  //----------------------------------------------------------------------------
  pub fn get_write_auxiliary_signals(&self) -> Option<bool>{
    self.config.write_auxiliary_signals.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_write_auxiliary_signals(&mut self, value: bool){
    self.config.write_auxiliary_signals = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_write_bath(&self) -> Option<bool>{
    self.config.write_bath.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_write_bath(&mut self, value: bool){
    self.config.write_bath = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_write_clusters(&self) -> Option<bool>{
    self.config.write_clusters.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_write_clusters(&mut self, value: bool){
    self.config.write_clusters = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_write_info(&self) -> Option<bool>{
    self.config.write_info.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_write_info(&mut self, value: bool){
    self.config.write_info = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_write_exchange_groups(&self) -> Option<bool>{
    self.config.write_exchange_groups.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_write_exchange_groups(&mut self, value: bool){
    self.config.write_exchange_groups = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_write_methyl_partitions(&self) -> Option<bool>{
    self.config.write_methyl_partitions.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_write_methyl_partitions(&mut self, value: bool){
    self.config.write_methyl_partitions = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_write_orientation_signals(&self) -> Option<bool>{
    self.config.write_orientation_signals.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_write_orientation_signals(&mut self, value: bool){
    self.config.write_orientation_signals = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_write_sans_spin_signals(&self) -> Option<bool>{
    self.config.write_sans_spin_signals.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_write_sans_spin_signals(&mut self, value: bool){
    self.config.write_sans_spin_signals = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_write_structure_pdb(&self) -> Option<bool>{
    self.config.write_structure_pdb.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_write_structure_pdb(&mut self, value: bool){
    self.config.write_structure_pdb = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_write_tensors(&self) -> Option<bool>{
    self.config.write_tensors.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_write_tensors(&mut self, value: bool){
    self.config.write_tensors = Some(value);
  }
  //----------------------------------------------------------------------------
}

