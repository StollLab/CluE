use pyo3::prelude::*;

use numpy::PyArray;

use ndarray::{Array1,Ix1};

use num_complex::Complex;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};


use clue_oxide::{
  config::{
    Config, 
    particle_config::ParticleProperties,
    DensityMatrixMethod,
    OrientationAveraging,
    command_line_input::CommandLineInput,
    config_toml::{
      DEFAULT_UNIT_ENERGY,
      DEFAULT_UNIT_DISTANCE,
      DEFAULT_UNIT_MAGNETIC_FIELD,
      DEFAULT_UNIT_TIME,
    }
    },
  info,
  physical_constants::*,
  space_3d::Vector3D,
};
use clue_oxide::io::FromTOMLString;
use crate::py_clue_errors::PyCluEError;

#[pyclass(name = "Config")]
#[derive(Debug)]
pub struct PyConfig{
  pub config: Config,
  pub rng: ChaCha20Rng,
  unit_of_distance: f64,
  unit_of_energy: f64,
  unit_of_magnetic_field: f64,
  unit_of_time: f64,
}

#[pymethods]
impl PyConfig{
 
  /*
  //----------------------------------------------------------------------------
  fn from_pydict(pydict: &PyDict) -> Result<Self,PyCluEError>
  {
    let toml = toml::to_string(pydict).unwrap();
    Self::from_input(toml)
  }
  */
  //----------------------------------------------------------------------------
  pub fn db_print(&self){
    println!("{:?}",self);
  }
  //----------------------------------------------------------------------------
  pub fn run(&self)// -> Result<(Vec::<f64>, Vec::<Complex::<f64>>),PyCluEError>{
      -> Result<
          (Py<PyArray<f64,Ix1>>,Py<PyArray<Complex::<f64>,Ix1>>),
          PyCluEError
         >
  {
    let (time_axis,signal) = clue_oxide::run(self.config.clone())?;

    let time_axis = Array1::from_vec(time_axis);

    let signal = Array1::from_vec(signal);

    let time_axis = Python::with_gil(|py|{
        PyArray::from_owned_array(py, time_axis).unbind()
    });
    let signal = Python::with_gil(|py|{
        PyArray::from_owned_array(py, signal).unbind()
    });
    Ok( (time_axis, signal) )
  }
  //----------------------------------------------------------------------------
  #[staticmethod]
  pub fn from_input(input: String) -> Result<Self,PyCluEError>
  {
  
    let config = if std::path::Path::new(&input).exists(){
      Config::from_toml_file(&input)
    }else{
      Config::from_toml_string(&input)
    }?;

    let rng = match config.rng_seed{
      Some(seed) => ChaCha20Rng::seed_from_u64(seed),
      None => ChaCha20Rng::from_entropy(),
    };

    Ok(Self{
      config,
      rng,
      unit_of_distance: distance_unit_to_meters(DEFAULT_UNIT_DISTANCE)?,
      unit_of_energy: energy_unit_to_hertz(DEFAULT_UNIT_ENERGY)?,
      unit_of_magnetic_field: magnetic_field_unit_to_tesla(DEFAULT_UNIT_MAGNETIC_FIELD)?,
      unit_of_time: time_unit_to_seconds(DEFAULT_UNIT_TIME)?,
    })
  }
  //----------------------------------------------------------------------------
  #[staticmethod]
  fn from_command_line_input(args: String) -> Result<Self,PyCluEError>{

    let mut input = Vec::<String>::new();
    input.push("clue_oxide".to_string());

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
      unit_of_distance: distance_unit_to_meters(DEFAULT_UNIT_DISTANCE)?,
      unit_of_energy: energy_unit_to_hertz(DEFAULT_UNIT_ENERGY)?,
      unit_of_magnetic_field: magnetic_field_unit_to_tesla(DEFAULT_UNIT_MAGNETIC_FIELD)?,
      unit_of_time: time_unit_to_seconds(DEFAULT_UNIT_TIME)?,
    })
  }
  //----------------------------------------------------------------------------
  pub fn get_do_replicate_unit_cell(&self) -> Option<bool>{
    self.config.do_replicate_unit_cell.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_do_replicate_unit_cell(&mut self, value:  bool){
    self.config.do_replicate_unit_cell = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_clash_distance(&self) -> Option<f64>{
    match self.config.clash_distance{
      Some(x) => Some(x/self.unit_of_distance),
      None => None
    }
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_clash_distance(&mut self, value: f64){
    self.config.clash_distance = Some(value*self.unit_of_distance);
  }
  //----------------------------------------------------------------------------
  pub fn get_clash_distance_pbc(&self) -> Option<f64>{
    match self.config.clash_distance_pbc{
      Some(x) => Some(x/self.unit_of_distance),
      None => None
    }
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_clash_distance_pbc(&mut self, value: f64 ){
    self.config.clash_distance_pbc = Some(value*self.unit_of_distance);
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
    let value = Vector3D::from([ 0.0, 0.0, value*self.unit_of_magnetic_field]);
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
  pub fn get_max_spins(&self) -> Option<usize>{
    self.config.max_spins.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_max_spins(&mut self, value: usize){
    self.config.max_spins = Some(value);
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
  pub fn get_neighbor_cutoff_coupling_xx_yy(&self) -> Option<f64>{
    match self.config.neighbor_cutoff_coupling_xx_yy{
      Some(x) => Some(x/self.unit_of_energy),
      None => None
    }
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_neighbor_cutoff_coupling_xx_yy(&mut self, value:  f64){
    self.config.neighbor_cutoff_coupling_xx_yy = Some(value*self.unit_of_energy);
  }
  //----------------------------------------------------------------------------
  pub fn get_neighbor_cutoff_delta_hyperfine_zz(&self) -> Option<f64>{
    match self.config.neighbor_cutoff_delta_hyperfine_zz{
      Some(x) => Some(x/self.unit_of_energy),
      None => None
    }
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_neighbor_cutoff_delta_hyperfine_zz(&mut self, value: f64){
    self.config.neighbor_cutoff_delta_hyperfine_zz = Some(value*self.unit_of_energy);
  }
  //----------------------------------------------------------------------------
  pub fn get_neighbor_cutoff_point_dipole_perpendicular(&self) -> Option<f64>{
    match self.config.neighbor_cutoff_point_dipole_perpendicular{
      Some(x) => Some(x/self.unit_of_energy),
      None => None
    }
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_neighbor_cutoff_point_dipole_perpendicular(&mut self, value: f64){
    self.config.neighbor_cutoff_point_dipole_perpendicular = Some(value*self.unit_of_energy);
  }
  //----------------------------------------------------------------------------
  pub fn get_neighbor_cutoff_distance(&self) -> Option<f64>{
    match self.config.neighbor_cutoff_distance{
      Some(x) => Some(x/self.unit_of_distance),
      None => None
    }
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_neighbor_cutoff_distance(&mut self, value: f64){
    self.config.neighbor_cutoff_distance = Some(value*self.unit_of_distance);
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
    let u = self.unit_of_energy*HZ_TO_RAD_PER_S;
    let u4 = u*u*u*u;
    match self.config.neighbor_cutoff_3_spin_hahn_taylor_4{
      Some(x) => Some(x/u4),
      None => None
    }
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_neighbor_cutoff_3_spin_hahn_taylor_4(&mut self, value: f64){
    let u = self.unit_of_energy*HZ_TO_RAD_PER_S;
    let u4 = u*u*u*u;
    self.config.neighbor_cutoff_3_spin_hahn_taylor_4 = Some(value*u4);
  }
  //----------------------------------------------------------------------------
  pub fn get_number_runs(&self) -> Option<usize>{
    self.config.number_runs.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_number_runs(&mut self, value:  usize){
    self.config.number_runs = Some(value);
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
  // TODO get/set partitioning
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
  //----------------------------------------------------------------------------
  pub fn get_radius(&self) -> Option<f64>{
    match self.config.radius{
      Some(x) => Some(x/self.unit_of_distance),
      None => None
    }
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_radius(&mut self, value: f64){
    self.config.radius = Some(value*self.unit_of_distance);
  }
  //----------------------------------------------------------------------------
  pub fn get_rng_seed(&self) -> Option<u64>{
    self.config.rng_seed.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_rng_seed(&mut self, value: u64){
    self.config.rng_seed = Some(value);
  }
  //----------------------------------------------------------------------------
  pub fn get_output_directory(&self) -> Option<String>{
    self.config.output_directory.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_output_directory(&mut self, value:  Option<String>){
    self.config.output_directory = value;
  }
  //----------------------------------------------------------------------------
  pub fn get_run_name(&self) -> Option<String>{
    self.config.run_name.clone()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_run_name(&mut self, value:  Option<String>){
    self.config.run_name = value;
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
  pub fn get_tau_increments(&self) -> Vec::<f64>{
    self.config.tau_increments.iter().map(|t| t/self.unit_of_time)
        .collect::<Vec::<f64>>()
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (value))]
  pub fn set_tau_increments(&mut self, value: Vec::<f64>){
    self.config.tau_increments = value.iter().map(|t| t*self.unit_of_time)
        .collect::<Vec::<f64>>();
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (group,isotope))]
  pub fn get_tunnel_splitting(&self,group: String,isotope: String) 
      -> Option<f64>
  {
    for particle_config in self.config.particles.iter(){
      if particle_config.label != group{
        continue;
      }

      let Some(properties) = &particle_config.properties else {
        return None;
      };
      
      let Some(spin_prop) = properties.isotope_properties.get(&isotope) else {
        return None;
      };

      if let Some(j) = spin_prop.exchange_coupling{
        return Some(j/C3_TUNNEL_SPLITTING_TO_EXCHANGE_COUPLING);
      }
      return None;
    }
    None
  }
  //----------------------------------------------------------------------------
  #[pyo3(signature = (group,isotope,value))]
  pub fn set_tunnel_splitting(&mut self,
      group: String, isotope: String, value: f64) 
      -> Result<(),PyCluEError>
  {
    let j =  C3_TUNNEL_SPLITTING_TO_EXCHANGE_COUPLING*value;

    for particle_config in self.config.particles.iter_mut(){
      if particle_config.label != group{
        continue;
      }

      if particle_config.properties.is_none(){
        particle_config.properties = Some(ParticleProperties::new());
      }

      let properties = particle_config.properties.as_mut().unwrap();
      
      let spin_prop = properties.isotope_properties
          .get_mut(&isotope).unwrap(); 

      spin_prop.exchange_coupling = Some(j)
    }
    Ok(())
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

