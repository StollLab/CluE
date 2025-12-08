//pub mod cif;
pub mod gro;
pub mod pdb;
pub mod extended_structure;
pub mod primary_structure;
pub mod exchange_groups;
pub mod particle;
pub mod particle_filter;

//use crate::config::particle_config;
use crate::config::{
  Config, 
  particle_config::{ParticleConfig,ParticleProperties,IsotopeProperties,
    TensorSpecifier},
};
use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::integration_grid::IntegrationGrid;
use crate::math;
use crate::space_3d::Vector3D;
use crate::structure::{
  exchange_groups::ExchangeGroupManager,
  particle::Particle,
  particle_filter::*
};
use crate::physical_constants::ANGSTROM;
use crate::isotopes::Isotope;

use crate::space_3d::SymmetricTensor3D;

use rand_chacha::ChaCha20Rng;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::collections::HashMap;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `DetectedSpin` holds the data about the detected spin.
/// `gamma_matrix` holds the effective gamma matrix.
/// `isotope` identifies the species of the spin.
/// `weighted_coordinates` holds the probability amplitude grid.
/// 'spin_multiplicity' determines the spin of the particle.
/// 'transition' selects which transition is being measured. 
#[derive(Debug,Clone)]
pub struct DetectedSpin{
  pub gamma_matrix: SymmetricTensor3D,
  pub isotope: Isotope,
  pub weighted_coordinates: IntegrationGrid,
  pub spin_multiplicity: usize,
  pub transition: [usize;2],
}
//------------------------------------------------------------------------------
impl DetectedSpin{
  /// This function returns a reference to the detected spin's gamma matrix.
  pub fn gyromagnetic_ratio_matrix(&self) 
    -> &SymmetricTensor3D
  {

    &self.gamma_matrix
  }
  //----------------------------------------------------------------------------
  /// This function returns the detected spin's multiplicity.
  pub fn spin_multiplicity(&self) -> usize{
    self.spin_multiplicity
  }
  //----------------------------------------------------------------------------
  pub fn get_average_spin_position(&self) -> Result<Vector3D,CluEError>{

    let mut r_ave = Vector3D::zeros();
    for ii in 0..self.weighted_coordinates.len(){
      let Ok(mut r) = self.weighted_coordinates.xyz(ii) else{
        return Err(CluEError::NoCentralSpinCoor);
      };

      let w = self.weighted_coordinates.weight(ii);
      
      r = r.scale(w);

      r_ave = &r_ave + &r;
    }

    Ok(r_ave)
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `Structure` contain info on the positions and identities of particles
/// within the system.
/// `detected_particl` is the particle giving rise to the detected signal.
/// `bath_particles` contains all particles that are not directly detected.
/// `bath_spins_indices` is a list of which bath_particles have spin.
/// `cell_indices` is a set of lists of all periodic copies of each particle 
/// from the primary cell 
/// `cell_offsets` is the set of displacements for periodic boundary conditions.
/// `connections` lists chemical bonds.
/// `molecules` is the set of molecules in the system.
/// `molecule_ids` is a list specifing which molecule each particle belongs to.
/// `cosubstitution_groups` is a collections of atoms that should be isotope 
/// substituted together
/// `exchange_groups` contains information on methyls.
/// `particle_config_ids` lists indices for which config.particle corresponds 
/// to each particle
/// `extracell_particle_config_ids` is a list of indices for which 
/// `config.extracell_particle` corresponds to each particle.
/// `pdb_origin` is the origin of the pdb frame.
/// `primary_cell_indices` lists indices indicating the particle that each 
/// particle is a periodic boundary condition copy of  
/// `nth_active_to_reference_index` lists indices for display and output files.
///
#[derive(Debug,Clone)]
pub struct Structure{
  pub detected_particle: Option<DetectedSpin>,
  pub bath_particles: Vec::<Particle>,
  //bath_spins_indices: Vec::<usize>,
  pub cell_indices: Vec::<Vec::<Option<usize>>>,
  pub cell_offsets: Vec::<Vector3D>,
  pub connections: AdjacencyList,
  molecules: Vec::<Vec::<usize>>,
  molecule_ids: Vec::<usize>,
  cosubstitution_groups: Vec::< Vec::<usize> >,
  pub exchange_groups: Option<ExchangeGroupManager>,
  particle_config_ids: Vec::<Option<usize>>,
  extracell_particle_config_ids: Vec::<Option<usize>>,
  pub pdb_origin: Vector3D,
  primary_cell_indices: Vec::<usize>,
  nth_active_to_reference_index: Vec::<usize>,
  pub active_indices: Vec::<usize>,
  pub bath_indices_to_active_indices: Vec::<Option<usize>>,
  pub bath_to_serial_indices: HashMap::<(u32,usize),usize>,
}

impl Structure{
  /// This function builds a new `Structure` from a list of `Particle`s
  /// their connectivity and cell info.
  pub fn new(
      bath_particles: Vec::<Particle>,
      connections: AdjacencyList,
      cell_offsets: Vec::<Vector3D>
      ) -> Self
  {
    Structure{
      detected_particle: None,
      bath_particles,
      //bath_spins_indices: Vec::<usize>::new(),
      cell_indices: Vec::<Vec::<Option<usize>>>::new(),
      cell_offsets,
      connections,
      extracell_particle_config_ids: Vec::<Option<usize>>::new(),
      molecules: Vec::<Vec::<usize>>::new(),
      molecule_ids: Vec::<usize>::new(),
      cosubstitution_groups: Vec::<Vec::<usize>>::new(),
      exchange_groups: None,
      particle_config_ids: Vec::<Option<usize>>::new(),
      pdb_origin: Vector3D::zeros(),
      primary_cell_indices: Vec::<usize>::new(),
      nth_active_to_reference_index: Vec::<usize>::new(),
      active_indices: Vec::<usize>::new(),
      bath_indices_to_active_indices: Vec::<Option<usize>>::new(),
      bath_to_serial_indices: HashMap::<(u32,usize),usize>::new(),
    }

  }
  //----------------------------------------------------------------------------
  /// This function builds a `Structure` from info in `Config`.
  pub fn build_structure(rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<Self,CluEError>
  {

    let Some(filename) = &config.input_structure_file else{
      return Err(CluEError::NoStructureFile);
    };

    let Some(model_idx) = config.pdb_model_index else{
      return Err(CluEError::NoModelIndex);
    };

    let mut structure = pdb::parse_pdb(filename,model_idx)?;

    structure.build_primary_structure(rng,config)?;

    structure.build_extended_structure(rng, config)?;

    structure.map_nth_active_to_reference_indices();

    structure.set_serials()?;

    Ok(structure)
  }
  //----------------------------------------------------------------------------
  /// This function returns the number of bath particles in the system.
  pub fn number(&self) -> usize{
    self.bath_particles.len()
  }
  //----------------------------------------------------------------------------
  /// This function returns the number of active bath particles in the system.
  pub fn number_active(&self) -> usize{
    let mut n_active = 0;
    for particle in self.bath_particles.iter(){
      if particle.active{ n_active += 1;}
    }
    n_active
  }
  //----------------------------------------------------------------------------
  /// This function finds the set of particles that pass `&ParticleFilter`
  /// and returns a list of references.
  pub fn find<'a>(&'a self, particle_filter: &ParticleFilter)
    -> Vec::<&'a Particle>
  {
    let indices = particle_filter.filter(self);
    let mut out = Vec::<&Particle>::with_capacity(indices.len());

    for idx in indices{
      out.push(&self.bath_particles[idx]);
    }

    out
  }
  //----------------------------------------------------------------------------
  /// This function identifies which pbc copy of the primary cell
  /// `particle_particle` with index `idx` resides in.
  pub fn cell_id(&self, idx: usize) -> Result<usize,CluEError>
  {
    if idx >= self.primary_cell_indices.len(){
      return Err(CluEError::CannotFindCellID(idx));
    }

    let idx0 = self.primary_cell_indices[idx];
    let cell_indices = &self.cell_indices[idx0];
    let Some(cell_id) = cell_indices.iter().position(|&x| {
        match x{ Some(id) => id == idx, None => false,} }) else{
      return Err(CluEError::CannotFindCellID(idx));
    };

    Ok(cell_id)
  }
  //----------------------------------------------------------------------------
  /// This function returns the ID number of the molecule that particle `idx`
  /// is in.
  pub fn molecule_id(&self,idx: usize) -> usize{

    let idx0 = self.primary_cell_indices[idx];
   
    self.molecule_ids[idx0]

  }
  //----------------------------------------------------------------------------
  /// This function retrieves the reference index of the bath particle with
  /// index `bath_index`.
  /// The reference index of the the detected spin in zero, and the first
  /// bath particle has index 1.
  pub fn get_reference_index_from_bath_index(&self, bath_index: usize) 
    -> Result<usize,CluEError>
  {
    if bath_index >= self.bath_particles.len(){
      return Err(CluEError::CannotFindRefIndexFromBathIndex(bath_index));
    }
    Ok(bath_index + 1)
  }
  //----------------------------------------------------------------------------
  /// This function retrieves the bath index of the bath particle with
  /// reference index `ref_index`.
  /// The reference index of the the detected spin in zero, and the first
  /// bath particle has index 1.
  pub fn get_bath_index_from_reference_index(&self, ref_index: usize) 
    -> Result<usize,CluEError>
  {
    let bath_index = ref_index - 1;
    if bath_index >= self.bath_particles.len() {
      return Err(CluEError::CannotFindBathIndexFromRefIndex(bath_index));
    }
    Ok(bath_index)
  }
  //----------------------------------------------------------------------------
  /// This function retrieves the bath index of `n`th active particle.
  pub fn get_bath_index_of_nth_active(&self, n: usize) 
    -> Result<usize,CluEError>
  {
    if n >= self.nth_active_to_reference_index.len(){
      return Err(CluEError::CannotFindRefIndexFromNthActive(n));
    }
    Ok(self.nth_active_to_reference_index[n] - 1)
  }
  //----------------------------------------------------------------------------
  pub fn get_positions(&self,indices: &[usize]) -> Vec::<Vector3D>{
    let mut positions = Vec::<Vector3D>::with_capacity(indices.len());

    for &index in indices.iter(){
      positions.push(self.bath_particles[index].coordinates.clone());
    }

    positions
  }
  //----------------------------------------------------------------------------
  /// This function retrieves the reference index of `n`th active particle.
  /// The reference index of the the detected spin in zero, and the first
  /// bath particle has index 1.
  /// The detected spin is always active and corresponds to `n=0`.
  pub fn get_reference_index_of_nth_active(&self, n: usize)
    -> Result<usize,CluEError>
  {
    if n >= self.nth_active_to_reference_index.len(){
      return Err(CluEError::CannotFindRefIndexFromNthActive(n));
    }
    Ok(self.nth_active_to_reference_index[n])

  }
  //----------------------------------------------------------------------------
  pub fn set_serials(&mut self) -> Result<(),CluEError>{
  
    self.bath_to_serial_indices = HashMap::<(u32,usize),usize>::with_capacity(
        self.bath_particles.len());

    for (bath_idx, particle) in self.bath_particles.iter().enumerate(){
      let cell_id = self.cell_id(bath_idx)?;
      if let Some(serial) = particle.serial{
        self.bath_to_serial_indices.insert((serial,cell_id),bath_idx);
      }
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function builds `nth_active_to_reference_index`, which maps
  /// the set of active particles to the set of all system particles.
  pub fn map_nth_active_to_reference_indices(&mut self){

    let number_active = self.number_active();
    self.active_indices = Vec::<usize>::with_capacity(number_active);
    self.bath_indices_to_active_indices = (0..self.bath_particles.len())
        .map(|_| None).collect::<Vec::<Option<usize>>>();

    let mut act_to_ref = Vec::<usize>::with_capacity(number_active + 1);

    act_to_ref.push(0);

    for (bath_index,particle) in self.bath_particles.iter().enumerate(){
      if particle.active{ 
        self.bath_indices_to_active_indices[bath_index] 
            = Some(act_to_ref.len());
        self.active_indices.push(bath_index);
        act_to_ref.push(bath_index + 1);
      }
    }

    self.nth_active_to_reference_index = act_to_ref;
  }
  //----------------------------------------------------------------------------
  /// This function takes a reference index and returns either `Ok(n)`,
  /// where when counting active particles in bath_particles, particle
  /// reference_index is the nth active particle.
  pub fn get_nth_active_from_reference_index(&self, reference_index: usize) 
    -> Result<usize,CluEError>
  { 
    // Check that the particle is a bath particle.
    if reference_index == 0{
      return Err(CluEError::DetectedSpinDoesNotHaveAnActiveIndex);
    }
    let bath_index = reference_index - 1;

    if bath_index >= self.bath_particles.len(){
      return Err(CluEError::CannotFindParticleForRefIndex(reference_index));
    }

    if !self.bath_particles[bath_index].active{
      return Err(CluEError::ParticleIsNotActive(reference_index));
    }

    // Count active particle upto and including the bath_index particle.
    let mut nth_active = 0;

    for (ii,particle) in self.bath_particles.iter().enumerate(){
      if particle.active{
        nth_active += 1;
      }
      if ii == bath_index{
        break;
      }
    }

    let ref_idx = self.get_reference_index_of_nth_active(nth_active)?;
    if ref_idx != reference_index {
      return Err(CluEError::CannotFindParticleForRefIndex(reference_index));
    }

    Ok(nth_active)
  }
  //----------------------------------------------------------------------------
  /// This function takes a list of PDB serial numbers an returns the 
  /// centroid position.
  pub fn centroid_over_serials(&self, serials: Vec::<u32>) 
    -> Result<Vector3D,CluEError>
  {
    let mut filter = ParticleFilter::new();
    filter.serials = serials;
    filter.indices = math::unique(self.primary_cell_indices.clone());

    self.centroid_over_groups(&[&filter])
  }
  //----------------------------------------------------------------------------
  /// This function finds the centroid over all particles that are in the
  /// group defined by 
  pub fn centroid_over_groups(&self, filters: &[&ParticleFilter])
    -> Result<Vector3D,CluEError>
  {
    let mut r_ave = Vector3D::zeros();
    let mut count = 0.0;

    for filter in filters.iter(){
      let indices = filter.filter(self);

      for &idx in indices.iter(){
        let r = &self.bath_particles[idx].coordinates;
        r_ave = &r_ave + r;
        count += 1.0;
      }
    }
    r_ave = r_ave.scale(1.0/count);

    Ok(r_ave)
  }
  //----------------------------------------------------------------------------
  /// This function returns the exchange coupling for particle `particle_index`.
  pub fn extract_exchange_coupling(&self, particle_index: usize,
      config: &Config)
    -> f64
  {
    let Some(isotope_properties) = self.extract_isotope_properties(
        particle_index,config) 
    else{
      return 0.0;
    };

    isotope_properties.exchange_coupling.unwrap_or(0.0)
  }
  //----------------------------------------------------------------------------
  /// This function returns the electric quadrupole coupling for particle 
  /// `particle_index`.
  pub fn extract_electric_quadrupole_specifier<'a>(&self, particle_index: usize,
      config: &'a Config)
    -> Option<&'a TensorSpecifier>
  {
    let isotope_properties = self.extract_isotope_properties(
        particle_index,config)?;

    isotope_properties.electric_quadrupole_coupling.as_ref()
  }
  //----------------------------------------------------------------------------
  /// This function returns the g-matrix for particle `particle_index`.
  pub fn extract_g_matrix_specifier<'a>(&self, particle_index: usize,
      config: &'a Config)
    -> Option<&'a TensorSpecifier>
  {
    let isotope_properties = self.extract_isotope_properties(
        particle_index,config)?;

    isotope_properties.g_matrix.as_ref()
  }
  //----------------------------------------------------------------------------
  /// This function returns the info to determine the hyperfine coupling for 
  /// particle `particle_index`.
  pub fn extract_hyperfine_specifier<'a>(&self, particle_index: usize,
      config: &'a Config)
    -> Option<&'a TensorSpecifier>
  {

    let isotope_properties = self.extract_isotope_properties(
        particle_index,config)?;

    isotope_properties.hyperfine_coupling.as_ref()
  }
  //----------------------------------------------------------------------------
  /// This function returns the `&IsotopeProperties` that appies to particle 
  /// `particle_index`.
  pub fn extract_isotope_properties<'a>(&self, particle_index: usize,
      config: &'a Config)
    -> Option<&'a IsotopeProperties>
  {
    let particle_index_0 = self.primary_cell_indices[particle_index];

    let p_cfg_id = self.particle_config_ids[particle_index_0]?;

    let properties = &config.particles[p_cfg_id].properties.as_ref()?;

    let key = self.bath_particles[particle_index].isotope.to_string();
    
    properties.isotope_properties.get(&key)

  }
  //----------------------------------------------------------------------------
  // This function searches the user specified particle configuration (if any),
  // and determines which one each particle goes with.
  // This function returns a error if multiple filter cover by the same 
  // particle. 
  fn pair_particle_configs(&self, particle_configs: &[ParticleConfig]) 
    -> Result< Vec::<Option<usize>>,CluEError>
  {
  
    // Initialize pairings.
    let n_particles = self.bath_particles.len();

    let mut particle_config_ids = (0..n_particles).map(|_| None).
      collect::<Vec::<Option<usize>>>();

    // Pair particles with configs.
    for (id,particle_config) in particle_configs.iter().enumerate(){
      if let Some(filter) = &particle_config.filter{
        let indices = filter.filter(self);

        for idx in indices{
          if let Some(prev_id) = particle_config_ids[idx]{
            return Err(CluEError::FiltersOverlap(
              particle_configs[prev_id].label.to_string(),
              particle_config.label.to_string()
            ));
          }
          particle_config_ids[idx] = Some(id);
        }
      
      } 
    }
    
    Ok(particle_config_ids)
  }
  //----------------------------------------------------------------------------
  /// This function saves the bath particles in csv format.
  pub fn bath_to_csv(&self, filename: &str, config: &Config) 
    -> Result<(),CluEError>
  {
    let Ok(file) = File::create(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()) );
    };

    let n_char_f64 = 16;
    let n_chars_per_line = 4*n_char_f64;
    let n_spins = self.bath_particles.len();

    let bytes_per_char = 32;

    let n_bytes = bytes_per_char*(n_spins+1)*n_chars_per_line;

    let mut stream = BufWriter::with_capacity(n_bytes,file);
    
    let line = "index,particle,x,y,z,active,cell_id,group\n";
    if stream.write(line.as_bytes()).is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }    

    for (ii, particle) in self.bath_particles.iter().enumerate(){

      let r = &particle.coordinates - &self.pdb_origin; 

      let idx = self.get_reference_index_from_bath_index(ii)?;

      let cell_id = self.cell_id(ii)?;

      let particle_index_0 = self.primary_cell_indices[ii];

      let group_label = if cell_id == 0{
        match self.particle_config_ids[particle_index_0]{
          Some(cfg_idx) => &config.particles[cfg_idx].label,
          None => "none",
        }
      }else{
        match self.extracell_particle_config_ids[particle_index_0]{
          Some(cfg_idx) => &config.extracell_particles[cfg_idx].label,
          None => "none",
        }
      };

      let line = format!("{},{},{},{},{},{},{},{}\n",
          idx, particle.isotope,
          r.x()/ANGSTROM, r.y()/ANGSTROM, r.z()/ANGSTROM,
          particle.active,cell_id,group_label);
    
      if stream.write(line.as_bytes()).is_err(){
        return Err(CluEError::CannotWriteFile(filename.to_string()) );
      }    
    }

    Ok(())
  }
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

  //use crate::config::lexer::get_tokens_from_line;
  use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};
  use crate::io::FromTOMLString;
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_structure(){
    let mut config = Config::from_toml_string(r##"
        input_structure_file = "assets/TEMPO.pdb"
        radius = 20
        detected_spin.position = [28,29]
        number_timepoints = [101]
        tau_increments = [1e-1]
        cluster_method = "cce"
        max_cluster_size = 2
        magnetic_field = 1.2
        replicate_unit_cell = true

        [[groups]]
        name = "tempo_h"
        selection.elements = ["H"]

        1H.c3_tunnel_splitting = 80e-3
        "##).unwrap();

    config.set_defaults().unwrap();

    let mut rng = ChaCha20Rng::from_entropy();

    let structure = Structure::build_structure(&mut rng,&config).unwrap();

    assert_eq!(structure.number_active(),19);
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
