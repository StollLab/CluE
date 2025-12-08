use crate::config::{Config,DetectedSpinCoordinates};
use crate::CluEError;
use crate::cluster::connected_subgraphs::separate_into_connected_subgraphs;
use crate::physical_constants::{HBAR,MUB,MUN};
use crate::elements::Element;
use crate::isotopes::Isotope;
use crate::quantum::tensors::construct_symmetric_tensor_from_tensor_specifier;
use crate::structure::ParticleFilter;
use crate::structure::{Structure, exchange_groups::*}; 
use crate::structure::ParticleProperties;
use crate::structure::SecondaryParticleFilter;
use crate::structure::IntegrationGrid;
use crate::structure::DetectedSpin;
use crate::space_3d;
use crate::space_3d::Vector3D;

use rand_chacha::ChaCha20Rng;

impl Structure{
  /// This method uses an input `Config` to set the structure's
  /// spins and exchange groups.  The number of bath particle is unchanged.
  pub fn build_primary_structure(&mut self, 
      rng: &mut ChaCha20Rng, config: &Config)
    -> Result<(),CluEError>
    {

    (self.molecules, self.molecule_ids)  =
      separate_into_connected_subgraphs(&self.connections);

    // TODO: one PBC on each side should be used to ensure reconect_bonds()
    // does not take a spin out of range.
    self.reconnect_bonds()?;
    
    self.set_spins(config)?;

    // Set methyl and primary amonium groups..
    self.set_exchange_groups();

    self.find_cosubstitution_groups(config)?;

    // Set the detected spin as the origin.
    self.set_detected_spin(rng,config)?;

    Ok(())
  }  
  //----------------------------------------------------------------------------
  // This method places the detected spin into the system.
  fn set_detected_spin(&mut self, rng: &mut ChaCha20Rng, config: &Config) 
      -> Result<(),CluEError>
  {

    let Some(detected_spin_position) = &config.detected_spin_position else{
      return Err(CluEError::NoCentralSpinCoor);
    };

    let mut weighted_coordinates = IntegrationGrid::new(3);
    match detected_spin_position{
      DetectedSpinCoordinates::CentroidOverGroup(labels) =>{
        
        let mut filters = Vec::<&ParticleFilter>::with_capacity(labels.len());
        
        for label in labels.iter(){

          let Some((_id, p_cfg)) = config.find_particle_config(label) else{
            return Err(CluEError::MissingFilter(label.to_string()))
          };

          let Some(filter) = &p_cfg.filter else{
            return Err(CluEError::MissingFilter(label.to_string()));
          };

          filters.push(filter);
        }
        
        let r_ave = self.centroid_over_groups(&filters)?;
        weighted_coordinates.push(vec![r_ave.x(),r_ave.y(),r_ave.z()],1.0)?;
      },
      DetectedSpinCoordinates::CentroidOverSerials(serials) => {
        let r_ave = self.centroid_over_serials(serials.clone())?;
        weighted_coordinates.push(vec![r_ave.x(),r_ave.y(),r_ave.z()],1.0)?;
      },
      DetectedSpinCoordinates::XYZ(r) => {
        weighted_coordinates.push(vec![r.x(),r.y(),r.z()],1.0)?;
      },
      DetectedSpinCoordinates::ProbabilityDistribution(grid) => 
        weighted_coordinates = grid.clone(),
    
    }

    let origin = weighted_coordinates.mean(); 
    self.pdb_origin = Vector3D::from([origin[0],origin[1],origin[2]])
      .scale(-1.0);

    weighted_coordinates.translate_mut(
        &[self.pdb_origin.x(),self.pdb_origin.y(),self.pdb_origin.z() ]); 

    for particle in self.bath_particles.iter_mut(){
      particle.coordinates = &particle.coordinates + &self.pdb_origin;
    }

    if let Some(exchange_group_manager) = &mut self.exchange_groups{
      for exchange_group in exchange_group_manager.exchange_groups
        .iter_mut(){
          exchange_group.translate(&self.pdb_origin);
      }
    }

    // TODO: move to config?
    let Some(isotope) = config.detected_spin_identity else{
      return Err(CluEError::NoCentralSpinIdentity);
    };
    
    let Some(transition) = config.detected_spin_transition else{
      return Err(CluEError::NoCentralSpinTransition);
    };

    let Some(spin_multiplicity) = config.detected_spin_multiplicity else{
      return Err(CluEError::NoDetectedSpinMultiplicity);
    };

    let Some(g_matrix_specifier) = &config.detected_spin_g_matrix else{
      return Err(CluEError::NoGMatrixSpecifier);
    };
    
    let g_matrix = construct_symmetric_tensor_from_tensor_specifier(rng,
        g_matrix_specifier, None, self,config)?;
    
    let mu: f64 = if isotope == Isotope::Electron{
      -MUB
    } else {
      MUN
    };
    let gamma_matrix = g_matrix.scale(mu/HBAR);


    self.detected_particle = Some(DetectedSpin{
        gamma_matrix,
        isotope,
        weighted_coordinates,
        spin_multiplicity,
        transition,
        });


    Ok(())

  }
  //----------------------------------------------------------------------------
  // This method sets self.bath_spins_indices so that each element corresponds
  // to an element of self_bath_particles that has the potential to have a spin.
  fn set_spins(&mut self,config: &Config) -> Result<(),CluEError>
  {
    

    let n_particles = self.bath_particles.len();
    self.primary_cell_indices = Vec::<usize>::with_capacity(n_particles);
    self.cell_indices = Vec::<Vec::<Option<usize>>>::with_capacity(n_particles);
    
    for idx in 0..n_particles{
      self.primary_cell_indices.push(idx);
      self.cell_indices.push(vec![Some(idx)]);
    }

    self.particle_config_ids = self.pair_particle_configs(&config.particles)?;
    self.extracell_particle_config_ids 
      = self.pair_particle_configs(&config.extracell_particles)?;

    // Find spins.
    //let n_spins = self.count_spins(config);

    //self.bath_spins_indices = Vec::<usize>::with_capacity(n_spins);


    for (idx, particle) in self.bath_particles.iter_mut().enumerate(){



      let mut spin_mult = particle.isotope.spin_multiplicity();
      if let Some(id) = self.particle_config_ids[idx]{
        if let Some(mult) 
          = config.max_spin_multiplicity_for_particle_config(id){
          spin_mult = mult;
        }
      }
      particle.active = spin_mult > 1;

    }

    Ok(()) 
  }
  //----------------------------------------------------------------------------
  // This method counts the to an elements of self_bath_particles that have
  // the potential to have a spin.
  /*
  fn count_spins(&self,config: &Config) -> usize{
    let mut n_spins = 0;
    for (idx, particle) in self.bath_particles.iter().enumerate(){

      let mut spin_mult = particle.isotope.spin_multiplicity();
      if let Some(id) = self.particle_config_ids[idx]{
        if let Some(mult) 
          = config.max_spin_multiplicity_for_particle_config(id){
          spin_mult = mult;
        }
      }

      if spin_mult <= 1 { continue; }

      n_spins += 1;
    }
    n_spins
  }
  */
  //----------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------
  // The method uses the PBCs to move atoms near atoms they are bonded to.
  fn reconnect_bonds(&mut self) -> Result<(),CluEError>
  {

    // TODO: Decide on better error handeling here.
    match self.cell_offsets.len(){
      0 => return Ok(()),
      3 => (),
      _ => return Err(CluEError::InorrectNumberOfCellOffsets(
            self.cell_offsets.len(),3)),
    }

    const MAXBOND: f64 = 3.0; // Angstroms.

    for idx0 in 0..self.bath_particles.len() {


      // Get the list of bonded indices, if there are any.
      let connections: &Vec::<usize>;
      if let Some(cnc)= self.connections.get_neighbors(idx0){
        connections = cnc;
      }else{
        continue;
      }

      // Get the coordinates of the atom of interest.
      let r0 = self.bath_particles[idx0].coordinates.clone();

      for idx in connections.iter(){

        // Get the coordinates of the bonded atom.
        let r = self.bath_particles[*idx].coordinates.clone();

        // Check if the atoms are near each other.
        let delta_r = &r -&r0;
        if delta_r.norm() < MAXBOND { continue; }

        // Loop through all 3 spatial dimensions.
        for ix in 0..3 {

          // Find the PBC copy of the neghbor that is closest.
          let r1 = space_3d::minimize_absolute_difference_for_vector3d_step(
            &r,&r0,&self.cell_offsets[ix]);

          self.bath_particles[*idx].coordinates = r1;
        }
      }

    }
    Ok(())
 }
 //-----------------------------------------------------------------------------
 // This function finds all the exchange groups and sets up the 
 // data structures for including exchange couplings in the spin Hamiltonian.
 // All the exchange couplings are initialized to zero.
 fn set_exchange_groups(&mut self){

   let exchange_groups = self.find_exchange_groups();
   
   // Initialize `exchange_group_ids`.
   let mut exchange_group_ids: Vec::<Option<usize>> 
     = (0..self.number()).map(|_| None).collect();

   // Loop through all exchange groups.
   for (ii,exchange_group) in exchange_groups.iter().enumerate(){
   
     //Loop through all hydrogens in the group.
     for h in exchange_group.indices(){

       // Set the group id for each hydrogen. 
       exchange_group_ids[h] = Some(ii);
     }
   }

   // Set default coupling to zero.
   let exchange_couplings: Vec::<f64> = (0..exchange_groups.len())
     .map(|_| 0.0).collect();

   self.exchange_groups = Some( ExchangeGroupManager{
     exchange_groups,
     exchange_group_ids,
     exchange_couplings,
   });
 }  
 //-----------------------------------------------------------------------------
 // This function searches through the structure and finds all methyl groups
 // and primary amonium groups, and returns the info as a 
 // `Vec::<ExchangeGroup>`.
 fn find_exchange_groups(&self) -> Vec::<ExchangeGroup> {
 
   // Five atoms are required to form a methyl group.
   let n_max: usize = self.number()/5;

   let mut exchange_groups = Vec::<ExchangeGroup>::with_capacity(n_max);

   for ii in 0..self.number(){

     let hydrogens = self.get_methyl_hydrogen_indices(ii);
      if let Some([h0,h1,h2]) = hydrogens {
        let r_carbon = self.bath_particles[ii].coordinates.clone();
        let r_h0 = self.bath_particles[h0].coordinates.clone();
        let r_h1 = self.bath_particles[h1].coordinates.clone();
        let r_h2 = self.bath_particles[h2].coordinates.clone();

        
        if self.bath_particles[ii].element == Element::Carbon{
          exchange_groups.push(ExchangeGroup::Methyl(
                C3Rotor::from(r_carbon, r_h0, r_h1, r_h2, [h0,h1,h2]) ));

        }else if self.bath_particles[ii].element == Element::Nitrogen{
          exchange_groups.push(ExchangeGroup::PrimaryAmonium(
                C3Rotor::from(r_carbon, r_h0, r_h1, r_h2, [h0,h1,h2]) ));
        }
      } 
   }

   exchange_groups
 }  
 //-----------------------------------------------------------------------------
 // This method takes an index of a methyl group carbon or of a primary amonium 
 // nitrogen and return the indices of the three bonded hyrogens.
 // If the index is not a C or N bonded to exacly three H, the the return value
 // is `None`; 
 // otherwise the return value is `Some([h[0], h[1], h[2] ])`,
 // where `h` is a list of the indices of the hydrogens in `bath_particles`.
 fn get_methyl_hydrogen_indices(&self, index: usize) ->Option<[usize;3] >{
   
   let particle = &self.bath_particles[index];

   if particle.element != Element::Carbon 
    && particle.element != Element::Nitrogen {
     return None;
   }

   // Get the list of bonded indices, if there are any.
   let connections: &Vec::<usize>;
   if let Some(cnc)= self.connections.get_neighbors(index){
     connections = cnc;
   }else{
     return None;
   }
   
   // Methyls and primary amoniums have four atoms bonded to the central
   // carbon and nitrogen respectively.
   if connections.len() != 4 { 
     return None;
   }

   let mut hydrons = Vec::<usize>::with_capacity(3);

   for idx in connections.iter() {

     if self.bath_particles[*idx].element == Element::Hydrogen{
       
       // Skip methane and amonium, which are not implemented. 
       if hydrons.len()==3{ return None;}

       hydrons.push(*idx);
     }
   }


   if hydrons.len()!=3{ return None;}

   Some([hydrons[0], hydrons[1], hydrons[2] ])
 }  
 //-----------------------------------------------------------------------------
   //----------------------------------------------------------------------------
  // TODO: find_cosubstitution_groups() is slow for large systems.
  // This function goes through each particle, applies secondary filters,
  // and records the cosubstitution_groups.
  fn find_cosubstitution_groups(&mut self, config: &Config)
    -> Result<(),CluEError>
  {

    let mut cosubstitution_group_ids
      = Vec::<Option<usize>>::with_capacity(self.bath_particles.len());
    for _ii in 0..self.bath_particles.len(){
      cosubstitution_group_ids.push(None);
    }

    let mut current_cosub_id = 0;


    // Loop through bath particles.
    for idx in 0..self.bath_particles.len(){

      let idx0 = self.primary_cell_indices[idx];
      // Check if there are custom properties for this particle.
      let Some(id) = self.particle_config_ids[idx0] else{continue;};
      let Some(properties) = &config.particles[id].properties else{continue;};
      if cosubstitution_group_ids[idx].is_some(){
        continue;
      }

      // Check if this particle has a filter.
      let filter: ParticleFilter;
      if let Some(fltr) = &config.particles[id].filter{
        filter= fltr.clone();
      }else{
        filter = ParticleFilter::new();
      }

      // Set cosubstitution group for this particle.
      update_cosubstitution_ids(&mut cosubstitution_group_ids,
        idx, current_cosub_id, filter.clone(),properties,self)?;
      current_cosub_id += 1;

    }
    for id in cosubstitution_group_ids.iter_mut(){
      if (*id).is_some() {continue;}
      *id = Some(current_cosub_id);
      current_cosub_id += 1;
    }
    let mut cosubstitution_groups
      = Vec::<Vec::<usize>>::with_capacity(current_cosub_id);
    for _ii in 0..current_cosub_id{
      cosubstitution_groups.push(Vec::<usize>::new());
    }
    for (idx,id_opt) in cosubstitution_group_ids.iter().enumerate(){
      let Some(id) = *id_opt else{
        return Err(CluEError::UnassignedCosubstitutionGroup(idx))
      };
      cosubstitution_groups[id].push(idx);
    }
    self.cosubstitution_groups = cosubstitution_groups;

    Ok(())
  }

 //-----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//------------------------------------------------------------------------------
// This function finds all the particles that cosubstitute with particle idx,
// and records them in cosubstitution_group_ids.
fn update_cosubstitution_ids(
    cosubstitution_group_ids: &mut [Option<usize>],
    idx: usize,
    current_cosub_id: usize,
    filter: ParticleFilter,
    properties: &ParticleProperties,
    structure: &Structure)
    -> Result<(),CluEError>
  {
      let idx0 = structure.primary_cell_indices[idx];
      let mut indices = Vec::<usize>::new();

      // Check if properties defines a cosubstitution set.
      match properties.cosubstitute{
        Some(SecondaryParticleFilter::Bonded) => {
          if let Some(bonded) = structure.connections.get_neighbors(idx0){
            indices = filter.filter_indices(structure,bonded);
          }
        },
        Some(SecondaryParticleFilter::Filter) 
          => indices = filter.filter(structure),
        Some(SecondaryParticleFilter::Particle) => (),
        Some(SecondaryParticleFilter::SameMolecule) => {
          let mol_id = structure.molecule_ids[idx0];
          indices = filter.filter_indices(structure,
              &structure.molecules[mol_id]);
        },
        None => (),
      }
      for index in indices.iter(){
        // Assign the particle to this cosubstitution group.
        cosubstitution_group_ids[*index] = Some(current_cosub_id)
      }

    Ok(())
  }
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



#[cfg(test)]
mod tests{
  use super::*;
  use crate::structure::pdb;
  use crate::config::Config;
  use crate::config::particle_config::ParticleConfig;
  use crate::config::DetectedSpinCoordinates;
  use rand::SeedableRng;

  //----------------------------------------------------------------------------
  #[test]
  fn test_find_cosubstitution_groups(){
    let mut rng = ChaCha20Rng::from_entropy();
    let filename = "./assets/a_TEMPO_a_water_a_glycerol.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    config.set_defaults().unwrap();
    structure.build_primary_structure(&mut rng, &config).unwrap();

    let mut filter_nx = ParticleFilter::new();
    filter_nx.elements = vec![Element::Hydrogen];
    filter_nx.not_bonded_elements = vec![Element::Oxygen];

    let mut properties = ParticleProperties::new();
    properties.cosubstitute = Some(SecondaryParticleFilter::SameMolecule);

    let mut particle_configs =vec![
      ParticleConfig::new("non-exchangeable".to_string())];
    particle_configs[0].properties = Some(properties);
    particle_configs[0].filter = Some(filter_nx);

    config.particles = particle_configs;

    structure.particle_config_ids 
      = structure.pair_particle_configs(&config.particles).unwrap();

    structure.find_cosubstitution_groups(&config).unwrap();
    let expected = vec![
        vec![2,3,4,6,7,8,10,11,13,14,16,17,20,21,22,24,25,26],
        vec![30,31,35,39,40],
        vec![0],vec![1],vec![5],vec![9],vec![12],
        vec![15],vec![18],vec![19],vec![23],vec![27],
        vec![28],vec![29],
        vec![32],vec![33],vec![34],vec![36],vec![37],
        vec![38],vec![41],vec![42],
        vec![43],vec![44],vec![45]
    ];

    assert_eq!(structure.cosubstitution_groups.len(), expected.len());
    assert_eq!(structure.cosubstitution_groups, expected);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_primary_structure_water(){
    let mut rng = ChaCha20Rng::from_entropy();
    let filename = "./assets/water.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    config.set_defaults().unwrap();
    structure.build_primary_structure(&mut rng, &config).unwrap();

  
    println!("DB: {:#?}",structure.connections);
    assert_eq!(structure.molecules.len(),1);
    assert_eq!(structure.molecules[0],vec![0,1,2]);
    assert_eq!(structure.molecule_ids,vec![0,0,0]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_primary_structure_tempo(){
    let mut rng = ChaCha20Rng::from_entropy();
    let filename = "./assets/TEMPO.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();

    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    config.set_defaults().unwrap();
    structure.build_primary_structure(&mut rng, &config).unwrap();

    assert_eq!(structure.bath_particles
        .iter().filter(|p| (*p).active).count(), 19);
    let exchange_group_manager = structure.exchange_groups
                                 .as_ref().unwrap();
    assert_eq!(exchange_group_manager.exchange_groups.len(), 4);

    let mut ids = vec![None; 29];
    ids[2] = Some(0);
    ids[3] = Some(0);
    ids[4] = Some(0);

    ids[6] = Some(1);
    ids[7] = Some(1);
    ids[8] = Some(1);

    ids[20] = Some(2);
    ids[21] = Some(2);
    ids[22] = Some(2);

    ids[24] = Some(3);
    ids[25] = Some(3);
    ids[26] = Some(3);

    for (idx, id) in exchange_group_manager.exchange_group_ids.iter()
      .enumerate(){
      assert_eq!(*id,ids[idx]);
    }

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_primary_structure_tempo_wat_gly_7nm(){
    let mut rng = ChaCha20Rng::from_entropy();
    // n    : Molecules    : Hydrons
    // 1    : TEMPO        :    18
    // 1500 : glycerols    : 12000
    // 7469 : waters       : 14938
    //
    // total hydrons =  26956.
    //
    let n_wat = 7469;
    let n_gly = 1500;

    let filename = "./assets/TEMPO_wat_gly_70A.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    config.set_defaults().unwrap();
    structure.build_primary_structure(&mut rng, &config).unwrap();
    
    assert_eq!(structure.molecules.len(), 1 + n_wat + n_gly);

    assert_eq!(structure.molecules[0].len(),29);
    for ii in 1..=1500{
      assert_eq!(structure.molecules[ii].len(),14);
    }
    for ii in 1501..structure.molecules.len(){
      assert_eq!(structure.molecules[ii].len(),3);
    }
  }
}



