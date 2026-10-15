use crate::config::{Config,DetectedSpinCoordinates};
use crate::CluEError;
use crate::cluster::connected_subgraphs::separate_into_connected_subgraphs;
use crate::physical_constants::{ANGSTROM,HBAR,MUB,MUN};
use crate::elements::Element;
use crate::isotopes::Isotope;
use crate::quantum::tensors::construct_symmetric_tensor_from_tensor_specifier;
use crate::structure::ParticleFilter;
use crate::structure::{Structure, exchange_groups::*}; 
use crate::structure::{ParticleConfig,ParticleProperties};
use crate::config::particle_config::Cosubstitute;
use crate::structure::IntegrationGrid;
use crate::structure::DetectedSpin;
use crate::structure::reciprocal_vectors;
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

    // Make molecules whole across the unit cell boundary.
    // Particles can end up outside the unit cell; `set_cell_shifts()` 
    // accounts for this when it chooses the number of unit cells.
    self.reconnect_bonds()?;

    // Particle filters need the primary cell indices.
    self.init_primary_cell_indices();

    // Set the detected spin as the origin, so that distance selections 
    // in the groups below are measured from the detected spin.
    let weighted_coordinates = self.set_detected_spin_position(config)?;
    
    self.set_spins(config)?;

    // Report particles that the groups cannot select for lack of a field.
    self.report_particles_missing_selection_fields(config);

    // Set methyl and primary amonium groups..
    self.set_exchange_groups();

    self.find_cosubstitution_groups(config)?;

    // Set the detected spin properties.
    self.set_detected_spin(rng,config,weighted_coordinates)?;

    Ok(())
  }  
  //----------------------------------------------------------------------------
  // This method sets up `primary_cell_indices`, `cell_indices`, and 
  // `cell_ids` for a structure that contains only the primary cell.
  fn init_primary_cell_indices(&mut self){
    let n_particles = self.bath_particles.len();
    self.primary_cell_indices = (0..n_particles).collect::<Vec::<usize>>();
    self.cell_ids = vec![0; n_particles];
    self.cell_indices = (0..n_particles).map(|idx| vec![Some(idx)])
      .collect::<Vec::<Vec::<Option<usize>>>>();
  }
  //----------------------------------------------------------------------------
  // This method finds the detected spin position, sets it as the origin by
  // translating the bath particles, and returns the translated
  // detected spin coordinates.
  // Distances are measured from the detected spin, so a `CentroidOverGroup`
  // group that uses a distance selection is an error.
  fn set_detected_spin_position(&mut self, config: &Config) 
      -> Result<IntegrationGrid,CluEError>
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

          if filter.within_distance.is_some() 
            || filter.not_within_distance.is_some(){
            return Err(CluEError::CentroidGroupHasDistanceSelection(
                  label.to_string()));
          }

          // Each group must contribute to the centroid.
          if filter.filter(self).is_empty(){
            return Err(CluEError::DetectedSpinSelectsNoParticles(format!(
                    "group \"{}\", but that group selects no particles",
                    label)));
          }

          filters.push(filter);
        }
        
        let r_ave = self.centroid_over_groups(&filters)?;
        weighted_coordinates.push(vec![r_ave.x(),r_ave.y(),r_ave.z()],1.0)?;
      },
      DetectedSpinCoordinates::CentroidOverSerials(serials) => {
        // Each serial must contribute to the centroid.
        let missing = serials.iter().filter(|&&serial| 
            !self.bath_particles.iter().any(|p| p.serial == Some(serial)))
          .map(|serial| serial.to_string()).collect::<Vec::<String>>();
        if !missing.is_empty(){
          return Err(CluEError::DetectedSpinSelectsNoParticles(format!(
                  "serials {:?}, but serial {} {} not in the structure",
                  serials, missing.join(", "), 
                  if missing.len() == 1 {"is"} else {"are"})));
        }

        let r_ave = self.centroid_over_serials(serials.clone())?;
        weighted_coordinates.push(vec![r_ave.x(),r_ave.y(),r_ave.z()],1.0)?;
      },
      DetectedSpinCoordinates::XYZ(r) => {
        weighted_coordinates.push(vec![r.x(),r.y(),r.z()],1.0)?;
      },
      DetectedSpinCoordinates::ProbabilityDistribution(grid) => 
        weighted_coordinates = grid.clone(),
    
    }

    // The origin is the weighted mean of the detected spin coordinates.
    let Some(origin) = weighted_coordinates.weighted_mean() else{
      return Err(CluEError::InvalidSpinGridWeights(
            weighted_coordinates.weights.iter().fold(0.0, |sum, w| sum + w)));
    };
    if origin.iter().any(|x| !x.is_finite()){
      return Err(CluEError::DetectedSpinPositionNotFinite(format!(
              "[{}, {}, {}]", origin[0]/ANGSTROM, origin[1]/ANGSTROM, 
              origin[2]/ANGSTROM)));
    }
    self.pdb_origin = Vector3D::from([origin[0],origin[1],origin[2]])
      .scale(-1.0);

    weighted_coordinates.translate_mut(
        &[self.pdb_origin.x(),self.pdb_origin.y(),self.pdb_origin.z() ]); 

    for particle in self.bath_particles.iter_mut(){
      particle.coordinates = &particle.coordinates + &self.pdb_origin;
    }

    Ok(weighted_coordinates)
  }
  //----------------------------------------------------------------------------
  // This method places the detected spin into the system.
  // The detected spin position should already be set by 
  // `set_detected_spin_position()`.
  fn set_detected_spin(&mut self, rng: &mut ChaCha20Rng, config: &Config,
      weighted_coordinates: IntegrationGrid) 
      -> Result<(),CluEError>
  {

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
  // This method prints a message for each group with a list of serials, 
  // atom names, residues, or residue sequence numbers, 
  // if some particles do not have that field and so cannot be selected by it.
  fn report_particles_missing_selection_fields(&self, config: &Config){
    let mut reported = std::collections::HashSet::<String>::new();
    for particle_config in config.particles.iter()
      .chain(config.extracell_particles.iter())
    {
      let Some(filter) = &particle_config.filter else{ continue; };
      if !reported.insert(particle_config.label.clone()){ continue; }

      for (key, field, n_missing) 
        in filter.count_particles_missing_fields(self)
      {
        let (particles_have, them) = match n_missing{
          0 => continue,
          1 => ("particle has", "it"),
          _ => ("particles have", "them"),
        };
        println!("Group \"{}\": {} {} no {}, \
so selection.{} cannot select {}.", 
            particle_config.label, n_missing, particles_have, field, key, them);
      }
    }
  }
  //----------------------------------------------------------------------------
  // This method sets self.bath_spins_indices so that each element corresponds
  // to an element of self_bath_particles that has the potential to have a spin.
  // The extracell groups are set in `set_extracell_groups()`, 
  // once the cell offsets are known.
  fn set_spins(&mut self,config: &Config) -> Result<(),CluEError>
  {

    self.particle_config_ids 
      = self.pair_particle_configs(&config.particles,None)?;

    // Find spins.
    //let n_spins = self.count_spins(config);

    //self.bath_spins_indices = Vec::<usize>::with_capacity(n_spins);


    // Only the primary cell groups apply here; 
    // the copies in other unit cells are set in `extend_structure()`.
    let actives = (0..self.bath_particles.len()).map(|idx| 
        self.max_spin_multiplicity_for_particle(idx, &self.particle_config_ids,
          &config.particles) > 1
      ).collect::<Vec::<bool>>();

    for (particle, active) in self.bath_particles.iter_mut().zip(actives){
      particle.active = active;
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
  // This method uses the PBCs to make molecules whole.
  // Starting from the lowest index particle of each molecule, each bonded
  // particle is moved by the lattice vector that brings it closest to the
  // particle it is bonded to, and the search continues from that particle.
  // Particles are only moved by lattice vectors, so the periodic system is
  // unchanged, and molecules that are already whole are not moved.
  fn reconnect_bonds(&mut self) -> Result<(),CluEError>
  {

    // TODO: Decide on better error handeling here.
    match self.cell_offsets.len(){
      0 => return Ok(()),
      3 => (),
      _ => return Err(CluEError::InorrectNumberOfCellOffsets(
            self.cell_offsets.len(),3)),
    }

    // Reciprocal vectors, which give the fractional coordinates 
    // along the cell edges a, b, and c.
    let Some(reciprocal_vectors) = reciprocal_vectors(&self.cell_offsets) 
    else{
      println!("The unit cell has no volume. \
Molecules will not be made whole across the cell boundary.");
      return Ok(());
    };
    let a = self.cell_offsets[0].clone();
    let b = self.cell_offsets[1].clone();
    let c = self.cell_offsets[2].clone();

    // Moving particles by the lattice vectors of a unit cell that does not 
    // describe the structure would scramble it.
    self.check_structure_fits_cell(&reciprocal_vectors)?;

    let mut is_placed = vec![false; self.bath_particles.len()];

    for molecule in self.molecules.iter(){
      let Some(&root) = molecule.iter().min() else{ continue; };
      is_placed[root] = true;
      let mut to_visit = vec![root];

      while let Some(idx0) = to_visit.pop(){
        let Some(neighbors) = self.connections.get_neighbors(idx0) else{
          continue;
        };

        for &idx in neighbors.iter(){
          if is_placed[idx]{ continue; }
          is_placed[idx] = true;
          to_visit.push(idx);

          // Find the lattice vector closest to the bond vector.
          let delta_r = &self.bath_particles[idx].coordinates
            - &self.bath_particles[idx0].coordinates;
          let n = reciprocal_vectors.iter()
            .map(|g| g.dot(&delta_r).round()).collect::<Vec::<f64>>();
          if n.iter().all(|&ni| ni == 0.0){ continue; }

          // Move the bonded particle next to particle idx0.
          let lattice_vector 
            = &(&a.scale(n[0]) + &b.scale(n[1])) + &c.scale(n[2]);
          self.bath_particles[idx].coordinates 
            = &self.bath_particles[idx].coordinates - &lattice_vector;
        }
      }
    }

    // Structures that are bonded to their own periodic copies, such as 
    // infinite polymers, cannot be made whole.
    const MAX_BOND: f64 = 3.0*ANGSTROM;
    let n_long = self.count_long_bonds(MAX_BOND);
    if n_long > 0{
      println!("{} bonds are longer than {} angstroms after making \
molecules whole across the cell boundary.", n_long, MAX_BOND/ANGSTROM);
    }

    Ok(())
 }
 //-----------------------------------------------------------------------------
 // This function counts the bonds that are longer than `max_bond`.
 fn count_long_bonds(&self, max_bond: f64) -> usize{
   let mut n_long = 0;
   for idx0 in 0..self.connections.len(){
     let Some(neighbors) = self.connections.get_neighbors(idx0) else{
       continue;
     };
     for &idx in neighbors.iter(){
       if idx <= idx0{ continue; }
       let r = &self.bath_particles[idx].coordinates 
         - &self.bath_particles[idx0].coordinates;
       if r.norm() > max_bond{
         n_long += 1;
       }
     }
   }
   n_long
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
    self.cosubstitution_groups = self.build_cosubstitution_groups(
        &self.particle_config_ids, &config.particles, None)?;

    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function returns the cosubstitution groups for the primary cell
  // particles paired to `particle_configs` by `config_ids`.
  // If `cell` is `Some(cell_id)`, the groups are for the copies of the 
  // primary cell particles in cell `cell_id`.
  pub(crate) fn build_cosubstitution_groups(&self, 
      config_ids: &[Option<usize>], particle_configs: &[ParticleConfig],
      cell: Option<usize>)
    -> Result<Vec::<Vec::<usize>>,CluEError>
  {

    let n_particles = config_ids.len();

    let mut cosubstitution_group_ids
      = Vec::<Option<usize>>::with_capacity(n_particles);
    for _ii in 0..n_particles{
      cosubstitution_group_ids.push(None);
    }

    let mut current_cosub_id = 0;


    // Loop through bath particles.
    for idx in 0..n_particles{

      let idx0 = self.primary_cell_indices[idx];
      // Check if there are custom properties for this particle.
      let Some(id) = config_ids[idx0] else{continue;};
      let Some(properties) = &particle_configs[id].properties else{continue;};
      if cosubstitution_group_ids[idx].is_some(){
        continue;
      }

      // Check if this particle has a filter.
      let filter: ParticleFilter;
      if let Some(fltr) = &particle_configs[id].filter{
        filter= fltr.clone();
      }else{
        filter = ParticleFilter::new();
      }

      // Set cosubstitution group for this particle.
      update_cosubstitution_ids(&mut cosubstitution_group_ids,
        idx, current_cosub_id, filter.clone(),properties,self,cell)?;
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

    Ok(cosubstitution_groups)
  }

 //-----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//------------------------------------------------------------------------------
// This function finds all the particles that cosubstitute with particle idx,
// and records them in cosubstitution_group_ids.
// If `cell` is `Some(cell_id)`, the filter is evaluated for the copies
// in cell `cell_id`.
fn update_cosubstitution_ids(
    cosubstitution_group_ids: &mut [Option<usize>],
    idx: usize,
    current_cosub_id: usize,
    filter: ParticleFilter,
    properties: &ParticleProperties,
    structure: &Structure,
    cell: Option<usize>)
    -> Result<(),CluEError>
  {
      let idx0 = structure.primary_cell_indices[idx];

      let filter_indices = |subset: &[usize]| -> Vec::<usize> {
        match cell{
          Some(cell_id) 
            => filter.filter_indices_in_cell(structure,subset,cell_id),
          None => filter.filter_indices(structure,subset),
        }
      };

      // Check if properties defines a cosubstitution set.
      // Each set contains particle `idx`, since a particle passes the filter
      // of its own group.
      let indices = match properties.cosubstitute{
        Some(Cosubstitute::Group) => match cell{
          Some(cell_id) => filter.filter_in_cell(structure,cell_id),
          None => filter.filter(structure),
        },
        Some(Cosubstitute::SameMolecule) => {
          let mol_id = structure.molecule_ids[idx0];
          filter_indices(&structure.molecules[mol_id])
        },
        Some(Cosubstitute::Particle) | None => vec![idx],
      };
      for index in indices.iter(){
        // A set is only built for an unassigned particle, and contains all
        // particles of its group that it shares a set with, so no particle 
        // is assigned twice.
        debug_assert!(cosubstitution_group_ids[*index].is_none());

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
  use crate::cluster::adjacency::AdjacencyList;
  use crate::structure::particle::Particle;
  use crate::physical_constants::ANGSTROM;
  use rand::SeedableRng;

  //----------------------------------------------------------------------------
  // This function moves the particles in `indices` by `shift`.
  fn shift_particles(structure: &mut Structure, indices: &[usize],
      shift: &Vector3D)
  {
    for &idx in indices.iter(){
      structure.bath_particles[idx].coordinates 
        = &structure.bath_particles[idx].coordinates + shift;
    }
  }
  //----------------------------------------------------------------------------
  // This function reads TEMPO, optionally replaces its cell, and splits it 
  // across the cell boundary: a methyl carbon and two of its hydrogens are 
  // moved by +a, and the oxygen is moved by -b + c.
  // It returns the split structure and the original coordinates.
  fn split_tempo(cell: Option<Vec::<Vector3D>>) 
    -> (Structure, Vec::<Vector3D>)
  {
    let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    if let Some(cell_edges) = cell{
      structure.cell_offsets = cell_edges;
    }
    assert_eq!(structure.cell_offsets.len(),3);
    let r_original = structure.bath_particles.iter()
      .map(|p| p.coordinates.clone()).collect::<Vec::<Vector3D>>();

    // Find a methyl group that does not contain the first particle.
    let (idx_c, hydrogens) = (1..structure.number()).find_map(|idx|
      structure.get_methyl_hydrogen_indices(idx).map(|h| (idx,h))
    ).unwrap();
    let idx_o = structure.bath_particles.iter()
      .position(|p| p.element == Element::Oxygen).unwrap();

    let a = structure.cell_offsets[0].clone();
    let b = structure.cell_offsets[1].clone();
    let c = structure.cell_offsets[2].clone();
    shift_particles(&mut structure,
        &[idx_c, hydrogens[0], hydrogens[1]], &a);
    shift_particles(&mut structure, &[idx_o], &(&c - &b));

    (structure, r_original)
  }
  //----------------------------------------------------------------------------
  // This function checks that `structure` has the coordinates `r_expected`.
  fn assert_coordinates(structure: &Structure, r_expected: &[Vector3D]){
    for (particle, r) in structure.bath_particles.iter().zip(r_expected){
      assert!((&particle.coordinates - r).norm() < 1e-6*ANGSTROM,
          "{} != {}",particle.coordinates, r);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_reconnect_bonds(){
    // Molecules split across the cell boundary are made whole again.
    let (mut structure, r_original) = split_tempo(None);
    (structure.molecules, structure.molecule_ids) 
      = separate_into_connected_subgraphs(&structure.connections);
    structure.reconnect_bonds().unwrap();
    assert_coordinates(&structure,&r_original);
    assert_eq!(structure.count_long_bonds(3.0*ANGSTROM),0);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_reconnect_bonds_hexagonal_cell(){
    // The lattice vectors of oblique cells are handled.
    let gamma = 120.0f64.to_radians();
    let cell = vec![
      Vector3D::from([20.0, 0.0, 0.0]).scale(ANGSTROM),
      Vector3D::from([20.0*gamma.cos(), 20.0*gamma.sin(), 0.0]).scale(ANGSTROM),
      Vector3D::from([0.0, 0.0, 25.0]).scale(ANGSTROM),
    ];
    let (mut structure, r_original) = split_tempo(Some(cell));
    (structure.molecules, structure.molecule_ids) 
      = separate_into_connected_subgraphs(&structure.connections);
    structure.reconnect_bonds().unwrap();
    assert_coordinates(&structure,&r_original);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_reconnect_bonds_whole_molecules_unchanged(){
    // Molecules that are already whole are not moved at all.
    let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    let r_original = structure.bath_particles.iter()
      .map(|p| p.coordinates.clone()).collect::<Vec::<Vector3D>>();
    (structure.molecules, structure.molecule_ids) 
      = separate_into_connected_subgraphs(&structure.connections);
    structure.reconnect_bonds().unwrap();
    for (particle, r) in structure.bath_particles.iter().zip(r_original){
      assert!(particle.coordinates.x() == r.x() 
          && particle.coordinates.y() == r.y()
          && particle.coordinates.z() == r.z());
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_reconnect_bonds_periodic_ring(){
    // A ring that is bonded to its own periodic copies cannot be made whole,
    // so one long bond remains.
    let cell_length = 10.0*ANGSTROM;
    let particles = (0..4).map(|ii| Particle::new(Element::Carbon,
          2.5*ANGSTROM*(ii as f64), 0.0, 0.0)).collect::<Vec::<Particle>>();
    let mut connections = AdjacencyList::with_capacity(4);
    for ii in 0..4{
      connections.connect(ii,(ii+1)%4);
    }
    let cell = vec![
      Vector3D::from([cell_length, 0.0, 0.0]),
      Vector3D::from([0.0, cell_length, 0.0]),
      Vector3D::from([0.0, 0.0, cell_length]),
    ];
    let mut structure = Structure::new(particles,connections,cell);
    (structure.molecules, structure.molecule_ids) 
      = separate_into_connected_subgraphs(&structure.connections);
    structure.reconnect_bonds().unwrap();
    assert_eq!(structure.count_long_bonds(3.0*ANGSTROM),1);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_primary_structure_split_tempo(){
    // The methyl rotors of a split molecule are built from bonded atoms.
    let (mut structure, _r_original) = split_tempo(None);
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    config.set_defaults().unwrap();
    structure.build_primary_structure(&mut rng, &config).unwrap();

    let exchange_group_manager = structure.exchange_groups.as_ref().unwrap();
    assert_eq!(exchange_group_manager.exchange_groups.len(), 4);
    for exchange_group in exchange_group_manager.exchange_groups.iter(){
      let h0 = exchange_group.indices()[0];
      let idx_c = structure.connections.get_neighbors(h0).unwrap()[0];
      let r_c = &structure.bath_particles[idx_c].coordinates;
      assert!((exchange_group.centroid() - r_c).norm() < 1.1*ANGSTROM);
    }
  }

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
    properties.cosubstitute = Some(Cosubstitute::SameMolecule);

    let mut particle_configs =vec![
      ParticleConfig::new("non-exchangeable".to_string())];
    particle_configs[0].properties = Some(properties);
    particle_configs[0].filter = Some(filter_nx);

    config.particles = particle_configs;

    structure.particle_config_ids 
      = structure.pair_particle_configs(&config.particles,None).unwrap();

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
  fn test_cosubstitution_groups_have_no_empty_groups(){
    // Every cosubstitution group has particles, and every particle is in 
    // exactly one group.
    use crate::io::FromTOMLString;
    for cosubstitute in ["", r#"cosubstitute = "particle""#, 
        r#"cosubstitute = "group""#, r#"cosubstitute = "filter""#,
        r#"cosubstitute = "same_molecule""#]
    {
      let config = Config::from_toml_string(&format!(r##"
          input_structure_file = "assets/TEMPO.pdb"
          radius = 20
          detected_spin.position = [28,29]
          replicate_unit_cell = false
          [[groups]]
          name = "h"
          selection.elements = ["H"]
          1H.abundance = 0.5
          2H.abundance = 0.5
          {}
          "##, cosubstitute)).unwrap();
      let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
      let mut rng = ChaCha20Rng::seed_from_u64(0);
      structure.build_primary_structure(&mut rng,&config).unwrap();

      let n_particles = structure.number();
      let groups = &structure.cosubstitution_groups;
      assert!(groups.iter().all(|group| !group.is_empty()), "{}", cosubstitute);
      let mut members = groups.concat();
      members.sort();
      assert_eq!(members, (0..n_particles).collect::<Vec::<usize>>());

      let n_h = 18;
      let expected_n_groups = match cosubstitute{
        "" | r#"cosubstitute = "particle""# => n_particles,
        _ => n_particles - n_h + 1,
      };
      assert_eq!(groups.len(), expected_n_groups, "{}", cosubstitute);
    }
  }
  //----------------------------------------------------------------------------
  // This function builds the primary structure of TEMPO with the detected 
  // spin position `position` and the groups `groups`.
  fn build_tempo_with_spin_position(position: &str, groups: &str)
    -> Result<Structure,CluEError>
  {
    use crate::io::FromTOMLString;
    let config = Config::from_toml_string(&format!(r##"
        input_structure_file = "assets/TEMPO.pdb"
        radius = 20
        detected_spin.position = {}
        replicate_unit_cell = false
        [[groups]]
        name = "tempo_n"
        selection.elements = ["N"]
        [[groups]]
        name = "nothing"
        selection.residues = ["XYZ"]
        {}
        "##, position, groups)).unwrap();
    let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    structure.build_primary_structure(&mut rng,&config)?;
    Ok(structure)
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_detected_spin_position_must_select_particles(){
    // Every group and serial that places the detected spin must match a 
    // particle, and the position must be finite.
    let empty_group = Some(CluEError::DetectedSpinSelectsNoParticles(
          "group \"nothing\", but that group selects no particles"
          .to_string()));
    assert_eq!(build_tempo_with_spin_position(r#"["nothing"]"#,"").err(),
        empty_group.clone());
    assert_eq!(build_tempo_with_spin_position(r#"["tempo_n","nothing"]"#,"")
        .err(), empty_group);
    assert_eq!(build_tempo_with_spin_position("[1000]","").err(),
        Some(CluEError::DetectedSpinSelectsNoParticles(
            "serials [1000], but serial 1000 is not in the structure"
            .to_string())));
    assert_eq!(build_tempo_with_spin_position("[28,1000,1001]","").err(),
        Some(CluEError::DetectedSpinSelectsNoParticles(
            "serials [28, 1000, 1001], but serial 1000, 1001 are not in the \
structure".to_string())));
    for position in ["[nan, 0.0, 0.0]", "[inf, 0.0, 0.0]"]{
      assert!(matches!(build_tempo_with_spin_position(position,""),
            Err(CluEError::DetectedSpinPositionNotFinite(_))), "{}", position);
    }

    // Positions that match particles still work.
    let structure = build_tempo_with_spin_position("[28,29]","").unwrap();
    let r_n = Vector3D::from([36.440, 36.900, 37.100]).scale(ANGSTROM);
    let r_o = Vector3D::from([35.290, 36.430, 37.810]).scale(ANGSTROM);
    assert_eq!(structure.pdb_origin, (&r_n + &r_o).scale(-0.5));
    let structure = build_tempo_with_spin_position(r#"["tempo_n"]"#,"")
      .unwrap();
    assert_eq!(structure.pdb_origin, r_n.scale(-1.0));
  }
  //----------------------------------------------------------------------------
  // This function builds TEMPO with a detected spin grid of two points: 
  // weight 0.9 at the N-O midpoint and weight 0.1 at 10 A from it along x.
  // It returns the structure, the N-O midpoint, and the original coordinates.
  fn build_tempo_with_spin_grid() -> (Structure, Vector3D, Vec::<Vector3D>){
    let r_n = Vector3D::from([36.440, 36.900, 37.100]).scale(ANGSTROM);
    let r_o = Vector3D::from([35.290, 36.430, 37.810]).scale(ANGSTROM);
    let r_no = (&r_n + &r_o).scale(0.5);
    let r_far = &r_no + &Vector3D::from([10.0*ANGSTROM, 0.0, 0.0]);

    let mut grid = IntegrationGrid::new(3);
    grid.push(vec![r_no.x(), r_no.y(), r_no.z()], 0.9).unwrap();
    grid.push(vec![r_far.x(), r_far.y(), r_far.z()], 0.1).unwrap();
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::ProbabilityDistribution(grid));
    config.set_defaults().unwrap();

    let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    let r_original = structure.bath_particles.iter()
      .map(|p| p.coordinates.clone()).collect::<Vec::<Vector3D>>();
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    structure.build_primary_structure(&mut rng,&config).unwrap();
    (structure, r_no, r_original)
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_detected_spin_grid_origin_is_weighted_mean(){
    let (structure, r_no, r_original) = build_tempo_with_spin_grid();

    // The weighted mean is 1 A from the N-O midpoint, 
    // and the unweighted mean would be 5 A from it.
    let r_mean = &r_no + &Vector3D::from([1.0*ANGSTROM, 0.0, 0.0]);
    assert!((&structure.pdb_origin + &r_mean).norm() < 1e-6*ANGSTROM);

    // The detected spin is at the origin.
    let detected_spin = structure.detected_particle.as_ref().unwrap();
    assert!(detected_spin.get_average_spin_position().unwrap().norm() 
        < 1e-6*ANGSTROM);

    // Vectors between the grid points and the particles are unchanged.
    let grid_original = [r_no.clone(), 
      &r_no + &Vector3D::from([10.0*ANGSTROM, 0.0, 0.0])];
    for (ii, r_grid) in grid_original.iter().enumerate(){
      let r_grid_new = detected_spin.weighted_coordinates.xyz(ii).unwrap();
      for (particle, r) in structure.bath_particles.iter().zip(&r_original){
        let delta = &(&r_grid_new - &particle.coordinates) - &(r_grid - r);
        assert!(delta.norm() < 1e-6*ANGSTROM);
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_detected_spin_position_in_nanometers(){
    // A position in nm gives the same origin as the serials it points to.
    use crate::io::FromTOMLString;
    let origin_of = |settings: &str| -> Vector3D{
      let config = Config::from_toml_string(&format!(r##"
          input_structure_file = "assets/TEMPO.pdb"
          radius = 2
          replicate_unit_cell = false
          {}
          "##, settings)).unwrap();
      let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
      let mut rng = ChaCha20Rng::seed_from_u64(0);
      structure.build_primary_structure(&mut rng,&config).unwrap();
      structure.pdb_origin
    };
    let r_serials = origin_of("detected_spin.position = [28, 29]");
    let r_nm = origin_of(r#"unit_of_distance = "nm"
          detected_spin.position = [3.5865, 3.6665, 3.7455]"#);
    assert!((&r_nm - &r_serials).norm() < 1e-6*ANGSTROM);

    // A one-point grid in nm agrees.
    let r_grid = origin_of(r#"unit_of_distance = "nm"
          detected_spin.position = [[3.5865, 3.6665, 3.7455, 1.0]]"#);
    assert!((&r_grid - &r_serials).norm() < 1e-6*ANGSTROM);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_detected_spin_position_empty_grid(){
    // The weights of an empty grid sum to zero, which set_defaults() rejects.
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::ProbabilityDistribution(
          IntegrationGrid::new(3)));
    assert_eq!(config.set_defaults(),
        Err(CluEError::SpinGridWeightsDoNotSumToOne(0.0, 1e-9)));

    // The structure also refuses grids whose weights do not have a positive
    // sum, in case set_defaults() was bypassed.
    let mut config = Config::new();
    config.detected_spin_position = Some(DetectedSpinCoordinates::XYZ(
          Vector3D::zeros()));
    config.set_defaults().unwrap();
    let mut grid = IntegrationGrid::new(3);
    grid.push(vec![0.0, 0.0, 0.0], 1.0).unwrap();
    grid.push(vec![1.0, 0.0, 0.0], -1.0).unwrap();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::ProbabilityDistribution(grid));
    let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    assert_eq!(structure.build_primary_structure(&mut rng,&config).err(),
        Some(CluEError::InvalidSpinGridWeights(0.0)));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_centroid_over_groups_of_no_particles(){
    let structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    let mut structure = structure;
    structure.init_primary_cell_indices();
    let mut filter = ParticleFilter::new();
    filter.residues = vec!["XYZ".to_string()];
    assert_eq!(structure.centroid_over_groups(&[&filter]),
        Err(CluEError::CentroidOfNoParticles));
    filter.residues = Vec::<String>::new();
    filter.serials = vec![28];
    assert_eq!(structure.centroid_over_groups(&[&filter]),
        Ok(Vector3D::from([36.440, 36.900, 37.100]).scale(ANGSTROM)));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_primary_structure_distance_selection(){
    // Distance selections are measured from the detected spin, not from
    // the origin of the input structure.
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    let mut config = Config::from_toml_file("./assets/solvent_box.toml")
      .unwrap();
    let filename = "./assets/water_gly_20pc_vol.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    let r_spin = Vector3D::from([25.0, 25.0, 25.0]).scale(ANGSTROM);
    config.detected_spin_position 
      = Some(DetectedSpinCoordinates::XYZ(r_spin.clone()));
    let r_pdb = structure.bath_particles.iter()
      .map(|p| p.coordinates.clone()).collect::<Vec::<Vector3D>>();

    structure.build_primary_structure(&mut rng, &config).unwrap();

    assert_eq!(structure.pdb_origin, r_spin.scale(-1.0));

    let (close_id,_) = config.find_particle_config("CLOSE_H").unwrap();
    let mut n_close = 0;
    for (idx,particle) in structure.bath_particles.iter().enumerate(){
      let r = particle.coordinates.norm();
      match structure.particle_config_ids[idx]{
        Some(id) if id == close_id => {
          assert!(r <= 4e-10);
          n_close += 1;
        },
        Some(_) => assert!(r >= 4e-10),
        None => (),
      }
      // The distance from the origin of the input structure is irrelevant.
      if r_pdb[idx].norm() <= 4e-10 && r > 4e-10{
        assert_ne!(structure.particle_config_ids[idx],Some(close_id));
      }
    }
    assert!(n_close > 0);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_centroid_over_group_with_distance_selection_errors(){
    // A group that places the detected spin cannot use a distance selection,
    // since distances are measured from the detected spin.
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    let filename = "./assets/TEMPO.pdb";
    let label = "nitroxide".to_string();

    let mut filter = ParticleFilter::new();
    filter.serials = vec![28,29];
    let mut config = Config::new();
    config.particles = vec![ParticleConfig::new(label.clone())];
    config.particles[0].filter = Some(filter.clone());
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverGroup(vec![label.clone()]) );
    config.set_defaults().unwrap();

    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    structure.build_primary_structure(&mut rng, &config).unwrap();

    for (r_max, r_min) in [(Some(10e-10),None), (None,Some(1e-10))]{
      let mut filter_r = filter.clone();
      filter_r.within_distance = r_max;
      filter_r.not_within_distance = r_min;
      config.particles[0].filter = Some(filter_r);

      let mut structure = pdb::parse_pdb(&filename,0).unwrap();
      let result = structure.build_primary_structure(&mut rng, &config);
      assert_eq!(result,
          Err(CluEError::CentroidGroupHasDistanceSelection(label.clone())));
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_primary_structure_water(){
    let mut rng = ChaCha20Rng::from_entropy();
    let filename = "./assets/water.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    let mut config = Config::new();
    // The detected spin is placed on the oxygen (serial 21030).
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![21030]) );
    config.set_defaults().unwrap();
    structure.build_primary_structure(&mut rng, &config).unwrap();

  
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



