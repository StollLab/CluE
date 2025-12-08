use crate::clue_errors::CluEError;
use crate::config::Config;
use crate::config::ReplicateUnitCell;
use crate::math;  
use crate::physical_constants::ANGSTROM;
use crate::space_3d::Vector3D;
use crate::structure::{Structure, exchange_groups::*};

use rand::seq::index::sample;
use rand_chacha::ChaCha20Rng;
use rand_distr::{Binomial, Distribution, Uniform};
use std::collections::HashMap;

impl Structure{
  /// This function takes a structure from build_primary_structure(), applies
  /// periodic boundary condition (if applicable), and sets isotope identities.
  pub fn build_extended_structure(&mut self,
      rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {

    // Construct offset vectors for each pbc application.
    self.set_cell_shifts(config)?;

    // Copy particles over to each PBC. 
    self.extend_structure(rng, config)?;

    // Set isotpic identities after adding voidable particles since the
    // non-void particles can potentially have multiple isotopic options. 
    self.set_isotopologue(rng, config)?;
  
    self.set_primary_cell_voidable_particles(rng, config)?; 

    self.update_exhange_groups(config)?;

    self.trim_system(config)?;

    self.check_primary_clashes(config)?;

    self.trim_pbc_clashes(config)?;

    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function applies periodic boudary conditions.
  fn extend_structure(&mut self, 
      rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {
    let n_particles0 = self.number();

    let n_to_reserve = n_particles0*(self.cell_offsets.len()-1);
    self.bath_particles.reserve( n_to_reserve );
    self.primary_cell_indices.reserve( n_to_reserve );

    for cindices in self.cell_indices.iter_mut(){
      cindices.reserve(self.cell_offsets.len()-1);
    }

    for particle_idx0 in 0..n_particles0{
      for _ii in 0..(self.cell_offsets.len()-1){
        self.cell_indices[particle_idx0].push(None); 
      }
    }
  

    // Check if there are any configurations.
    let particle_configs = &config.extracell_particles;

    // Initialize the HashMap for indices with the same void probability
    // The key is the void probability in parts per trillion.
    let mut property_map = HashMap::<u64, Vec::<usize>>::with_capacity(
        self.cosubstitution_groups.len());

    // Find all the cosubstitution_groups with the same substitution chance.
    for (cosub_idx,cosubstitution_group) in self.cosubstitution_groups.iter()
      .enumerate()
    {
      if cosubstitution_group.is_empty(){ continue; }

      let particle_idx0 = cosubstitution_group[0];
      if !self.bath_particles[particle_idx0].active {continue;}
      
      let mut p_keep = 1.0;

      // Check if this particle has a custom config.
      if let Some(config_id) = self.extracell_particle_config_ids[particle_idx0] 
      {

        // Check if this particle has any custom properties.
        if let Some(properties) = &particle_configs[config_id].properties 
        {

          // Get void probability.
          if let Some(p_remove) = properties.isotopic_distribution
            .void_probability{ 
          
              p_keep = 1.0 - p_remove;
          }
        }
      }
      

      // Determine the chance of not void, 
      // and convert it into parts per trillion (as a u64). 
      let ppt = (p_keep*1e12) as u64;

      // Add the indices to the property_map.
      if let Some(indices) = &mut property_map.get_mut(&ppt){
        indices.push(cosub_idx);
      }else{
        let mut indices = Vec::<usize>::with_capacity(n_particles0);
        indices.push(cosub_idx);
        property_map.insert(ppt,indices);
      }

    }

    // Loop through each unit cell other than the primary cell.
    for (cell_id,offset) in self.cell_offsets.iter().enumerate().skip(1){

      // Loop through all void probability, indices pairs.
      for (ppt, indices) in property_map.iter(){ 

        // Convert the keep probability to a floating point number.
        let p_keep = (*ppt as f64)/1e12; 
        let p_remove = 1.0 - p_keep;
        let p = if p_keep <= p_remove{
          p_keep
        }else{
          p_remove
        };
        
        // Generate a binomial distribution with n_group elements
        // a success chance of p.
        let n_groups = indices.len();
        let Ok(bin) = Binomial::new(n_groups as u64, p) else {
          return Err(CluEError::CannotSampleBinomialDistribution(
                n_groups,p));
        };


        // Sample from the binomial distribution for how many elements to keep. 
        let n_indices = bin.sample(rng) as usize;

        // Sample n_indices out of n_groups elements without repeats to keep.
        let random_indices = sample(rng,n_groups,n_indices);
        let mut selected_indices: Vec::<Option<usize>>;
        if p_keep <= p_remove{
          // random_indices contains indices to keep.
          selected_indices = (0..n_groups).map(|_n| None).collect();
          for idx in random_indices.iter(){
            selected_indices[idx] = Some(idx);
          }
        }else{
          // random_indices contains indices to remove.
          selected_indices = (0..n_groups).map(Some).collect();
          for idx in random_indices.iter(){
            selected_indices[idx] = None;
          }
        };


        // Loop over particles to keep.
        for cosub_idx_opt in selected_indices.iter(){
          let Some(cosub_idx) = cosub_idx_opt else{ continue;};
          let cosub_idx = indices[*cosub_idx];

          // Loop over all coorelated particles.
          for particle_idx0 in self.cosubstitution_groups[cosub_idx].iter(){

            // Copy the equivalent particle from the primary cell.
            let mut new_particle = self.bath_particles[*particle_idx0].clone();

            // Move the new particle to the correct cell.
            new_particle.coordinates = &new_particle.coordinates + offset;

            // Save the new particle.
            self.bath_particles.push(new_particle);
        
            // Update metadata.
            self.primary_cell_indices.push(*particle_idx0);
            let particle_idx = self.bath_particles.len() - 1;
            self.cell_indices[*particle_idx0][cell_id] = Some(particle_idx);
          }
        }

      }
    } 

    Ok(())

  }
  //----------------------------------------------------------------------------
  // This function runs through the primary cell and randomly removes
  // particles that the user has set to be randomly removable.
  fn set_primary_cell_voidable_particles(&mut self, 
      rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {

    // Check if there are any configurations.
    let particle_configs = &config.particles;

    // Initialize the HashMap for indices with the same void probability
    // The key is the void probability in parts per trillion.
    let mut property_map = HashMap::<u64, Vec::<usize>>::with_capacity(
        self.cosubstitution_groups.len());

    // Find all the cosubstitution_groups with the same substitution chance.
    for (cosub_idx,cosubstitution_group) in self.cosubstitution_groups.iter()
      .enumerate()
    {
      if cosubstitution_group.is_empty(){ continue; }

      let particle_idx0 = cosubstitution_group[0];
      if !self.bath_particles[particle_idx0].active {continue;}
      
      let mut p_keep = 1.0;

      // Check if this particle has a custom config.
      if let Some(config_id) = self.particle_config_ids[particle_idx0] 
      {

        // Check if this particle has any custom properties.
        if let Some(properties) = &particle_configs[config_id].properties 
        {

          // Get void probability.
          if let Some(p_remove) = properties.isotopic_distribution
            .void_probability{ 
          
              p_keep = 1.0 - p_remove;
          }
        }
      }
      

      // Determine the chance of not void, 
      // and convert it into parts per trillion (as a u64). 
      let ppt = (p_keep*1e12) as u64;

      // Add the indices to the property_map.
      if let Some(indices) = &mut property_map.get_mut(&ppt){
        indices.push(cosub_idx);
      }else{
        let n_particles0 = self.cell_indices.len();
        let mut indices = Vec::<usize>::with_capacity(n_particles0);
        indices.push(cosub_idx);
        property_map.insert(ppt,indices);
      }

    }

    // Loop through all void probability, indices pairs.
    for (ppt, indices) in property_map.iter(){ 

      // Convert the keep probability to a floating point number.
      let p_keep = (*ppt as f64)/1e12; 
      let p_remove = 1.0 - p_keep;
      let p = if p_keep <= p_remove{
        p_keep
      }else{
        p_remove
      };
      
      // Generate a binomial distribution with n_group elements
      // a success chance of p.
      let n_groups = indices.len();
      let Ok(bin) = Binomial::new(n_groups as u64, p) else {
        return Err(CluEError::CannotSampleBinomialDistribution(
              n_groups,p));
      };


      // Sample from the binomial distribution for how many elements to keep. 
      let n_indices = bin.sample(rng) as usize;

      // Sample n_indices out of n_groups elements without repeats to keep.
      let random_indices = sample(rng,n_groups,n_indices);
      let mut selected_indices: Vec::<Option<usize>>;
      if p_remove <= p_keep{
        // random_indices contains indices to remove.
        selected_indices = (0..n_groups).map(|_n| None).collect();
        for idx in random_indices.iter(){
          selected_indices[idx] = Some(idx);
        }
      }else{
        // random_indices contains indices to keep.
        selected_indices = (0..n_groups).map(Some).collect();
        for idx in random_indices.iter(){
          selected_indices[idx] = None;
        }
      };


      // Loop over particles to keep.
      for cosub_idx_opt in selected_indices.iter(){
        let Some(cosub_idx) = cosub_idx_opt else{ continue;};
        let cosub_idx = indices[*cosub_idx];

        // Loop over all coorelated particles.
        for &particle_idx0 in self.cosubstitution_groups[cosub_idx].iter(){
          self.bath_particles[particle_idx0].active = false;
        }
      }

    }

    Ok(())

  }
  //----------------------------------------------------------------------------
  // This function trims the extended system down to the user specified shape. 
  fn trim_system(&mut self, config: &Config) -> Result<(),CluEError>{
    let Some(radius) = config.radius else {
      return Err(CluEError::NoRadius);
    };

    let n_particles = self.bath_particles.len();
    for idx in 0..n_particles{

      let indices: Vec::<usize>;
      let r: &Vector3D;
      if let Some(exchange_group_manager) = &self.exchange_groups{
        match exchange_group_manager.exchange_group_ids[idx]{
          Some(id) 
            => {
              r = exchange_group_manager.exchange_groups[id].centroid();
              indices = exchange_group_manager.exchange_groups[id]
                .indices().clone();
            },
          None => {
            r = &self.bath_particles[idx].coordinates;
            indices = vec![idx];
          },
        }
      }else{
        r = &self.bath_particles[idx].coordinates;
        indices = vec![idx];
      }

      if r.norm() > radius{
        for &index in indices.iter(){
          self.bath_particles[index].active = false;
        }
      }
    }

    Ok(())
  }
  
  //----------------------------------------------------------------------------
  // This function checks for clashes between particles in the primary cell
  // and PBC particles.
  fn check_primary_clashes(&mut self, config: &Config) -> Result<(),CluEError>{

    let Some(clash_distance) = &config.clash_distance else{
      return Ok(());
    };

    for (idx0, particle0) in self.bath_particles.iter().enumerate(){
      if !particle0.active { continue; }
      let cell_id0 = self.cell_id(idx0)?;
      if cell_id0 != 0 { break; } 
      let r0 = &self.bath_particles[idx0].coordinates;

      for (idx1, particle1) in self.bath_particles.iter().enumerate()
        .skip(idx0+1){
        if !particle1.active { continue; }
        let cell_id1 = self.cell_id(idx1)?;
        if cell_id1 != 0 { break; } 
        let r1 = &self.bath_particles[idx1].coordinates;

        let delta_r = (r1 -r0).norm();

        if delta_r < *clash_distance{
          return
            Err(CluEError::ParticlesClash(idx0,particle0.element.to_string(),
                  idx1,particle1.element.to_string(),
                  delta_r/ANGSTROM, clash_distance/ANGSTROM));
        }
      }
    }

    Ok(())
  }
  //----------------------------------------------------------------------------
  // This removes particles in in PBC copies that closer than a user specified
  // distance.
  fn trim_pbc_clashes(&mut self, config: &Config) -> Result<(),CluEError>{

    let Some(clash_distance) = &config.clash_distance_pbc else {
      return Ok(());
    };

    let mut to_remove = (0..self.bath_particles.len()).map(|_| false)
      .collect::<Vec::<bool>>();
    for (idx0, particle0) in self.bath_particles.iter().enumerate(){
      if !particle0.active { continue; }
      let cell_id0 = self.cell_id(idx0)?;
      let r0 = &self.bath_particles[idx0].coordinates;

      for (idx1, particle1) in self.bath_particles.iter().enumerate()
        .skip(idx0){
        if !particle1.active { continue; }
        if to_remove[idx1] { continue; }

        let cell_id1 = self.cell_id(idx1)?;

        if cell_id1 <= cell_id0 { continue; }

        let r1 = &self.bath_particles[idx1].coordinates;

        let delta_r = (r1 -r0).norm();

        to_remove[idx1] = delta_r < *clash_distance;
      } 
    }

    for (idx, remove) in to_remove.iter().enumerate(){
      if *remove{

        self.bath_particles[idx].active = false;
      }
    }

    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function updates the exchange groups to account for the extended
  // system.
  fn update_exhange_groups(&mut self, config: &Config) -> Result<(),CluEError>
  {
    let Some(exchange_group_manager0) = &self.exchange_groups else {
      return Ok(());
    };

    let n_ex = exchange_group_manager0.exchange_groups.len()
      *self.cell_offsets.len();
    
    let mut exchange_groups = Vec::<ExchangeGroup>::with_capacity(n_ex);
    let mut exchange_group_ids = self.bath_particles.iter().map(|_p| None)
      .collect::<Vec::<Option<usize>>>();
    let mut exchange_couplings = Vec::<f64>::with_capacity(n_ex);


    // Loop over unit cells.
    for icell in 0..self.cell_offsets.len(){ 
      for exchange_group_0 in exchange_group_manager0.exchange_groups.iter()
      {
        let indices_0 = exchange_group_0.indices();
        let mut indices = Vec::<usize>::with_capacity(indices_0.len());

        for &idx0 in indices_0.iter(){
          let Some(idx) = self.cell_indices[idx0][icell] else{
            break;
          };
          indices.push(idx);
        }
        if indices.len() != indices_0.len(){
          continue;
        }

        let mut all_same_isotope = true;
        for &idx in indices.iter(){
          if self.bath_particles[indices[0]].isotope 
            != self.bath_particles[idx].isotope{
            all_same_isotope = false;
          }
        }
        if !all_same_isotope{
          continue;
        }

        let center = exchange_group_0.centroid() + &self.cell_offsets[icell];
        let exchange_group: ExchangeGroup = match exchange_group_0{
          ExchangeGroup::Methyl(rotor) 
          => ExchangeGroup::Methyl(C3Rotor{center,normal: rotor.normal.clone(), 
              indices: [indices[0],indices[1],indices[2]]}),
          ExchangeGroup::PrimaryAmonium(rotor)
          => ExchangeGroup::PrimaryAmonium(C3Rotor{ 
              center, normal: rotor.normal.clone(),
              indices: [indices[0],indices[1],indices[2]]}),
        };

        exchange_groups.push(exchange_group);


        for &idx in indices.iter(){
          exchange_group_ids[idx] = Some(exchange_groups.len()-1);
        }


        let exchange_coupling: f64 
          = self.extract_exchange_coupling(indices[0],config);
        exchange_couplings.push(exchange_coupling);

     
      }
    }
    self.exchange_groups = Some(ExchangeGroupManager{
      exchange_groups,
      exchange_group_ids,
      exchange_couplings,
    });
    
    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function draws a random isotopologue from the user defined
  // distribution.
  fn set_isotopologue(&mut self, rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {
    // Initialize rng range.
    let range = Uniform::new(0.0f64, 1.0);

    // Check if there are any configurations.
    let particle_configs = &config.particles;

    for unit_cell_id in 0..self.cell_offsets.len(){

      // Loop over bath particle.
      for cosubstitution_group in self.cosubstitution_groups.iter(){

        // Generate a random number.
        let random_number = range.sample(rng);

        for particle_idx0 in cosubstitution_group.iter(){
          if !self.bath_particles[*particle_idx0].active {continue;}
        
          let Some(particle_idx) 
            = self.cell_indices[*particle_idx0][unit_cell_id] else{continue;};

          let particle = &mut self.bath_particles[particle_idx];
      
          // Check if this particle has a custom config.
          let config_id = match unit_cell_id {
            0 => {
              let Some(cid) = self.particle_config_ids[*particle_idx0] 
                else { continue};
              cid
            },
            _ => {
              let Some(cid) = self.extracell_particle_config_ids[*particle_idx0] 
              else { continue};
              cid
            },
          };

          // Check if this particle has any custom properties.
          let Some(properties) = &particle_configs[config_id].properties else{
            continue;
          };


          let mut cdf = 0.0;
     
          let isotopic_distribution = &properties.isotopic_distribution;
          

          for iso in isotopic_distribution.isotope_abundances.iter(){
            cdf += iso.abundance;
            if cdf >= random_number{
              particle.isotope = iso.isotope;
              particle.active = particle.isotope.spin_multiplicity() > 1;
              break;
            } 
          }

          if let Some(isotope_properties) = properties.isotope_properties
              .get(&particle.isotope.to_string()){
            
            if let Some(active) = isotope_properties.active{
              particle.active = active;
            } 
          }

        }
      }
    }

    Ok(())

  }
  //----------------------------------------------------------------------------
  // This function sets the offset vector between PBC copies.
  fn set_cell_shifts(&mut self, config: &Config) -> Result<(),CluEError>{
    let cell_edges = self.cell_offsets.clone();

    match config.replicate_unit_cell{
      Some(ReplicateUnitCell::No) => {
        self.cell_offsets = vec![Vector3D::zeros()];
        return Ok(());
      },
      Some(_) => (),
      None => {return Err(CluEError::NoApplyPBC) }
    }

    if cell_edges.is_empty(){
      self.cell_offsets = vec![Vector3D::zeros()];
      println!("No unit cell information found. \
Periodic boudary conditions will not be applied.");
      return Ok(());
    }else if cell_edges.len() != 3 {
      return Err(CluEError::IncorrectNumberOfAxes(cell_edges.len(),3));
    }

    self.cell_offsets = build_cell_shifts(cell_edges, config)?;
    Ok(())
  }
}
  //----------------------------------------------------------------------------
  // This function builds the offset vector between PBC copies.
  fn build_cell_shifts(cell_edges: Vec::<Vector3D>, config: &Config) 
    -> Result<Vec::<Vector3D>,CluEError>
  {
    
    let n_cells_per_dim = count_extra_cells_per_dim(&cell_edges,config)?;

    let mut n_cells = 1;
    for ix in 0..3{
      n_cells *= 1 + 2*(n_cells_per_dim[ix] as usize);
    }

    let mut cell_offsets = Vec::<Vector3D>::with_capacity(n_cells);
  
    for ix in -n_cells_per_dim[0]..=n_cells_per_dim[0]{
      for iy in -n_cells_per_dim[1]..=n_cells_per_dim[1]{
        for iz in -n_cells_per_dim[2]..=n_cells_per_dim[2]{

          let offset = &(&cell_edges[0].scale(ix as f64)
            + &cell_edges[1].scale(iy as f64)) 
            + &cell_edges[2].scale(iz as f64);

          cell_offsets.push(offset);
    }}}
    cell_offsets.sort();
    Ok(cell_offsets)
  }
  //----------------------------------------------------------------------------
  fn count_extra_cells_per_dim(cell_edges: &[Vector3D], config: &Config) 
      -> Result<[i32;3],CluEError>
  {

    match config.replicate_unit_cell{
      Some(ReplicateUnitCell::Auto) => (),
      Some(ReplicateUnitCell::Specified(nc)) => {
        return Ok([nc[0] as i32, nc[1] as i32, nc[2] as i32]);
      },
      Some(ReplicateUnitCell::No) => return Ok([0,0,0]),
      None => return Err(CluEError::NoApplyPBC),
    }

    let mut n_cells_per_dim = [1,1,1];

    let Some(radius) = config.radius else {
      return Err(CluEError::NoRadius);
    };


    for ix in 0..3{
      n_cells_per_dim[ix] = std::cmp::max(1,
       math::ceil(radius/cell_edges[ix].norm()) as i32
          );

    }
    Ok(n_cells_per_dim)
  
  }
  //----------------------------------------------------------------------------




#[cfg(test)]
mod tests{
  use super::*;
  use crate::structure::pdb;
  use crate::structure::particle_filter::*;
  use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};
  use crate::config::DetectedSpinCoordinates;

  use crate::config::particle_config::{ParticleConfig,
    ParticleProperties,IsotopeAbundance};
  use crate::elements::Element;
  use crate::isotopes::Isotope;

  use crate::io::FromTOMLString;
  
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_cell_shifts(){
    let cell_edges = vec![
      Vector3D::from([1.0, 0.0, 0.0]),
      Vector3D::from([0.0, 0.5, 0.0]),
      Vector3D::from([0.0, 0.0, 2.0]),
    ];

    let config = Config::from_toml_string(r##"
        replicate_unit_cell = [1,2,1]
        "##).unwrap();

    let cell_shifts = build_cell_shifts(cell_edges, &config).unwrap(); 
    assert_eq!(cell_shifts.len(), 3*5*3);
    assert_eq!(cell_shifts, vec![
      Vector3D::from([0.0, 0.0, 0.0]),
      Vector3D::from([0.0, -0.5, 0.0]),
      Vector3D::from([0.0, 0.5, 0.0]),
      Vector3D::from([-1.0, 0.0, 0.0]),
      Vector3D::from([0.0, -1.0, 0.0]),
      Vector3D::from([0.0, 1.0, 0.0]),
      Vector3D::from([1.0, 0.0, 0.0]),
      Vector3D::from([-1.0, -0.5, 0.0]),
      Vector3D::from([-1.0, 0.5, 0.0]),
      Vector3D::from([1.0, -0.5, 0.0]),
      Vector3D::from([1.0, 0.5, 0.0]),
      Vector3D::from([-1.0, -1.0, 0.0]),
      Vector3D::from([-1.0, 1.0, 0.0]),
      Vector3D::from([1.0, -1.0, 0.0]),
      Vector3D::from([1.0, 1.0, 0.0]),
      Vector3D::from([0.0, 0.0, -2.0]),
      Vector3D::from([0.0, 0.0, 2.0]),
      Vector3D::from([0.0, -0.5, -2.0]),
      Vector3D::from([0.0, -0.5, 2.0]),
      Vector3D::from([0.0, 0.5, -2.0]),
      Vector3D::from([0.0, 0.5, 2.0]),
      Vector3D::from([-1.0, 0.0, -2.0]),
      Vector3D::from([-1.0, 0.0, 2.0]),
      Vector3D::from([0.0, -1.0, -2.0]),
      Vector3D::from([0.0, -1.0, 2.0]),
      Vector3D::from([0.0, 1.0, -2.0]),
      Vector3D::from([0.0, 1.0, 2.0]),
      Vector3D::from([1.0, 0.0, -2.0]),
      Vector3D::from([1.0, 0.0, 2.0]),
      Vector3D::from([-1.0, -0.5, -2.0]),
      Vector3D::from([-1.0, -0.5, 2.0]),
      Vector3D::from([-1.0, 0.5, -2.0]),
      Vector3D::from([-1.0, 0.5, 2.0]),
      Vector3D::from([1.0, -0.5, -2.0]),
      Vector3D::from([1.0, -0.5, 2.0]),
      Vector3D::from([1.0, 0.5, -2.0]),
      Vector3D::from([1.0, 0.5, 2.0]),
      Vector3D::from([-1.0, -1.0, -2.0]),
      Vector3D::from([-1.0, -1.0, 2.0]),
      Vector3D::from([-1.0, 1.0, -2.0]),
      Vector3D::from([-1.0, 1.0, 2.0]),
      Vector3D::from([1.0, -1.0, -2.0]),
      Vector3D::from([1.0, -1.0, 2.0]),
      Vector3D::from([1.0, 1.0, -2.0]),
      Vector3D::from([1.0, 1.0, 2.0]),
    ]);

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_extended_structure_tempo_wat_gly_7nm(){
    // n    : Molecules    : Hydrons 
    // 1    : TEMPO        :    18
    // 1500 : glycerols    : 12000
    // 7469 : waters       : 14938
    //
    // total hydrons =  26956.
    //
    let n_wat = 7469;
    let n_gly = 1500;
    let n_ex = 2*n_wat + 3*n_gly; 
    let n_nx = 5*n_gly + 18;

    let filename = "./assets/TEMPO_wat_gly_70A.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();


    let mut particle_configs =vec![
      ParticleConfig::new("exchangeable_hydrogens".to_string()),
      ParticleConfig::new("nonexchangeable_hydrogens".to_string()),
    ];

    let mut filter_ex = ParticleFilter::new();
    filter_ex.elements = vec![Element::Hydrogen];
    filter_ex.bonded_elements = vec![Element::Oxygen];
    particle_configs[0].filter = Some(filter_ex.clone());

    let mut filter_nx = ParticleFilter::new();
    filter_nx.elements = vec![Element::Hydrogen];
    filter_nx.not_bonded_elements = vec![Element::Oxygen];
    particle_configs[1].filter = Some(filter_nx.clone());
    
    let mut properties = ParticleProperties::new();
    properties.isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen1, abundance: 0.5});
    properties.isotopic_distribution.isotope_abundances.push(
        IsotopeAbundance{isotope: Isotope::Hydrogen2, abundance: 0.5});


    particle_configs[0].properties = Some(properties.clone());

    properties.cosubstitute = Some(SecondaryParticleFilter::SameMolecule);
    particle_configs[1].properties = Some(properties);
    
    let mut config = Config::new();
    config.particles = particle_configs;
    config.radius = Some(73.5676e-10);
    config.replicate_unit_cell = Some(ReplicateUnitCell::Auto);
    let n_uc = 125;
    let r_nitrogen=Vector3D::from([36.440*1e-10, 36.900*1e-10,  37.100*1e-10]);
    let r_oxygen = Vector3D::from([35.290*1e-10, 36.430*1e-10, 37.810*1e-10]);
    let r_e = (&r_nitrogen + &r_oxygen).scale(0.5);
    config.detected_spin_position = Some( 
        DetectedSpinCoordinates::XYZ(r_e) );


    config.set_defaults().unwrap();

    let mut rng =  ChaCha20Rng::from_entropy();
    structure.build_primary_structure(&mut rng,&config).unwrap();

    structure.build_extended_structure(&mut rng, &config).unwrap();

    assert_eq!(structure.cell_offsets.len(),n_uc);



    // Check for for the correct exchanges.
    let [n_h, n_h1, n_h2,n_h1_h1, n_h1_h2, n_h2_h1, n_h2_h2,_n_mol] =
      get_conditional_hydron_stats(&structure,&filter_ex);

    assert_eq!(n_h,n_uc*n_ex);
    let approx_eq = |m: usize,n: usize| -> bool{
      (1.0 - (m as f64)/(n as f64) ).abs() <1e3 
    };
    assert!(n_h1_h1>0);
    assert!(n_h2_h2>0);
    assert!(n_h2_h1>0);
    assert!(n_h1_h2>0);
    assert!( approx_eq(n_h1,n_h2) );
    assert!( approx_eq(n_h1_h1,n_h2_h2) );
    assert!( approx_eq(n_h1_h1,n_h2_h1) );
    assert!( approx_eq(n_h1_h1,n_h1_h2) );

    let [n_h, n_h1, n_h2, n_h1_h1, n_h1_h2, n_h2_h1, n_h2_h2,_n_mol] =
      get_conditional_hydron_stats(&structure,&filter_nx);

    assert_eq!(n_h,n_uc*n_nx);
    assert!(n_h1_h1>0);
    assert!(n_h2_h2>0);
    assert_eq!(n_h2_h1,0);
    assert_eq!(n_h1_h2,0);
    assert!( approx_eq(n_h1,n_h2) );
    assert!( approx_eq(n_h1_h1,n_h2_h2) );


   
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  fn get_conditional_hydron_stats(structure: &Structure,filter: &ParticleFilter)
    -> [usize;8]
  {

    let mut n_h = 0;
    let mut n_h1 = 0;
    let mut n_h2 = 0;
    let mut n_h1_h1 = 0;
    let mut n_h1_h2 = 0;
    let mut n_h2_h1 = 0;
    let mut n_h2_h2 = 0;
    let mut n_mol = 0;

    let mut is_ref_proton = true;
    let mut ref_mol_id = 1000000;
    let mut ref_cell_id = 1000000;

    let indices = filter.filter(structure);

    for idx in indices{
      let particle = &structure.bath_particles[idx];
      if (*particle).element == Element::Hydrogen{
        n_h += 1;
      }else{
        assert!((*particle).element != Element::Oxygen);
      }
      let mol_id = structure.molecule_id(idx);
      let idx0 = structure.primary_cell_indices[idx];
      let cell_indices = &structure.cell_indices[idx0];
      let cell_id = cell_indices.iter().position(|&x| {
            match x{ Some(id) => id == idx, None => false,}
          }).unwrap();
      let is_proton = (*particle).isotope == Isotope::Hydrogen1;
      let is_deuteron = (*particle).isotope == Isotope::Hydrogen2;
      
      if is_proton{
        n_h1 += 1;
      }else if is_deuteron{
        n_h2 += 1;
      }

      if !is_proton && !is_deuteron{ continue; }

      if mol_id != ref_mol_id || cell_id != ref_cell_id{
        ref_mol_id = mol_id;
        ref_cell_id = cell_id;
        is_ref_proton = is_proton;
        n_mol += 1;
        continue;
      }

      match (is_proton,is_ref_proton){
        (true,true) => n_h1_h1 +=1, 
        (true,false) => n_h1_h2 +=1, 
        (false,true) => n_h2_h1 +=1, 
        (false,false) => n_h2_h2 +=1, 
      }

    }
    [n_h, n_h1, n_h2, n_h1_h1, n_h1_h2, n_h2_h1, n_h2_h2, n_mol]
  }

  //----------------------------------------------------------------------------
  #[test]
  fn test_build_extended_structure_tempo_3gly_1npr(){

    let config = Config::from_toml_string(r##"
      input_structure_file = "./assets/TEMPO_3gly_1npr_50A.pdb"
      detected_spin.position = [28,29]
      radius = 73.5676
      replicate_unit_cell = true
      
      pulse_sequence = "hahn"
      tau_increments = [1e-2]
      number_timepoints = [101]
      
      [[groups]]
      name = "exchangeable_hydrogens"
      selection.elements = ["H"]
      selection.bonded_elements = ["O"]

      1H.abundance = 0.5
      2H.abundance = 0.5
      
      [[groups]]
      name = "nonexchangeable_hydrogens"
      selection.elements = ["H"]
      selection.not_bonded_elements = ["O"]
      
      1H.abundance = 0.5
      2H.abundance = 0.5
      
      cosubstitute = "same_molecule"
        "##).unwrap();
    
    let mut rng = ChaCha20Rng::from_entropy();
    let structure = Structure::build_structure(&mut rng,&config).unwrap();
    
    let n_uc = 125;

    // Test: set_cell_shifts().
    assert_eq!(structure.cell_offsets.len(),n_uc);

    // Test: set_isotopologue().
    let mut group_gly_oh = ParticleFilter::new();
    group_gly_oh.elements = vec![Element::Hydrogen];
    group_gly_oh.residues = vec!["GLY".to_string()];
    group_gly_oh.bonded_elements = vec![Element::Oxygen];
    let indices_gly_oh = group_gly_oh.filter(&structure);

    let mut group_gly_ch = ParticleFilter::new();
    group_gly_ch.elements = vec![Element::Hydrogen];
    group_gly_ch.residues = vec!["GLY".to_string()];
    group_gly_ch.bonded_elements = vec![Element::Carbon];
    let indices_gly_ch = group_gly_ch.filter(&structure);

    let mut group_npr_ch = ParticleFilter::new();
    group_npr_ch.residues = vec!["NPR".to_string()];
    group_npr_ch.elements = vec![Element::Hydrogen];
    group_npr_ch.bonded_elements = vec![Element::Carbon];
    let indices_npr_ch = group_npr_ch.filter(&structure);

    let mut group_npr_oh = ParticleFilter::new();
    group_npr_oh.residues = vec!["NPR".to_string()];
    group_npr_oh.elements = vec![Element::Hydrogen];
    group_npr_oh.bonded_elements = vec![Element::Oxygen];
    let indices_npr_oh = group_npr_oh.filter(&structure);

    let mut group_tempo_h = ParticleFilter::new();
    group_tempo_h.residues = vec!["TEM".to_string()];
    group_tempo_h.elements = vec![Element::Hydrogen];
    let indices_tempo_h = group_tempo_h.filter(&structure);


    assert_eq!(indices_gly_oh.len(), n_uc*775*3);
    assert_eq!(indices_gly_ch.len(), n_uc*775*5);
    assert_eq!(indices_npr_oh.len(), n_uc*251*1);
    assert_eq!(indices_npr_ch.len(), n_uc*251*7);
    assert_eq!(indices_tempo_h.len(), n_uc*18);

    let assert_corr_approx_eq = |left: &[f64],right: &[f64], name: &str|{
      assert_eq!(left.len(),right.len());

      let tol = 0.05;
      for (ii,&l) in left.iter().enumerate(){
        let r = right[ii];
        if (l-r).abs() >= tol{
          panic!("{}: MAD(left({:?}), right({:?})) > {}.", name,left,right,tol);
        } 
      }
    };

    let corr_gly_oh = get_1h_2h_correlations(&structure,&indices_gly_oh);
    assert_corr_approx_eq(&corr_gly_oh,&vec![0.5,0.5],"gly_oh");

    let corr_gly_ch = get_1h_2h_correlations(&structure,&indices_gly_ch);
    assert_corr_approx_eq(&corr_gly_ch,&vec![0.9,0.5],"gly_ch");

    let corr_npr_oh = get_1h_2h_correlations(&structure,&indices_npr_oh);
    assert_corr_approx_eq(&corr_npr_oh,&vec![0.5,0.5],"npr_oh");

    let corr_npr_ch = get_1h_2h_correlations(&structure,&indices_npr_ch);
    assert_corr_approx_eq(&corr_npr_ch,&vec![0.93,0.5],"npr_h");

  }
  //----------------------------------------------------------------------------
  fn get_1h_2h_correlations(structure: &Structure, indices: &[usize])
    -> Vec::<f64>
  {
    let mut corr_1 = 0.0;
    let mut corr_9 = 0.0;

    for (ii,idx) in indices.iter().enumerate().skip(1){
      let idx0 = indices[ii-1];
      let particle0 = &structure.bath_particles[idx0];
      let particle = &structure.bath_particles[*idx];
      if (*particle).isotope == (*particle0).isotope{
        corr_1 += 1.0;
      }

      if ii < 9 { continue; }
      
      let idx0 = indices[ii-9];
      let particle0 = &structure.bath_particles[idx0];
      let particle = &structure.bath_particles[*idx];
      if (*particle).isotope == (*particle0).isotope{
        corr_9 += 1.0;
      }
    }
    corr_1 /= indices.len() as f64 - 1.0;
    corr_9 /= indices.len() as f64 - 9.0;

    vec![corr_1, corr_9]
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_add_voidable_particles(){
    let config0 = Config::from_toml_string(r##"
      input_structure_file = "./assets/TEMPO_3gly_1npr_50A.pdb"
      detected_spin.position = [28,29]
      radius = 70.0
      replicate_unit_cell = true
      
      pulse_sequence = "hahn"
      tau_increments = [1e-2]
      number_timepoints = [101]
      
      [[groups]]
      name = "h"
      selection.elements = ["H"]

      [[groups]]
      name = "non_h"
      selection.not_elements = ["H"]
      drop_probability = 1.0

    "##).unwrap();

    let mut rng = ChaCha20Rng::from_entropy();
    let structure0 = Structure::build_structure(&mut rng,&config0).unwrap();

    let config1 = Config::from_toml_string(r##"
      input_structure_file = "./assets/TEMPO_3gly_1npr_50A.pdb"
      detected_spin.position = [28,29]
      radius = 70.0
      replicate_unit_cell = true
      
      pulse_sequence = "hahn"
      tau_increments = [1e-2]
      number_timepoints = [101]
      
      [[groups]]
      name = "h"
      selection.elements = ["H"]
      drop_probability = 0.5

      [[groups]]
      name = "non_h"
      selection.not_elements = ["H"]
      drop_probability = 1.0


    "##).unwrap();

    let structure1 = Structure::build_structure(&mut rng,&config1).unwrap();


    let n1 = structure1.number_active() as f64;
    let n0 = structure0.number_active() as f64;
    let std = (0.25*n0).sqrt();
    let n_min = 0.5*n0 - 5.0*std;
    let n_max = 0.5*n0 + 5.0*std;
    assert!(n1 > n_min);
    assert!(n1 < n_max);

  }
  //----------------------------------------------------------------------------

}
