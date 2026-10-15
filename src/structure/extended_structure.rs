use crate::clue_errors::CluEError;
use crate::config::Config;
use crate::config::ReplicateUnitCell;
use crate::config::particle_config::ParticleConfig;
use crate::math;  
use crate::physical_constants::ANGSTROM;
use crate::space_3d::Vector3D;
use crate::structure::{Structure, exchange_groups::*, reciprocal_vectors};
use crate::structure::particle::Particle;

use rand::seq::index::sample;
use rand_chacha::ChaCha20Rng;
use rand_distr::{Binomial, Distribution, Uniform};
use std::collections::{BTreeMap, HashMap};

impl Structure{
  /// This function takes a structure from build_primary_structure(), applies
  /// periodic boundary condition (if applicable), and sets isotope identities.
  pub fn build_extended_structure(&mut self,
      rng: &mut ChaCha20Rng, config: &Config) 
    -> Result<(),CluEError>
  {

    // Construct offset vectors for each pbc application.
    self.set_cell_shifts(config)?;

    // Assign groups and cosubstitution groups outside the primary cell.
    self.set_extracell_groups(config)?;

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
  // This function pairs the copies of the primary cell particles in the other
  // unit cells with `config.extracell_particles`, and finds their 
  // cosubstitution groups.
  // If no extracell filter depends on the unit cell, one assignment is 
  // shared by all unit cells; otherwise, each unit cell is assigned 
  // separately, since for example distances from the detected spin differ.
  // Unit cells with the same assignment share it, and since the 
  // cosubstitution groups only depend on the assignment, they share 
  // those too.
  fn set_extracell_groups(&mut self, config: &Config) -> Result<(),CluEError>
  {
    let particle_configs = &config.extracell_particles;

    let depends_on_cell = particle_configs.iter().any(|p_cfg| {
      match &p_cfg.filter{
        Some(filter) => filter.depends_on_cell(),
        None => false,
      }
    });

    let n_cells = self.cell_offsets.len();

    if !depends_on_cell || n_cells == 1{
      let config_ids = self.pair_particle_configs(particle_configs,None)?;
      let cosubstitution_groups = self.build_cosubstitution_groups(
          &config_ids, particle_configs, None)?;

      self.extracell_particle_config_ids = vec![config_ids];
      self.extracell_cosubstitution_groups = vec![cosubstitution_groups];
      self.extracell_group_indices = vec![0; n_cells];
      return Ok(());
    }

    // Map from each distinct assignment to its index.
    let mut group_index_map = HashMap::<Vec::<Option<usize>>,usize>::new();
    let mut extracell_cosubstitution_groups = Vec::<Vec::<Vec::<usize>>>::new();

    // The primary cell uses `particle_config_ids` and 
    // `cosubstitution_groups` instead, so its index is not used.
    let mut extracell_group_indices = Vec::<usize>::with_capacity(n_cells);
    extracell_group_indices.push(0);

    for cell in 1..n_cells{
      let config_ids = self.pair_particle_configs(particle_configs,
          Some(cell))?;

      if let Some(&group_idx) = group_index_map.get(&config_ids){
        extracell_group_indices.push(group_idx);
        continue;
      }

      let cosubstitution_groups = self.build_cosubstitution_groups(
          &config_ids, particle_configs, Some(cell))?;
      let group_idx = extracell_cosubstitution_groups.len();
      extracell_cosubstitution_groups.push(cosubstitution_groups);
      group_index_map.insert(config_ids, group_idx);
      extracell_group_indices.push(group_idx);
    }

    let mut extracell_particle_config_ids 
      = vec![Vec::<Option<usize>>::new(); group_index_map.len()];
    for (config_ids, group_idx) in group_index_map.into_iter(){
      extracell_particle_config_ids[group_idx] = config_ids;
    }

    self.extracell_particle_config_ids = extracell_particle_config_ids;
    self.extracell_cosubstitution_groups = extracell_cosubstitution_groups;
    self.extracell_group_indices = extracell_group_indices;

    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function returns a map from the chance to keep a cosubstitution 
  // group, in parts per trillion, to the indices of the groups 
  // in `cosubstitution_groups` with that chance.
  // The map is ordered, so that random numbers are drawn in the same order
  // for the same seed.
  fn build_keep_probability_map(&self, 
      cosubstitution_groups: &[Vec::<usize>],
      config_ids: &[Option<usize>],
      particle_configs: &[ParticleConfig])
    -> BTreeMap::<u64, Vec::<usize>>
  {
    let n_particles0 = self.cell_indices.len();

    // Initialize the map for indices with the same void probability
    // The key is the void probability in parts per trillion.
    let mut property_map = BTreeMap::<u64, Vec::<usize>>::new();

    // Find all the cosubstitution_groups with the same substitution chance.
    for (cosub_idx,cosubstitution_group) in cosubstitution_groups.iter()
      .enumerate()
    {
      if cosubstitution_group.is_empty(){ continue; }

      // Skip groups in which no particle can have a spin 
      // under `particle_configs`.
      if !self.group_can_have_spin(cosubstitution_group, config_ids,
          particle_configs) {continue;}

      // All particles in a cosubstitution group are in the same group.
      let particle_idx0 = cosubstitution_group[0];
      
      let mut p_keep = 1.0;

      // Check if this particle has a custom config.
      if let Some(config_id) = config_ids[particle_idx0] 
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

    property_map
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
    self.cell_ids.reserve( n_to_reserve );

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

    // Unit cells that share their groups also share the property map.
    let property_maps = (0..self.extracell_cosubstitution_groups.len())
      .map(|group_idx| self.build_keep_probability_map(
            &self.extracell_cosubstitution_groups[group_idx],
            &self.extracell_particle_config_ids[group_idx],
            particle_configs)
      ).collect::<Vec::<BTreeMap::<u64, Vec::<usize>>>>();

    // Copies that cannot be within the radius are not stored
    // (see `copy_can_be_active()`).
    let Some(radius) = config.radius else {
      return Err(CluEError::NoRadius);
    };
    let delta = self.max_distance_from_exchange_group_center();
    let store_all = store_all_copies();

    // Loop through each unit cell other than the primary cell.
    for cell_id in 1..self.cell_offsets.len(){
      let offset = self.cell_offsets[cell_id].clone();

      // Get the groups of this unit cell.
      let group_idx = self.extracell_group_index(cell_id);
      let property_map = &property_maps[group_idx];

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
          for particle_idx0 
            in self.extracell_cosubstitution_groups[group_idx][cosub_idx]
            .iter()
          {

            // Spinless particles are not copied.
            if self.max_spin_multiplicity_for_particle(*particle_idx0,
                &self.extracell_particle_config_ids[group_idx],
                particle_configs) < 2 { continue; }

            // Copies that can never be active are not stored.
            // The random numbers above are drawn per cosubstitution group,
            // so skipping a copy does not change them.
            if !store_all && !self.copy_can_be_active(*particle_idx0, &offset,
                radius, delta) { continue; }

            // Copy the equivalent particle from the primary cell.
            let mut new_particle = self.bath_particles[*particle_idx0].clone();

            // Move the new particle to the correct cell.
            new_particle.coordinates = &new_particle.coordinates + &offset;

            // The copy can have a spin from its own group.
            new_particle.active = true;

            // Save the new particle.
            self.bath_particles.push(new_particle);
        
            // Update metadata.
            self.primary_cell_indices.push(*particle_idx0);
            self.cell_ids.push(cell_id);
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

    // Initialize the map for indices with the same void probability
    // The key is the void probability in parts per trillion.
    // The map is ordered, so that random numbers are drawn in the same order
    // for the same seed.
    let mut property_map = BTreeMap::<u64, Vec::<usize>>::new();

    // Find all the cosubstitution_groups with the same substitution chance.
    for (cosub_idx,cosubstitution_group) in self.cosubstitution_groups.iter()
      .enumerate()
    {
      if cosubstitution_group.is_empty(){ continue; }

      // Skip groups in which no particle can have a spin. 
      // This is decided from the groups rather than from `active`, 
      // so that it does not depend on the sampled isotopes.
      if !self.group_can_have_spin(cosubstitution_group, 
          &self.particle_config_ids, particle_configs) {continue;}

      // All particles in a cosubstitution group are in the same group.
      let particle_idx0 = cosubstitution_group[0];
      
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
    
    // Loop through particles by index.
    for idx in 0..n_particles{

      let indices: Vec::<usize>;
      let r: &Vector3D;

      // Check if the particle is part of an exchange group.
      if let Some(exchange_group_manager) = &self.exchange_groups{
        
        // Get center of group.
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
        // The particle is a singleton 
        // and the center is the particle coordinates
        r = &self.bath_particles[idx].coordinates;
        indices = vec![idx];
      }

      if r.norm() > radius{
        for &index in indices.iter(){
          self.bath_particles[index].active = false;
        }
      }
    }
    for particle in self.bath_particles.iter_mut(){
      if particle.isotope.spin_multiplicity() < 2{
        particle.active = false;
      }
    }
    Ok(())
  }
  
  //----------------------------------------------------------------------------
  // This function checks for clashes between active particles in the 
  // primary cell, and returns an error for the first clashing pair 
  // (idx0, idx1), ordered by idx0 and then by idx1.
  fn check_primary_clashes(&mut self, config: &Config) -> Result<(),CluEError>{

    let Some(clash_distance) = config.clash_distance else{
      return Ok(());
    };
    if !(clash_distance > 0.0){
      return Ok(());
    }

    let n_particles0 = self.cell_indices.len();
    let indices = (0..n_particles0)
      .filter(|&idx| self.bath_particles[idx].active)
      .collect::<Vec::<usize>>();
    let bins = bin_particles(&self.bath_particles, &indices, clash_distance);

    let mut first_clash: Option<(usize,usize,f64)> = None;
    for &idx0 in indices.iter(){
      let r0 = &self.bath_particles[idx0].coordinates;
      for idx1 in neighbor_bin_particles(&bins, r0, clash_distance){
        if idx1 <= idx0 { continue; }
        let delta_r = (&self.bath_particles[idx1].coordinates - r0).norm();
        if delta_r < clash_distance{
          let is_first = match first_clash{
            Some((i0,i1,_)) => (idx0,idx1) < (i0,i1),
            None => true,
          };
          if is_first{
            first_clash = Some((idx0,idx1,delta_r));
          }
        }
      }
    }

    if let Some((idx0,idx1,delta_r)) = first_clash{
      return Err(CluEError::ParticlesClash(
            idx0, self.bath_particles[idx0].element.to_string(),
            idx1, self.bath_particles[idx1].element.to_string(),
            delta_r/ANGSTROM, clash_distance/ANGSTROM));
    }

    Ok(())
  }
  //----------------------------------------------------------------------------
  // This removes particles in PBC copies that are closer than a user 
  // specified distance: an active particle is made inactive if an active 
  // particle in a unit cell with a lower ID is closer than the clash 
  // distance.
  fn trim_pbc_clashes(&mut self, config: &Config) -> Result<(),CluEError>{

    match config.replicate_unit_cell{
      Some(ReplicateUnitCell::No) => return Ok(()),
      None => return Ok(()),
      _ => (),
    }

    let Some(clash_distance) = config.clash_distance_pbc else {
      return Ok(());
    };
    let clash_distance = clash_distance.abs();
    if clash_distance == 0.0 || !clash_distance.is_finite(){
      return Ok(());
    }
    let clash_distance_squared = clash_distance*clash_distance;

    let indices = (0..self.bath_particles.len())
      .filter(|&idx| self.bath_particles[idx].active)
      .collect::<Vec::<usize>>();
    let bins = bin_particles(&self.bath_particles, &indices, clash_distance);

    let mut to_remove = vec![false; self.bath_particles.len()];

    for &idx1 in indices.iter(){
      let cell_id1 = self.cell_id(idx1)?;
      let r1 = &self.bath_particles[idx1].coordinates;

      for idx0 in neighbor_bin_particles(&bins, r1, clash_distance){
        if self.cell_id(idx0)? >= cell_id1 { continue; }

        let delta_r = r1 - &self.bath_particles[idx0].coordinates;
        if delta_r.dot(&delta_r) < clash_distance_squared{
          to_remove[idx1] = true;
          break;
        }
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

        // A group is only copied if all of its members are stored;
        // members are dropped, spinless, or beyond the radius 
        // (`copy_can_be_active()`) together.
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

    for unit_cell_id in 0..self.cell_offsets.len(){

      // Get the cosubstitution groups of this unit cell.
      let group_idx = self.extracell_group_index(unit_cell_id);
      let cosubstitution_groups = match unit_cell_id{
        0 => &self.cosubstitution_groups,
        _ => &self.extracell_cosubstitution_groups[group_idx],
      };

      // Loop over bath particle.
      for cosubstitution_group in cosubstitution_groups.iter(){

        // Generate a random number.
        let random_number = range.sample(rng);

        // Loop through particle that should be changed together.
        for particle_idx0 in cosubstitution_group.iter(){
        
          let Some(particle_idx) 
            = self.cell_indices[*particle_idx0][unit_cell_id] else{continue;};

          let particle = &mut self.bath_particles[particle_idx];
      
          // Check if this particle has a custom config.
          let particle_config = match unit_cell_id {
            0 => {
              let Some(cid) = self.particle_config_ids[*particle_idx0] 
                else { continue};
              &config.particles[cid]
            },
            _ => {
              let Some(cid) 
                = self.extracell_particle_config_ids[group_idx][*particle_idx0]
              else { continue};
              &config.extracell_particles[cid]
            },
          };

          // Check if this particle has any custom properties.
          let Some(properties) = &particle_config.properties else{
            continue;
          };


          let mut cdf = 0.0;
     
          let isotopic_distribution = &properties.isotopic_distribution;
          
          // Loop through all potential substitutes.
          // The abundances sum to one (see `Config::check_particle_configs`).
          // For example, if 1H has a 50% abundance, 2H 40%, and 3H 10%,
          // and random_number = 0.95, the group will be set as 3H,
          // since cdf = 0.5 < 0.95 on the first pass, cdf = 0.9 < 0.95 on 
          // the second pass, and cdf = 1.0 >= 0.95 on the third pass.
          // If floating point rounding leaves cdf below random_number
          // after the last isotope, the last isotope is used.
          // Groups without abundances keep the isotope from the structure 
          // file.
          let abundances = &isotopic_distribution.isotope_abundances;
          let isotope = abundances.iter().find(|iso| {
            cdf += iso.abundance;
            cdf >= random_number
          }).or(abundances.last()).map(|iso| iso.isotope);

          if let Some(isotope) = isotope{
            particle.isotope = isotope;
            particle.active = particle.isotope.spin_multiplicity() > 1;
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

    // Only lattice translations are applied, so other symmetry operations
    // of the space group would be missing.
    if let Some(space_group) = &self.space_group{
      if space_group.split_whitespace().collect::<String>() != "P1"{
        return Err(CluEError::SpaceGroupNotP1(space_group.to_string()));
      }
    }

    let Some(reciprocal_vectors) = reciprocal_vectors(&cell_edges) else{
      self.cell_offsets = vec![Vector3D::zeros()];
      println!("The unit cell has no volume. \
Periodic boudary conditions will not be applied.");
      return Ok(());
    };

    let n_cells_per_dim = self.count_extra_cells_per_dim(&cell_edges,
        &reciprocal_vectors, config)?;

    self.cell_offsets = build_cell_shifts(cell_edges, n_cells_per_dim);
    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function returns the number of unit cells to add on each side of 
  // the primary cell along each cell edge.
  // It should be called with the detected spin at the origin.
  fn count_extra_cells_per_dim(&self, cell_edges: &[Vector3D], 
      reciprocal_vectors: &[Vector3D;3], config: &Config) 
      -> Result<[i32;3],CluEError>
  {
    let positions = self.bath_particles.iter()
      .map(|particle| particle.coordinates.clone())
      .collect::<Vec::<Vector3D>>();

    // Exchange groups are kept if their center is within the radius,
    // so their particles can be slightly further away.
    let delta = self.max_distance_from_exchange_group_center();

    match config.replicate_unit_cell{
      Some(ReplicateUnitCell::Auto) => {
        let Some(radius) = config.radius else {
          return Err(CluEError::NoRadius);
        };
        Ok(auto_extra_cells_per_dim(cell_edges, reciprocal_vectors,
              &positions, radius, delta))
      },
      Some(ReplicateUnitCell::Specified(nc)) => {
        let n_cells_per_dim = [nc[0] as i32, nc[1] as i32, nc[2] as i32];
        if let Some(radius) = config.radius{
          let n_required = required_extra_cells_per_dim(reciprocal_vectors,
              &positions, radius + delta);
          if (0..3).any(|ii| n_cells_per_dim[ii] < n_required[ii]){
            println!("replicate_unit_cell = [{}, {}, {}] is less than the \
[{}, {}, {}] unit cells needed to fill the radius.",
                n_cells_per_dim[0], n_cells_per_dim[1], n_cells_per_dim[2],
                n_required[0], n_required[1], n_required[2]);
          }
        }
        Ok(n_cells_per_dim)
      },
      Some(ReplicateUnitCell::No) => Ok([0,0,0]),
      None => Err(CluEError::NoApplyPBC),
    }
  }
  //----------------------------------------------------------------------------
  // This function returns the largest distance between a particle in an 
  // exchange group and the center of the exchange group.
  fn max_distance_from_exchange_group_center(&self) -> f64{
    let Some(exchange_group_manager) = &self.exchange_groups else{
      return 0.0;
    };
    let mut max_distance = 0.0;
    for exchange_group in exchange_group_manager.exchange_groups.iter(){
      for idx in exchange_group.indices(){
        let r = &self.bath_particles[idx].coordinates 
          - exchange_group.centroid();
        max_distance = f64::max(max_distance, r.norm());
      }
    }
    max_distance
  }
  //----------------------------------------------------------------------------
  // This function returns `false` if `trim_system()` would deactivate the
  // copy of primary cell particle `idx0` displaced by `offset`, whatever
  // its isotope, and `true` otherwise.
  // Coordinates are relative to the detected spin.
  // `trim_system()` measures a particle in an exchange group from the
  // centroid of the group, and any other particle from its own position.
  // A periodic copy of an exchange group only forms if all of its members
  // are stored and have the same isotope (`update_exhange_groups()`),
  // so the members of an exchange group are kept or skipped together:
  // they are kept if the centroid is within `radius + delta`, where `delta`
  // is the largest distance of a member from the centroid of its group
  // (`max_distance_from_exchange_group_center()`).
  // This keeps every group that is trimmed by its centroid, and every member
  // within `radius` of a group that does not form because of mixed isotopes.
  // A small margin makes sure rounding errors can only keep more copies.
  // Other particles use the same coordinates and comparison as
  // `trim_system()`.
  fn copy_can_be_active(&self, idx0: usize, offset: &Vector3D, radius: f64,
      delta: f64) -> bool
  {
    if let Some(exchange_group_manager) = &self.exchange_groups{
      if let Some(id) = exchange_group_manager.exchange_group_ids[idx0]{
        let r = exchange_group_manager.exchange_groups[id].centroid()
          + offset;
        return r.norm() <= radius + delta + 1e-6*ANGSTROM;
      }
    }

    let r = &self.bath_particles[idx0].coordinates + offset;
    r.norm() <= radius
  }
}
//------------------------------------------------------------------------------
#[cfg(test)]
thread_local!{
  // Tests can set this to store every periodic copy, as `extend_structure()`
  // did before it skipped copies that cannot be within the radius.
  static STORE_ALL_COPIES: std::cell::Cell<bool>
    = const { std::cell::Cell::new(false) };
}
//------------------------------------------------------------------------------
// This function returns `true` if `extend_structure()` should store every
// periodic copy; this is only possible in tests.
fn store_all_copies() -> bool{
  #[cfg(test)]
  { STORE_ALL_COPIES.with(|store_all| store_all.get()) }
  #[cfg(not(test))]
  { false }
}
//------------------------------------------------------------------------------
// This function sets whether `extend_structure()` stores every periodic copy
// in the calling thread.
#[cfg(test)]
pub(crate) fn set_store_all_copies(store_all: bool){
  STORE_ALL_COPIES.with(|store_all_copies| store_all_copies.set(store_all));
}
  //----------------------------------------------------------------------------
  // This function returns the cubic bin, with edge `bin_size`, 
  // that contains `r`.
  fn bin_of(r: &Vector3D, bin_size: f64) -> [i64;3]{
    [ (r.x()/bin_size).floor() as i64,
      (r.y()/bin_size).floor() as i64,
      (r.z()/bin_size).floor() as i64 ]
  }
  //----------------------------------------------------------------------------
  // This function sorts the particles with indices `indices` into cubic bins
  // with edge `bin_size`, so that particles closer than `bin_size` are 
  // in the same or neighboring bins.
  fn bin_particles(particles: &[Particle], indices: &[usize], bin_size: f64)
    -> HashMap::<[i64;3],Vec::<usize>>
  {
    let mut bins = HashMap::<[i64;3],Vec::<usize>>::new();
    for &idx in indices.iter(){
      bins.entry(bin_of(&particles[idx].coordinates, bin_size))
        .or_default().push(idx);
    }
    bins
  }
  //----------------------------------------------------------------------------
  // This function returns the indices of the particles in `bins` that are
  // in the bin containing `r` or in its 26 neighboring bins.
  fn neighbor_bin_particles<'a>(bins: &'a HashMap::<[i64;3],Vec::<usize>>,
      r: &Vector3D, bin_size: f64) -> impl Iterator<Item = usize> + 'a
  {
    let [ix,iy,iz] = bin_of(r, bin_size);
    (-1..=1).flat_map(move |dx| (-1..=1).flat_map(move |dy| 
        (-1..=1).map(move |dz| [ix + dx, iy + dy, iz + dz])))
      .filter_map(move |bin| bins.get(&bin))
      .flat_map(|indices| indices.iter().copied())
  }
  //----------------------------------------------------------------------------
  // This function builds the offset vector between PBC copies.
  fn build_cell_shifts(cell_edges: Vec::<Vector3D>, n_cells_per_dim: [i32;3]) 
    -> Vec::<Vector3D>
  {
    let mut n_cells = 1;
    for &n in n_cells_per_dim.iter(){  
      n_cells *= 1 + 2*(n as usize);
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
    cell_offsets
  }
  //----------------------------------------------------------------------------
  // This function returns the number of unit cells needed on each side of 
  // the primary cell, so that every periodic image within `radius` of the
  // origin of every particle at `positions` is included.
  // With reciprocal vectors g_i, the image r + sum_j k_j a_j has the 
  // fractional coordinate g_i·r + k_i along a_i, and 
  // |g_i·r + k_i| <= |g_i| radius if the image is within the radius,
  // so |k_i| <= floor(|g_i| radius + |g_i·r|).
  fn required_extra_cells_per_dim(reciprocal_vectors: &[Vector3D;3],
      positions: &[Vector3D], radius: f64) -> [i32;3]
  {
    let mut n_cells_per_dim = [0,0,0];
    for (ii, g) in reciprocal_vectors.iter().enumerate(){
      let f_max = positions.iter().map(|r| g.dot(r).abs())
        .fold(0.0, f64::max);
      n_cells_per_dim[ii] = (g.norm()*radius + f_max).floor() as i32;
    }
    n_cells_per_dim
  }
  //----------------------------------------------------------------------------
  // This function returns the number of unit cells on each side of 
  // the primary cell for replicate_unit_cell = "auto".
  // This is the number required to fill the radius, `radius` + `delta`, 
  // but at least one and at least ceil(radius/|a_i|), as in earlier versions.
  fn auto_extra_cells_per_dim(cell_edges: &[Vector3D], 
      reciprocal_vectors: &[Vector3D;3], positions: &[Vector3D], 
      radius: f64, delta: f64) -> [i32;3]
  {
    let n_required = required_extra_cells_per_dim(reciprocal_vectors,
        positions, radius + delta);

    let mut n_cells_per_dim = [1,1,1];
    for ii in 0..3{
      let n_edge = math::ceil(radius/cell_edges[ii].norm()) as i32;
      n_cells_per_dim[ii] = *[1, n_edge, n_required[ii]].iter().max()
        .unwrap();
    }
    n_cells_per_dim
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
    ParticleProperties,IsotopeAbundance,Cosubstitute};
  use crate::elements::Element;
  use crate::isotopes::Isotope;

  use crate::io::FromTOMLString;
  use crate::cluster::adjacency::AdjacencyList;
  use crate::structure::particle::Particle;
  use std::collections::HashMap;
  
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

    let structure = Structure::new(Vec::<Particle>::new(), 
        AdjacencyList::with_capacity(0), cell_edges.clone());
    let reciprocal_vectors = reciprocal_vectors(&cell_edges).unwrap();
    let n_cells_per_dim = structure.count_extra_cells_per_dim(&cell_edges,
        &reciprocal_vectors, &config).unwrap();
    assert_eq!(n_cells_per_dim,[1,2,1]);

    let cell_shifts = build_cell_shifts(cell_edges, n_cells_per_dim); 
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
    // This test checks copies in every unit cell, also beyond the radius.
    set_store_all_copies(true);
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

    properties.cosubstitute = Some(Cosubstitute::SameMolecule);
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
    // This test checks copies in every unit cell, also beyond the radius.
    set_store_all_copies(true);

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
  #[test]
  fn test_build_extended_structure_solvent_box_isotopes(){
    // This test checks copies in every unit cell, also beyond the radius.
    set_store_all_copies(true);
    // 3224 waters and 250 glycerols in a 50.146 A cubic box.
    let mut config = Config::from_toml_file("./assets/solvent_box.toml")
      .unwrap();
    config.input_structure_file
      = Some("./assets/water_gly_20pc_vol.pdb".to_string());
    // The toml position is in angstroms, so put the spin at the box center.
    config.detected_spin_position = Some(DetectedSpinCoordinates::XYZ(
          Vector3D::from([25.0, 25.0, 25.0]).scale(ANGSTROM)));
    let radius = config.radius.unwrap();

    let mut rng = ChaCha20Rng::seed_from_u64(20261012);
    let structure = Structure::build_structure(&mut rng,&config).unwrap();

    let n_cells = structure.cell_offsets.len();
    assert_eq!(n_cells,125);

    let label_of = |idx0: usize| -> Option<&str>{
      structure.particle_config_ids[idx0]
        .map(|id| config.particles[id].label.as_str())
    };
    let extracell_label_of = |idx0: usize, cell: usize| -> Option<&str>{
      let group_idx = structure.extracell_group_index(cell);
      structure.extracell_particle_config_ids[group_idx][idx0]
        .map(|id| config.extracell_particles[id].label.as_str())
    };
    let p_keep_of = |label: &str| -> f64{
      let (_id,p_cfg) = config.find_particle_config(label).unwrap();
      let p_drop = p_cfg.properties.as_ref().unwrap()
        .isotopic_distribution.void_probability.unwrap_or(0.0);
      1.0 - p_drop
    };

    // Group sanity.
    let n_particles0 = structure.cell_indices.len();
    let mut group_indices = HashMap::<String,Vec::<usize>>::new();
    for idx0 in 0..n_particles0{
      let particle = &structure.bath_particles[idx0];
      let label = label_of(idx0);
      if particle.element == Element::Hydrogen{
        assert!(label.is_some(), "hydrogen {} has no group", idx0);
      }
      // Distance selections are measured from the detected spin.
      let r = particle.coordinates.norm();
      match label{
        Some("CLOSE_H") => assert!(r <= 4e-10),
        Some(_) => assert!(r >= 4e-10),
        None => (),
      }
      if let Some(label) = label{
        group_indices.entry(label.to_string()).or_default().push(idx0);
      }
    }
    assert!(!group_indices["CLOSE_H"].is_empty());
    // The groups of each cell are checked in 
    // `test_set_extracell_groups_shared_by_cells()`.

    // Outside the primary cell, no copy is within 4 A of the detected spin,
    // so each hydrogen copy is in a solvent group.
    for cell in 1..n_cells{
      for idx0 in 0..n_particles0{
        let label = extracell_label_of(idx0,cell);
        assert_ne!(label,Some("CLOSE_H"));
        let particle = &structure.bath_particles[idx0];
        if particle.element != Element::Hydrogen{ continue; }
        let expected = match particle.residue.as_deref(){
          Some("SOL") => vec!["WATER_H"],
          _ => vec!["GLYCEROL_OH","GLYCEROL_CH"],
        };
        assert!(expected.contains(&label.unwrap()),
            "hydrogen {} in cell {} has group {:?}",idx0,cell,label);
      }
      // Glycerol C-H copies cosubstitute as whole molecules.
      let group_idx = structure.extracell_group_index(cell);
      for group in structure.extracell_cosubstitution_groups[group_idx].iter(){
        if group.is_empty() 
          || extracell_label_of(group[0],cell) != Some("GLYCEROL_CH"){
          continue;
        }
        assert_eq!(group.len(),5);
      }
    }

    let water_h = &group_indices["WATER_H"];
    let gly_oh = &group_indices["GLYCEROL_OH"];
    let gly_ch = &group_indices["GLYCEROL_CH"];
    assert!(water_h.len() > 6400);
    assert!(gly_oh.len() > 740);
    assert!(gly_ch.len() > 1240);

    let gly_ch_groups = structure.cosubstitution_groups.iter()
      .filter(|g| !g.is_empty() && label_of(g[0]) == Some("GLYCEROL_CH"))
      .collect::<Vec::<&Vec::<usize>>>();
    let mut n_full_groups = 0;
    for group in gly_ch_groups.iter(){
      assert!(group.len() <= 5);
      if group.len() == 5{ n_full_groups += 1; }
      let mol_id = structure.molecule_ids[group[0]];
      for &idx0 in group.iter(){
        assert_eq!(label_of(idx0),Some("GLYCEROL_CH"));
        assert_eq!(structure.molecule_ids[idx0],mol_id);
      }
    }
    assert!(n_full_groups >= 245);
    assert_eq!(gly_ch_groups.iter().map(|g| g.len()).sum::<usize>(),
        gly_ch.len());

    // Isotopes: every copy of every group hydrogen is 1H.
    // (CLOSE_H has no element selection, so it can include non-hydrogens.)
    for indices in group_indices.values(){
      for &idx0 in indices.iter(){
        if structure.bath_particles[idx0].element != Element::Hydrogen{
          continue;
        }
        for idx in structure.cell_indices[idx0].iter().flatten(){
          assert_eq!(structure.bath_particles[*idx].isotope,
              Isotope::Hydrogen1);
        }
      }
    }

    // Drop ratios per cell and over all cells.
    let count_kept = |indices: &[usize], cell: usize| -> (usize,usize){
      let mut k = 0;
      let mut n = 0;
      for &idx0 in indices.iter(){
        let Some(kept) = is_kept(&structure,idx0,cell,radius) else{
          continue;
        };
        n += 1;
        if kept { k += 1; }
      }
      (k,n)
    };
    for (label,indices) in [("WATER_H",water_h), ("GLYCEROL_OH",gly_oh)]{
      let p = p_keep_of(label);
      let (mut k_tot, mut n_tot) = (0,0);
      for cell in 0..n_cells{
        let (k,n) = count_kept(indices,cell);
        assert_binomial(k,n,p,&format!("{}, cell {}",label,cell));
        k_tot += k;
        n_tot += n;
      }
      assert_binomial(k_tot,n_tot,p,&format!("{}, all cells",label));
    }

    // Glycerol C-H: every C-H on a glycerol is the same in each cell.
    let p = p_keep_of("GLYCEROL_CH");
    let (mut k_tot, mut n_tot) = (0,0);
    for cell in 0..n_cells{
      let (mut k, mut n) = (0,0);
      for group in gly_ch_groups.iter(){
        let states = group.iter()
          .filter_map(|&idx0| is_kept(&structure,idx0,cell,radius))
          .collect::<Vec::<bool>>();
        if states.is_empty(){ continue; }
        assert!(states.iter().all(|&s| s == states[0]),
            "glycerol C-H {:?} differ in cell {}: {:?}",group,cell,states);
        let isotopes = group.iter()
          .filter_map(|&idx0| structure.cell_indices[idx0][cell])
          .map(|idx| structure.bath_particles[idx].isotope)
          .collect::<Vec::<Isotope>>();
        assert!(isotopes.iter().all(|&iso| iso == isotopes[0]));
        n += 1;
        if states[0] { k += 1; }
      }
      assert_binomial(k,n,p,&format!("GLYCEROL_CH, cell {}",cell));
      k_tot += k;
      n_tot += n;
    }
    assert!(k_tot > 0);
    assert_binomial(k_tot,n_tot,p,"GLYCEROL_CH, all cells");

    // No cosubstitution: hydrons on the same molecule are independent.
    for (label,indices) in [("WATER_H",water_h), ("GLYCEROL_OH",gly_oh)]{
      let p = p_keep_of(label);
      let mut by_molecule = HashMap::<usize,Vec::<usize>>::new();
      for &idx0 in indices.iter(){
        by_molecule.entry(structure.molecule_ids[idx0]).or_default()
          .push(idx0);
      }
      let mut n_multiple = 0;
      let mut expected = 0.0;
      for cell in 0..n_cells{
        for hydrons in by_molecule.values(){
          let (k,n) = count_kept(hydrons,cell);
          if k >= 2 { n_multiple += 1; }
          let nf = n as f64;
          expected += nf*(nf-1.0)/2.0*p*p;
        }
      }
      assert!((n_multiple as f64) < 4.0*expected + 5.0,
          "{}: {} molecules with multiple kept hydrons, expected {}.",
          label,n_multiple,expected);
    }

    // CLOSE_H is never dropped and is inactive in the primary cell.
    let close_h = &group_indices["CLOSE_H"];
    for &idx0 in close_h.iter(){
      assert!(!structure.bath_particles[idx0].active);
    }
    // Copies of CLOSE_H hydrogens in other cells follow their solvent group:
    // they are dropped at the group rate and are not set inactive.
    let close_h_hydrogens = close_h.iter()
      .filter(|&&idx0| structure.bath_particles[idx0].element 
          == Element::Hydrogen)
      .copied().collect::<Vec::<usize>>();
    assert!(!close_h_hydrogens.is_empty());
    let (mut k_tot, mut n_tot) = (0,0);
    let mut p_ave = 0.0;
    for cell in 1..n_cells{
      for &idx0 in close_h_hydrogens.iter(){
        let label = extracell_label_of(idx0,cell).unwrap();
        p_ave += p_keep_of(label);
        n_tot += 1;
        if let Some(idx) = structure.cell_indices[idx0][cell]{
          k_tot += 1;
          if structure.bath_particles[idx].coordinates.norm() <= radius{
            assert!(structure.bath_particles[idx].active);
          }
        }
      }
    }
    assert_binomial(k_tot,n_tot,p_ave/(n_tot as f64),"CLOSE_H copies");
  }
  //----------------------------------------------------------------------------
  // This function returns whether the copy of primary cell particle `idx0`
  // in cell `cell` survived random dropping, or `None` if this cannot be
  // determined because the primary cell particle was trimmed.
  fn is_kept(structure: &Structure, idx0: usize, cell: usize, radius: f64)
    -> Option<bool>
  {
    if cell == 0{
      let particle = &structure.bath_particles[idx0];
      if particle.coordinates.norm() > radius{
        return None;
      }
      return Some(particle.active);
    }
    Some(structure.cell_indices[idx0][cell].is_some())
  }
  //----------------------------------------------------------------------------
  fn assert_binomial(k: usize, n: usize, p: f64, label: &str){
    let mean = (n as f64)*p;
    let std = (mean*(1.0-p)).sqrt();
    assert!(((k as f64) - mean).abs() <= 5.0*std + 1.0,
        "{}: kept {} of {}, expected {} +/- {}.",label,k,n,mean,std);
  }
  //----------------------------------------------------------------------------
  // This function builds TEMPO with periodic copies of the unit cell 
  // and the specified groups.
  fn build_periodic_tempo(groups: &str, seed: u64) -> Structure{
    let config = Config::from_toml_string(&format!(r##"
        input_structure_file = "assets/TEMPO.pdb"
        radius = 20
        detected_spin.position = [28,29]
        number_timepoints = [101]
        tau_increments = [1e-1]
        cluster_method = "cce"
        max_cluster_size = 2
        magnetic_field = 1.2
        replicate_unit_cell = true
        {}
        "##,groups)).unwrap();

    let mut rng = ChaCha20Rng::seed_from_u64(seed);
    let structure = Structure::build_structure(&mut rng,&config).unwrap();
    assert!(structure.cell_offsets.len() > 1);
    structure
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_structure_is_reproducible_with_drop_probabilities(){
    // With several drop probabilities, the random draws must happen 
    // in the same order on every build with the same seed.
    let groups = r##"
        [[groups]]
        name = "methyl_h"
        selection.serials = [3,4,5,7,8,9,21,22,23,25,26,27]
        drop_probability = 0.3

        [[groups]]
        name = "other_h"
        selection.elements = ["H"]
        selection.not_serials = [3,4,5,7,8,9,21,22,23,25,26,27]
        drop_probability = 0.6

        [[groups]]
        name = "nitrogen"
        selection.elements = ["N"]
        drop_probability = 0.5
        "##;

    let reference = build_periodic_tempo(groups,7);
    for _ii in 0..8{
      let structure = build_periodic_tempo(groups,7);
      assert_eq!(structure.number(),reference.number());
      assert_eq!(structure.cell_indices,reference.cell_indices);
      for (p, p_ref) in structure.bath_particles.iter()
        .zip(reference.bath_particles.iter())
      {
        assert_eq!(p.coordinates,p_ref.coordinates);
        assert_eq!(p.isotope,p_ref.isotope);
        assert_eq!(p.active,p_ref.active);
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_structure_extracell_only_group(){
    // This test checks copies in every unit cell, also beyond the radius.
    set_store_all_copies(true);
    // A group that only applies outside the primary cell must be used
    // regardless of the order of the groups.
    let c_primary = r##"
        [[groups]]
        name = "c_primary"
        selection.elements = ["C"]
        selection.cell_type = "primary_cell"
        "##;
    let h_primary = r##"
        [[groups]]
        name = "h_primary"
        selection.elements = ["H"]
        selection.cell_type = "primary_cell"
        "##;
    let c_extra = r##"
        [[groups]]
        name = "c_extra"
        selection.elements = ["C"]
        selection.cell_type = "extracells"
        13C.abundance = 1.0
        "##;
    let n_c0 = 9;

    for groups in [format!("{}{}{}",h_primary,c_primary,c_extra),
                   format!("{}{}",c_extra,c_primary)]
    {
      let structure = build_periodic_tempo(&groups,0);
      let n_cells = structure.cell_offsets.len();

      let mut n_c = [0,0];
      let mut n_13c = [0,0];
      for (idx,particle) in structure.bath_particles.iter().enumerate(){
        if particle.element != Element::Carbon{ continue; }
        let k = if structure.cell_id(idx) == Ok(0){ 0 }else{ 1 };
        n_c[k] += 1;
        if particle.isotope == Isotope::Carbon13{
          n_13c[k] += 1;
        }else{
          assert_eq!(particle.isotope,Isotope::Carbon12);
          assert!(!particle.active);
        }
      }
      assert_eq!(n_c,[n_c0, n_c0*(n_cells-1)]);
      assert_eq!(n_13c,[0, n_c0*(n_cells-1)]);
    }
  }
  //----------------------------------------------------------------------------
  // This function builds TEMPO with the unit cell `cell_edges`, and checks 
  // that every periodic image of every hydrogen within the radius exists.
  fn assert_hydrogen_images_within_radius_exist(cell_edges: Vec::<Vector3D>,
      toml: &str)
  {
    let config = Config::from_toml_string(toml).unwrap();
    let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    structure.cell_offsets = cell_edges.clone();
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    structure.build_primary_structure(&mut rng,&config).unwrap();
    structure.build_extended_structure(&mut rng,&config).unwrap();
    let radius = config.radius.unwrap();

    let n_particles0 = structure.cell_indices.len();
    let mut n_checked = 0;
    for idx0 in 0..n_particles0{
      if structure.bath_particles[idx0].element != Element::Hydrogen{ 
        continue;
      }
      let r0 = structure.bath_particles[idx0].coordinates.clone();
      for ii in -8..=8{ for jj in -8..=8{ for kk in -8..=8{
        let offset = &(&cell_edges[0].scale(ii as f64)
            + &cell_edges[1].scale(jj as f64)) 
            + &cell_edges[2].scale(kk as f64);
        if (&r0 + &offset).norm() > radius - 1e-3*ANGSTROM{ continue; }
        n_checked += 1;
        let cell = structure.cell_offsets.iter().position(|r| *r == offset);
        let Some(cell) = cell else{
          panic!("no unit cell for offset {} of hydrogen {}",offset,idx0);
        };
        assert!(structure.cell_indices[idx0][cell].is_some());
      }}}
    }
    assert!(n_checked > n_particles0);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_auto_replication_covers_radius_oblique_cell(){
    // The cell thickness normal to the a-b plane is 10 A, 
    // half the edge length.
    let gamma = 150.0f64.to_radians();
    let cell_edges = vec![
      Vector3D::from([20.0, 0.0, 0.0]).scale(ANGSTROM),
      Vector3D::from([20.0*gamma.cos(), 20.0*gamma.sin(), 0.0]).scale(ANGSTROM),
      Vector3D::from([0.0, 0.0, 25.0]).scale(ANGSTROM),
    ];
    assert_hydrogen_images_within_radius_exist(cell_edges, r##"
        input_structure_file = "assets/TEMPO.pdb"
        radius = 38
        detected_spin.position = [28,29]
        replicate_unit_cell = true
        clash_distance_pbc = 0.0
        "##);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_auto_replication_covers_radius_distant_atoms(){
    // The detected spin is 1.5 unit cells from the molecule.
    let cell_edges = vec![
      Vector3D::from([20.0, 0.0, 0.0]).scale(ANGSTROM),
      Vector3D::from([0.0, 20.0, 0.0]).scale(ANGSTROM),
      Vector3D::from([0.0, 0.0, 20.0]).scale(ANGSTROM),
    ];
    assert_hydrogen_images_within_radius_exist(cell_edges, r##"
        input_structure_file = "assets/TEMPO.pdb"
        radius = 15
        detected_spin.position = [66.44, 36.90, 37.10]
        replicate_unit_cell = true
        clash_distance_pbc = 0.0
        "##);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_auto_extra_cells_per_dim_cubic_cells_unchanged(){
    // For cubic cells with the detected spin at the center, the number of 
    // unit cells is the same as in earlier versions, max(1,ceil(R/L)).
    // The cases are those of the structures in assets/.
    let delta = 1.1*ANGSTROM;
    for (radius, length) in [(73.5676, 72.5676), (60.0, 50.146), 
        (20.0, 72.5676), (70.0, 49.3802), (80.0, 72.5676)]
    {
      let (radius, length) = (radius*ANGSTROM, length*ANGSTROM);
      let cell_edges = vec![
        Vector3D::from([length, 0.0, 0.0]),
        Vector3D::from([0.0, length, 0.0]),
        Vector3D::from([0.0, 0.0, length]),
      ];
      let reciprocal_vectors = reciprocal_vectors(&cell_edges).unwrap();
      let mut positions = Vec::<Vector3D>::new();
      for x in [-0.499, 0.499]{ for y in [-0.499, 0.499]{ for z in [-0.499, 0.499]{
        positions.push(Vector3D::from([x, y, z]).scale(length));
      }}}
      let n = std::cmp::max(1, math::ceil(radius/length) as i32);
      assert_eq!(auto_extra_cells_per_dim(&cell_edges, &reciprocal_vectors,
            &positions, radius, delta), [n,n,n]);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_set_cell_shifts_space_group(){
    // Only lattice translations are applied, so other space groups 
    // cannot be replicated.
    for (replicate, space_group, expected) in [
      ("true", Some("P 21 21 21"), 
       Err(CluEError::SpaceGroupNotP1("P 21 21 21".to_string()))),
      ("false", Some("P 21 21 21"), Ok(())),
      ("true", Some("P 1"), Ok(())),
      ("true", Some("P1"), Ok(())),
      ("true", None, Ok(())),
    ]{
      let config = Config::from_toml_string(&format!(r##"
          input_structure_file = "assets/TEMPO.pdb"
          radius = 20
          detected_spin.position = [28,29]
          replicate_unit_cell = {}
          "##,replicate)).unwrap();
      let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
      assert_eq!(structure.space_group, Some("P 1".to_string()));
      structure.space_group = space_group.map(|s| s.to_string());
      let mut rng = ChaCha20Rng::seed_from_u64(0);
      structure.build_primary_structure(&mut rng,&config).unwrap();
      let result = structure.build_extended_structure(&mut rng,&config);
      assert_eq!(result,expected);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_unit_cell_smaller_than_structure(){
    // A unit cell that is much smaller than the structure does not 
    // describe it, so it is not used.
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    config.set_defaults().unwrap();

    let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    let length = 2.0*ANGSTROM;
    structure.cell_offsets = vec![
      Vector3D::from([length, 0.0, 0.0]),
      Vector3D::from([0.0, length, 0.0]),
      Vector3D::from([0.0, 0.0, length]),
    ];
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    let result = structure.build_primary_structure(&mut rng, &config);
    assert!(matches!(result, 
          Err(CluEError::StructureLargerThanCell(_,span)) if span > 2.0),
        "{:?}",result);
  }
  //----------------------------------------------------------------------------
  // This function builds TEMPO in a cubic unit cell with edge `length` 
  // and the specified settings.
  fn build_tempo_in_cell(length: f64, settings: &str) 
    -> Result<Structure,CluEError>
  {
    let config = Config::from_toml_string(&format!(r##"
        input_structure_file = "assets/TEMPO.pdb"
        radius = 20
        detected_spin.position = [28,29]
        replicate_unit_cell = true
        {}
        "##,settings)).unwrap();
    let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    structure.cell_offsets = vec![
      Vector3D::from([length, 0.0, 0.0]),
      Vector3D::from([0.0, length, 0.0]),
      Vector3D::from([0.0, 0.0, length]),
    ];
    let mut rng = ChaCha20Rng::seed_from_u64(5);
    structure.build_primary_structure(&mut rng,&config)?;
    structure.build_extended_structure(&mut rng,&config)?;
    Ok(structure)
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_cell_id(){
    let structure = build_tempo_in_cell(10.0*ANGSTROM,"").unwrap();
    assert!(structure.cell_offsets.len() > 1);
    for idx in 0..structure.number(){
      let idx0 = structure.primary_cell_indices[idx];
      let expected = structure.cell_indices[idx0].iter()
        .position(|&x| x == Some(idx)).unwrap();
      assert_eq!(structure.cell_id(idx), Ok(expected));
    }
    assert!(structure.cell_id(structure.number()).is_err());
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_check_primary_clashes(){
    // The reported clash is the first pair (idx0, idx1) of active primary 
    // cell particles closer than the clash distance.
    let clash_distance = 2.0*ANGSTROM;
    let reference = build_tempo_in_cell(10.0*ANGSTROM,"").unwrap();
    let n_particles0 = reference.cell_indices.len();
    let mut expected = None;
    'outer: for idx0 in 0..n_particles0{
      if !reference.bath_particles[idx0].active{ continue; }
      for idx1 in (idx0+1)..n_particles0{
        if !reference.bath_particles[idx1].active{ continue; }
        let r = (&reference.bath_particles[idx1].coordinates
            - &reference.bath_particles[idx0].coordinates).norm();
        if r < clash_distance{
          expected = Some(CluEError::ParticlesClash(
              idx0, reference.bath_particles[idx0].element.to_string(),
              idx1, reference.bath_particles[idx1].element.to_string(),
              r/ANGSTROM, clash_distance/ANGSTROM));
          break 'outer;
        }
      }
    }
    let expected = expected.unwrap();

    let result = build_tempo_in_cell(10.0*ANGSTROM,"clash_distance = 2.0");
    assert_eq!(result.err(), Some(expected));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_trim_pbc_clashes(){
    // An active particle is removed if an active particle in a unit cell 
    // with a lower ID is closer than the clash distance.
    let clash_distance = 4.0*ANGSTROM;
    let reference = build_tempo_in_cell(10.0*ANGSTROM,"").unwrap();
    let structure = build_tempo_in_cell(10.0*ANGSTROM,
        "clash_distance_pbc = 4.0").unwrap();
    assert_eq!(structure.number(), reference.number());

    let mut n_removed = 0;
    for idx1 in 0..reference.number(){
      let particle1 = &reference.bath_particles[idx1];
      let cell_id1 = reference.cell_id(idx1).unwrap();
      let is_clashing = particle1.active 
        && (0..reference.number()).any(|idx0| {
          let particle0 = &reference.bath_particles[idx0];
          particle0.active && reference.cell_id(idx0).unwrap() < cell_id1
            && (&particle1.coordinates - &particle0.coordinates).norm() 
            < clash_distance
        });
      if is_clashing{ n_removed += 1; }
      assert_eq!(structure.bath_particles[idx1].active,
          particle1.active && !is_clashing);
    }
    assert!(n_removed > 0);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_set_extracell_groups_shared_by_cells(){
    // Each unit cell has the groups and cosubstitution groups of its own,
    // and unit cells with the same groups share them.
    let mut config = Config::from_toml_file("./assets/solvent_box.toml")
      .unwrap();
    config.input_structure_file
      = Some("./assets/water_gly_20pc_vol.pdb".to_string());
    config.detected_spin_position = Some(DetectedSpinCoordinates::XYZ(
          Vector3D::from([25.0, 25.0, 25.0]).scale(ANGSTROM)));
    let mut rng = ChaCha20Rng::seed_from_u64(20261012);
    let structure = Structure::build_structure(&mut rng,&config).unwrap();
    let n_cells = structure.cell_offsets.len();
    assert_eq!(n_cells,125);

    let particle_configs = &config.extracell_particles;
    for cell in 1..n_cells{
      let group_idx = structure.extracell_group_index(cell);
      let config_ids = structure.pair_particle_configs(particle_configs,
          Some(cell)).unwrap();
      let cosubstitution_groups = structure.build_cosubstitution_groups(
          &config_ids, particle_configs, Some(cell)).unwrap();
      assert_eq!(structure.extracell_particle_config_ids[group_idx],
          config_ids);
      assert_eq!(structure.extracell_cosubstitution_groups[group_idx],
          cosubstitution_groups);
    }

    // No copy is within 4 A of the detected spin, 
    // so all other unit cells have the same groups.
    let group_indices = (1..n_cells)
      .map(|cell| structure.extracell_group_index(cell))
      .collect::<std::collections::HashSet::<usize>>();
    assert_eq!(group_indices.len(),1);
    assert_eq!(structure.extracell_particle_config_ids.len(),1);

    // Every cosubstitution group has particles, and every particle is in 
    // exactly one group.
    let n_particles0 = structure.cell_indices.len();
    for cosubstitution_groups in structure.extracell_cosubstitution_groups
      .iter()
    {
      assert!(cosubstitution_groups.iter().all(|group| !group.is_empty()));
      let mut members = cosubstitution_groups.concat();
      members.sort();
      assert_eq!(members, (0..n_particles0).collect::<Vec::<usize>>());
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_isotopes_drawn_from_listed_abundances(){
    // Every particle in a group with abundances gets a listed isotope.
    let build = |groups: &str, seed: u64| -> Structure{
      let config = Config::from_toml_string(&format!(r##"
          input_structure_file = "assets/TEMPO.pdb"
          radius = 20
          detected_spin.position = [28,29]
          replicate_unit_cell = false
          {}
          "##, groups)).unwrap();
      let mut rng = ChaCha20Rng::seed_from_u64(seed);
      Structure::build_structure(&mut rng,&config).unwrap()
    };

    let (mut n_h1, mut n_h2) = (0, 0);
    for seed in 0..200{
      let structure = build(r##"
          [[groups]]
          name = "h"
          selection.elements = ["H"]
          1H.abundance = 0.5
          2H.abundance = 0.5
          "##, seed);
      for particle in structure.bath_particles.iter()
        .filter(|p| p.element == Element::Hydrogen)
      {
        match particle.isotope{
          Isotope::Hydrogen1 => n_h1 += 1,
          Isotope::Hydrogen2 => n_h2 += 1,
          isotope => panic!("unexpected isotope {}", isotope),
        }
      }
    }
    assert_binomial(n_h2, n_h1 + n_h2, 0.5, "2H");

    let structure = build(r##"
        [[groups]]
        name = "c"
        selection.elements = ["C"]
        13C.abundance = 1.0
        "##, 0);
    let carbons = structure.bath_particles.iter()
      .filter(|p| p.element == Element::Carbon).collect::<Vec::<_>>();
    assert_eq!(carbons.len(), 9);
    assert!(carbons.iter().all(|p| p.isotope == Isotope::Carbon13));

    // Without abundances, the most common isotope is kept.
    let structure = build(r##"
        [[groups]]
        name = "h"
        selection.elements = ["H"]
        "##, 0);
    assert!(structure.bath_particles.iter()
        .filter(|p| p.element == Element::Hydrogen)
        .all(|p| p.isotope == Isotope::Hydrogen1));
  }
  //----------------------------------------------------------------------------
  // This function counts the hydrogens and carbons of TEMPO in each unit cell.
  fn count_h_and_c_per_cell(structure: &Structure) -> Vec::<[usize;2]>{
    let mut counts = vec![[0,0]; structure.cell_offsets.len()];
    for (idx,particle) in structure.bath_particles.iter().enumerate(){
      let cell = structure.cell_id(idx).unwrap();
      match particle.element{
        Element::Hydrogen => counts[cell][0] += 1,
        Element::Carbon => counts[cell][1] += 1,
        _ => (),
      }
    }
    counts
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_build_structure_same_molecule_group_with_spinless_first_member(){
    // This test checks copies in every unit cell, also beyond the radius.
    set_store_all_copies(true);
    // The first particle of TEMPO is a carbon without spin; the molecule
    // must still be copied, but only the particles that can have a spin.
    let structure = build_periodic_tempo(r##"
        [[groups]]
        name = "tempo"
        selection.residues = ["TEM"]
        cosubstitute = "same_molecule"
        "##,0);
    assert_eq!(structure.bath_particles[0].element, Element::Carbon);

    let counts = count_h_and_c_per_cell(&structure);
    assert_eq!(counts[0],[18,9]);
    for count in counts.iter().skip(1){
      assert_eq!(*count,[18,0]);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_drop_same_molecule_group_with_spinless_first_member(){
    // This test checks copies in every unit cell, also beyond the radius.
    set_store_all_copies(true);
    // Whole molecules are dropped, also when the first particle 
    // of the molecule cannot have a spin.
    let structure = build_periodic_tempo(r##"
        [[groups]]
        name = "tempo"
        selection.residues = ["TEM"]
        cosubstitute = "same_molecule"
        drop_probability = 0.5
        "##,3);

    let counts = count_h_and_c_per_cell(&structure);
    let mut n_kept = 0;
    for count in counts.iter().skip(1){
      assert!(*count == [18,0] || *count == [0,0], "{:?}",count);
      if count[0] == 18 { n_kept += 1; }
    }
    assert!(n_kept > 0 && n_kept < counts.len()-1);
    assert_binomial(n_kept, counts.len()-1, 0.5, "TEMPO copies");

    // In the primary cell, all hydrogens are kept or dropped together.
    let n_active_h = structure.bath_particles.iter().enumerate()
      .filter(|(idx,p)| structure.cell_id(*idx) == Ok(0)
          && p.element == Element::Hydrogen && p.active).count();
    assert!(n_active_h == 18 || n_active_h == 0);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_drop_waters_in_primary_cell(){
    // The first particle of each water is an oxygen without spin; 
    // waters in the primary cell must still be dropped as whole molecules.
    let config = Config::from_toml_string(r##"
        input_structure_file = "./assets/water_gly_20pc_vol.pdb"
        detected_spin.position = [25.0, 25.0, 25.0]
        radius = 60
        replicate_unit_cell = false
        number_timepoints = [101]
        tau_increments = [1e-1]
        magnetic_field = 1.2

        [[groups]]
        name = "water"
        selection.residues = ["SOL"]
        cosubstitute = "same_molecule"
        drop_probability = 0.5
        "##).unwrap();
    let mut rng = ChaCha20Rng::seed_from_u64(11);
    let structure = Structure::build_structure(&mut rng,&config).unwrap();
    let radius = config.radius.unwrap();

    let mut kept_by_molecule = HashMap::<usize,Vec::<bool>>::new();
    for (idx0,particle) in structure.bath_particles.iter().enumerate(){
      if particle.element != Element::Hydrogen 
        || particle.residue.as_deref() != Some("SOL"){ continue; }
      let kept = is_kept(&structure,idx0,0,radius).unwrap();
      kept_by_molecule.entry(structure.molecule_id(idx0)).or_default()
        .push(kept);
    }
    let n_waters = kept_by_molecule.len();
    assert_eq!(n_waters,3224);

    let mut n_kept = 0;
    for kept in kept_by_molecule.values(){
      assert_eq!(kept.len(),2);
      assert_eq!(kept[0],kept[1]);
      if kept[0]{ n_kept += 1; }
    }
    assert_binomial(n_kept, n_waters, 0.5, "waters");
  }
  //----------------------------------------------------------------------------
  // This function lists the active particles in bath order, and the exchange
  // groups that have an active member, by primary cell index and unit cell,
  // so that structures with different bath indices can be compared.
  fn active_set(structure: &Structure)
    -> (Vec::<(usize,usize,Isotope,Vector3D)>,
        Vec::<(Vec::<(usize,usize)>,Vector3D,f64)>)
  {
    let label = |idx: usize|
      (structure.primary_cell_indices[idx], structure.cell_id(idx).unwrap());

    let particles = structure.bath_particles.iter().enumerate()
      .filter(|(_idx,p)| p.active)
      .map(|(idx,p)| {
        let (idx0,cell) = label(idx);
        (idx0, cell, p.isotope, p.coordinates.clone())
      }).collect::<Vec::<_>>();

    let mut exchange_groups = Vec::<(Vec::<(usize,usize)>,Vector3D,f64)>::new();
    if let Some(manager) = &structure.exchange_groups{
      for (id,group) in manager.exchange_groups.iter().enumerate(){
        let indices = group.indices();
        if !indices.iter().any(|&idx| structure.bath_particles[idx].active){
          continue;
        }
        exchange_groups.push((indices.iter().map(|&idx| label(idx)).collect(),
            group.centroid().clone(), manager.exchange_couplings[id]));
      }
    }
    (particles, exchange_groups)
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_skipped_copies_do_not_change_active_set(){
    // Skipping periodic copies that cannot be within the radius must not
    // change the active particles, their order, or the exchange groups.
    // The radii cut through the methyls of the nearest copies of TEMPO
    // (72.57 Å away), and the random hydrogen isotopes give methyls
    // with mixed isotopes, which do not form exchange groups.
    let groups = r##"
        [[groups]]
        name = "methyl_h"
        selection.serials = [3,4,5,7,8,9,21,22,23,25,26,27]
        1H.abundance = 0.5
        2H.abundance = 0.5

        [[groups]]
        name = "other_h"
        selection.elements = ["H"]
        selection.not_serials = [3,4,5,7,8,9,21,22,23,25,26,27]
        drop_probability = 0.3
        "##;

    let mut n_cut_groups = 0;
    let mut n_mixed_groups = 0;
    for radius in [68,70,71,72,73,74,75,76,78]{
      let config = Config::from_toml_string(&format!(r##"
          input_structure_file = "assets/TEMPO.pdb"
          radius = {}
          detected_spin.position = [28,29]
          replicate_unit_cell = true
          {}
          "##,radius,groups)).unwrap();
      let radius = config.radius.unwrap();

      for seed in 0..4{
        set_store_all_copies(true);
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let reference = Structure::build_structure(&mut rng,&config).unwrap();

        set_store_all_copies(false);
        let mut rng = ChaCha20Rng::seed_from_u64(seed);
        let structure = Structure::build_structure(&mut rng,&config).unwrap();

        assert!(structure.number() < reference.number());
        assert_eq!(active_set(&structure), active_set(&reference),
            "radius {} Å, seed {}",radius/ANGSTROM,seed);

        // Count the methyl copies that the radius cuts through, and the
        // ones that do not form exchange groups because of mixed isotopes.
        let manager = reference.exchange_groups.as_ref().unwrap();
        let delta = 1.1*ANGSTROM;
        for group in manager.exchange_groups.iter(){
          if (group.centroid().norm() - radius).abs() < delta{
            n_cut_groups += 1;
          }
        }
        for cell in 1..reference.cell_offsets.len(){
          for &idx0 in [2,6,20,24].iter(){
            let indices = (0..3).filter_map(|ii|
                reference.cell_indices[idx0+ii][cell])
              .collect::<Vec::<usize>>();
            if indices.len() != 3 { continue; }
            let isotopes = indices.iter()
              .map(|&idx| reference.bath_particles[idx].isotope)
              .collect::<Vec::<Isotope>>();
            let r = reference.bath_particles[indices[0]].coordinates.norm();
            if (r - radius).abs() < delta
              && isotopes.iter().any(|iso| *iso != isotopes[0]){
              n_mixed_groups += 1;
            }
          }
        }
      }
    }
    assert!(n_cut_groups > 0);
    assert!(n_mixed_groups > 0);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_copy_can_be_active(){
    // Members of an exchange group are kept if the centroid of the group
    // is within the radius plus the largest member distance from the
    // centroid, and other particles if they are within the radius.
    let config = Config::from_toml_string(r##"
        input_structure_file = "assets/TEMPO.pdb"
        radius = 20
        detected_spin.position = [28,29]
        replicate_unit_cell = false
        "##).unwrap();
    let mut structure = pdb::parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    let mut rng = ChaCha20Rng::seed_from_u64(0);
    structure.build_primary_structure(&mut rng,&config).unwrap();
    let radius = config.radius.unwrap();
    let delta = structure.max_distance_from_exchange_group_center();
    assert!(delta > 0.9*ANGSTROM && delta < 1.2*ANGSTROM);

    // Returns the offset that moves `r0` to distance `distance`
    // along `r0`.
    let offset_to = |r0: &Vector3D, distance: f64| -> Vector3D{
      &r0.scale(distance/r0.norm()) - r0
    };

    // A methyl: all members get the answer of the centroid.
    let manager = structure.exchange_groups.as_ref().unwrap();
    let group = &manager.exchange_groups[0];
    let centroid = group.centroid().clone();
    for (distance, expected) in [(radius + 0.99*delta, true),
        (radius + 1.01*delta, false)]
    {
      let offset = offset_to(&centroid, distance);
      for idx0 in group.indices(){
        assert_eq!(structure.copy_can_be_active(idx0,&offset,radius,delta),
            expected);
      }
    }

    // The nitrogen is not in an exchange group.
    let idx_n = (0..structure.number()).find(|&idx|
        structure.bath_particles[idx].element == Element::Nitrogen).unwrap();
    assert!(manager.exchange_group_ids[idx_n].is_none());
    let r_n = structure.bath_particles[idx_n].coordinates.clone();
    for (distance, expected) in [(radius*(1.0 - 1e-6), true),
        (radius*(1.0 + 1e-6), false)]
    {
      let offset = offset_to(&r_n, distance);
      assert_eq!(structure.copy_can_be_active(idx_n,&offset,radius,delta),
          expected);
    }
  }
  //----------------------------------------------------------------------------

}
