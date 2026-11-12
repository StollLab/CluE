
use crate::build_adjacency_list;
use crate::find_clusters;
use crate::config::{
  Config,
  ClusterMethod,
  OrientationAveraging,
  ClusterSource,
  StructureInput,
  SAVE_FILE_BATH, 
  SAVE_FILE_CLUSTERS, 
  SAVE_FILE_DET_SPIN_COOR,
  SAVE_DIR_INFO, 
  SAVE_FILE_EXCHANGE_GROUPS,
  SAVE_FILE_METHYL_PARTITIONS,
  SAVE_DIR_ORIENTATION_SIGNALS,
  SAVE_FILE_ORIENTATIONS,
  SAVE_FILE_SANS_SPIN_SIGNALS,
  SAVE_FILE_STRUCTURE_PDB,
  SAVE_FILE_TENSORS,
};
use crate::CluEError;
use crate::cluster::methyl_clusters::partition_cluster_set_by_exchange_groups;
use crate::cluster::{
  cluster_set::ClusterSet, 
  read_clusters::read_cluster_file,
};
use crate::cluster::partition::{
  expand_block_clusters,
  get_partition_table,
  PartitioningMethod,
  partition_system,
};
use crate::HamiltonianTensors;
use crate::integration_grid::IntegrationGrid;
use crate::physical_constants::{ANGSTROM,ONE,PI};
use crate::signal::{
  Signal,
  load_csv_to_vec_signals,
  write_vec_signals,
  cluster_signals::calculate_cluster_signals,
};
use crate::output_manifest::{
  OutputState,
  SIGNAL_FILE,
  draw_seeds,
  instance_seeds,
  structure_seed,
};
use crate::Structure;
use crate::quantum::{
  general_spin_hamiltonian::build_detected_spin_hamiltonian,
  cluster_operators::ClusterSpinOperators,
  spin_states::SpinStates,
};
use crate::math;
use crate::space_3d::UnitSpherePoint;

use num_complex::Complex64;
use ndarray::{Array1,Array2};
type CxMat = Array2::<Complex64>;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

//------------------------------------------------------------------------------
/// This function calculates the signal of each input structure, averaged 
/// over `number_runs` instances of building the system, and returns the
/// average over the structures, weighted by their weights.
/// Every random draw derives from the main seed `run.rng_seed`: each 
/// structure has a seed (`output_manifest::structure_seed()`), each 
/// instance a seed drawn from the structure seed, and each orientation a 
/// seed drawn from the instance generator after the structure, the 
/// coupling tensors, and the orientation grid are built.
/// `output` records the progress; complete instances and orientations
/// whose files are present are loaded instead of computed.
pub fn calculate_signals(config: &Config, output: &mut OutputState) 
  -> Result<Vec::<Signal>,CluEError>
{
  let Some(max_cluster_size) = config.max_cluster_size else {
    return Err(CluEError::NoMaxClusterSize);
  };
  let Some(number_runs) = config.number_runs else {
    return Err(CluEError::NoNumberSystemInstances);
  }; 
  let Some(main_seed) = config.rng_seed else{
    return Err(CluEError::NoRngSeed);
  };

  let n_tot = config.get_total_number_timesteps();

  let structures = config.structure_inputs()?;
  let total_weight = structures.iter().map(|s| s.weight).sum::<f64>();

  let mut signals = (0..=max_cluster_size).map(|_ii| Signal::zeros(n_tot))
    .collect::<Vec::<Signal>>();

  for structure_input in structures.iter(){

    let seed = structure_seed(main_seed, structure_input.content_hash);
    let new_seeds = instance_seeds(seed, number_runs);
    let s_idx = output.begin_structure(structure_input, seed, &new_seeds)?;

    let mut structure_signals = (0..=max_cluster_size)
      .map(|_ii| Signal::zeros(n_tot))
      .collect::<Vec::<Signal>>();

    for (ii,&seed) in new_seeds.iter().enumerate(){

      println!("\nSystem instance: {}/{}.", ii+1, number_runs);

      let saved_signals = if output.instance_is_complete(s_idx, ii){
        output.instance_dir(s_idx, ii).and_then(|dir| 
            load_saved_signals(&format!("{}/{}", dir, SIGNAL_FILE),
              max_cluster_size, n_tot))
      }else{
        None
      };

      let inst_signals = match saved_signals{
        Some(signals) => {
          println!("Loaded the saved signal of this instance.");
          signals
        },
        None => calculate_structure_signal(seed, config, structure_input,
            s_idx, ii, output)?,
      };

      for (order_idx, sig) in inst_signals.iter().enumerate(){
        structure_signals[order_idx] = &structure_signals[order_idx] + sig;
      }
    }

    for sig in structure_signals.iter_mut(){
      sig.scale_mut(ONE/(number_runs as f64));
    }

    if let Some(structure_dir) = output.structure_dir(s_idx){
      write_signals(&structure_signals, 
          &format!("{}/{}", structure_dir, SIGNAL_FILE))?;
    }

    // With one structure, the weight is 1, and the scaling is exact.
    let weight = Complex64{re: structure_input.weight/total_weight, im: 0.0};
    for (order_idx, sig) in structure_signals.iter().enumerate(){
      signals[order_idx] = &signals[order_idx] + &sig.scale(weight);
    }
  }

  output.write_structures_index()?;

  Ok(signals)
}
//------------------------------------------------------------------------------
// This function loads saved signals of orders 0, ..., `max_cluster_size` 
// from `filename`. It returns `None`, with a message, if the file cannot be
// read or does not have that many signals with `n_tot` time points each;
// the signals are then computed again.
fn load_saved_signals(filename: &str, max_cluster_size: usize, n_tot: usize)
  -> Option<Vec::<Signal>>
{
  match load_csv_to_vec_signals(filename){
    Ok(mut signals) => {
      if signals.len() > max_cluster_size 
        && signals.iter().all(|s| s.len() == n_tot)
      {
        signals.truncate(max_cluster_size + 1);
        Some(signals)
      }else{
        println!("\"{}\" does not have {} signals of {} time points; \
it is computed again.", filename, max_cluster_size + 1, n_tot);
        None
      }
    },
    Err(err) => {
      println!("{}; it is computed again.", err);
      None
    },
  }
}
//------------------------------------------------------------------------------
// This function writes the signals of orders 0, 1, ..., to `filename`.
fn write_signals(signals: &[Signal], filename: &str) -> Result<(),CluEError>{
  let headers = (0..signals.len())
    .map(|ii| format!("signal_{}",ii))
    .collect::<Vec::<String>>();
  write_vec_signals(signals, headers, filename)
}
//------------------------------------------------------------------------------
// This function builds instance `k` (from 0) of structure `structure_input`
// (record `s_idx` of `output`) from the instance seed `seed`, and 
// calculates the signal averaged over the orientations.
// Orientations that `output` records as complete are loaded.
fn calculate_structure_signal(seed: u64, config: &Config, 
    structure_input: &StructureInput, s_idx: usize, k: usize,
    output: &mut OutputState) -> Result<Vec::<Signal>,CluEError>
{
  let mut rng = ChaCha20Rng::seed_from_u64(seed);
  output.begin_instance(s_idx, k)?;

  // Determine maximum cluster size.
  let Some(max_cluster_size) = config.max_cluster_size else{
    return Err(CluEError::NoMaxClusterSize);
  };

  
  // Initialize output.
  let mut order_n_signals = {
    // Get number of data points per trace.
    let n_tot = config.get_total_number_timesteps();

    (0..=max_cluster_size).map(|_| Signal::zeros(n_tot))
    .collect::<Vec::<Signal>>()
  };


  // Build input structure.
  let structure = Structure::build_structure_from_file(&mut rng, config,
      &structure_input.file)?;

  // Generate coupling tensors.
  let tensors = HamiltonianTensors::generate(&mut rng, &structure, config)?;


  // Determine if/where to save results.
  let save_dir_opt = output.instance_dir(s_idx, k);
  if let Some(save_dir) = &save_dir_opt{
    write_instance_info(save_dir, &structure, config)?;
  }

  optionally_save_tensors(&save_dir_opt,&tensors,&structure,config)?;


  let n_electron = if config.cluster_method == Some(ClusterMethod::GCCE){
    1
  }else{
    0
  };

  let spin_ops = {

    let spin_multiplicity_set =
        math::unique(tensors.spin_multiplicities.clone());

   let max_spins_per_cluster_unit = match config.partitioning{
     Some(PartitioningMethod::Particles) => 1,
     Some(PartitioningMethod::ExchangeGroupsAndParticles) => 3,
     Some(PartitioningMethod::KMeans(kmeans_size)) => kmeans_size,
     Some(PartitioningMethod::RestrictedKMeans(kmeans_size)) => kmeans_size,
     None => return Err(CluEError::NoPartitioningMethod),
   };


   let potential_max_order = n_electron 
       + max_cluster_size*max_spins_per_cluster_unit;
   
   let max_spins = match config.max_spins{
     Some(s) => s + n_electron,
     None => potential_max_order,
   };

   let det_hamiltonian = if config.cluster_method == Some(ClusterMethod::GCCE){
     build_detected_spin_hamiltonian(&tensors)?
   }else{
     (Array1::<f64>::zeros(1) ,CxMat::eye(1))
   };

   ClusterSpinOperators::new(det_hamiltonian,&spin_multiplicity_set,
     max_spins,config)?
  };


  // Determin orientation averaging method.
  let integration_grid = match &config.orientation_grid{
    Some(OrientationAveraging::Grid(grid)) => grid.clone(),
    Some(OrientationAveraging::Lebedev(n_ori)) 
      => IntegrationGrid::lebedev(*n_ori)?.remove_3d_hemisphere(),
    Some(OrientationAveraging::Random(n_ori)) 
      => IntegrationGrid::random_unit_sphere(*n_ori,&mut rng),
    None => IntegrationGrid::z_3d(),
  };

  if config.write_orientation_grid == Some(true){
    if let Some(save_dir) = &save_dir_opt{
      
      let info_path = format!("{}/{}",save_dir,SAVE_DIR_INFO);

      if let Err(err) = std::fs::create_dir_all(info_path.clone()){
        return Err(CluEError::CannotCreateDir(info_path,
              err.to_string()));
      }
      let ori_file = format!("{}/{}.xyzw",info_path,SAVE_FILE_ORIENTATIONS);

      integration_grid.to_csv(&ori_file)?;

      let ori_file = format!("{}/{}_polar.csv",info_path,SAVE_FILE_ORIENTATIONS);

      integration_grid.to_polar_csv(&ori_file)?;
    }
  }

  // One seed per orientation, so that each orientation can be computed,
  // and resumed, on its own.
  let orientation_seeds = draw_seeds(&mut rng, integration_grid.len());
  output.set_orientation_seeds(s_idx, k, &orientation_seeds)?;

  let n_tot = config.get_total_number_timesteps();

  // Loop over orientations.
  for iori in 0..integration_grid.len(){
    
    let rot_dir = UnitSpherePoint::from(integration_grid.xyz(iori)?);

    let (ori_name, ori_path_opt) 
      = get_orientation_save_dir(&rot_dir,config,&save_dir_opt)?;

    let saved_signals = if output.orientation_is_complete(s_idx, k, iori){
      ori_path_opt.as_ref().and_then(|dir| 
          load_saved_signals(&format!("{}/{}", dir, SIGNAL_FILE),
            max_cluster_size, n_tot))
    }else{
      None
    };
    
    let mut ori_sigs = match saved_signals{
      Some(signals) => {
        println!("Loaded the saved signal of this orientation.");
        signals
      },
      None => {
        output.begin_orientation(s_idx, k, iori)?;
        let mut ori_rng = ChaCha20Rng::seed_from_u64(orientation_seeds[iori]);
        let signals = calculate_signal_at_orientation(&mut ori_rng, rot_dir,
            &spin_ops,tensors.clone(), &structure,config, &ori_name,
            &ori_path_opt)?;
        output.finish_orientation(s_idx, k, iori, &ori_path_opt)?;
        signals
      },
    };

    let weight = Complex64{ re: integration_grid.weight(iori), im: 0.0};

    for size_idx in 0..=max_cluster_size{
      ori_sigs[size_idx].scale_mut(weight);
      order_n_signals[size_idx] 
        = &order_n_signals[size_idx] + &ori_sigs[size_idx];
    }

  }


  // Optionally save signals.
  if let Some(save_dir) =  &save_dir_opt{
    write_signals(&order_n_signals, &format!("{}/{}", save_dir, SIGNAL_FILE))?;
  }
  output.finish_instance(s_idx, k)?;

  Ok(order_n_signals)
}
//------------------------------------------------------------------------------
// This function calculates the signal for the input structure at the input
// orientation.
// `save_dir` is the name of the orientation directory, and `save_path_opt`
// its path, if the orientation's files are saved.
fn calculate_signal_at_orientation(rng: &mut ChaCha20Rng, 
    rot_dir: UnitSpherePoint,
    spin_ops: &ClusterSpinOperators,
    mut tensors: HamiltonianTensors, structure: &Structure,
    config: &Config, save_dir: &str, save_path_opt: &Option<String>)
  -> Result<Vec::<Signal>,CluEError>
{

  // Rotate the coupling tensor to the specified orientation.
  tensors.rotate_pasive(&rot_dir)?;

  let states = SpinStates::generate(rng, &tensors,config)?;
  
  let Some(ensemble_cce) = config.ensemble_cce else{
    return Err(CluEError::NoEnsembleCCE);
  };
  let Some(do_mean_field_averaging) = config.mean_fields else{
    return Err(CluEError::NoMeanFields);
  };
  
  if ensemble_cce && do_mean_field_averaging{
    return Err(CluEError::MeanFieldECCENotImplemented);
  }

  if do_mean_field_averaging{
    tensors.set_mean_field_couplings(&states.states)?;
  }

  // Determine maximum cluster size.
  let Some(max_cluster_size) = config.max_cluster_size else{
    return Err(CluEError::NoMaxClusterSize);
  };


  // Find clusters.
  let mut cluster_set = match &config.cluster_source{
    Some(ClusterSource::File(clusters_file)) 
      => read_cluster_file(clusters_file, structure)?,
    Some(ClusterSource::OrientationsDir(ori_dir)) => {
      let clusters_file = format!("{}/{}/clusters.toml",ori_dir,save_dir);
      read_cluster_file(&clusters_file, structure)?
    },
    Some(ClusterSource::Structure) => {

      // Determine spin adjacencies.
      let spin_adjacency_list 
          = build_adjacency_list(&tensors, structure, config)?;

      // Partition the system into blocks.
      // The partition table references spins by their tensor/active index.
      let partition_table = get_partition_table(rng,
         &spin_adjacency_list,
         &tensors, structure, config)?;

      if let Some(path) = &save_path_opt{
        if config.write_partition_table == Some(true){
          partition_table.to_csv(path,structure)?;
        }
      } 
      // Determine block adjacencies for the partitioned system.
      let block_adjacency_list 
          = partition_system(&spin_adjacency_list,&partition_table)?;

      // Find clusters of blocks.
      let block_cluster_set 
          = find_clusters(&block_adjacency_list, max_cluster_size)?;

      // Expand block clusters out to spin clusters.
      let mut clu_set = expand_block_clusters(
            block_cluster_set,&partition_table)?;

      if let Some(max_spins) = config.max_spins{
        clu_set.prune_large_clusters(max_spins)?;
      } 

     clu_set
    },
    None => return Err(CluEError::NoClusterSource),
  };

  for (size_idx,clusters_of_size) in cluster_set.clusters
      .iter().enumerate().skip(1)
  {
    println!("Found {} clusters of size {}.",clusters_of_size.len(),size_idx);
  }

  if let Some(path) = &save_path_opt{
    if config.write_clusters == Some(true){
      let cluster_save_path = format!("{}/{}.toml",path,SAVE_FILE_CLUSTERS);
      cluster_set.save(&cluster_save_path,structure)?;
    }
  }


  // Calculate cluster signals.
  let order_n_signals = match &config.cluster_method{
    //Some(ClusterMethod::AnalyticRestricted2CCE) => 
    //  calculate_analytic_restricted_2cluster_signals(&mut cluster_set, &tensors,
    //      config,&save_path_opt,structure)?,
    Some(_cce) => {
      calculate_cluster_signals(&mut cluster_set, &states, spin_ops, &tensors, 
          config,&save_path_opt,structure)?
    },
    None => return Err(CluEError::NoClusterMethod)
  };

  if let Some(save_dir) =  &save_path_opt{
    let save_path = format!("{}/signal.csv", save_dir);
    let headers = (0..=max_cluster_size)
      .map(|ii| format!("signal_{}",ii))
      .collect::<Vec::<String>>();

    write_vec_signals(&order_n_signals, headers, &save_path)?;

    if config.write_sans_spin_signals == Some(true){
      let save_path = format!("{}/{}",save_dir,SAVE_FILE_SANS_SPIN_SIGNALS);

      caculate_sans_spin_signals(&order_n_signals[max_cluster_size - 1],
        &cluster_set, structure, config, &save_path)?;
    }


    if config.write_methyl_partitions == Some(true){
      let save_path = format!("{}/{}",save_dir,SAVE_FILE_METHYL_PARTITIONS);
      
      calculate_methyl_partition_cce(cluster_set, structure, config, 
          &save_path)?;
    }
  }

  Ok(order_n_signals)
} 
//------------------------------------------------------------------------------
// For each active spin in structure, this function calculates the signal,
// where the spin is dropped. 
fn caculate_sans_spin_signals(ref_signal: &Signal,
    cluster_set: &ClusterSet, structure: &Structure, config: &Config, 
    save_path: &str)
  -> Result<(),CluEError>
{
  
  let mut spin_signals 
    = caculate_bath_spin_contributions(cluster_set, structure, config)?;

  for sig in spin_signals.iter_mut(){
    *sig = ref_signal / sig;
  }

  let headers = structure.bath_particles.iter().enumerate()
    .filter_map(|(idx,particle)| 
        if particle.active{
          Some( format!("sans_{}",idx) )
        }else{
          None
        }
    ).collect::< Vec::<String> >();

  write_vec_signals(&spin_signals,headers,save_path)
}
//------------------------------------------------------------------------------
// For each active spin in structure, this function calculates the product of
// auxiliary signals for clusters that contain the spin.
// Note that since clusters contain multiple spins, the contributions from
// multiple spins are not independent of each other, 
// meaning that the product of contributions will not in general equal the 
// calculated signal.
fn caculate_bath_spin_contributions(
    cluster_set: &ClusterSet, structure: &Structure, config: &Config)
  -> Result<Vec::<Signal>,CluEError>
{

  // Get number of data points per trace.
  let n_tot = config.get_total_number_timesteps();

  let mut spin_contributions = (0..structure.number_active())
    .map(|_| Signal::ones(n_tot)).collect::<Vec::<Signal>>();

  let mut index = 0;
  let vertex_to_index = structure.bath_particles.iter()
    .filter_map(|particle| 
        if particle.active{
          let idx = index;
          index += 1;
          Some(idx)
        }else{
          None
        }  
    ).collect::< Vec::<usize> >();

  for clusters_of_size in cluster_set.clusters.iter(){

    for cluster in clusters_of_size{

      for &vertex in cluster.vertices(){

        let idx = vertex_to_index[vertex-1];

        match &cluster.signal{
          Ok(Some(signal)) => 
            spin_contributions[idx] = &spin_contributions[idx] * signal,
          Ok(None) => (),
          Err(err) => return Err(err.clone()),
        }
      } 
    }

  }

  Ok(spin_contributions)
}
//------------------------------------------------------------------------------
// This function partition the auxiliary signals by which methyl groups
// contribute to them and saves the partitions to csv files.
fn calculate_methyl_partition_cce(
    cluster_set: ClusterSet, structure: &Structure, config: &Config, 
    save_path: &str)
  -> Result<(),CluEError>
{
  match std::fs::create_dir_all(save_path){
    Ok(_) => (),
    Err(err) => return Err(CluEError::CannotCreateDir(save_path.to_string(),
          err.to_string())),
  }

  if let Some(exchange_group_manager) = &structure.exchange_groups{

    let cluster_partitions = partition_cluster_set_by_exchange_groups(
        cluster_set, exchange_group_manager, structure)?;

    let n_tot = config.get_total_number_timesteps();

    for (key_name, value_cluster_set) in cluster_partitions.iter(){

      let mut part_signal = Signal::ones(n_tot);

      for size_idx in 0..(*value_cluster_set).len() {
        for cluster in value_cluster_set.clusters[size_idx].iter(){
          if let Ok(Some(aux_signal)) = &cluster.signal{
            part_signal = &part_signal * aux_signal;
          }
        }
      }

      let part_save_path = format!("{}/{}.csv",
          save_path,key_name);
      part_signal.write_to_csv(&part_save_path)?;
    }
  }

  //let signal = do_cluster_correlation_expansion_product(&cluster_set,config)?;
  Ok(())
} 
//------------------------------------------------------------------------------
// This function creates the directory `save_dir` of an instance and
// writes its info files (`output.info` and the switches below it).
fn write_instance_info(
    save_dir: &str,
    structure: &Structure,
    config: &Config,
    ) -> Result<(),CluEError>
{
  match std::fs::create_dir_all(save_dir){
    Ok(_) => (),
    Err(err) => return Err(CluEError::CannotCreateDir(save_dir.to_string(),
          err.to_string())),
  }

  if config.write_info == Some(true){

    let info_path = format!("{}/{}",save_dir,SAVE_DIR_INFO);

    if let Err(err) = std::fs::create_dir_all(info_path.clone()){
      return Err(CluEError::CannotCreateDir(info_path,
            err.to_string()));
    }

    if config.write_detected_spin == Some(true){
      if let Some(det) = &structure.detected_particle{
        let mut coors = det.weighted_coordinates.clone();
        let pdb_origin = &structure.pdb_origin;
        coors.translate_mut(
            &[-pdb_origin.x(),-pdb_origin.y(),-pdb_origin.z()]
        );
        coors = coors.scale(1.0/ANGSTROM);
        coors.to_csv(
          &format!("{}/{}.xyzw", info_path, SAVE_FILE_DET_SPIN_COOR),
        )?;
      }
    }

    if config.write_bath == Some(true){
      structure.bath_to_csv(
          &format!("{}/{}.csv", info_path, SAVE_FILE_BATH),
          config)?;
    }

    if config.write_structure_pdb == Some(true){
      structure.write_pdb(&format!("{}/{}.pdb", 
            info_path, SAVE_FILE_STRUCTURE_PDB))?;
    }

    if let Some(exchange_group_manager) = &structure.exchange_groups{
      if config.write_exchange_groups == Some(true){
        let csv_file = format!("{}/{}.csv",
            info_path,SAVE_FILE_EXCHANGE_GROUPS);
        exchange_group_manager.to_csv(&csv_file,structure)?;
      }
    }
  }

  Ok(())
}
//------------------------------------------------------------------------------
fn optionally_save_tensors(
      save_dir_opt: &Option<String>,
      tensors: &HamiltonianTensors,
      structure: &Structure,
      config: &Config
      ) -> Result<(),CluEError>
{
  if let Some(save_dir) = &save_dir_opt{
    if config.write_info == Some(true){
      let info_path = format!("{}/{}",save_dir,SAVE_DIR_INFO);
      
      if config.write_tensors == Some(true){
        if let Err(err) = std::fs::create_dir_all(info_path.clone()){
          return Err(CluEError::CannotCreateDir(info_path,
                err.to_string()));
        }
        let tensor_path = format!("{}/{}.txt",info_path,SAVE_FILE_TENSORS);

        tensors.save(&tensor_path,structure)?;
      }
    }
  }
  Ok(())
}
//------------------------------------------------------------------------------
fn get_orientation_save_dir(
      rot_dir: &UnitSpherePoint,
      config: &Config,
      path_opt: &Option<String>,
      ) -> Result<(String,Option<String>),CluEError>
{
  let theta_degrees = rot_dir.theta()*180.0/PI;
  let phi_degrees = rot_dir.phi()*180.0/PI;
  let save_dir = format!("theta_{}deg_phi_{}deg",theta_degrees, phi_degrees);
  println!(
      "\nMagnetic field orientation: θ = {} degrees; φ = {} degrees.", 
      theta_degrees, phi_degrees);
  // Determine if/where to save results.
  let save_path_opt = match path_opt{
    Some(path) => {

      if config.write_orientation_signals == Some(true){ 
      
        let save_path = format!("{}/{}/{}",
            path, SAVE_DIR_ORIENTATION_SIGNALS,save_dir);
      
        match std::fs::create_dir_all(save_path.clone()){
          Ok(_) => (),
          Err(err) => return Err(CluEError::CannotCreateDir(save_path,
                err.to_string())),
        }

        Some(save_path)
      }else{
        None
      }
    },
    None => None,
  };

  Ok((save_dir,save_path_opt))
}
//------------------------------------------------------------------------------
  



