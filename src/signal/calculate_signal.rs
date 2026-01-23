
use crate::build_adjacency_list;
use crate::find_clusters;
use crate::config::{
  Config,
  ClusterMethod,
  OrientationAveraging,
  ClusterSource,
  SAVE_FILE_BATH, 
  SAVE_FILE_CLUSTERS, 
  SAVE_FILE_DET_SPIN_COOR,
  SAVE_DIR_INFO, 
  SAVE_FILE_EXCHANGE_GROUPS,
  SAVE_FILE_METHYL_PARTITIONS,
  SAVE_DIR_ORIENTATION_SIGNALS,
  SAVE_FILE_ORIENTATIONS,
  SAVE_FILE_PARTITION_TABLE,
  SAVE_FILE_SANS_SPIN_SIGNALS,
  SAVE_FILE_STRUCTURE_PDB,
  SAVE_FILE_TENSORS,
};
use crate::CluEError;
use crate::cluster::methyl_clusters::partition_cluster_set_by_exchange_groups;
use crate::cluster::{
  cluster_set::ClusterSet, 
  read_clusters::read_cluster_file,
  unit_of_clustering::UnitOfClustering,
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
use crate::signal::Signal;
use crate::signal::write_vec_signals;
use crate::signal::cluster_correlation_expansion::*;
use crate::signal::calculate_analytic_restricted_2cluster_signals::{
  calculate_analytic_restricted_2cluster_signals};
use crate::Structure;
use crate::quantum::spin_hamiltonian::*;
use crate::math;
use crate::space_3d::UnitSpherePoint;

use num_complex::Complex;
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};
use rand_distr::{Distribution,Uniform};

//------------------------------------------------------------------------------
/// This function averages of `number_runs` instances of building
/// a system and calculating the signal.
pub fn calculate_signals(rng: &mut ChaCha20Rng, config: &Config,
    path_opt: &Option<String>) -> Result<Vec::<Signal>,CluEError>
{
  let Some(max_cluster_size) = config.max_cluster_size else {
    return Err(CluEError::NoMaxClusterSize);
  };
  let Some(number_runs) = config.number_runs else {
    return Err(CluEError::NoNumberSystemInstances);
  }; 

  let range = Uniform::from(0..u64::MAX);
  let new_seeds = (0..number_runs).map(|_ii|
      range.sample(rng) )
    .collect::<Vec::<u64>>();

  let n_tot = config.number_timepoints.iter().sum::<usize>();

  let mut signals = (0..max_cluster_size).map(|_ii| Signal::zeros(n_tot))
    .collect::<Vec::<Signal>>();

  for (ii,&seed) in new_seeds.iter().enumerate(){

    println!("\nSystem instance: {}/{}.", ii+1, number_runs);

    let mut inst_rng = ChaCha20Rng::seed_from_u64(seed);
    let inst_signals = calculate_structure_signal(&mut inst_rng, 
        config, path_opt)?;

    for (order_idx, sig) in inst_signals.iter().enumerate(){
      signals[order_idx] = &signals[order_idx] + sig;
    }
  }

  for sig in signals.iter_mut(){
    sig.mut_scale(ONE/(number_runs as f64));
  }
  Ok(signals)
}
//------------------------------------------------------------------------------
// This function builds the spin structure and calculates the signal over
// one or multiple orientations.
fn calculate_structure_signal(rng: &mut ChaCha20Rng, config: &Config,
    path_opt: &Option<String>) -> Result<Vec::<Signal>,CluEError>
{

  // Determine maximum cluster size.
  let Some(max_cluster_size) = config.max_cluster_size else{
    return Err(CluEError::NoMaxClusterSize);
  };

  
  // Initialize output.
  let mut order_n_signals = {
    // Get number of data points per trace.
    let n_tot = config.number_timepoints.iter().sum::<usize>();

    (0..max_cluster_size).map(|_| Signal::zeros(n_tot))
    .collect::<Vec::<Signal>>()
  };


  // Build input structure.
  let structure = Structure::build_structure(rng,config)?;


  // Generate coupling tensors.
  let tensors = HamiltonianTensors::generate(rng, &structure, config)?;


  // Determine if/where to save results.
  let save_dir_opt = get_system_save_dir_opt(path_opt,&structure,config)?;

  optionally_save_tensors(&save_dir_opt,&tensors,&structure,config)?;


  let spin_ops = {

    let spin_multiplicity_set =
        math::unique(tensors.spin_multiplicities.clone());

    match config.unit_of_clustering{
    Some(UnitOfClustering::Spin) 
        => ClusterSpinOperators::new(&spin_multiplicity_set,max_cluster_size)?,

    Some(UnitOfClustering::Set) => {
      let max_spins_per_cluster_unit = match config.partitioning{
        Some(PartitioningMethod::Particles) => 1,
        Some(PartitioningMethod::ExchangeGroupsAndParticles) => 3,
        Some(PartitioningMethod::KMeans(_)) => usize::MAX,
        None => return Err(CluEError::NoPartitioningMethod),
      };


      let potential_max_order = max_cluster_size*max_spins_per_cluster_unit;
      
      let max_spins = match config.max_spins{
        Some(s) => s,
        None => potential_max_order,
      };

      println!("DB: s = {}; p ={}; c = {:?}",
          max_spins,potential_max_order,config.max_spins);
      ClusterSpinOperators::new(&spin_multiplicity_set,
        max_spins)?
    },
    None => return Err(CluEError::NoUnitOfClustering),
    }
  };


  // Determin orientation averaging method.
  let integration_grid = match &config.orientation_grid{
    Some(OrientationAveraging::Grid(grid)) => grid.clone(),
    Some(OrientationAveraging::Lebedev(n_ori)) 
      => IntegrationGrid::lebedev(*n_ori)?.remove_3d_hemisphere(),
    Some(OrientationAveraging::Random(n_ori)) 
      => IntegrationGrid::random_unit_sphere(*n_ori,rng),
    None => IntegrationGrid::z_3d(),
  };

  if config.write_orientation_grid == Some(true){
    if let Some(save_dir) = &save_dir_opt{
      
      let info_path = format!("{}/{}",save_dir,SAVE_DIR_INFO);

      if std::fs::create_dir_all(info_path.clone()).is_err(){
        return Err(CluEError::CannotCreateDir(info_path));
      }
      let ori_file = format!("{}/{}.xyzw",info_path,SAVE_FILE_ORIENTATIONS);

      integration_grid.to_csv(&ori_file)?;

      let ori_file = format!("{}/{}_polar.csv",info_path,SAVE_FILE_ORIENTATIONS);

      integration_grid.to_polar_csv(&ori_file)?;
    }
  }


  // Loop over orientations.
  for iori in 0..integration_grid.len(){
    
    let rot_dir = UnitSpherePoint::from(integration_grid.xyz(iori)?);
    
    let mut ori_sigs = calculate_signal_at_orientation(rng,rot_dir,
        &spin_ops,tensors.clone(), &structure,config, &save_dir_opt )?;

    let weight = Complex::<f64>{ re: integration_grid.weight(iori), im: 0.0};

    for size_idx in 0..max_cluster_size{
      ori_sigs[size_idx].mut_scale(weight);
      order_n_signals[size_idx] 
        = &order_n_signals[size_idx] + &ori_sigs[size_idx];
    }

  }


  // Optionally save signals.
  if let Some(save_dir) =  &save_dir_opt{
    let save_path = format!("{}/signal.csv", save_dir);
    let headers = (1..=max_cluster_size)
      .map(|ii| format!("signal_{}",ii))
      .collect::<Vec::<String>>();
    write_vec_signals(&order_n_signals, headers, &save_path)?;
  }

  Ok(order_n_signals)
}
//------------------------------------------------------------------------------
// This function calculates the signal for the input structure at the input
// orientation.
fn calculate_signal_at_orientation(rng: &mut ChaCha20Rng, 
    rot_dir: UnitSpherePoint,
    spin_ops: &ClusterSpinOperators,
    mut tensors: HamiltonianTensors, structure: &Structure,
    config: &Config, path_opt: &Option<String>)
  -> Result<Vec::<Signal>,CluEError>
{
  
  let (save_dir,save_path_opt) 
    = get_orientation_save_dir(&rot_dir,config,path_opt)?;

  // Rotate the coupling tensor to the specified orientation.
  tensors.rotate_pasive(&rot_dir);


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

      // Pull out unit_of_clustering early, before more expensive calculations.
      let Some(unit_of_clustering) = &config.unit_of_clustering else{
        return Err(CluEError::NoUnitOfClustering);
      };

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
          let save_path = format!("{}/{}",path,save_dir);
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
      let mut clu_set = expand_block_clusters(block_cluster_set,&partition_table,
          unit_of_clustering)?;

      if let Some(max_spins) = config.max_spins{
        clu_set.prune_large_clusters(max_spins)?;
      } 

     clu_set
    },
    None => return Err(CluEError::NoClusterSource),
  };

  for (size_idx,clusters_of_size) in cluster_set.clusters.iter().enumerate(){
    println!("Found {} clusters of size {}.",clusters_of_size.len(),size_idx+1);
  }

  if let Some(path) = &save_path_opt{
    if config.write_clusters == Some(true){
      let cluster_save_path = format!("{}/{}.toml",path,SAVE_FILE_CLUSTERS);
      cluster_set.save(&cluster_save_path,structure)?;
    }
  }


  // Calculate cluster signals.
  let order_n_signals = match &config.cluster_method{
    Some(ClusterMethod::AnalyticRestricted2CCE) => 
      calculate_analytic_restricted_2cluster_signals(&mut cluster_set, &tensors,
          config,&save_path_opt,structure)?,
    Some(_cce) => {
      do_cluster_correlation_expansion(&mut cluster_set, spin_ops, &tensors, 
          config,&save_path_opt,structure)?
    },
    None => return Err(CluEError::NoClusterMethod)
  };


  if let Some(save_dir) =  &save_path_opt{
    let save_path = format!("{}/signal.csv", save_dir);
    let headers = (1..=max_cluster_size)
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
  let n_tot = config.number_timepoints.iter().sum::<usize>();

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
    Err(_) => return Err(CluEError::CannotCreateDir(save_path.to_string())),
  }

  if let Some(exchange_group_manager) = &structure.exchange_groups{

    let cluster_partitions = partition_cluster_set_by_exchange_groups(
        cluster_set, exchange_group_manager, structure)?;

    let n_tot = config.number_timepoints.iter().sum::<usize>();

    for (key_name, value_cluster_set) in cluster_partitions.iter(){

      if key_name == "methyl_1718_1719_1720_methyl_4448_4449_4450"{
        println!("DB: {}",key_name);
      }

      let mut part_signal = Signal::ones(n_tot);

      for size_idx in 0..(*value_cluster_set).len() {
        for cluster in value_cluster_set.clusters[size_idx].iter(){
          if key_name == "methyl_1718_1719_1720_methyl_4448_4449_4450"{
            println!("DB: {:?}",cluster.vertices);
          }
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
fn get_system_save_dir_opt(
    path_opt: &Option<String>,
    structure: &Structure,
    config: &Config,
    ) -> Result<Option<String>,CluEError>
{
  let save_dir_opt = match path_opt{
    Some(path) => {
      
      
      let save_dir = match &config.run_name{
        Some(run_name) => {
          format!("{}/{}",path,run_name)
        },
        None => {
          let structure_hash = math::str_hash(&structure);
          format!("{}/run-{}",path,structure_hash)
        }
      };

      match std::fs::create_dir_all(save_dir.clone()){
        Ok(_) => (),
        Err(_) => return Err(CluEError::CannotCreateDir(save_dir)),
      }
  
      if config.write_info == Some(true){

        let info_path = format!("{}/{}",save_dir,SAVE_DIR_INFO);

        if std::fs::create_dir_all(info_path.clone()).is_err(){
          return Err(CluEError::CannotCreateDir(info_path));
        }

        if config.write_detected_spin == Some(true){
          if let Some(det) = &structure.detected_particle{
            let mut coors = det.weighted_coordinates.clone();
            let pdb_origin = &structure.pdb_origin;
            coors.translate_mut(
                &vec![-pdb_origin.x(),-pdb_origin.y(),-pdb_origin.z()]
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
      Some(save_dir)
    },
    None => None,
  };

  Ok(save_dir_opt)
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
        if std::fs::create_dir_all(info_path.clone()).is_err(){
          return Err(CluEError::CannotCreateDir(info_path));
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
          Err(_) => return Err(CluEError::CannotCreateDir(save_path)),
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
  



