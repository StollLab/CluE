//! The record of an output directory (`clue_manifest.toml`), the seeds of
//! the structures, instances, and orientations of a run, and the decision
//! whether an existing output directory is new, complete, can be resumed,
//! or belongs to a different simulation.
//!
//! Layout of an output directory:
//! ```text
//! <output>/
//!   clue_manifest.toml   this record
//!   rng_seed.txt         main seed
//!   tau_axis.csv
//!   signal.csv           weighted average over the structures
//!   structures.csv       id, file, content_hash, weight, number_runs, directory
//!   structures/<id>/
//!     signal.csv         average over the instances of the structure
//!     instance_<k>/      k = 1, ..., number_runs
//!       signal.csv
//!       info/
//!       orientations/theta_<θ>deg_phi_<φ>deg/
//! ```
use crate::clue_errors::CluEError;
use crate::config::{
  ClusterSource, Config, PartitionTableConfig, StructureInput,
  SAVE_DIR_ORIENTATION_SIGNALS,
};
use crate::io::write_atomically;
use crate::math::fnv1a_64;

use rand_chacha::{rand_core::{RngCore,SeedableRng}, ChaCha20Rng};
use rand_distr::{Distribution,Uniform};
use serde::{Serialize,Deserialize};

use std::path::Path;

/// The name of the record of an output directory.
pub const MANIFEST_FILE: &str = "clue_manifest.toml";
/// The name of the file with the main seed.
pub const RNG_SEED_FILE: &str = "rng_seed.txt";
/// The name of the index of the structures of an output directory.
pub const STRUCTURES_FILE: &str = "structures.csv";
/// The directory of the structures in an output directory.
pub const SAVE_DIR_STRUCTURES: &str = "structures";
/// The name of the signal files.
pub const SIGNAL_FILE: &str = "signal.csv";
/// The name of the time axis file.
pub const TAU_AXIS_FILE: &str = "tau_axis.csv";

// The version of the layout of `clue_manifest.toml`.
const MANIFEST_FORMAT: u32 = 1;

/// The version of CluE that writes the output.
pub const CLUE_VERSION: &str = env!("CARGO_PKG_VERSION");
/// The version of the Rust compiler that compiled CluE (set by build.rs).
pub const RUSTC_VERSION: &str = env!("CLUE_RUSTC_VERSION");

// The names of the output switches, as recorded in the manifest.
const OUT_INFO_BATH: &str = "info.bath";
const OUT_INFO_DETECTED_SPIN: &str = "info.detected_spin";
const OUT_INFO_EXCHANGE_GROUPS: &str = "info.exchange_groups";
const OUT_INFO_ORIENTATION_GRID: &str = "info.orientation_grid";
const OUT_INFO_STRUCTURE_PDB: &str = "info.structure_pdb";
const OUT_INFO_TENSORS: &str = "info.tensors";
const OUT_ORIENTATION_SIGNALS: &str = "orientation_signals";
const OUT_AUXILIARY_SIGNALS: &str = "auxiliary_signals";
const OUT_CLUSTERS: &str = "clusters";
const OUT_METHYL_PARTITIONS: &str = "methyl_partitions";
const OUT_PARTITION_TABLE: &str = "partition_table";
const OUT_SANS_SPIN_SIGNALS: &str = "sans_spin_signals";

// The most differences of the configuration listed in an error message.
const MAX_LISTED_DIFFERENCES: usize = 40;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// Seeds.
//------------------------------------------------------------------------------
/// This function draws a main seed from the operating system's entropy.
/// The seed is at most `i64::MAX`, so that it can be written as
/// `rng_seed` in the input file, whose integers are signed 64-bit integers.
pub fn draw_main_seed() -> u64{
  rand::random::<u64>() >> 1
}
//------------------------------------------------------------------------------
/// This function returns the seed of the structure whose file contents have
/// the hash `content_hash`: the first number of the stream `content_hash`
/// of the generator seeded with `main_seed`.
/// The seed does not depend on the position of the structure in a list of
/// structures, so adding structures leaves the seeds of the others
/// unchanged.
pub fn structure_seed(main_seed: u64, content_hash: u64) -> u64{
  let mut rng = ChaCha20Rng::seed_from_u64(main_seed);
  rng.set_stream(content_hash);
  rng.next_u64()
}
//------------------------------------------------------------------------------
/// This function draws `n` seeds from `rng`.
pub fn draw_seeds(rng: &mut ChaCha20Rng, n: usize) -> Vec::<u64>{
  let range = Uniform::from(0..u64::MAX);
  (0..n).map(|_| range.sample(rng)).collect::<Vec::<u64>>()
}
//------------------------------------------------------------------------------
/// This function returns the seeds of the first `n` instances of the
/// structure with seed `structure_seed`.
/// The first k seeds are the same for every n >= k.
pub fn instance_seeds(structure_seed: u64, n: usize) -> Vec::<u64>{
  let mut rng = ChaCha20Rng::seed_from_u64(structure_seed);
  draw_seeds(&mut rng, n)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// The record of an output directory.

// Seeds are written as decimal strings, since TOML integers are signed 
// 64-bit integers and the derived seeds use the full range of u64.
mod seed_text{
  use serde::{Deserialize,Deserializer,Serializer};
  pub fn serialize<S: Serializer>(seed: &u64, serializer: S) 
    -> Result<S::Ok,S::Error>
  {
    serializer.serialize_str(&seed.to_string())
  }
  pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) 
    -> Result<u64,D::Error>
  {
    let text = String::deserialize(deserializer)?;
    text.parse::<u64>().map_err(serde::de::Error::custom)
  }
}
mod seed_texts{
  use serde::{Deserialize,Deserializer,Serializer,ser::SerializeSeq};
  pub fn serialize<S: Serializer>(seeds: &[u64], serializer: S) 
    -> Result<S::Ok,S::Error>
  {
    let mut seq = serializer.serialize_seq(Some(seeds.len()))?;
    for seed in seeds.iter(){
      seq.serialize_element(&seed.to_string())?;
    }
    seq.end()
  }
  pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) 
    -> Result<Vec::<u64>,D::Error>
  {
    let texts = Vec::<String>::deserialize(deserializer)?;
    texts.iter().map(|text| text.parse::<u64>()
        .map_err(serde::de::Error::custom))
      .collect()
  }
}

/// `Manifest` is the record of an output directory (`clue_manifest.toml`).
/// `config_fingerprint` is `Config::physics_fingerprint()`;
/// `input_files` are the files read during the run, with their content 
/// hashes; `structures` holds the progress of each structure.
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
pub struct Manifest{
  pub manifest_format: u32,
  pub clue_version: String,
  pub rustc_version: String,
  #[serde(with = "seed_text")]
  pub rng_seed: u64,
  pub cluster_batch_size: usize,
  pub config_fingerprint: String,
  #[serde(default)]
  pub input_files: Vec::<InputFileRecord>,
  #[serde(default)]
  pub structures: Vec::<StructureRecord>,
}

/// `InputFileRecord` is a file read during a run and the FNV-1a hash of 
/// its contents (16 hexadecimal digits).
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
pub struct InputFileRecord{
  pub path: String,
  pub content_hash: String,
}

/// `StructureRecord` is the progress of one input structure: its id 
/// (the name of its directory under `structures/`), file, content hash,
/// weight, seed, number of instances, and the instances.
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
pub struct StructureRecord{
  pub id: String,
  pub file: String,
  pub content_hash: String,
  pub weight: f64,
  #[serde(with = "seed_text")]
  pub seed: u64,
  pub number_runs: usize,
  #[serde(default)]
  pub instances: Vec::<InstanceRecord>,
}

/// `InstanceRecord` is the progress of one instance of a structure.
/// `index` counts from 1. `orientation_seeds` are known once the instance 
/// has built its orientation grid.
/// `complete` is `true` once the instance's signal is saved; `outputs` 
/// are the output switches that were on (see `instance_outputs()`), and 
/// `files` the files of the instance outside its orientation directories,
/// relative to the instance directory 
/// (`structures/<id>/instance_<index>`).
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
pub struct InstanceRecord{
  pub index: usize,
  #[serde(with = "seed_text")]
  pub seed: u64,
  #[serde(with = "seed_texts", default)]
  pub orientation_seeds: Vec::<u64>,
  pub complete: bool,
  #[serde(default)]
  pub outputs: Vec::<String>,
  #[serde(default)]
  pub files: Vec::<String>,
  #[serde(default)]
  pub orientations: Vec::<OrientationRecord>,
}

/// `OrientationRecord` is the progress of one orientation of an instance.
/// An orientation is recorded as complete only when its directory is 
/// written (`output.orientation_signals = true`).
/// `directory` is the orientation directory, relative to the instance 
/// directory; `outputs` are the output switches that were on (see 
/// `orientation_outputs()`), and `files` the files of the orientation 
/// directory, relative to it.
#[derive(Debug,Clone,PartialEq,Default,Serialize,Deserialize)]
pub struct OrientationRecord{
  pub complete: bool,
  #[serde(default)]
  pub directory: String,
  #[serde(default)]
  pub outputs: Vec::<String>,
  #[serde(default)]
  pub files: Vec::<String>,
}

impl Manifest{
  //----------------------------------------------------------------------------
  // This function returns the record of a new simulation of `config`, 
  // which must have a seed.
  fn new(config: &Config) -> Result<Self,CluEError>{
    let Some(rng_seed) = config.rng_seed else{
      return Err(CluEError::NoRngSeed);
    };
    let Some(cluster_batch_size) = config.cluster_batch_size else{
      return Err(CluEError::NoClusterBatchSize);
    };
    Ok(Manifest{
      manifest_format: MANIFEST_FORMAT,
      clue_version: CLUE_VERSION.to_string(),
      rustc_version: RUSTC_VERSION.to_string(),
      rng_seed,
      cluster_batch_size,
      config_fingerprint: config.physics_fingerprint(),
      input_files: input_files(config)?,
      structures: Vec::<StructureRecord>::new(),
    })
  }
  //----------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
/// This function reads the record of the output directory `dir`.
/// It returns `None` if `dir` has no record.
pub fn read_manifest(dir: &str) -> Result<Option<Manifest>,CluEError>{
  let filename = manifest_path(dir);
  if !Path::new(&filename).exists(){
    return Ok(None);
  }
  let text = match std::fs::read_to_string(&filename){
    Ok(text) => text,
    Err(err) => return Err(CluEError::CannotParseManifest(filename,
          err.to_string())),
  };
  let manifest: Manifest = match toml::from_str(&text){
    Ok(manifest) => manifest,
    Err(err) => return Err(CluEError::CannotParseManifest(filename,
          err.to_string())),
  };
  if manifest.manifest_format != MANIFEST_FORMAT{
    return Err(CluEError::CannotParseManifest(filename, format!(
            "manifest_format = {} is not supported; this version of CluE \
reads manifest_format = {}", manifest.manifest_format, MANIFEST_FORMAT)));
  }
  Ok(Some(manifest))
}
//------------------------------------------------------------------------------
// This function writes `manifest` to the output directory `dir`.
fn write_manifest(dir: &str, manifest: &Manifest) -> Result<(),CluEError>{
  let filename = manifest_path(dir);
  let text = match toml::to_string(manifest){
    Ok(text) => text,
    Err(err) => return Err(CluEError::CannotWriteFile(filename,
          err.to_string())),
  };
  write_text(&filename, &text)
}
//------------------------------------------------------------------------------
// This function writes `text` to `filename` through a temporary file.
fn write_text(filename: &str, text: &str) -> Result<(),CluEError>{
  write_atomically(filename, |tmp_file| 
    match std::fs::write(tmp_file, text){
      Ok(()) => Ok(()),
      Err(err) => Err(CluEError::CannotWriteFile(filename.to_string(),
            err.to_string())),
    }
  )
}
//------------------------------------------------------------------------------
// This function returns the path of the record of the output directory `dir`.
fn manifest_path(dir: &str) -> String{
  format!("{}/{}", dir, MANIFEST_FILE)
}
//------------------------------------------------------------------------------
/// This function returns, for each structure recorded in the output 
/// directory `output_dir`, its file, its id, and its directory, which holds
/// everything computed for that structure.
pub fn structure_output_dirs(output_dir: &str) 
  -> Result<Vec::<(String,String,String)>,CluEError>
{
  let Some(manifest) = read_manifest(output_dir)? else{
    return Err(CluEError::CannotParseManifest(manifest_path(output_dir),
          "the file does not exist".to_string()));
  };
  Ok(manifest.structures.iter()
    .map(|s| (s.file.clone(), s.id.clone(), 
          format!("{}/{}/{}", output_dir, SAVE_DIR_STRUCTURES, s.id)))
    .collect())
}
//------------------------------------------------------------------------------
/// This function returns the files that a run of `config` reads besides the
/// structure file, with the hashes of their contents: the cluster file, or
/// every `clusters.toml` of a directory of orientation directories 
/// (`cluster_expansion.clustering.file`), and a partition-table csv file
/// (`cluster_expansion.clustering.partition_table`), which is read only
/// when clusters are built from the structure.
/// Other input files (`orientations.file`, a csv `detected_spin.position`)
/// are read when the input is parsed, so their contents are part of the
/// config.
pub fn input_files(config: &Config) 
  -> Result<Vec::<InputFileRecord>,CluEError>
{
  let mut paths = Vec::<String>::new();
  match &config.cluster_source{
    Some(ClusterSource::File(file)) => paths.push(file.clone()),
    Some(ClusterSource::OrientationsDir(dir)) => {
      let entries = match std::fs::read_dir(dir){
        Ok(entries) => entries,
        Err(err) => return Err(CluEError::CannotOpenFile(dir.clone(),
              err.to_string())),
      };
      let mut files = Vec::<String>::new();
      for entry in entries.flatten(){
        let clusters_file = entry.path().join("clusters.toml");
        if clusters_file.is_file(){
          files.push(clusters_file.to_string_lossy().to_string());
        }
      }
      files.sort();
      paths.append(&mut files);
    },
    Some(ClusterSource::Structure) | None => {
      if let Some(PartitionTableConfig::CSV(file)) = &config.partition_table{
        paths.push(file.clone());
      }
    },
  }

  let mut records = Vec::<InputFileRecord>::with_capacity(paths.len());
  for path in paths{
    let contents = match std::fs::read(&path){
      Ok(contents) => contents,
      Err(err) => return Err(CluEError::CannotOpenFile(path,
            err.to_string())),
    };
    records.push(InputFileRecord{
      content_hash: format!("{:016x}", fnv1a_64(&contents)),
      path,
    });
  }
  Ok(records)
}
//------------------------------------------------------------------------------
/// This function returns the output switches of `config` that write files
/// of an instance outside its orientation directories, and whether 
/// orientation directories are written.
pub fn instance_outputs(config: &Config) -> Vec::<String>{
  let on = |switch: Option<bool>| switch == Some(true);
  let mut outputs = Vec::<String>::new();
  // As in calculate_signal::write_instance_info() and 
  // calculate_signal::optionally_save_tensors().
  if on(config.write_info){
    for (switch, name) in [(config.write_bath, OUT_INFO_BATH),
      (config.write_detected_spin, OUT_INFO_DETECTED_SPIN),
      (config.write_exchange_groups, OUT_INFO_EXCHANGE_GROUPS),
      (config.write_structure_pdb, OUT_INFO_STRUCTURE_PDB),
      (config.write_tensors, OUT_INFO_TENSORS)]
    {
      if on(switch){
        outputs.push(name.to_string());
      }
    }
  }
  // The orientation grid is written without output.info.
  if on(config.write_orientation_grid){
    outputs.push(OUT_INFO_ORIENTATION_GRID.to_string());
  }
  if on(config.write_orientation_signals){
    outputs.push(OUT_ORIENTATION_SIGNALS.to_string());
  }
  outputs.sort();
  outputs
}
//------------------------------------------------------------------------------
/// This function returns the output switches of `config` that write files
/// in an orientation directory. 
/// Orientation directories exist only with 
/// `output.orientation_signals = true`; otherwise the list is empty.
pub fn orientation_outputs(config: &Config) -> Vec::<String>{
  let on = |switch: Option<bool>| switch == Some(true);
  let mut outputs = Vec::<String>::new();
  if !on(config.write_orientation_signals){
    return outputs;
  }
  // As in calculate_signal::calculate_signal_at_orientation() and
  // cluster_signals::calculate_auxiliary_signals().
  let builds_clusters = matches!(config.cluster_source, 
      Some(ClusterSource::Structure) | None);
  for (switch, name) in [(config.write_auxiliary_signals,OUT_AUXILIARY_SIGNALS),
    (config.write_clusters, OUT_CLUSTERS),
    (config.write_methyl_partitions, OUT_METHYL_PARTITIONS),
    (config.write_sans_spin_signals, OUT_SANS_SPIN_SIGNALS)]
  {
    if on(switch){
      outputs.push(name.to_string());
    }
  }
  if builds_clusters && on(config.write_partition_table){
    outputs.push(OUT_PARTITION_TABLE.to_string());
  }
  outputs.sort();
  outputs
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// The state of an existing output directory.

/// `OutputStatus` is what a run would do with its output directory.
/// * `NotSaved`: `output.directory = ""`; nothing is saved.
/// * `NewWithoutSeed`: `output.directory` and `run.rng_seed` are not set; 
///   the default directory is named after the seed that the run draws, 
///   so it is new.
/// * `New(dir)`: `dir` does not exist or is empty; a new simulation.
/// * `Complete(dir)`: `dir` holds the complete output with all requested
///   files; nothing needs to be calculated.
/// * `Resume(dir, summaries)`: `dir` holds a compatible, incomplete output;
///   `summaries` has one entry per input structure.
#[derive(Debug,Clone,PartialEq)]
pub enum OutputStatus{
  NotSaved,
  NewWithoutSeed,
  New(String),
  Complete(String),
  Resume(String,Vec::<ResumeSummary>),
}

/// `ResumeSummary` describes how a run resumes one structure.
/// `new_structure` is `true` if the output has no record of the structure.
/// Of the `number_runs` instances, `instances_complete` are complete and
/// are loaded, and `instances_added` are not recorded yet (`number_runs`
/// was increased); `orientations_complete` counts the complete 
/// orientations of the recorded instances that are not complete.
#[derive(Debug,Clone,PartialEq)]
pub struct ResumeSummary{
  pub file: String,
  pub id: String,
  pub new_structure: bool,
  pub number_runs: usize,
  pub instances_complete: usize,
  pub instances_added: usize,
  pub orientations_complete: usize,
}

// Where the output of a config goes.
enum OutputDirectory{
  NotSaved,
  Default,
  Given(String),
}
//------------------------------------------------------------------------------
// This function returns where the output of `config` goes.
fn output_directory(config: &Config) -> OutputDirectory{
  match &config.output_directory{
    Some(dir) if dir.is_empty() => OutputDirectory::NotSaved,
    Some(dir) => OutputDirectory::Given(dir.clone()),
    None => OutputDirectory::Default,
  }
}
//------------------------------------------------------------------------------
/// This function returns the default output directory of `config` with the
/// main seed `main_seed`: `CluE-<FNV-1a hash of the fingerprint and seed>`.
/// It does not depend on the output settings, `run.number_runs`, the 
/// weights, or the structure files, so that the allowed changes of an 
/// output (more instances, other output files) stay in one directory.
pub fn default_directory(config: &Config, main_seed: u64) -> String{
  let text = format!("{}\nrng_seed = {}", config.physics_fingerprint(), 
      main_seed);
  format!("CluE-{:016x}", fnv1a_64(text.as_bytes()))
}
//------------------------------------------------------------------------------
/// This function returns what a run of `config` would do with its output
/// directory, or the error that stops the run because the directory holds
/// the output of a different simulation.
/// It only reads files: it never creates, writes, renames, or deletes a file
/// or directory, so `clue_oxide --check` can call it.
/// `config` must have its time axis set (`Config::set_tau_axis()`).
pub fn inspect_output(config: &Config) -> Result<OutputStatus,CluEError>{
  let dir = match output_directory(config){
    OutputDirectory::NotSaved => return Ok(OutputStatus::NotSaved),
    OutputDirectory::Default => match config.rng_seed{
      Some(seed) => default_directory(config, seed),
      None => return Ok(OutputStatus::NewWithoutSeed),
    },
    OutputDirectory::Given(dir) => dir,
  };

  let path = Path::new(&dir);
  if !path.exists(){
    return Ok(OutputStatus::New(dir));
  }
  if !path.is_dir(){
    return Err(CluEError::CannotCreateDir(dir, 
          "a file with this name exists".to_string()));
  }
  let entries = list_directory(&dir)?;
  if entries.is_empty(){
    return Ok(OutputStatus::New(dir));
  }

  let Some(manifest) = read_manifest(&dir)? else{
    const MAX_LISTED: usize = 10;
    let mut listing = entries.iter().take(MAX_LISTED).cloned()
      .collect::<Vec::<String>>();
    if entries.len() > MAX_LISTED{
      listing.push(format!("and {} more", entries.len() - MAX_LISTED));
    }
    return Err(CluEError::UnrecognizedOutput(dir, listing));
  };

  check_compatible(&manifest, config, &dir)?;

  let (complete, summaries) = completeness(&manifest, config, &dir)?;
  if complete{
    return Ok(OutputStatus::Complete(dir));
  }

  // Saved auxiliary signals are loaded by batch position.
  let Some(batch_size) = config.cluster_batch_size else{
    return Err(CluEError::NoClusterBatchSize);
  };
  if batch_size != manifest.cluster_batch_size{
    return Err(CluEError::OutputBatchSizeMismatch(dir, 
          manifest.cluster_batch_size, batch_size));
  }

  Ok(OutputStatus::Resume(dir, summaries))
}
//------------------------------------------------------------------------------
// This function returns the error that stops a run of `config` into the 
// output directory `dir` with the record `manifest`, if the output belongs
// to a different simulation.
// The checks follow the table in report_plan_checkpoints.md (§E).
fn check_compatible(manifest: &Manifest, config: &Config, dir: &str)
  -> Result<(),CluEError>
{
  if manifest.clue_version != CLUE_VERSION{
    return Err(CluEError::OutputVersionMismatch(dir.to_string(),
          manifest.clue_version.clone(), CLUE_VERSION.to_string()));
  }

  let mut differences = fingerprint_differences(&manifest.config_fingerprint,
      &config.physics_fingerprint());
  differences.append(&mut input_file_differences(&manifest.input_files,
        &input_files(config)?));
  if !differences.is_empty(){
    let rustc = if manifest.rustc_version != RUSTC_VERSION{
      Some((manifest.rustc_version.clone(), RUSTC_VERSION.to_string()))
    }else{
      None
    };
    return Err(CluEError::OutputConfigMismatch(dir.to_string(), differences,
          rustc));
  }

  let Some(main_seed) = config.rng_seed else{
    return Err(CluEError::OutputNeedsSeed(dir.to_string(), 
          manifest.rng_seed));
  };
  if main_seed != manifest.rng_seed{
    return Err(CluEError::OutputSeedMismatch(dir.to_string(), 
          manifest.rng_seed, main_seed));
  }

  let structures = config.structure_inputs()?;
  let Some(number_runs) = config.number_runs else{
    return Err(CluEError::NoNumberSystemInstances);
  };
  for record in manifest.structures.iter(){
    let Some(structure) = structures.iter().find(|s| s.id == record.id) else{
      return Err(CluEError::OutputStructureRemoved(dir.to_string(),
            record.file.clone(), record.id.clone()));
    };
    if number_runs < record.number_runs{
      return Err(CluEError::OutputFewerRuns(dir.to_string(),
            record.file.clone(), record.number_runs, number_runs));
    }
    check_recorded_seeds(record, structure, main_seed, dir)?;
  }
  Ok(())
}
//------------------------------------------------------------------------------
// This function checks that the seeds recorded for `structure` are the 
// seeds derived from `main_seed`.
// The orientation seeds are drawn only after the structure is built, so
// they are checked by the run (`OutputState::set_orientation_seeds()`).
fn check_recorded_seeds(record: &StructureRecord, structure: &StructureInput,
    main_seed: u64, dir: &str) -> Result<(),CluEError>
{
  let seed = structure_seed(main_seed, structure.content_hash);
  if seed != record.seed{
    return Err(CluEError::OutputSeedSequenceChanged(dir.to_string(),
          format!("structure \"{}\": recorded seed {}, derived seed {}",
            record.file, record.seed, seed)));
  }
  let seeds = instance_seeds(seed, record.instances.len());
  for (instance, derived) in record.instances.iter().zip(seeds.iter()){
    if instance.seed != *derived{
      return Err(CluEError::OutputSeedSequenceChanged(dir.to_string(),
            format!("structure \"{}\", instance {}: recorded seed {}, \
derived seed {}", record.file, instance.index, instance.seed, derived)));
    }
  }
  Ok(())
}
//------------------------------------------------------------------------------
// This function decides whether the output in `dir` with the record 
// `manifest` is complete for `config`, and summarizes, for each structure,
// what a run would load and compute.
// `manifest` must be compatible with `config` (`check_compatible()`).
fn completeness(manifest: &Manifest, config: &Config, dir: &str) 
  -> Result<(bool,Vec::<ResumeSummary>),CluEError>
{
  let Some(number_runs) = config.number_runs else{
    return Err(CluEError::NoNumberSystemInstances);
  };
  let inst_outputs = instance_outputs(config);
  let ori_outputs = orientation_outputs(config);

  let mut complete = [RNG_SEED_FILE, TAU_AXIS_FILE, SIGNAL_FILE, 
    STRUCTURES_FILE].iter()
    .all(|file| Path::new(dir).join(file).is_file());

  let mut summaries = Vec::<ResumeSummary>::new();
  for structure in config.structure_inputs()?.iter(){
    let Some(record) = manifest.structures.iter()
      .find(|r| r.id == structure.id) else
    {
      complete = false;
      summaries.push(ResumeSummary{
        file: structure.file.clone(),
        id: structure.id.clone(),
        new_structure: true,
        number_runs,
        instances_complete: 0,
        instances_added: number_runs,
        orientations_complete: 0,
      });
      continue;
    };

    let mut instances_complete = 0;
    let mut orientations_complete = 0;
    for instance in record.instances.iter().take(number_runs){
      let instance_dir = Path::new(dir).join(SAVE_DIR_STRUCTURES)
        .join(&record.id).join(format!("instance_{}", instance.index));
      if instance_is_complete(instance, &inst_outputs, &ori_outputs, 
          &instance_dir)
      {
        instances_complete += 1;
      }else{
        orientations_complete += instance.orientations.iter()
          .filter(|o| orientation_is_complete(o, &ori_outputs, &instance_dir))
          .count();
      }
    }
    let instances_added = number_runs.saturating_sub(record.instances.len());

    if instances_complete < number_runs 
      || record.weight != structure.weight
      || !Path::new(dir).join(SAVE_DIR_STRUCTURES).join(&record.id)
        .join(SIGNAL_FILE).is_file()
    {
      complete = false;
    }

    summaries.push(ResumeSummary{
      file: structure.file.clone(),
      id: structure.id.clone(),
      new_structure: false,
      number_runs,
      instances_complete,
      instances_added,
      orientations_complete,
    });
  }
  Ok((complete, summaries))
}
//------------------------------------------------------------------------------
// This function returns `true` if `instance`, in `instance_dir`, is 
// complete, was written with at least the outputs `inst_outputs`, its 
// files are present, and, if orientation directories are requested, each
// orientation is complete with at least the outputs `ori_outputs`.
fn instance_is_complete(instance: &InstanceRecord, inst_outputs: &[String],
    ori_outputs: &[String], instance_dir: &Path) -> bool
{
  instance.complete
    && is_subset(inst_outputs, &instance.outputs)
    && files_exist(&instance.files, instance_dir)
    && (!inst_outputs.iter().any(|o| o == OUT_ORIENTATION_SIGNALS)
        || (instance.orientations.len() == instance.orientation_seeds.len()
          && instance.orientations.iter()
          .all(|o| orientation_is_complete(o, ori_outputs, instance_dir))))
}
//------------------------------------------------------------------------------
// This function returns `true` if `orientation`, of the instance in 
// `instance_dir`, is complete, was written with at least the outputs 
// `ori_outputs`, and its files are present.
fn orientation_is_complete(orientation: &OrientationRecord, 
    ori_outputs: &[String], instance_dir: &Path) -> bool
{
  orientation.complete
    && is_subset(ori_outputs, &orientation.outputs)
    && files_exist(&orientation.files, &instance_dir.join(
          &orientation.directory))
}
//------------------------------------------------------------------------------
// This function returns `true` if every element of `subset` is in `set`.
fn is_subset(subset: &[String], set: &[String]) -> bool{
  subset.iter().all(|item| set.contains(item))
}
//------------------------------------------------------------------------------
// This function returns `true` if every file of `files`, relative to `dir`,
// exists.
fn files_exist(files: &[String], dir: &Path) -> bool{
  files.iter().all(|file| dir.join(file).is_file())
}
//------------------------------------------------------------------------------
// This function returns the sorted names of the entries of the directory 
// `dir`.
fn list_directory(dir: &str) -> Result<Vec::<String>,CluEError>{
  let entries = match std::fs::read_dir(dir){
    Ok(entries) => entries,
    Err(err) => return Err(CluEError::CannotOpenFile(dir.to_string(),
          err.to_string())),
  };
  let mut names = entries.flatten()
    .map(|entry| entry.file_name().to_string_lossy().to_string())
    .collect::<Vec::<String>>();
  names.sort();
  Ok(names)
}
//------------------------------------------------------------------------------
// This function returns the files under the directory `dir`, as sorted 
// paths relative to `dir`, without the subdirectory `exclude` of `dir` and
// without temporary `.tmp` files.
fn list_files(dir: &str, exclude: Option<&str>)
  -> Result<Vec::<String>,CluEError>
{
  let root = Path::new(dir);
  let mut files = Vec::<String>::new();
  let mut stack = vec![root.to_path_buf()];
  let excluded = exclude.map(|e| root.join(e));
  while let Some(path) = stack.pop(){
    let entries = match std::fs::read_dir(&path){
      Ok(entries) => entries,
      Err(err) => return Err(CluEError::CannotOpenFile(
            path.to_string_lossy().to_string(), err.to_string())),
    };
    for entry in entries.flatten(){
      let entry_path = entry.path();
      if entry_path.is_dir(){
        if Some(&entry_path) != excluded.as_ref(){
          stack.push(entry_path);
        }
      }else if entry_path.extension().map_or(true, |ext| ext != "tmp"){
        if let Ok(relative) = entry_path.strip_prefix(root){
          files.push(relative.to_string_lossy().replace('\\', "/"));
        }
      }
    }
  }
  files.sort();
  Ok(files)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// Differences between configurations.

//------------------------------------------------------------------------------
// This function lists the differences between the recorded fingerprint 
// `old` and the fingerprint `new` (`Config::physics_fingerprint()`), 
// line by line, each group of changed lines under the path of the field 
// that contains it.
fn fingerprint_differences(old: &str, new: &str) -> Vec::<String>{
  if old == new{
    return Vec::<String>::new();
  }
  let a = old.lines().collect::<Vec::<&str>>();
  let b = new.lines().collect::<Vec::<&str>>();

  // Lines shared at the start and at the end.
  let mut start = 0;
  while start < a.len() && start < b.len() && a[start] == b[start]{
    start += 1;
  }
  let mut end = 0;
  while end < a.len() - start && end < b.len() - start 
    && a[a.len() - 1 - end] == b[b.len() - 1 - end]
  {
    end += 1;
  }
  let a_mid = &a[start..a.len() - end];
  let b_mid = &b[start..b.len() - end];

  // Edit script of the middle: (old line, new line) indices, with `None`
  // for a removed or added line.
  let script = edit_script(a_mid, b_mid);

  let mut differences = Vec::<String>::new();
  let mut ii = 0;
  while ii < script.len(){
    if let (Some(_), Some(_)) = script[ii]{
      ii += 1;
      continue;
    }
    // A group of changed lines, named after the field of its first line.
    let path = match script[ii]{
      (Some(i), _) => field_path(&a, start + i),
      (None, Some(j)) => field_path(&b, start + j),
      (None, None) => String::new(),
    };
    if !path.is_empty(){
      differences.push(format!("{}:", path));
    }
    let mut removed = Vec::<String>::new();
    let mut added = Vec::<String>::new();
    while ii < script.len(){
      match script[ii]{
        (Some(_), Some(_)) => break,
        (Some(i), None) => removed.push(format!("  - {}", a_mid[i].trim())),
        (None, Some(j)) => added.push(format!("  + {}", b_mid[j].trim())),
        (None, None) => (),
      }
      ii += 1;
    }
    differences.append(&mut removed);
    differences.append(&mut added);
  }

  if differences.len() > MAX_LISTED_DIFFERENCES{
    let n_more = differences.len() - MAX_LISTED_DIFFERENCES;
    differences.truncate(MAX_LISTED_DIFFERENCES);
    differences.push(format!("... and {} more lines", n_more));
  }
  differences
}
//------------------------------------------------------------------------------
// This function returns a shortest edit script from the lines `a` to the
// lines `b`, from their longest common subsequence: `(Some(i), Some(j))` 
// keeps `a[i] == b[j]`, `(Some(i), None)` removes `a[i]`, and 
// `(None, Some(j))` adds `b[j]`.
// For long inputs, all of `a` is replaced by all of `b`.
fn edit_script(a: &[&str], b: &[&str]) 
  -> Vec::<(Option<usize>,Option<usize>)>
{
  const MAX_TABLE_SIZE: usize = 4_000_000;
  let (n, m) = (a.len(), b.len());
  if (n + 1)*(m + 1) > MAX_TABLE_SIZE{
    return (0..n).map(|i| (Some(i), None))
      .chain((0..m).map(|j| (None, Some(j))))
      .collect();
  }

  // lcs[i][j] is the length of the longest common subsequence of a[i..]
  // and b[j..].
  let mut lcs = vec![vec![0usize; m + 1]; n + 1];
  for i in (0..n).rev(){
    for j in (0..m).rev(){
      lcs[i][j] = if a[i] == b[j]{
        lcs[i + 1][j + 1] + 1
      }else{
        std::cmp::max(lcs[i + 1][j], lcs[i][j + 1])
      };
    }
  }

  let mut script = Vec::<(Option<usize>,Option<usize>)>::new();
  let (mut i, mut j) = (0, 0);
  while i < n && j < m{
    if a[i] == b[j]{
      script.push((Some(i), Some(j)));
      i += 1;
      j += 1;
    }else if lcs[i + 1][j] >= lcs[i][j + 1]{
      script.push((Some(i), None));
      i += 1;
    }else{
      script.push((None, Some(j)));
      j += 1;
    }
  }
  script.extend((i..n).map(|i| (Some(i), None)));
  script.extend((j..m).map(|j| (None, Some(j))));
  script
}
//------------------------------------------------------------------------------
// This function returns the path of the field that contains line `index` 
// of a pretty-printed `Debug` text, such as 
// `particles > properties > isotope_properties`: the names of the 
// enclosing blocks that are fields (`name: ...`).
fn field_path(lines: &[&str], index: usize) -> String{
  let indent = |line: &str| line.len() - line.trim_start().len();
  let mut path = Vec::<String>::new();
  let mut level = indent(lines[index]);
  for line in lines[..index].iter().rev(){
    let d = indent(line);
    if d < level{
      level = d;
      let text = line.trim();
      if let Some((name, _)) = text.split_once(": "){
        path.push(name.to_string());
      }
      if d == 0{
        break;
      }
    }
  }
  path.reverse();
  path.join(" > ")
}
//------------------------------------------------------------------------------
// This function lists the differences between the recorded input files 
// `old` and the input files `new`.
fn input_file_differences(old: &[InputFileRecord], new: &[InputFileRecord])
  -> Vec::<String>
{
  let mut differences = Vec::<String>::new();
  for record in old.iter(){
    match new.iter().find(|r| r.path == record.path){
      None => differences.push(format!("input file \"{}\": no longer read",
            record.path)),
      Some(r) if r.content_hash != record.content_hash 
        => differences.push(format!("input file \"{}\": contents changed",
              record.path)),
      Some(_) => (),
    }
  }
  for record in new.iter(){
    if !old.iter().any(|r| r.path == record.path){
      differences.push(format!("input file \"{}\": not read by the recorded \
run", record.path));
    }
  }
  differences
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `OutputState` is the output directory of a run and its record, which
/// the run updates as structures, instances, and orientations finish.
/// Without an output directory (`output.directory = ""`), the record is 
/// kept in memory only.
#[derive(Debug,Clone)]
pub struct OutputState{
  dir: Option<String>,
  manifest: Manifest,
  instance_outputs: Vec::<String>,
  orientation_outputs: Vec::<String>,
}

impl OutputState{
  //----------------------------------------------------------------------------
  /// This function starts the output of a run of `config` with the status
  /// `status` (from `inspect_output()`). `config` must have a seed.
  /// For a new output, it creates the directory; for a new or resumed 
  /// output, it writes `rng_seed.txt` and the record.
  pub fn start(config: &Config, status: &OutputStatus) 
    -> Result<Self,CluEError>
  {
    let Some(main_seed) = config.rng_seed else{
      return Err(CluEError::NoRngSeed);
    };
    let (dir, manifest) = match status{
      OutputStatus::NotSaved => (None, Manifest::new(config)?),
      // inspect_output() names the directory once the run has its seed.
      OutputStatus::NewWithoutSeed => (Some(default_directory(config, 
              main_seed)), Manifest::new(config)?),
      OutputStatus::New(dir) => (Some(dir.clone()), Manifest::new(config)?),
      OutputStatus::Complete(dir) | OutputStatus::Resume(dir,_) => {
        let Some(mut manifest) = read_manifest(dir)? else{
          return Err(CluEError::CannotParseManifest(manifest_path(dir),
                "the file does not exist".to_string()));
        };
        // The configuration is unchanged (checked by inspect_output()); 
        // the record names the latest compiler and the current input files.
        manifest.rustc_version = RUSTC_VERSION.to_string();
        manifest.input_files = input_files(config)?;
        (Some(dir.clone()), manifest)
      },
    };

    let state = OutputState{
      dir,
      manifest,
      instance_outputs: instance_outputs(config),
      orientation_outputs: orientation_outputs(config),
    };

    if let Some(dir) = &state.dir{
      create_dir(dir)?;
      write_text(&format!("{}/{}", dir, RNG_SEED_FILE), 
          &format!("{}\n", main_seed))?;
      state.save()?;
    }
    Ok(state)
  }
  //----------------------------------------------------------------------------
  /// This function returns the output directory, if the output is saved.
  pub fn directory(&self) -> Option<&str>{
    self.dir.as_deref()
  }
  //----------------------------------------------------------------------------
  /// This function returns the record.
  pub fn manifest(&self) -> &Manifest{
    &self.manifest
  }
  //----------------------------------------------------------------------------
  // This function writes the record, if the output is saved.
  fn save(&self) -> Result<(),CluEError>{
    match &self.dir{
      Some(dir) => write_manifest(dir, &self.manifest),
      None => Ok(()),
    }
  }
  //----------------------------------------------------------------------------
  /// This function starts `structure`, with seed `seed` and instance seeds
  /// `instance_seeds`, and returns the index of its record.
  /// A new structure, and instances added by a larger `number_runs`, are
  /// added to the record; recorded seeds must equal the given seeds.
  /// The structure directory is created.
  pub fn begin_structure(&mut self, structure: &StructureInput, seed: u64,
      instance_seeds: &[u64]) -> Result<usize,CluEError>
  {
    let s_idx = match self.manifest.structures.iter()
      .position(|r| r.id == structure.id)
    {
      Some(s_idx) => s_idx,
      None => {
        self.manifest.structures.push(StructureRecord{
          id: structure.id.clone(),
          file: structure.file.clone(),
          content_hash: structure.content_hash_hex(),
          weight: structure.weight,
          seed,
          number_runs: 0,
          instances: Vec::<InstanceRecord>::new(),
        });
        self.manifest.structures.len() - 1
      },
    };

    let dir = self.dir.clone().unwrap_or_default();
    let record = &mut self.manifest.structures[s_idx];
    if record.seed != seed{
      return Err(CluEError::OutputSeedSequenceChanged(dir,
            format!("structure \"{}\": recorded seed {}, derived seed {}",
              record.file, record.seed, seed)));
    }
    for (instance, &derived) in record.instances.iter()
      .zip(instance_seeds.iter())
    {
      if instance.seed != derived{
        return Err(CluEError::OutputSeedSequenceChanged(dir,
              format!("structure \"{}\", instance {}: recorded seed {}, \
derived seed {}", record.file, instance.index, instance.seed, derived)));
      }
    }
    for (k, &instance_seed) in instance_seeds.iter().enumerate()
      .skip(record.instances.len())
    {
      record.instances.push(InstanceRecord{
        index: k + 1,
        seed: instance_seed,
        orientation_seeds: Vec::<u64>::new(),
        complete: false,
        outputs: Vec::<String>::new(),
        files: Vec::<String>::new(),
        orientations: Vec::<OrientationRecord>::new(),
      });
    }
    record.file = structure.file.clone();
    record.weight = structure.weight;
    record.number_runs = std::cmp::max(record.number_runs, 
        instance_seeds.len());

    if let Some(structure_dir) = self.structure_dir(s_idx){
      create_dir(&structure_dir)?;
    }
    self.save()?;
    Ok(s_idx)
  }
  //----------------------------------------------------------------------------
  /// This function returns the directory of structure `s_idx`, if the 
  /// output is saved.
  pub fn structure_dir(&self, s_idx: usize) -> Option<String>{
    let dir = self.dir.as_ref()?;
    Some(format!("{}/{}/{}", dir, SAVE_DIR_STRUCTURES, 
          self.manifest.structures[s_idx].id))
  }
  //----------------------------------------------------------------------------
  /// This function returns the directory of instance `k` (from 0) of
  /// structure `s_idx`, if the output is saved.
  pub fn instance_dir(&self, s_idx: usize, k: usize) -> Option<String>{
    Some(format!("{}/instance_{}", self.structure_dir(s_idx)?, k + 1))
  }
  //----------------------------------------------------------------------------
  /// This function returns `true` if instance `k` of structure `s_idx` is
  /// complete with all requested files, so that its signal can be loaded.
  pub fn instance_is_complete(&self, s_idx: usize, k: usize) -> bool{
    let Some(instance_dir) = self.instance_dir(s_idx, k) else{
      return false;
    };
    match self.manifest.structures[s_idx].instances.get(k){
      Some(instance) => instance_is_complete(instance, 
          &self.instance_outputs, &self.orientation_outputs, 
          Path::new(&instance_dir)),
      None => false,
    }
  }
  //----------------------------------------------------------------------------
  /// This function marks instance `k` of structure `s_idx` as not complete,
  /// before it is computed.
  pub fn begin_instance(&mut self, s_idx: usize, k: usize) 
    -> Result<(),CluEError>
  {
    let instance = &mut self.manifest.structures[s_idx].instances[k];
    if instance.complete{
      instance.complete = false;
      self.save()?;
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function records the orientation seeds of instance `k` of 
  /// structure `s_idx`. Recorded seeds must equal `seeds`.
  pub fn set_orientation_seeds(&mut self, s_idx: usize, k: usize, 
      seeds: &[u64]) -> Result<(),CluEError>
  {
    let dir = self.dir.clone().unwrap_or_default();
    let record = &mut self.manifest.structures[s_idx];
    let instance = &mut record.instances[k];
    if !instance.orientation_seeds.is_empty() 
      && instance.orientation_seeds != seeds
    {
      return Err(CluEError::OutputSeedSequenceChanged(dir,
            format!("structure \"{}\", instance {}: the orientation seeds \
differ", record.file, instance.index)));
    }
    instance.orientation_seeds = seeds.to_vec();
    instance.orientations.resize(seeds.len(), OrientationRecord::default());
    self.save()
  }
  //----------------------------------------------------------------------------
  /// This function returns `true` if orientation `i` of instance `k` of
  /// structure `s_idx` is complete with all requested files, so that its 
  /// signal can be loaded.
  pub fn orientation_is_complete(&self, s_idx: usize, k: usize, i: usize)
    -> bool
  {
    let Some(instance_dir) = self.instance_dir(s_idx, k) else{
      return false;
    };
    match self.manifest.structures[s_idx].instances[k].orientations.get(i){
      Some(orientation) => orientation_is_complete(orientation,
          &self.orientation_outputs, Path::new(&instance_dir)),
      None => false,
    }
  }
  //----------------------------------------------------------------------------
  /// This function marks orientation `i` of instance `k` of structure 
  /// `s_idx` as not complete, before it is computed.
  pub fn begin_orientation(&mut self, s_idx: usize, k: usize, i: usize) 
    -> Result<(),CluEError>
  {
    let orientation 
      = &mut self.manifest.structures[s_idx].instances[k].orientations[i];
    if orientation.complete{
      orientation.complete = false;
      self.save()?;
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function records orientation `i` of instance `k` of structure 
  /// `s_idx` as complete with the files in `orientation_dir`. 
  /// Without an orientation directory, nothing is recorded, and the 
  /// orientation is computed again if the instance is resumed.
  pub fn finish_orientation(&mut self, s_idx: usize, k: usize, i: usize,
      orientation_dir: &Option<String>) -> Result<(),CluEError>
  {
    let (Some(instance_dir), Some(orientation_dir)) 
      = (self.instance_dir(s_idx, k), orientation_dir) 
    else{
      return Ok(());
    };
    let files = list_files(orientation_dir, None)?;
    let orientation 
      = &mut self.manifest.structures[s_idx].instances[k].orientations[i];
    orientation.directory = relative_to(orientation_dir, &instance_dir);
    orientation.files = files;
    orientation.outputs = self.orientation_outputs.clone();
    orientation.complete = true;
    self.save()
  }
  //----------------------------------------------------------------------------
  /// This function records instance `k` of structure `s_idx` as complete,
  /// with the files of its directory outside its orientation directories.
  pub fn finish_instance(&mut self, s_idx: usize, k: usize) 
    -> Result<(),CluEError>
  {
    let files = match self.instance_dir(s_idx, k){
      Some(instance_dir) => list_files(&instance_dir, 
          Some(SAVE_DIR_ORIENTATION_SIGNALS))?,
      None => Vec::<String>::new(),
    };
    let instance = &mut self.manifest.structures[s_idx].instances[k];
    instance.files = files;
    instance.outputs = self.instance_outputs.clone();
    instance.complete = true;
    self.save()
  }
  //----------------------------------------------------------------------------
  /// This function writes the index of the structures (`structures.csv`),
  /// if the output is saved.
  pub fn write_structures_index(&self) -> Result<(),CluEError>{
    let Some(dir) = &self.dir else{
      return Ok(());
    };
    let records = &self.manifest.structures;
    let column = |f: &dyn Fn(&StructureRecord) -> String| 
      records.iter().map(f).collect::<Vec::<String>>();
    let data = vec![
      column(&|r| r.id.clone()),
      column(&|r| r.file.clone()),
      column(&|r| r.content_hash.clone()),
      column(&|r| r.weight.to_string()),
      column(&|r| r.number_runs.to_string()),
      column(&|r| format!("{}/{}", SAVE_DIR_STRUCTURES, r.id)),
    ];
    let headers = ["id", "file", "content_hash", "weight", "number_runs",
      "directory"].iter().map(|h| h.to_string()).collect::<Vec::<String>>();
    crate::io::write_data(&data, &format!("{}/{}", dir, STRUCTURES_FILE), 
        headers)
  }
  //----------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
// This function creates the directory `dir` and its parents.
fn create_dir(dir: &str) -> Result<(),CluEError>{
  match std::fs::create_dir_all(dir){
    Ok(()) => Ok(()),
    Err(err) => Err(CluEError::CannotCreateDir(dir.to_string(),
          err.to_string())),
  }
}
//------------------------------------------------------------------------------
// This function returns `path`, which starts with `dir/`, relative to `dir`.
fn relative_to(path: &str, dir: &str) -> String{
  path.strip_prefix(dir).unwrap_or(path).trim_start_matches('/').to_string()
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#[cfg(test)]
mod tests{
  use super::*;
  use crate::config::command_line_input::CommandLineInput;
  use crate::io::FromTOMLString;
  use crate::physical_constants::ONE;
  use crate::signal::{Signal,load_csv_to_vec_signals};

  const ASSET: &str = "assets/1omp_K26R1_0.003228966616048703.toml";

  // Random bath populations differ between instances and orientations.
  const RANDOM: [(&str,&str); 2] = [("ensemble = true","ensemble = false"),
    ("populations = \"thermal\"","populations = \"random\"")];
  const AUX: (&str,&str) 
    = ("auxiliary_signals = false","auxiliary_signals = true");
  // Two orientations with equal weights.
  const TWO_ORIENTATIONS: [(&str,&str); 2] = [
    ("grid = \"vector\"","grid = \"vector_grid\""),
    ("vector = [0.67436015, 0.67436015, 0.3007936 ]",
     "vector_grid = [[0.67436015, 0.67436015, 0.3007936, 0.5], \
[0.0, 0.0, 1.0, 0.5]]"),
  ];
  //----------------------------------------------------------------------------
  // This function returns a path for test `name` under the system temp 
  // directory; nothing exists at the path.
  fn temp_dir(name: &str) -> String{
    let dir = std::env::temp_dir().join(format!(
          "clue_oxide_test_output_{}_{}", name, std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    dir.to_str().unwrap().to_string()
  }
  //----------------------------------------------------------------------------
  // This function returns the 1omp test input with a radius of 8 Å and
  // 11 time points, writing to `dir`, with the replacements `changes` and
  // the keys `run` in the [run] table.
  fn small_toml(dir: &str, run: &str, changes: &[(&str,&str)]) -> String{
    let mut toml = std::fs::read_to_string(ASSET).unwrap()
      .replace("radius = 300","radius = 8")
      .replace("number_timepoints = [301]","number_timepoints = [11]")
      .replace("[output]", &format!("[output]\n  directory = '{}'", dir));
    for (old, new) in changes.iter(){
      assert!(toml.contains(old), "{}", old);
      toml = toml.replace(old, new);
    }
    format!("{}\n[run]\n{}\n", toml, run)
  }
  //----------------------------------------------------------------------------
  fn small_config(dir: &str, run: &str, changes: &[(&str,&str)]) -> Config{
    let mut config = Config::from_toml_string(&small_toml(dir, run, changes))
      .unwrap();
    config.set_tau_axis().unwrap();
    config
  }
  //----------------------------------------------------------------------------
  fn read_signals(filename: &str) -> Vec::<Signal>{
    load_csv_to_vec_signals(filename).unwrap()
  }
  //----------------------------------------------------------------------------
  // This function returns the directory of instance `k` (from 1) of the 
  // only structure of the output in `out`.
  fn instance_dir(out: &str, k: usize) -> String{
    let structures = structure_output_dirs(out).unwrap();
    assert_eq!(structures.len(), 1);
    format!("{}/instance_{}", structures[0].2, k)
  }
  //----------------------------------------------------------------------------
  // This function returns the sorted orientation directories of an 
  // instance directory.
  fn orientation_dirs(instance_dir: &str) -> Vec::<String>{
    let mut dirs = std::fs::read_dir(format!("{}/{}", instance_dir, 
          SAVE_DIR_ORIENTATION_SIGNALS)).unwrap()
      .map(|e| e.unwrap().path().to_str().unwrap().to_string())
      .collect::<Vec::<String>>();
    dirs.sort();
    dirs
  }
  //----------------------------------------------------------------------------
  // This function returns every file under `dir` with its contents.
  fn tree_snapshot(dir: &str) -> Vec::<(String,Vec::<u8>)>{
    let mut files = Vec::<(String,Vec::<u8>)>::new();
    if !Path::new(dir).exists(){
      return files;
    }
    for file in list_files(dir, None).unwrap(){
      let contents = std::fs::read(Path::new(dir).join(&file)).unwrap();
      files.push((file, contents));
    }
    files
  }
  //----------------------------------------------------------------------------
  // This function edits the record of the output in `out`.
  fn edit_manifest<F: FnOnce(&mut Manifest)>(out: &str, edit: F){
    let mut manifest = read_manifest(out).unwrap().unwrap();
    edit(&mut manifest);
    write_manifest(out, &manifest).unwrap();
  }
  //----------------------------------------------------------------------------
  // Plan tests 1, 2, and 3: separate instance directories, a complete 
  // rerun, and resuming a deleted orientation.
  #[test]
  fn test_instances_complete_rerun_and_resume(){
    let out = temp_dir("instances");
    let config = small_config(&out, 
        "number_runs = 2\nrng_seed = 11\ncluster_batch_size = 4", 
        &[RANDOM[0], RANDOM[1], AUX]);
    let (_, signal) = crate::run(config.clone()).unwrap();

    // 1. Each instance has its own directory and checkpoints, and the 
    // average is over both instances (N04).
    let (dir_1, dir_2) = (instance_dir(&out, 1), instance_dir(&out, 2));
    let s1 = read_signals(&format!("{}/{}", dir_1, SIGNAL_FILE));
    let s2 = read_signals(&format!("{}/{}", dir_2, SIGNAL_FILE));
    assert_ne!(s1, s2);
    let batch = |dir: &str| std::fs::read_to_string(format!(
          "{}/auxiliary_signals/cluster_size_1_batch_0.csv", 
          orientation_dirs(dir)[0])).unwrap();
    assert_ne!(batch(&dir_1), batch(&dir_2));

    let top = read_signals(&format!("{}/{}", out, SIGNAL_FILE));
    for (order, top_signal) in top.iter().enumerate(){
      let n = s1[order].len();
      let expected = &(&Signal::zeros(n) + &s1[order]) + &s2[order];
      assert_eq!(*top_signal, expected.scale(ONE/(2_f64)));
    }
    assert_ne!(top.last(), s1.last());
    assert_eq!(signal, top.last().unwrap().data);
    assert!(Path::new(&format!("{}/{}", out, RNG_SEED_FILE)).is_file());
    assert_eq!(std::fs::read_to_string(format!("{}/{}", out, RNG_SEED_FILE))
        .unwrap(), "11\n");

    // 2. A complete output is loaded, without calculations.
    assert_eq!(inspect_output(&config), Ok(OutputStatus::Complete(out.clone())));
    let before = tree_snapshot(&out);
    let (_, signal_again) = crate::run(config.clone()).unwrap();
    assert_eq!(signal, signal_again);
    assert_eq!(before, tree_snapshot(&out));

    // 3. Without the orientation signal of instance 2, instance 1 is loaded
    // and instance 2 is resumed from its checkpoints, with the same result.
    let ori_signal = format!("{}/{}", orientation_dirs(&dir_2)[0], SIGNAL_FILE);
    std::fs::remove_file(&ori_signal).unwrap();
    match inspect_output(&config){
      Ok(OutputStatus::Resume(dir, summaries)) => {
        assert_eq!(dir, out);
        assert_eq!(summaries.len(), 1);
        assert_eq!(summaries[0].instances_complete, 1);
        assert_eq!(summaries[0].instances_added, 0);
        assert_eq!(summaries[0].orientations_complete, 0);
        assert!(!summaries[0].new_structure);
      },
      other => panic!("expected Resume, got {:?}", other),
    }
    let (_, signal_resumed) = crate::run(config.clone()).unwrap();
    assert_eq!(signal, signal_resumed);
    assert!(Path::new(&ori_signal).is_file());
    assert_eq!(top, read_signals(&format!("{}/{}", out, SIGNAL_FILE)));

    let _ = std::fs::remove_dir_all(&out);
  }
  //----------------------------------------------------------------------------
  // Plan test 4: an entropy seed is saved, a rerun without the seed stops,
  // and the saved seed resumes.
  #[test]
  fn test_entropy_seed_is_saved_and_needed_to_resume(){
    let out = temp_dir("entropy_seed");
    let mut config = small_config(&out, "number_runs = 1", &RANDOM);
    assert_eq!(config.rng_seed, None);
    assert_eq!(inspect_output(&config), Ok(OutputStatus::New(out.clone())));

    let (_, signal) = crate::run(config.clone()).unwrap();
    let seed = std::fs::read_to_string(format!("{}/{}", out, RNG_SEED_FILE))
      .unwrap().trim().parse::<u64>().unwrap();
    assert!(seed <= i64::MAX as u64);

    let err = inspect_output(&config).unwrap_err();
    assert_eq!(err, CluEError::OutputNeedsSeed(out.clone(), seed));
    assert_eq!(crate::run(config.clone()).unwrap_err(), err);
    assert!(err.to_string().contains(&format!("set rng_seed = {}", seed)));

    config.rng_seed = Some(seed);
    assert_eq!(inspect_output(&config), Ok(OutputStatus::Complete(out.clone())));
    assert_eq!(crate::run(config).unwrap().1, signal);

    let _ = std::fs::remove_dir_all(&out);
  }
  //----------------------------------------------------------------------------
  // Plan test 5: increasing number_runs from 1 to 2 gives the result of a
  // fresh run with number_runs = 2.
  #[test]
  fn test_more_runs_equal_a_fresh_run(){
    let out_a = temp_dir("runs_1_2");
    let out_b = temp_dir("runs_2");

    let config_1 = small_config(&out_a, "number_runs = 1\nrng_seed = 3", 
        &RANDOM);
    crate::run(config_1).unwrap();

    let config_2 = small_config(&out_a, "number_runs = 2\nrng_seed = 3", 
        &RANDOM);
    match inspect_output(&config_2){
      Ok(OutputStatus::Resume(_, summaries)) => {
        assert_eq!(summaries[0].instances_complete, 1);
        assert_eq!(summaries[0].instances_added, 1);
      },
      other => panic!("expected Resume, got {:?}", other),
    }
    let instance_1 = tree_snapshot(&instance_dir(&out_a, 1));
    let (_, resumed) = crate::run(config_2).unwrap();
    assert_eq!(instance_1, tree_snapshot(&instance_dir(&out_a, 1)));

    let fresh_config = small_config(&out_b, "number_runs = 2\nrng_seed = 3", 
        &RANDOM);
    let (_, fresh) = crate::run(fresh_config).unwrap();
    assert_eq!(resumed, fresh);

    let _ = std::fs::remove_dir_all(&out_a);
    let _ = std::fs::remove_dir_all(&out_b);
  }
  //----------------------------------------------------------------------------
  // Plan test 6: outputs of different simulations stop the run.
  #[test]
  fn test_different_simulations_stop_the_run(){
    let out = temp_dir("aborts");
    let config = small_config(&out, "number_runs = 2\nrng_seed = 5", &[]);
    crate::run(config.clone()).unwrap();

    // Fewer instances.
    let mut fewer = config.clone();
    fewer.number_runs = Some(1);
    assert_eq!(inspect_output(&fewer), Err(CluEError::OutputFewerRuns(
            out.clone(), config.input_structure_file.clone().unwrap(), 2, 1)));
    assert_eq!(crate::run(fewer).unwrap_err(), CluEError::OutputFewerRuns(
            out.clone(), config.input_structure_file.clone().unwrap(), 2, 1));

    // A different seed.
    let mut reseeded = config.clone();
    reseeded.rng_seed = Some(6);
    assert_eq!(inspect_output(&reseeded), 
        Err(CluEError::OutputSeedMismatch(out.clone(), 5, 6)));

    // A changed physics setting, with and without a rustc note.
    let changed = small_config(&out, "number_runs = 2\nrng_seed = 5", 
        &[("magnetic_field = 1.2","magnetic_field = 1.3")]);
    match inspect_output(&changed){
      Err(CluEError::OutputConfigMismatch(dir, differences, None)) => {
        assert_eq!(dir, out);
        assert_eq!(differences, vec![
            "magnetic_field > elements:".to_string(),
            "  - 1.2,".to_string(),
            "  + 1.3,".to_string()]);
      },
      other => panic!("expected OutputConfigMismatch, got {:?}", other),
    }
    edit_manifest(&out, |m| m.rustc_version = "rustc 0.0.0 (old)".to_string());
    let err = inspect_output(&changed).unwrap_err();
    assert!(matches!(&err, CluEError::OutputConfigMismatch(_, _, 
            Some((old, new))) if old == "rustc 0.0.0 (old)" 
          && new == RUSTC_VERSION));
    assert!(err.to_string().contains("compiled with rustc 0.0.0 (old); this \
build uses"));
    // Unchanged settings are not affected by the compiler version.
    assert_eq!(inspect_output(&config), Ok(OutputStatus::Complete(out.clone())));

    // Another CluE version.
    edit_manifest(&out, |m| m.clue_version = "0.0.1".to_string());
    assert_eq!(inspect_output(&config), Err(CluEError::OutputVersionMismatch(
            out.clone(), "0.0.1".to_string(), CLUE_VERSION.to_string())));
    edit_manifest(&out, |m| m.clue_version = CLUE_VERSION.to_string());

    // A structure file with other contents is a removed structure.
    let structure_copy = format!("{}_structure.pdb", out);
    let mut contents = std::fs::read_to_string(
        config.input_structure_file.as_ref().unwrap()).unwrap();
    contents.insert_str(0, "REMARK changed copy\n");
    std::fs::write(&structure_copy, contents).unwrap();
    let mut other_structure = config.clone();
    other_structure.input_structure_file = Some(structure_copy.clone());
    let record = &read_manifest(&out).unwrap().unwrap().structures[0];
    assert_eq!(inspect_output(&other_structure), 
        Err(CluEError::OutputStructureRemoved(out.clone(), record.file.clone(),
            record.id.clone())));

    // A recorded instance seed that this build does not derive.
    edit_manifest(&out, |m| m.structures[0].instances[1].seed += 1);
    assert!(matches!(inspect_output(&config), 
          Err(CluEError::OutputSeedSequenceChanged(..))));

    let _ = std::fs::remove_file(&structure_copy);
    let _ = std::fs::remove_dir_all(&out);
  }
  //----------------------------------------------------------------------------
  // Plan test 7: requesting methyl partitions after a run without them 
  // computes them, with the same signal.
  #[test]
  fn test_missing_outputs_are_computed_again(){
    let out = temp_dir("methyl");
    let without = small_config(&out, "rng_seed = 2", 
        &[("methyl_partitions = true","methyl_partitions = false")]);
    let (_, signal) = crate::run(without).unwrap();
    let ori_dir = orientation_dirs(&instance_dir(&out, 1))[0].clone();
    let methyl_dir = format!("{}/{}", ori_dir, 
        crate::config::SAVE_FILE_METHYL_PARTITIONS);
    assert!(!Path::new(&methyl_dir).exists());

    let with = small_config(&out, "rng_seed = 2", &[]);
    assert_eq!(with.write_methyl_partitions, Some(true));
    assert!(matches!(inspect_output(&with), 
          Ok(OutputStatus::Resume(..))));
    let (_, signal_with) = crate::run(with.clone()).unwrap();
    assert_eq!(signal, signal_with);
    assert!(std::fs::read_dir(&methyl_dir).unwrap().count() > 0);
    assert_eq!(inspect_output(&with), Ok(OutputStatus::Complete(out.clone())));

    let _ = std::fs::remove_dir_all(&out);
  }
  //----------------------------------------------------------------------------
  // Plan test 8: run_name is a removed key, also for clue_oxide --check.
  #[test]
  fn test_run_name_is_removed(){
    let out = temp_dir("run_name");
    let toml = small_toml(&out, "rng_seed = 1", 
        &[("[output]","[output]\n  run_name = \"system\"")]);
    let err = Config::from_toml_string(&toml).unwrap_err();
    let CluEError::InvalidTOMLKeys(problems) = &err else{
      panic!("expected InvalidTOMLKeys, got {:?}", err);
    };
    assert_eq!(problems.len(), 1);
    assert!(problems[0].starts_with("output.run_name: removed; "));

    let input_file = format!("{}.toml", out);
    std::fs::write(&input_file, &toml).unwrap();
    let report = crate::check(CommandLineInput::new(vec!["clue_oxide".into(),
          "-c".into(), input_file.clone()]).unwrap());
    assert_eq!(report.input, Err(err));
    assert_eq!(report.output, None);
    assert!(!report.is_ok());
    assert!(!Path::new(&out).exists());
    let _ = std::fs::remove_file(&input_file);
  }
  //----------------------------------------------------------------------------
  // Plan test 9: an orientation computed on its own gives the same signal
  // as in the full run.
  #[test]
  fn test_orientations_are_independent(){
    let out = temp_dir("orientations");
    let changes = [RANDOM[0], RANDOM[1], TWO_ORIENTATIONS[0], 
      TWO_ORIENTATIONS[1]];
    let config = small_config(&out, "rng_seed = 8", &changes);
    let (_, signal) = crate::run(config.clone()).unwrap();

    let dirs = orientation_dirs(&instance_dir(&out, 1));
    assert_eq!(dirs.len(), 2);
    let signals = dirs.iter()
      .map(|d| read_signals(&format!("{}/{}", d, SIGNAL_FILE)))
      .collect::<Vec::<_>>();
    assert_ne!(signals[0], signals[1]);

    // Recompute the second orientation only.
    let ori_signal = format!("{}/{}", dirs[1], SIGNAL_FILE);
    std::fs::remove_file(&ori_signal).unwrap();
    let first = tree_snapshot(&dirs[0]);
    let (_, signal_again) = crate::run(config).unwrap();
    assert_eq!(read_signals(&ori_signal), signals[1]);
    assert_eq!(first, tree_snapshot(&dirs[0]));
    assert_eq!(signal, signal_again);

    let _ = std::fs::remove_dir_all(&out);
  }
  //----------------------------------------------------------------------------
  // Plan test 10: the seeds of a structure and its instances.
  #[test]
  fn test_seeds(){
    let (h1, h2) = (0x0123456789abcdef_u64, 0xfedcba9876543210_u64);
    // The structure seed depends on the main seed and the content hash,
    // not on the position of the structure in a list.
    let seeds_12 = [h1, h2].iter().map(|&h| structure_seed(7, h))
      .collect::<Vec::<u64>>();
    let seeds_21 = [h2, h1].iter().map(|&h| structure_seed(7, h))
      .collect::<Vec::<u64>>();
    assert_eq!(seeds_12[0], seeds_21[1]);
    assert_eq!(seeds_12[1], seeds_21[0]);
    assert_ne!(seeds_12[0], seeds_12[1]);
    assert_ne!(structure_seed(7, h1), structure_seed(8, h1));

    // The first instance seeds do not depend on the number of instances.
    let five = instance_seeds(seeds_12[0], 5);
    assert_eq!(instance_seeds(seeds_12[0], 2), five[..2].to_vec());
    assert_eq!(five.len(), 5);

    // Main seeds can be written as TOML integers.
    for _ in 0..100{
      assert!(draw_main_seed() <= i64::MAX as u64);
    }
  }
  //----------------------------------------------------------------------------
  // Plan test 11: without orientation directories, an unfinished instance
  // is computed again from its seed.
  #[test]
  fn test_without_orientation_directories(){
    let out = temp_dir("no_orientations");
    let config = small_config(&out, "rng_seed = 4", &[RANDOM[0], RANDOM[1],
        ("orientation_signals = true","orientation_signals = false")]);
    let (_, signal) = crate::run(config.clone()).unwrap();
    let dir_1 = instance_dir(&out, 1);
    assert!(!Path::new(&format!("{}/{}", dir_1, 
            SAVE_DIR_ORIENTATION_SIGNALS)).exists());
    let manifest = read_manifest(&out).unwrap().unwrap();
    assert!(manifest.structures[0].instances[0].orientations.iter()
        .all(|o| !o.complete));

    // An interrupted instance: not complete, without its signal.
    edit_manifest(&out, |m| m.structures[0].instances[0].complete = false);
    std::fs::remove_file(format!("{}/{}", dir_1, SIGNAL_FILE)).unwrap();
    assert!(matches!(inspect_output(&config), Ok(OutputStatus::Resume(..))));
    let (_, signal_again) = crate::run(config.clone()).unwrap();
    assert_eq!(signal, signal_again);
    assert_eq!(inspect_output(&config), Ok(OutputStatus::Complete(out.clone())));

    let _ = std::fs::remove_dir_all(&out);
  }
  //----------------------------------------------------------------------------
  // Plan test 12: clue_oxide --check reports the input and the output, 
  // and writes nothing.
  #[test]
  fn test_check_reports_input_and_output(){
    let out = temp_dir("check");
    let input_file = format!("{}.toml", out);
    let check = |toml: &str| -> CheckReport_{
      std::fs::write(&input_file, toml).unwrap();
      let before = tree_snapshot(&out);
      let report = crate::check(CommandLineInput::new(vec![
            "clue_oxide".into(), "-c".into(), input_file.clone()]).unwrap());
      assert_eq!(before, tree_snapshot(&out), "--check wrote files");
      report
    };
    type CheckReport_ = crate::CheckReport;
    let output_line = |report: &CheckReport_| report.messages()[1].clone();

    // Not saved.
    let toml = small_toml(&out, "rng_seed = 1", &[]).replace(
        &format!("directory = '{}'", out), "directory = ''");
    let report = check(&toml);
    assert_eq!(report.messages(), vec![(false, "input OK".to_string()),
        (false, "output not saved (output.directory = \"\"); nothing to \
check".to_string())]);
    assert!(report.is_ok());

    // No directory and no seed.
    let toml = small_toml(&out, "number_runs = 1", &[]).replace(
        &format!("directory = '{}'", out), "");
    let report = check(&toml);
    assert_eq!(report.output, Some(Ok(OutputStatus::NewWithoutSeed)));
    assert!(report.is_ok());

    // Missing and empty directory.
    let toml = small_toml(&out, "number_runs = 2\nrng_seed = 1", &[]);
    assert_eq!(output_line(&check(&toml)), (false, format!("output {}: does \
not exist or is empty; a run starts a new simulation", out)));
    std::fs::create_dir_all(&out).unwrap();
    assert_eq!(check(&toml).output, Some(Ok(OutputStatus::New(out.clone()))));

    // Complete.
    let config = Config::from_toml_string(&toml).unwrap();
    crate::run(config.clone()).unwrap();
    assert_eq!(output_line(&check(&toml)), (false, format!("output {}: \
complete; a run will do no calculations", out)));

    // Resume: one complete instance and one added instance.
    let toml_3 = small_toml(&out, "number_runs = 3\nrng_seed = 1", &[]);
    std::fs::remove_file(format!("{}/{}", instance_dir(&out, 2), SIGNAL_FILE))
      .unwrap();
    let report = check(&toml_3);
    assert_eq!(output_line(&report), (false, format!("output {}: a run \
resumes: 1 of 3 instances complete (1 orientations complete in unfinished \
instances), 1 instances added; the orientation seeds are checked by the run",
        out)));
    assert!(report.is_ok());

    // An output that stops the run: the message of the run.
    let toml_1 = small_toml(&out, "number_runs = 1\nrng_seed = 1", &[]);
    let report = check(&toml_1);
    let run_err = crate::run(Config::from_toml_string(&toml_1).unwrap())
      .unwrap_err();
    assert_eq!(output_line(&report), (true, format!("output: CluE Error: {}.",
            run_err)));
    assert!(!report.is_ok());

    // Input with errors: the output is not checked.
    let report = check(&toml.replace("radius = 8", "radius = \"8\""));
    assert!(report.input.is_err());
    assert_eq!(output_line(&report), (true, "output: not checked, because the \
input has errors".to_string()));
    assert!(!report.is_ok());

    let _ = std::fs::remove_file(&input_file);
    let _ = std::fs::remove_dir_all(&out);
  }
  //----------------------------------------------------------------------------
  // Checkpoints: a damaged batch file is computed again; a different batch
  // size stops a resume; a directory without a record stops the run.
  #[test]
  fn test_checkpoint_files(){
    let out = temp_dir("batches");
    let config = small_config(&out, "rng_seed = 9\ncluster_batch_size = 4", 
        &[AUX]);
    let (_, signal) = crate::run(config.clone()).unwrap();
    let ori_dir = orientation_dirs(&instance_dir(&out, 1))[0].clone();
    let batch = format!("{}/auxiliary_signals/cluster_size_2_batch_1.csv", 
        ori_dir);
    let full_batch = std::fs::read_to_string(&batch).unwrap();

    // A truncated batch file and a missing orientation signal.
    let truncated = full_batch.lines().take(3).collect::<Vec::<_>>()
      .join("\n");
    std::fs::write(&batch, truncated).unwrap();
    std::fs::remove_file(format!("{}/{}", ori_dir, SIGNAL_FILE)).unwrap();

    // Another batch size cannot resume.
    let mut other_batch_size = config.clone();
    other_batch_size.cluster_batch_size = Some(5);
    assert_eq!(inspect_output(&other_batch_size), 
        Err(CluEError::OutputBatchSizeMismatch(out.clone(), 4, 5)));

    let (_, signal_again) = crate::run(config.clone()).unwrap();
    assert_eq!(signal, signal_again);
    assert_eq!(std::fs::read_to_string(&batch).unwrap(), full_batch);

    // A complete output does not depend on the batch size.
    assert_eq!(inspect_output(&other_batch_size), 
        Ok(OutputStatus::Complete(out.clone())));

    let _ = std::fs::remove_dir_all(&out);

    // A directory with files but without a record.
    std::fs::create_dir_all(&out).unwrap();
    std::fs::write(format!("{}/notes.txt", out), "x").unwrap();
    assert_eq!(inspect_output(&config), Err(CluEError::UnrecognizedOutput(
            out.clone(), vec!["notes.txt".to_string()])));
    let _ = std::fs::remove_dir_all(&out);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_fingerprint_differences(){
    let old = "Config {\n    radius: Some(\n        8.0,\n    ),\n    \
number_runs: None,\n    particles: [\n        ParticleConfig {\n            \
label: \"h\",\n        },\n    ],\n}";
    let new = old.replace("8.0", "9.0").replace("\"h\"", "\"d\"");
    assert_eq!(fingerprint_differences(old, &new), vec![
        "radius:", "  - 8.0,", "  + 9.0,",
        "particles:", "  - label: \"h\",", "  + label: \"d\","]
        .iter().map(|s| s.to_string()).collect::<Vec::<String>>());
    assert!(fingerprint_differences(old, old).is_empty());

    // An added line.
    let added = old.replace("    number_runs: None,\n", 
        "    number_runs: None,\n    rng_seed: None,\n");
    assert_eq!(fingerprint_differences(old, &added), 
        vec!["  + rng_seed: None,".to_string()]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_manifest_round_trip(){
    let manifest = Manifest{
      manifest_format: MANIFEST_FORMAT,
      clue_version: CLUE_VERSION.to_string(),
      rustc_version: RUSTC_VERSION.to_string(),
      rng_seed: u64::MAX,
      cluster_batch_size: 10,
      config_fingerprint: "Config {\n    a: \"'''\\\"\",\n}".to_string(),
      input_files: vec![InputFileRecord{path: "c.toml".to_string(), 
        content_hash: "0011223344556677".to_string()}],
      structures: vec![StructureRecord{
        id: "s_1".to_string(),
        file: "a.pdb".to_string(),
        content_hash: "1".to_string(),
        weight: 0.25,
        seed: u64::MAX - 1,
        number_runs: 1,
        instances: vec![InstanceRecord{
          index: 1,
          seed: 0,
          orientation_seeds: vec![u64::MAX, 3],
          complete: true,
          outputs: vec!["info.bath".to_string()],
          files: vec!["signal.csv".to_string()],
          orientations: vec![OrientationRecord{complete: true, 
            directory: "orientations/x".to_string(),
            outputs: vec![], files: vec!["signal.csv".to_string()]},
            OrientationRecord::default()],
        }],
      }],
    };
    let text = toml::to_string(&manifest).unwrap();
    let read: Manifest = toml::from_str(&text).unwrap();
    assert_eq!(read, manifest);

    // An unsupported format is reported.
    let out = temp_dir("manifest_format");
    std::fs::create_dir_all(&out).unwrap();
    let mut future = manifest.clone();
    future.manifest_format = MANIFEST_FORMAT + 1;
    write_manifest(&out, &future).unwrap();
    assert!(matches!(read_manifest(&out), 
          Err(CluEError::CannotParseManifest(..))));
    let _ = std::fs::remove_dir_all(&out);
  }
  //----------------------------------------------------------------------------
}
