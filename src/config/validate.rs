//! Checks of the input before a run: the key check of the parsed TOML
//! (`check_toml_keys`), and the completeness check of a `Config`
//! (`Config::check_completeness`).
use crate::clue_errors::CluEError;
use crate::config::{
  ClusterMethod, ClusterPopulations, ClusterSource, Config,
  DetectedPopulation, DetectionOp, PartitionTableConfig,
};
use crate::config::toml_keys::*;
use crate::isotopes::Isotope;

use std::path::Path;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// Keys of the tables of the input file.
// The key names are the ones read by `ConfigTOML` and the structs of its
// tables (src/config/config_toml.rs), `ParticleConfig`, 
// `ParticleProperties`, `IsotopeProperties`, `EigSpecifier`
// (src/config/particle_config.rs), and `ParticleFilter` and
// `VectorSpecifier` (src/structure/particle_filter.rs).

// The top level holds only tables.
const TOP_LEVEL_KEYS: [&str; 10] = [
  "units",
  "structure",
  "detected_spin",
  "groups",
  "conditions",
  "orientations",
  "pulse_sequence",
  "cluster_expansion",
  "run",
  "output",
];

const UNITS_KEYS: [&str; 4] = [
  "distance",
  "energy",
  "magnetic_field",
  "time",
];

const STRUCTURE_KEYS: [&str; 6] = [
  "file",
  "pdb_model_index",
  "radius",
  "replicate_unit_cell",
  "clash_distance",
  "clash_distance_pbc",
];

const DETECTED_SPIN_KEYS: [&str; 11] = [
  "density_matrix",
  "detection_operator",
  "frame",
  "g_matrix",
  "identity",
  "max_spherical_tensor_rank",
  "multiplicity",
  "position",
  "transition",
  "weight_tolerance",
  "zerofield",
];

const CONDITIONS_KEYS: [&str; 2] = [
  "magnetic_field",
  "temperature",
];

const ORIENTATIONS_KEYS: [&str; 6] = [
  "file",
  "grid",
  "number",
  "vector",
  "vector_grid",
  "weight_tolerance",
];

const PULSE_SEQUENCE_KEYS: [&str; 5] = [
  "sequence",
  "number_timepoints",
  "tau_increments",
  "number_timepoints2",
  "tau2_increments",
];

const CLUSTER_EXPANSION_KEYS: [&str; 7] = [
  "method",
  "max_cluster_size",
  "ensemble",
  "populations",
  "mean_fields",
  "thermalize_cluster_populations",
  "clustering",
];

const CLUSTERING_KEYS: [&str; 8] = [
  "file",
  "max_spins",
  "partitioning",
  "kmeans_size",
  "kmeans_iter_cap",
  "partition_table",
  "connect_exchange_groups",
  "pair_cutoffs",
];

const PAIR_CUTOFFS_KEYS: [&str; 7] = [
  KEY_CUTOFF_COUPLING,
  KEY_CUTOFF_DELTA_HF,
  KEY_CUTOFF_DIPOLE_PERP,
  KEY_CUTOFF_DISTANCE,
  KEY_CUTOFF_HAHN_MOD_DEPTH,
  KEY_CUTOFF_HAHN_TAYLOR_4,
  KEY_CUTOFF_DELTA_ZEEMAN,
];

const RUN_KEYS: [&str; 4] = [
  "rng_seed",
  "number_runs",
  "parallel",
  "cluster_batch_size",
];

const OUTPUT_KEYS: [&str; 15] = [
  "directory",
  "run_name",
  KEY_OUT_AUX_SIGS,
  KEY_OUT_BATH,
  KEY_OUT_CLUSTERS,
  KEY_OUT_DET_SPIN,
  KEY_OUT_EXCHANGE_GROUPS,
  KEY_OUT_INFO,
  KEY_OUT_METHYL_PARTITIONS,
  KEY_OUT_ORI_GRID,
  KEY_OUT_ORI_SIGS,
  KEY_OUT_PART_TAB,
  KEY_OUT_SANS_SPIN_SIGS,
  KEY_OUT_STRUC_PDB,
  KEY_OUT_TENSORS,
];

// The isotope tables of a group are checked by `is_isotope_key()`.
const GROUP_KEYS: [&str; 7] = [
  KEY_NAME,
  KEY_ABUNDANCE_TOL,
  KEY_DROP_PROB,
  KEY_ISO_COSUBSTITUTE,
  KEY_CELLS,
  KEY_INCLUDE,
  KEY_EXCLUDE,
];

const INCLUDE_KEYS: [&str; 9] = [
  KEY_SELE_CHAIN_IDS,
  KEY_SELE_ELEMENTS,
  KEY_SELE_INDICES,
  KEY_SELE_NAMES,
  KEY_SELE_RESIDUES,
  KEY_SELE_RES_NUMS,
  KEY_SELE_SERIALS,
  KEY_SELE_DISTANCE,
  KEY_BONDED_TO,
];

const EXCLUDE_KEYS: [&str; 8] = [
  KEY_SELE_CHAIN_IDS,
  KEY_SELE_ELEMENTS,
  KEY_SELE_INDICES,
  KEY_SELE_NAMES,
  KEY_SELE_RESIDUES,
  KEY_SELE_RES_NUMS,
  KEY_SELE_SERIALS,
  KEY_BONDED_TO,
];

// The seven atom attributes.
const BONDED_TO_KEYS: [&str; 7] = [
  KEY_SELE_CHAIN_IDS,
  KEY_SELE_ELEMENTS,
  KEY_SELE_INDICES,
  KEY_SELE_NAMES,
  KEY_SELE_RESIDUES,
  KEY_SELE_RES_NUMS,
  KEY_SELE_SERIALS,
];

const ISOTOPE_KEYS: [&str; 8] = [
  KEY_ISO_ABUNDACE,
  KEY_ISO_ACTIVE,
  KEY_ISO_C3_TUNNEL_SPLITTING,
  KEY_ISO_ELEC_QUADRUPOLE,
  KEY_ISO_EXCHANGE_COUPLING,
  KEY_ISO_G_MATRIX,
  KEY_ISO_HYPERFINE,
  KEY_ISO_ZEROFIELD,
];

const TENSOR_KEYS: [&str; 4] = [
  KEY_EIG_VALUES,
  KEY_EIG_AXES,
  KEY_IST_COEF,
  KEY_STEVENS_COEF,
];

const AXES_KEYS: [&str; 3] = [KEY_EIG_X_AXIS, KEY_EIG_Y_AXIS, KEY_EIG_Z_AXIS];

const AXIS_KEYS: [&str; 6] = [
  KEY_VEC_SPECIFIER_FROM,
  KEY_VEC_SPECIFIER_FROM_BONDED_TO,
  KEY_VEC_SPECIFIER_FROM_SAME_MOLECULE_AS,
  KEY_VEC_SPECIFIER_TO,
  KEY_VEC_SPECIFIER_TO_BONDED_TO,
  KEY_VEC_SPECIFIER_TO_SAME_MOLECULE_AS,
];
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// Keys of the input format before 1.0.0, which are rejected with a hint to
// what replaced them.

// The top-level keys and tables of the old format, and where they moved.
const MOVED_TOP_LEVEL_KEYS: [(&str,&str); 37] = [
  ("clash_distance", "structure.clash_distance"),
  ("clash_distance_pbc", "structure.clash_distance_pbc"),
  ("cluster_batch_size", "run.cluster_batch_size"),
  ("cluster_method", "cluster_expansion.method"),
  ("cluster_source", "cluster_expansion.clustering.file"),
  ("connect_exchange_groups", 
   "cluster_expansion.clustering.connect_exchange_groups"),
  ("ensemble_cce", "cluster_expansion.ensemble"),
  ("input_structure_file", "structure.file"),
  ("kmeans_iter_cap", "cluster_expansion.clustering.kmeans_iter_cap"),
  ("kmeans_size", "cluster_expansion.clustering.kmeans_size"),
  ("magnetic_field", "conditions.magnetic_field"),
  ("max_cluster_size", "cluster_expansion.max_cluster_size"),
  ("max_spherical_tensor_rank", "detected_spin.max_spherical_tensor_rank"),
  ("max_spins", "cluster_expansion.clustering.max_spins"),
  ("mean_fields", "cluster_expansion.mean_fields"),
  ("number_runs", "run.number_runs"),
  ("number_timepoints", "pulse_sequence.number_timepoints"),
  ("number_timepoints2", "pulse_sequence.number_timepoints2"),
  ("output_directory", "output.directory"),
  ("partition_table", "cluster_expansion.clustering.partition_table"),
  ("partitioning", "cluster_expansion.clustering.partitioning"),
  ("pdb_model_index", "structure.pdb_model_index"),
  ("populations", "cluster_expansion.populations"),
  ("radius", "structure.radius"),
  ("replicate_unit_cell", "structure.replicate_unit_cell"),
  ("rng_seed", "run.rng_seed"),
  ("run_in_parallel", "run.parallel"),
  ("run_name", "output.run_name"),
  ("tau_increments", "pulse_sequence.tau_increments"),
  ("tau2_increments", "pulse_sequence.tau2_increments"),
  ("temperature", "conditions.temperature"),
  ("thermalize_cluster_populations", 
   "cluster_expansion.thermalize_cluster_populations"),
  ("unit_of_distance", "units.distance"),
  ("unit_of_energy", "units.energy"),
  ("unit_of_magnetic_field", "units.magnetic_field"),
  ("unit_of_time", "units.time"),
  ("pair_cutoffs", "cluster_expansion.clustering.pair_cutoffs"),
];

// The old key pulse_sequence is now the table [pulse_sequence]; 
// an old value, which is not a table, moved to pulse_sequence.sequence.
const OLD_KEY_PULSE_SEQUENCE: &str = "pulse_sequence";

// The keys that were removed because the code did not use them, 
// by the table they were in, with the reason.
const REMOVED_KEYS: [(Schema,&str,&str); 7] = [
  (Schema::TopLevel, "min_cell_size", "it was read but never used"),
  (Schema::TopLevel, "max_cell_size", "it was read but never used"),
  (Schema::TopLevel, "pulses", "[[pulses]] was read but never used; the \
pulses are built from detected_spin.transition"),
  (Schema::DetectedSpin, "electric_quadrupole", "it was read but never \
used; the isotope key electric_quadrupole of [[groups]] is unchanged"),
  (Schema::PairCutoffs, "coupling", "it was never read; the coupling cutoff \
is coupling_xx_yy"),
  (Schema::Output, "cluster_signals", "it was read but never used; cluster \
signals are written with auxiliary_signals"),
  (Schema::Output, "config", "it was read but never used"),
];

// The old group table that the include and exclude tables and the group
// key cells replace.
const OLD_KEY_SELECTION: &str = "selection";

// The keys of the old selection table that the group key cells replaces.
const OLD_CELL_KEYS: [&str; 2] = ["cell_type", "primary_cell"];

// The atom attributes of the old selection table, and their new names.
const OLD_SELECTION_ATTRIBUTES: [(&str,&str); 7] = [
  ("chain_ids", KEY_SELE_CHAIN_IDS),
  ("elements", KEY_SELE_ELEMENTS),
  ("indices", KEY_SELE_INDICES),
  ("names", KEY_SELE_NAMES),
  ("residues", KEY_SELE_RESIDUES),
  ("residue_sequence_numbers", KEY_SELE_RES_NUMS),
  ("serials", KEY_SELE_SERIALS),
];
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `Schema` lists the kinds of tables in the input file.
#[derive(Debug,Clone,Copy,PartialEq)]
enum Schema{
  TopLevel,
  Units,
  Structure,
  DetectedSpin,
  Group,
  Include,
  Exclude,
  BondedTo,
  Isotope,
  Tensor,
  Axes,
  Axis,
  Conditions,
  Orientations,
  PulseSequence,
  ClusterExpansion,
  Clustering,
  PairCutoffs,
  Run,
  Output,
}

const ALL_SCHEMAS: [Schema; 20] = [
  Schema::TopLevel,
  Schema::Units,
  Schema::Structure,
  Schema::DetectedSpin,
  Schema::Group,
  Schema::Include,
  Schema::Exclude,
  Schema::BondedTo,
  Schema::Isotope,
  Schema::Tensor,
  Schema::Axes,
  Schema::Axis,
  Schema::Conditions,
  Schema::Orientations,
  Schema::PulseSequence,
  Schema::ClusterExpansion,
  Schema::Clustering,
  Schema::PairCutoffs,
  Schema::Run,
  Schema::Output,
];

impl Schema{
  //----------------------------------------------------------------------------
  // This function returns the fixed keys of the table.
  fn keys(&self) -> &'static [&'static str]{
    match self{
      Schema::TopLevel => &TOP_LEVEL_KEYS,
      Schema::Units => &UNITS_KEYS,
      Schema::Structure => &STRUCTURE_KEYS,
      Schema::DetectedSpin => &DETECTED_SPIN_KEYS,
      Schema::Group => &GROUP_KEYS,
      Schema::Include => &INCLUDE_KEYS,
      Schema::Exclude => &EXCLUDE_KEYS,
      Schema::BondedTo => &BONDED_TO_KEYS,
      Schema::Isotope => &ISOTOPE_KEYS,
      Schema::Tensor => &TENSOR_KEYS,
      Schema::Axes => &AXES_KEYS,
      Schema::Axis => &AXIS_KEYS,
      Schema::Conditions => &CONDITIONS_KEYS,
      Schema::Orientations => &ORIENTATIONS_KEYS,
      Schema::PulseSequence => &PULSE_SEQUENCE_KEYS,
      Schema::ClusterExpansion => &CLUSTER_EXPANSION_KEYS,
      Schema::Clustering => &CLUSTERING_KEYS,
      Schema::PairCutoffs => &PAIR_CUTOFFS_KEYS,
      Schema::Run => &RUN_KEYS,
      Schema::Output => &OUTPUT_KEYS,
    }
  }
  //----------------------------------------------------------------------------
  // This function returns `true` if `key` is allowed in the table.
  fn allows(&self, key: &str) -> bool{
    self.keys().contains(&key)
      || (*self == Schema::Group && is_isotope_key(key))
  }
  //----------------------------------------------------------------------------
  // This function returns the kind of table that the value of `key` is,
  // if the value is a table (`false`) or an array of tables (`true`).
  // Keys whose values are not tables give `None`.
  // Keys whose values may be a table or an array, such as tensors and
  // axes, are tables only when given as tables.
  fn child(&self, key: &str) -> Option<(Schema,bool)>{
    match self{
      Schema::TopLevel => match key{
        "units" => Some((Schema::Units,false)),
        "structure" => Some((Schema::Structure,false)),
        "detected_spin" => Some((Schema::DetectedSpin,false)),
        "groups" => Some((Schema::Group,true)),
        "conditions" => Some((Schema::Conditions,false)),
        "orientations" => Some((Schema::Orientations,false)),
        "pulse_sequence" => Some((Schema::PulseSequence,false)),
        "cluster_expansion" => Some((Schema::ClusterExpansion,false)),
        "run" => Some((Schema::Run,false)),
        "output" => Some((Schema::Output,false)),
        _ => None,
      },
      Schema::DetectedSpin => match key{
        "g_matrix" | "zerofield" => Some((Schema::Tensor,false)),
        _ => None,
      },
      Schema::Group => {
        if key == KEY_INCLUDE{
          Some((Schema::Include,false))
        }else if key == KEY_EXCLUDE{
          Some((Schema::Exclude,false))
        }else if is_isotope_key(key){
          Some((Schema::Isotope,false))
        }else{
          None
        }
      },
      Schema::Include | Schema::Exclude => match key{
        KEY_BONDED_TO => Some((Schema::BondedTo,false)),
        _ => None,
      },
      Schema::Isotope => match key{
        KEY_ISO_G_MATRIX | KEY_ISO_HYPERFINE | KEY_ISO_ELEC_QUADRUPOLE
          | KEY_ISO_ZEROFIELD => Some((Schema::Tensor,false)),
        _ => None,
      },
      Schema::Tensor => match key{
        KEY_EIG_AXES => Some((Schema::Axes,false)),
        _ => None,
      },
      Schema::Axes => Some((Schema::Axis,false)),
      Schema::ClusterExpansion => match key{
        "clustering" => Some((Schema::Clustering,false)),
        _ => None,
      },
      Schema::Clustering => match key{
        "pair_cutoffs" => Some((Schema::PairCutoffs,false)),
        _ => None,
      },
      Schema::Units | Schema::Structure | Schema::BondedTo | Schema::Axis
        | Schema::Conditions | Schema::Orientations | Schema::PulseSequence
        | Schema::PairCutoffs | Schema::Run | Schema::Output => None,
    }
  }
  //----------------------------------------------------------------------------
  // This function names the table in messages.
  fn description(&self) -> &'static str{
    match self{
      Schema::TopLevel => "the top level",
      Schema::Units => "[units]",
      Schema::Structure => "[structure]",
      Schema::DetectedSpin => "[detected_spin]",
      Schema::Group => "[[groups]]",
      Schema::Include => "the include table of [[groups]]",
      Schema::Exclude => "the exclude table of [[groups]]",
      Schema::BondedTo => "a bonded_to table of [[groups]] \
(include.bonded_to, exclude.bonded_to)",
      Schema::Isotope => "an isotope table of [[groups]]",
      Schema::Tensor => "a tensor table (values, axes, \
spherical_coefficients, stevens_coefficients)",
      Schema::Axes => "the axes table of a tensor (x, y, z)",
      Schema::Axis => "an axis table of a tensor (from, from_bonded_to, \
from_same_molecule_as, to, to_bonded_to, to_same_molecule_as)",
      Schema::Conditions => "[conditions]",
      Schema::Orientations => "[orientations]",
      Schema::PulseSequence => "[pulse_sequence]",
      Schema::ClusterExpansion => "[cluster_expansion]",
      Schema::Clustering => "[cluster_expansion.clustering]",
      Schema::PairCutoffs => "[cluster_expansion.clustering.pair_cutoffs]",
      Schema::Run => "[run]",
      Schema::Output => "[output]",
    }
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//------------------------------------------------------------------------------
// This function returns `true` if `key` names an isotope table of a group.
// The groups read their isotope tables by `Isotope::to_string()`
// (`ParticleProperties::set_from_toml_table()`), so a key is accepted only
// if it is the string of the isotope it parses to.
fn is_isotope_key(key: &str) -> bool{
  match Isotope::from(key){
    Ok(isotope) => isotope.to_string() == key,
    Err(_) => false,
  }
}
//------------------------------------------------------------------------------
/// This function checks that every key of the parsed input `table` is
/// allowed in the table that contains it.
/// All keys that are not allowed are reported together, each with its full
/// path, for example `groups[2].include.elemnts`;
/// `groups[i]` counts the [[groups]] entries from 0.
/// Keys of the input format before 1.0.0 are reported with what replaced 
/// them, for example `radius: moved to structure.radius`.
/// Types and values are not checked here.
pub fn check_toml_keys(table: &toml::Table) -> Result<(),CluEError>{
  let mut errors = Vec::<String>::new();
  check_table(table, Schema::TopLevel, "", None, &mut errors);
  if errors.is_empty(){
    Ok(())
  }else{
    Err(CluEError::InvalidTOMLKeys(errors))
  }
}
//------------------------------------------------------------------------------
// This function checks the keys of `table`, which is of kind `schema` and
// at `path`, and of its subtables, and adds a message to `errors` for every
// key that is not allowed.
// `group` is the name of the group that contains the table, if any.
fn check_table(table: &toml::Table, schema: Schema, path: &str,
    group: Option<&str>, errors: &mut Vec::<String>)
{
  for (key, value) in table.iter(){
    let key_path = join_path(path, key);

    if !schema.allows(key){
      // Each key of an old selection table gets its own hint.
      if let (Schema::Group, OLD_KEY_SELECTION, Some(selection)) 
        = (schema, key.as_str(), value.as_table())
      {
        if !selection.is_empty(){
          // Without other keys, the old selection selected every atom of
          // its cells.
          let only_cells = selection.keys()
            .all(|k| OLD_CELL_KEYS.contains(&k.as_str()));
          for old_key in selection.keys(){
            errors.push(with_group(
                  &format!("{}.{}", key_path, old_key), group,
                  &old_selection_key_hint(path, old_key, only_cells)));
          }
          continue;
        }
      }
      errors.push(not_allowed_message(path, key, schema, group));
      continue;
    }

    if schema == Schema::TopLevel && key == OLD_KEY_PULSE_SEQUENCE
      && !value.is_table()
    {
      errors.push(with_group(&key_path, group, &format!("moved to {}.{}",
              OLD_KEY_PULSE_SEQUENCE, "sequence")));
      continue;
    }

    let Some((child, is_array)) = schema.child(key) else{
      continue;
    };

    // Values of the wrong type are reported by the parser.
    match (value, is_array){
      (toml::Value::Table(subtable), false)
        => check_table(subtable, child, &key_path, group, errors),
      (toml::Value::Array(array), true) => {
        for (ii, element) in array.iter().enumerate(){
          let toml::Value::Table(subtable) = element else{
            continue;
          };
          let element_group = if child == Schema::Group{
            subtable.get(KEY_NAME).and_then(|v| v.as_str())
          }else{
            group
          };
          check_table(subtable, child, &format!("{}[{}]", key_path, ii),
              element_group, errors);
        }
      },
      _ => (),
    }
  }
}
//------------------------------------------------------------------------------
// This function returns `key` at `path`.
fn join_path(path: &str, key: &str) -> String{
  if path.is_empty(){
    key.to_string()
  }else{
    format!("{}.{}", path, key)
  }
}
//------------------------------------------------------------------------------
// This function writes the message `problem` for the key at `key_path`,
// in the group named `group`, if any.
fn with_group(key_path: &str, group: Option<&str>, problem: &str) -> String{
  match group{
    Some(name) => format!("{} (group \"{}\"): {}", key_path, name, problem),
    None => format!("{}: {}", key_path, problem),
  }
}
//------------------------------------------------------------------------------
// This function writes the message for `key` in the table at `path`, 
// which is of kind `schema`, where `key` is not allowed.
// A key of the old input format names what replaced it; a key that is 
// allowed in other tables names them.
fn not_allowed_message(path: &str, key: &str, schema: Schema,
    group: Option<&str>) -> String
{
  let key_path = join_path(path, key);

  if schema == Schema::TopLevel{
    if let Some((_, new_path)) = MOVED_TOP_LEVEL_KEYS.iter()
      .find(|(old, _)| *old == key)
    {
      return with_group(&key_path, group, &format!("moved to {}", new_path));
    }
  }

  if let Some((_, _, reason)) = REMOVED_KEYS.iter()
    .find(|(s, old, _)| *s == schema && *old == key)
  {
    return with_group(&key_path, group, &format!("removed; {}", reason));
  }

  if schema == Schema::Group && key == OLD_KEY_SELECTION{
    return with_group(&key_path, group, &format!("replaced by {}, {}, and \
{}", join_path(path, KEY_INCLUDE), join_path(path, KEY_EXCLUDE), 
          join_path(path, KEY_CELLS)));
  }

  let places = ALL_SCHEMAS.iter()
    .filter(|s| s.allows(key))
    .map(|s| s.description())
    .collect::<Vec::<&str>>();

  if places.is_empty(){
    with_group(&key_path, group, "unknown key")
  }else{
    with_group(&key_path, group, &format!("not allowed in {}; this key \
belongs in {}", schema.description(), places.join(" or ")))
  }
}
//------------------------------------------------------------------------------
// This function returns what replaced the key `old_key` of the old
// selection table of the group at `group_path`.
// `only_cells` is `true` if the selection table had no keys other than 
// `OLD_CELL_KEYS`, so that it selected every atom of its cells; 
// the group then also needs `include = {}`.
fn old_selection_key_hint(group_path: &str, old_key: &str, only_cells: bool)
  -> String
{
  let at = |key: &str| join_path(group_path, key);
  let include_all = if only_cells{
    format!(", together with {} = {{}} to select all atoms as before", 
        at(KEY_INCLUDE))
  }else{
    String::new()
  };

  match old_key{
    "cell_type" => return format!("replaced by {} = \"all\", \"primary\", \
or \"extra\" (was \"all_cells\", \"primary_cell\", \"extracells\"){}", 
        at(KEY_CELLS), include_all),
    "primary_cell" => return format!("replaced by {} = \"primary\" (true) \
or \"extra\" (false){}", at(KEY_CELLS), include_all),
    "cell_ids" | "not_cell_ids" => return "removed; selection by unit-cell \
index is no longer supported".to_string(),
    "within_distance" => return format!("replaced by {} = [0.0, max]",
        at("include.distance")),
    "not_within_distance" => return format!("replaced by {} = [min, inf]",
        at("include.distance")),
    _ => (),
  }

  // X, not_X, bonded_X, and not_bonded_X.
  let (table, attribute) = if let Some(a) 
    = old_key.strip_prefix("not_bonded_")
  {
    ("exclude.bonded_to", a)
  }else if let Some(a) = old_key.strip_prefix("bonded_"){
    ("include.bonded_to", a)
  }else if let Some(a) = old_key.strip_prefix("not_"){
    ("exclude", a)
  }else{
    ("include", old_key)
  };

  match OLD_SELECTION_ATTRIBUTES.iter().find(|(old, _)| *old == attribute){
    Some((_, new)) => format!("moved to {}.{}", at(table), new),
    None => format!("unknown key; the selection table is replaced by {}, \
{}, and {}", at(KEY_INCLUDE), at(KEY_EXCLUDE), at(KEY_CELLS)),
  }
}
//------------------------------------------------------------------------------


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl Config{
  //----------------------------------------------------------------------------
  /// This function checks that the settings a run needs are set,
  /// that settings required by other settings are set, and that
  /// referenced files and directories exist.
  /// All problems are reported together, named by their TOML key.
  /// Relative paths are resolved against the current directory,
  /// as in a run.
  /// Settings that are read when the input is parsed
  /// (`cluster_expansion.clustering.partitioning` with its `kmeans_size`, 
  /// `orientations.file`, and a csv `detected_spin.position`) are checked
  /// by the parser.
  /// Errors that depend on the structure (group names, overlapping groups,
  /// empty selections) are found only by a run.
  pub fn check_completeness(&self) -> Result<(),CluEError>{
    let mut problems = Vec::<String>::new();

    // Settings without a default.
    let required = [
      (self.input_structure_file.is_some(), "structure.file"),
      (self.radius.is_some(), "structure.radius"),
      (self.magnetic_field.is_some(), "conditions.magnetic_field"),
      (self.detected_spin_position.is_some(), "detected_spin.position"),
      (self.pulse_sequence.is_some(), "pulse_sequence.sequence"),
      (!self.number_timepoints.is_empty(), 
       "pulse_sequence.number_timepoints"),
      (!self.tau_increments.is_empty(), "pulse_sequence.tau_increments"),
      (self.cluster_method.is_some(), "cluster_expansion.method"),
      (self.max_cluster_size.is_some(), "cluster_expansion.max_cluster_size"),
      // Set from units.time or its default when the input is parsed;
      // it can be unset only in a `Config` built in code.
      (self.unit_of_time_to_seconds.is_some(), "units.time"),
    ];
    for (is_set, key) in required.iter(){
      if !is_set{
        problems.push(format!("{} is not set", key));
      }
    }

    // Settings required by other settings.
    let has_temperature = self.temperature.is_some();

    // The bath populations are generated for every cluster method
    // (`SpinStates::generate()`).
    if self.cluster_populations == Some(ClusterPopulations::Thermal)
      && !has_temperature
    {
      problems.push(format!("conditions.temperature is required for \
cluster_expansion.populations = \"{}\"", KEY_DENSITY_MATRIX_THERMAL));
    }

    // Mean-field couplings are not implemented for ensemble CCE, for any
    // cluster method (`calculate_signal_at_orientation()`).
    if self.mean_fields == Some(true) && self.ensemble_cce == Some(true){
      problems.push(CluEError::MeanFieldECCENotImplemented.to_string());
    }

    // The detected spin needs at least two levels to have a coherence
    // that a pulse sequence can create and detect.
    // A multiplicity of 1 also follows from a spin-0 detected_spin.identity.
    let multiplicity = self.detected_spin_multiplicity;
    if let Some(m) = multiplicity{
      if m < 2{
        problems.push(format!("detected_spin.multiplicity must be at least \
2, but is {}", m));
      }
    }
    let valid_multiplicity = multiplicity.filter(|&m| m >= 2);

    // CCE (`build_block_diag_hamiltonian()`) and gCCE (the pulses of
    // `build_detection_operators()`) act on the transition between the two
    // levels of detected_spin.transition. It has a default only for
    // multiplicity 2 (`Config::set_detected_spin()`); for S > 1/2 the user
    // must choose the transition. APPA and LCE do not use it.
    let is_cce = self.cluster_method == Some(ClusterMethod::CCE);
    let is_gcce = self.cluster_method == Some(ClusterMethod::GCCE);
    if let (true, Some(m)) = (is_cce || is_gcce, valid_multiplicity){
      match self.detected_spin_transition{
        None => problems.push(format!("detected_spin.transition is required \
for detected_spin.multiplicity = {}", m)),
        Some([level_0, level_1]) => {
          if level_0 == level_1 || level_0 >= m || level_1 >= m{
            problems.push(format!("detected_spin.transition = [{}, {}] must \
name two different levels from 0 to {}", level_0, level_1, m - 1));
          }
        },
      }
    }

    // The detected spin operators are built only for gCCE, and only if
    // the detected spin has more than one state
    // (`build_detection_operators()`).
    if let (true, Some(m)) = (is_gcce, valid_multiplicity){
      match self.detection_operator{
        None => problems.push("detected_spin.detection_operator is required \
for cluster_expansion.method = \"gcce\"".to_string()),
        // [i, i] is the projector |i><i|, a valid observable, so only the
        // range is checked.
        Some(DetectionOp::Transition(level_0, level_1)) => {
          if level_0 >= m || level_1 >= m{
            problems.push(format!("detected_spin.detection_operator = \
[{}, {}] must name levels from 0 to {}", level_0, level_1, m - 1));
          }
        },
        Some(_) => (),
      }
      match self.detected_population{
        None => problems.push("detected_spin.density_matrix is required \
for cluster_expansion.method = \"gcce\"".to_string()),
        Some(DetectedPopulation::Thermal) if !has_temperature =>
          problems.push(format!("conditions.temperature is required for \
detected_spin.density_matrix = \"{}\"", KEY_DENSITY_MATRIX_THERMAL)),
        Some(_) => (),
      }
    }

    // thermalize_cluster_populations is used only by gCCE (`gcce()`).
    if is_gcce && self.thermalize_cluster_populations == Some(true){
      if self.ensemble_cce != Some(true){
        problems.push(CluEError::ThermalizationRequiresEnsembleCCE
            .to_string());
      }
      if !has_temperature{
        problems.push("conditions.temperature is required for \
cluster_expansion.thermalize_cluster_populations = true".to_string());
      }
    }

    // Files and directories.
    if let Some(filename) = &self.input_structure_file{
      if let Some(problem) = check_path(filename, false){
        problems.push(format!("structure.file: {}", problem));
      }
    }

    match &self.cluster_source{
      Some(ClusterSource::File(filename)) => {
        if let Some(problem) = check_path(filename, false){
          problems.push(format!("cluster_expansion.clustering.file: {}", 
                problem));
        }
      },
      Some(ClusterSource::OrientationsDir(dirname)) => {
        if let Some(problem) = check_path(dirname, true){
          problems.push(format!("cluster_expansion.clustering.file: {}", 
                problem));
        }
      },
      // The partition table is read only when clusters are built from
      // the structure.
      Some(ClusterSource::Structure) | None => {
        if let Some(PartitionTableConfig::CSV(filename))
          = &self.partition_table
        {
          if let Some(problem) = check_path(filename, false){
            problems.push(format!(
                  "cluster_expansion.clustering.partition_table: {}", 
                  problem));
          }
        }
      },
    }

    if problems.is_empty(){
      Ok(())
    }else{
      Err(CluEError::IncompleteConfig(problems))
    }
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//------------------------------------------------------------------------------
// This function returns a description of the problem if `path` is not an
// existing file (`is_dir == false`) or directory (`is_dir == true`).
fn check_path(path: &str, is_dir: bool) -> Option<String>{
  let p = Path::new(path);
  if !p.exists(){
    return Some(format!("\"{}\" does not exist", path));
  }
  if is_dir && !p.is_dir(){
    return Some(format!("\"{}\" is not a directory", path));
  }
  if !is_dir && !p.is_file(){
    return Some(format!("\"{}\" is not a file", path));
  }
  None
}
//------------------------------------------------------------------------------


//==============================================================================
#[cfg(test)]
mod tests{
  use super::*;
  use crate::io::FromTOMLString;
  use crate::config::config_toml::ConfigTOML;
  use strum::IntoEnumIterator;

  //----------------------------------------------------------------------------
  // This function returns the messages of the key check of `toml_str`.
  fn key_errors(toml_str: &str) -> Vec::<String>{
    let table: toml::Table = toml::from_str(toml_str).unwrap();
    match check_toml_keys(&table){
      Ok(()) => Vec::<String>::new(),
      Err(CluEError::InvalidTOMLKeys(errors)) => errors,
      Err(err) => panic!("expected InvalidTOMLKeys, got {:?}", err),
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_assets_pass_key_check(){
    for filename in ["assets/1omp_K26R1_0.003228966616048703.toml",
        "assets/1omp_K26R1_0.003228966616048703_gcce.toml",
        "assets/solvent_box.toml"]
    {
      let toml_str = std::fs::read_to_string(filename).unwrap();
      assert_eq!(key_errors(&toml_str), Vec::<String>::new(), "{}", filename);
      assert!(ConfigTOML::from_toml_string(&toml_str).is_ok(), "{}", filename);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_all_table_forms_pass_key_check(){
    // Every table kind, with tensors and axes given as tables.
    let errors = key_errors(r##"
        [units]
        distance = "Å"
        energy = "MHz"
        magnetic_field = "T"
        time = "μs"
        [structure]
        file = "a.pdb"
        pdb_model_index = 0
        radius = 10
        replicate_unit_cell = false
        clash_distance = 0.01
        clash_distance_pbc = 0.001
        [detected_spin]
        position = [0.0, 0.0, 0.0]
        max_spherical_tensor_rank = 2
        g_matrix.values = [2.0, 2.0, 2.0]
        g_matrix.axes.x = {from = "a", to_bonded_to = "b"}
        g_matrix.axes.y = "random"
        zerofield.values = [1.0, 1.0, -2.0]
        zerofield.axes.x = [1.0, 0.0, 0.0]
        zerofield.axes.z.from_same_molecule_as = "a"
        zerofield.axes.z.to_same_molecule_as = "b"
        zerofield.spherical_coefficients = [[0.0, 0.0, 0.0, 0.0, 0.0]]
        [conditions]
        magnetic_field = 1.2
        temperature = 50
        [orientations]
        grid = "lebedev"
        number = 170
        [pulse_sequence]
        sequence = "hahn"
        number_timepoints = [10]
        tau_increments = [0.1]
        number_timepoints2 = [1]
        tau2_increments = [0.0]
        [cluster_expansion]
        method = "cce"
        max_cluster_size = 2
        ensemble = false
        populations = "uniform"
        mean_fields = false
        thermalize_cluster_populations = false
        [cluster_expansion.clustering]
        file = "clusters.toml"
        max_spins = 4
        partitioning = "kmeans"
        kmeans_size = 3
        kmeans_iter_cap = 10
        partition_table = [1, 1]
        connect_exchange_groups = false
        [cluster_expansion.clustering.pair_cutoffs]
        distance = 10
        coupling_xx_yy = 1.0
        [run]
        rng_seed = 0
        number_runs = 1
        parallel = true
        cluster_batch_size = 10
        [output]
        directory = "out"
        run_name = "system"
        bath = true
        partition_table = true
        detected_spin = true
        [[groups]]
        name = "h"
        drop_probability = 0.1
        cells = "all"
        include = {distance = [0, 5], bonded_to = {residues = "SOL"}}
        exclude.residue_numbers = 3
        exclude.bonded_to.names = ["O1"]
        1H.abundance = 1.0
        1H.hyperfine.values = [1.0, 1.0, -2.0]
        1H.hyperfine.axes.x = {from_bonded_to = "a", to = "b"}
        1H.hyperfine.stevens_coefficients = [[0.0, 0.0, 0.0, 0.0, 0.0]]
        e.g_matrix = [2.0]
        [[groups]]
        name = "o"
        exclude.bonded_to.residues = ["SOL"]
        "##);
    assert_eq!(errors, Vec::<String>::new());
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_unknown_keys_are_reported_together(){
    let errors = key_errors(r##"
        lebedev = 1
        [structure]
        radius = 10
        radiuss = 10
        [detected_spin]
        positon = [0.0, 0.0, 0.0]
        g_matrix.valeus = [2.0, 2.0, 2.0]
        g_matrix.axes.w = [1.0, 0.0, 0.0]
        g_matrix.axes.x = {form = "a"}
        [[groups]]
        name = "a"
        include.elements = ["H"]
        [[groups]]
        name = "b"
        exclude.elemnts = ["H"]
        include.bonded_to.distance = [0, 1]
        1H.abundanse = 1.0
        1H.hyperfine.axes.q = [1.0, 0.0, 0.0]
        1h.abundance = 1.0
        "##);
    assert_eq!(errors, vec![
        "detected_spin.g_matrix.axes.w: unknown key".to_string(),
        "detected_spin.g_matrix.axes.x.form: unknown key".to_string(),
        "detected_spin.g_matrix.valeus: unknown key".to_string(),
        "detected_spin.positon: unknown key".to_string(),
        "groups[1].1H.abundanse (group \"b\"): unknown key".to_string(),
        "groups[1].1H.hyperfine.axes.q (group \"b\"): unknown key".to_string(),
        "groups[1].1h (group \"b\"): unknown key".to_string(),
        "groups[1].exclude.elemnts (group \"b\"): unknown key".to_string(),
        "groups[1].include.bonded_to.distance (group \"b\"): not allowed in \
a bonded_to table of [[groups]] (include.bonded_to, exclude.bonded_to); this \
key belongs in [units] or the include table of [[groups]] or \
[cluster_expansion.clustering.pair_cutoffs]".to_string(),
        "lebedev: unknown key".to_string(),
        "structure.radiuss: unknown key".to_string(),
    ]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_misplaced_keys_name_their_table(){
    let errors = key_errors(r##"
        1H.abundance = 1.0
        [orientations]
        temperature = 50
        [detected_spin]
        radius = 10
        [output]
        number_runs = 2
        [[groups]]
        name = "a"
        elements = ["H"]
        exclude.distance = [0, 4]
        "##);
    assert_eq!(errors, vec![
        "1H: not allowed in the top level; this key belongs in [[groups]]"
          .to_string(),
        "detected_spin.radius: not allowed in [detected_spin]; this key \
belongs in [structure]".to_string(),
        "groups[0].elements (group \"a\"): not allowed in [[groups]]; this \
key belongs in the include table of [[groups]] or the exclude table of \
[[groups]] or a bonded_to table of [[groups]] (include.bonded_to, \
exclude.bonded_to)".to_string(),
        "groups[0].exclude.distance (group \"a\"): not allowed in the \
exclude table of [[groups]]; this key belongs in [units] or the include \
table of [[groups]] or [cluster_expansion.clustering.pair_cutoffs]".to_string(),
        "orientations.temperature: not allowed in [orientations]; this key \
belongs in [conditions]".to_string(),
        "output.number_runs: not allowed in [output]; this key belongs in \
[run]".to_string(),
    ]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_old_keys_get_migration_hints(){
    // The top level of assets/1omp_K26R1_0.003228966616048703.toml before 
    // the input redesign, with the removed keys.
    let errors = key_errors(r##"
        input_structure_file = "a.pdb"
        populations = "thermal"
        cluster_method = "cce"
        ensemble_cce = true
        magnetic_field = 1.2
        max_cluster_size = 2
        partitioning = "exchange_groups"
        pulse_sequence = "hahn"
        run_name = "system"
        unit_of_time = "ns"
        min_cell_size = 1
        [pair_cutoffs]
        hahn_mod_depth = 2e-7
        [detected_spin]
        electric_quadrupole = [1.0, 1.0, -2.0]
        [output]
        cluster_signals = true
        config = true
        [[pulses]]
        name = "pi"
        [cluster_expansion.clustering.pair_cutoffs]
        coupling = 1.0
        "##);
    assert_eq!(errors, vec![
        "cluster_expansion.clustering.pair_cutoffs.coupling: removed; it was \
never read; the coupling cutoff is coupling_xx_yy".to_string(),
        "cluster_method: moved to cluster_expansion.method".to_string(),
        "detected_spin.electric_quadrupole: removed; it was read but never \
used; the isotope key electric_quadrupole of [[groups]] is unchanged"
          .to_string(),
        "ensemble_cce: moved to cluster_expansion.ensemble".to_string(),
        "input_structure_file: moved to structure.file".to_string(),
        "magnetic_field: moved to conditions.magnetic_field".to_string(),
        "max_cluster_size: moved to cluster_expansion.max_cluster_size"
          .to_string(),
        "min_cell_size: removed; it was read but never used".to_string(),
        "output.cluster_signals: removed; it was read but never used; \
cluster signals are written with auxiliary_signals".to_string(),
        "output.config: removed; it was read but never used".to_string(),
        "pair_cutoffs: moved to cluster_expansion.clustering.pair_cutoffs"
          .to_string(),
        "partitioning: moved to cluster_expansion.clustering.partitioning"
          .to_string(),
        "populations: moved to cluster_expansion.populations".to_string(),
        "pulse_sequence: moved to pulse_sequence.sequence".to_string(),
        "pulses: removed; [[pulses]] was read but never used; the pulses are \
built from detected_spin.transition".to_string(),
        "run_name: moved to output.run_name".to_string(),
        "unit_of_time: moved to units.time".to_string(),
    ]);

    // Every moved key names a key that the new layout accepts.
    for (old, new) in MOVED_TOP_LEVEL_KEYS.iter(){
      let mut schema = Schema::TopLevel;
      let parts = new.split('.').collect::<Vec::<&str>>();
      for (ii, part) in parts.iter().enumerate(){
        assert!(schema.allows(part), "{} -> {}", old, new);
        if ii + 1 < parts.len(){
          schema = schema.child(part).unwrap().0;
        }
      }
      // The old key is not a key of the new top level.
      assert!(!Schema::TopLevel.allows(old), "{}", old);
    }
    // Removed keys are not accepted where they were.
    for (schema, key, _) in REMOVED_KEYS.iter(){
      assert!(!schema.allows(key), "{:?}: {}", schema, key);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_old_selection_keys_get_migration_hints(){
    let errors = key_errors(r##"
        [[groups]]
        name = "pro_ex"
        selection.not_residues = ["R1C"]
        selection.elements = ["H"]
        selection.bonded_elements = ["O","N"]
        selection.not_bonded_elements = ["C"]
        selection.residue_sequence_numbers = [1]
        selection.within_distance = 4.0
        selection.not_within_distance = 1.0
        selection.cell_type = "primary_cell"
        selection.primary_cell = true
        selection.cell_ids = [0]
        selection.not_cell_ids = [1]
        selection.elemnts = ["H"]
        [[groups]]
        name = "empty"
        selection = {}
        [[groups]]
        name = "primary_only"
        selection.cell_type = "primary_cell"
        [[groups]]
        name = "extra_only"
        selection.primary_cell = false
        "##);
    let g = "groups[0].selection";
    let n = "(group \"pro_ex\")";
    assert_eq!(errors, vec![
        format!("{}.bonded_elements {}: moved to \
groups[0].include.bonded_to.elements", g, n),
        format!("{}.cell_ids {}: removed; selection by unit-cell index is no \
longer supported", g, n),
        format!("{}.cell_type {}: replaced by groups[0].cells = \"all\", \
\"primary\", or \"extra\" (was \"all_cells\", \"primary_cell\", \
\"extracells\")", g, n),
        format!("{}.elements {}: moved to groups[0].include.elements", g, n),
        format!("{}.elemnts {}: unknown key; the selection table is replaced \
by groups[0].include, groups[0].exclude, and groups[0].cells", g, n),
        format!("{}.not_bonded_elements {}: moved to \
groups[0].exclude.bonded_to.elements", g, n),
        format!("{}.not_cell_ids {}: removed; selection by unit-cell index is \
no longer supported", g, n),
        format!("{}.not_residues {}: moved to groups[0].exclude.residues", 
            g, n),
        format!("{}.not_within_distance {}: replaced by \
groups[0].include.distance = [min, inf]", g, n),
        format!("{}.primary_cell {}: replaced by groups[0].cells = \
\"primary\" (true) or \"extra\" (false)", g, n),
        format!("{}.residue_sequence_numbers {}: moved to \
groups[0].include.residue_numbers", g, n),
        format!("{}.within_distance {}: replaced by \
groups[0].include.distance = [0.0, max]", g, n),
        "groups[1].selection (group \"empty\"): replaced by \
groups[1].include, groups[1].exclude, and groups[1].cells".to_string(),
        // A selection with only cell keys selected every atom of its cells.
        "groups[2].selection.cell_type (group \"primary_only\"): replaced by \
groups[2].cells = \"all\", \"primary\", or \"extra\" (was \"all_cells\", \
\"primary_cell\", \"extracells\"), together with groups[2].include = {} \
to select all atoms as before".to_string(),
        "groups[3].selection.primary_cell (group \"extra_only\"): replaced \
by groups[3].cells = \"primary\" (true) or \"extra\" (false), together \
with groups[3].include = {} to select all atoms as before".to_string(),
    ]);

    // Every attribute key of the old selection table, in its four forms, 
    // moves to a key that the new tables accept.
    for (old, new) in OLD_SELECTION_ATTRIBUTES.iter(){
      for (prefix, table, schema) in [("", "include", Schema::Include),
          ("not_", "exclude", Schema::Exclude),
          ("bonded_", "include.bonded_to", Schema::BondedTo),
          ("not_bonded_", "exclude.bonded_to", Schema::BondedTo)]
      {
        assert!(schema.allows(new), "{}", new);
        assert_eq!(old_selection_key_hint("groups[0]", 
              &format!("{}{}", prefix, old), false),
            format!("moved to groups[0].{}.{}", table, new));
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_key_check_runs_when_toml_is_parsed(){
    let result = Config::from_toml_string(r##"
        [orientations]
        grid = "lebedev"
        numbr = 170
        "##);
    assert_eq!(result.err(), Some(CluEError::InvalidTOMLKeys(vec![
          "orientations.numbr: unknown key".to_string()])));

    // A TOML syntax error is reported before the key check.
    let result = Config::from_toml_string("structure.radius = ");
    assert!(matches!(result, Err(CluEError::CannotReadTOMLFile(_))));

    // The message lists every key on its own line.
    let error = CluEError::InvalidTOMLKeys(vec!["a: unknown key".to_string(),
        "b: unknown key".to_string()]);
    assert_eq!(error.to_string(), "the input contains keys that are not \
allowed where they appear:\n  a: unknown key\n  b: unknown key");
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_values_of_other_types_are_not_descended(){
    // Tensors and axes as arrays or strings, and tables of the wrong type,
    // are left to the parser.
    let errors = key_errors(r##"
        groups = {name = "a"}
        detected_spin = 1
        [orientations]
        vector = [0.0, 0.0, 1.0]
        "##);
    assert_eq!(errors, Vec::<String>::new());
    let errors = key_errors(r##"
        [[groups]]
        name = "b"
        include = 1
        exclude.bonded_to = "C"
        "##);
    assert_eq!(errors, Vec::<String>::new());
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_is_isotope_key(){
    // Every isotope that a group reads is accepted.
    for isotope in Isotope::iter(){
      assert!(is_isotope_key(&isotope.to_string()), "{}", isotope);
    }
    assert!(is_isotope_key("e"));
    for key in ["1h", "H", "H1", "", "include", "lebedev"]{
      assert!(!is_isotope_key(key), "{}", key);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_schema_keys_are_unique(){
    for schema in ALL_SCHEMAS.iter(){
      let keys = schema.keys();
      for (ii, key) in keys.iter().enumerate(){
        assert!(!keys[ii+1..].contains(key), "{:?}: {}", schema, key);
        // Fixed keys are never isotope names.
        assert!(!is_isotope_key(key), "{:?}: {}", schema, key);
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_every_schema_is_reachable(){
    // Every table kind in ALL_SCHEMAS is the child of an allowed key, so 
    // that its keys are checked; this also catches a missing entry in 
    // ALL_SCHEMAS, since a child that is not in it fails here.
    let mut reached = vec![Schema::TopLevel];
    let mut ii = 0;
    while ii < reached.len(){
      let schema = reached[ii];
      let keys = schema.keys().iter().map(|k| k.to_string())
        .chain(std::iter::once("1H".to_string()))
        .collect::<Vec::<String>>();
      for key in keys.iter(){
        if !schema.allows(key){
          continue;
        }
        if let Some((child, _)) = schema.child(key){
          assert!(ALL_SCHEMAS.contains(&child), "{:?}", child);
          if !reached.contains(&child){
            reached.push(child);
          }
        }
      }
      ii += 1;
    }
    for schema in ALL_SCHEMAS.iter(){
      assert!(reached.contains(schema), "{:?} is not reachable", schema);
    }
  }
  //----------------------------------------------------------------------------
  // This function returns the problems that the completeness check reports
  // for `config`.
  fn problems_of(config: &Config) -> Vec::<String>{
    match config.check_completeness(){
      Ok(()) => Vec::<String>::new(),
      Err(CluEError::IncompleteConfig(problems)) => problems,
      Err(err) => panic!("expected IncompleteConfig, got {:?}", err),
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_assets_are_complete(){
    for filename in ["assets/1omp_K26R1_0.003228966616048703.toml",
        "assets/1omp_K26R1_0.003228966616048703_gcce.toml"]
    {
      let config = Config::from_toml_file(filename).unwrap();
      assert_eq!(problems_of(&config), Vec::<String>::new(), "{}", filename);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_missing_settings_are_reported_together(){
    let config = Config::from_toml_string("").unwrap();
    assert_eq!(problems_of(&config), vec![
        "structure.file is not set".to_string(),
        "structure.radius is not set".to_string(),
        "conditions.magnetic_field is not set".to_string(),
        "detected_spin.position is not set".to_string(),
        "pulse_sequence.sequence is not set".to_string(),
        "pulse_sequence.number_timepoints is not set".to_string(),
        "pulse_sequence.tau_increments is not set".to_string(),
        "cluster_expansion.method is not set".to_string(),
        "cluster_expansion.max_cluster_size is not set".to_string(),
    ]);

    // units.time is set when the input is parsed; 
    // a Config built in code can lack it.
    let mut config = Config::from_toml_file(
        "assets/1omp_K26R1_0.003228966616048703.toml").unwrap();
    config.unit_of_time_to_seconds = None;
    assert_eq!(problems_of(&config), vec![
        "units.time is not set".to_string()]);

    let message = CluEError::IncompleteConfig(vec!["a is not set".to_string(),
        "b is not set".to_string()]).to_string();
    assert_eq!(message, 
        "the input is incomplete:\n  a is not set\n  b is not set");
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_conditionally_required_settings(){
    let cce_file = "assets/1omp_K26R1_0.003228966616048703.toml";
    let gcce_file = "assets/1omp_K26R1_0.003228966616048703_gcce.toml";

    // populations = "thermal" needs a temperature.
    let mut config = Config::from_toml_file(cce_file).unwrap();
    config.temperature = None;
    assert_eq!(problems_of(&config), vec![
        "conditions.temperature is required for \
cluster_expansion.populations = \"thermal\"".to_string()]);

    // A field set after parsing: gCCE needs the detected spin operators.
    let mut config = Config::from_toml_file(cce_file).unwrap();
    config.cluster_method = Some(ClusterMethod::GCCE);
    assert_eq!(problems_of(&config), vec![
        "detected_spin.detection_operator is required for \
cluster_expansion.method = \"gcce\"".to_string(),
        "detected_spin.density_matrix is required for \
cluster_expansion.method = \"gcce\"".to_string(),
    ]);
    // A detected spin with one state has no coherence; gCCE would panic
    // on the missing detected spin operators.
    config.detected_spin_multiplicity = Some(1);
    assert_eq!(problems_of(&config), vec![
        "detected_spin.multiplicity must be at least 2, but is 1".to_string()]);

    // The gCCE asset uses a thermal detected spin density matrix and 
    // thermalize_cluster_populations = true.
    let mut config = Config::from_toml_file(gcce_file).unwrap();
    config.temperature = None;
    config.ensemble_cce = Some(false);
    assert_eq!(problems_of(&config), vec![
        "conditions.temperature is required for \
cluster_expansion.populations = \"thermal\"".to_string(),
        "conditions.temperature is required for detected_spin.density_matrix \
= \"thermal\"".to_string(),
        "cluster_expansion.thermalize_cluster_populations = true requires \
cluster_expansion.ensemble = true".to_string(),
        "conditions.temperature is required for \
cluster_expansion.thermalize_cluster_populations = true".to_string(),
    ]);

    // These settings are used only by gCCE.
    config.cluster_method = Some(ClusterMethod::CCE);
    config.cluster_populations = Some(ClusterPopulations::Uniform);
    assert_eq!(problems_of(&config), Vec::<String>::new());
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_detected_spin_levels(){
    let cce_str = std::fs::read_to_string(
        "assets/1omp_K26R1_0.003228966616048703.toml").unwrap();
    let gcce_file = "assets/1omp_K26R1_0.003228966616048703_gcce.toml";

    // Without a transition, S = 1 has no default transition, since it is 
    // not clear which pair of levels the pulses drive.
    let triplet = cce_str.replace("[detected_spin]",
        "[detected_spin]\n  multiplicity = 3");
    let mut config = Config::from_toml_string(&triplet).unwrap();
    assert_eq!(config.detected_spin_transition, None);
    assert_eq!(problems_of(&config), vec![
        "detected_spin.transition is required for \
detected_spin.multiplicity = 3".to_string()]);

    // gCCE also builds its pulses from the transition.
    config.cluster_method = Some(ClusterMethod::GCCE);
    config.detection_operator = Some(DetectionOp::Transition(0,1));
    config.detected_population = Some(DetectedPopulation::Thermal);
    assert_eq!(problems_of(&config), vec![
        "detected_spin.transition is required for \
detected_spin.multiplicity = 3".to_string()]);

    // APPA and LCE do not use the transition.
    for method in [ClusterMethod::APPA, ClusterMethod::LCE]{
      config.cluster_method = Some(method);
      assert_eq!(problems_of(&config), Vec::<String>::new());
    }

    // A transition between two different levels 0..=2.
    config.cluster_method = Some(ClusterMethod::CCE);
    for levels in [[0,1], [1,2], [2,0]]{
      config.detected_spin_transition = Some(levels);
      assert_eq!(problems_of(&config), Vec::<String>::new(), "{:?}", levels);
    }
    for levels in [[1,1], [0,3], [3,0]]{
      config.detected_spin_transition = Some(levels);
      assert_eq!(problems_of(&config), vec![format!(
            "detected_spin.transition = [{}, {}] must name two different \
levels from 0 to 2", levels[0], levels[1])], "{:?}", levels);
    }

    // A detected spin with one level has no coherence, for every method;
    // the transition is then not checked.
    for method in [ClusterMethod::CCE, ClusterMethod::GCCE, 
        ClusterMethod::APPA, ClusterMethod::LCE]
    {
      config.cluster_method = Some(method);
      config.detected_spin_multiplicity = Some(1);
      assert_eq!(problems_of(&config), vec![
          "detected_spin.multiplicity must be at least 2, but is 1"
            .to_string()]);
    }

    // A spin-0 identity gives multiplicity 1 by default.
    let spin_0 = cce_str.replace("[detected_spin]",
        "[detected_spin]\n  identity = \"12C\"");
    let config = Config::from_toml_string(&spin_0).unwrap();
    assert_eq!(config.detected_spin_multiplicity, Some(1));
    assert_eq!(problems_of(&config), vec![
        "detected_spin.multiplicity must be at least 2, but is 1".to_string()]);

    // The levels of a detection operator [i, j] must exist; [i, i] is the
    // projector onto level i and is allowed.
    let mut config = Config::from_toml_file(gcce_file).unwrap();
    config.detection_operator = Some(DetectionOp::Transition(1,1));
    assert_eq!(problems_of(&config), Vec::<String>::new());
    config.detection_operator = Some(DetectionOp::Transition(0,2));
    assert_eq!(problems_of(&config), vec![
        "detected_spin.detection_operator = [0, 2] must name levels from 0 \
to 1".to_string()]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_mean_fields_with_ensemble_cce(){
    // The asset uses cluster_expansion.ensemble = true and mean_fields = false.
    let mut config = Config::from_toml_file(
        "assets/1omp_K26R1_0.003228966616048703.toml").unwrap();
    config.mean_fields = Some(true);
    assert_eq!(problems_of(&config), vec![
        "mean field averaging is not implemented for ensemble CCE \
(cluster_expansion.mean_fields = true with cluster_expansion.ensemble = true)"
          .to_string()]);
    config.ensemble_cce = Some(false);
    assert_eq!(problems_of(&config), Vec::<String>::new());
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_missing_files_are_reported(){
    let mut config = Config::from_toml_file(
        "assets/1omp_K26R1_0.003228966616048703.toml").unwrap();
    config.input_structure_file = Some("assets/no_such_file.pdb".to_string());
    config.partition_table = Some(PartitionTableConfig::CSV(
          "assets".to_string()));
    assert_eq!(problems_of(&config), vec![
        "structure.file: \"assets/no_such_file.pdb\" does not exist"
          .to_string(),
        "cluster_expansion.clustering.partition_table: \"assets\" is not a \
file".to_string(),
    ]);

    // The partition table is not read for clusters from a file.
    config.input_structure_file = Some(
        "assets/1omp_K26R1_0.003228966616048703.pdb".to_string());
    config.cluster_source = Some(ClusterSource::File(
          "assets/no_clusters.toml".to_string()));
    assert_eq!(problems_of(&config), vec![
        "cluster_expansion.clustering.file: \"assets/no_clusters.toml\" does \
not exist".to_string()]);

    config.cluster_source = Some(ClusterSource::OrientationsDir(
          "assets/TEMPO.pdb".to_string()));
    assert_eq!(problems_of(&config), vec![
        "cluster_expansion.clustering.file: \"assets/TEMPO.pdb\" is not a \
directory".to_string()]);

    config.cluster_source = Some(ClusterSource::OrientationsDir(
          "assets".to_string()));
    assert_eq!(problems_of(&config), Vec::<String>::new());
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_run_checks_completeness_first(){
    let mut config = Config::from_toml_file(
        "assets/1omp_K26R1_0.003228966616048703.toml").unwrap();
    config.radius = None;
    // No output directory is created for an incomplete config.
    config.output_directory = Some(
        "CluE-test_run_checks_completeness_first".to_string());
    assert_eq!(crate::run(config).err(), Some(CluEError::IncompleteConfig(
            vec!["structure.radius is not set".to_string()])));
    assert!(!Path::new("CluE-test_run_checks_completeness_first").exists());
  }
  //----------------------------------------------------------------------------
}
