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
// The key names are the ones read by `ConfigTOML` (src/config/config_toml.rs),
// `Config::set_from_config_toml` (src/config.rs), `ParticleConfig`,
// `ParticleProperties`, `IsotopeProperties`, `EigSpecifier`
// (src/config/particle_config.rs), and `ParticleFilter` and
// `VectorSpecifier` (src/structure/particle_filter.rs).

const TOP_LEVEL_KEYS: [&str; 45] = [
  "clash_distance",
  "clash_distance_pbc",
  "cluster_batch_size",
  "cluster_method",
  "cluster_source",
  "connect_exchange_groups",
  "ensemble_cce",
  "input_structure_file",
  "kmeans_iter_cap",
  "kmeans_size",
  "magnetic_field",
  "max_cell_size",
  "max_cluster_size",
  "max_spherical_tensor_rank",
  "max_spins",
  "mean_fields",
  "min_cell_size",
  "number_runs",
  "number_timepoints",
  "number_timepoints2",
  "output_directory",
  "partition_table",
  "partitioning",
  "pdb_model_index",
  "populations",
  "pulse_sequence",
  "radius",
  "replicate_unit_cell",
  "rng_seed",
  "run_in_parallel",
  "run_name",
  "tau_increments",
  "tau2_increments",
  "temperature",
  "thermalize_cluster_populations",
  "unit_of_distance",
  "unit_of_energy",
  "unit_of_magnetic_field",
  "unit_of_time",
  // Tables.
  "detected_spin",
  "groups",
  "orientations",
  "output",
  "pair_cutoffs",
  "pulses",
];

const DETECTED_SPIN_KEYS: [&str; 11] = [
  "density_matrix",
  "detection_operator",
  "electric_quadrupole",
  "frame",
  "g_matrix",
  "identity",
  "multiplicity",
  "position",
  "transition",
  "weight_tolerance",
  "zerofield",
];

const ORIENTATIONS_KEYS: [&str; 6] = [
  "file",
  "grid",
  "number",
  "vector",
  "vector_grid",
  "weight_tolerance",
];

const OUTPUT_KEYS: [&str; 15] = [
  KEY_OUT_AUX_SIGS,
  KEY_OUT_BATH,
  KEY_OUT_CLU_SIGS,
  KEY_OUT_CLUSTERS,
  KEY_OUT_CONFIG,
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

// "coupling" is accepted but not read (see CluE_input_redesign.md,
// "Removed keys"); it is removed with the new input format.
const PAIR_CUTOFFS_KEYS: [&str; 8] = [
  "coupling",
  KEY_CUTOFF_COUPLING,
  KEY_CUTOFF_DELTA_HF,
  KEY_CUTOFF_DIPOLE_PERP,
  KEY_CUTOFF_DISTANCE,
  KEY_CUTOFF_HAHN_MOD_DEPTH,
  KEY_CUTOFF_HAHN_TAYLOR_4,
  KEY_CUTOFF_DELTA_ZEEMAN,
];

// [[pulses]] is accepted but not read (see CluE_input_redesign.md,
// "Removed keys"); it is removed with the new input format.
const PULSE_KEYS: [&str; 2] = [KEY_NAME, "matrix"];

// The isotope tables of a group are checked by `is_isotope_key()`.
const GROUP_KEYS: [&str; 5] = [
  KEY_NAME,
  KEY_ABUNDANCE_TOL,
  KEY_DROP_PROB,
  KEY_ISO_COSUBSTITUTE,
  KEY_SELECTION,
];

const SELECTION_KEYS: [&str; 34] = [
  KEY_SELE_CHAIN_IDS,
  KEY_SELE_ELEMENTS,
  KEY_SELE_INDICES,
  KEY_SELE_NAMES,
  KEY_SELE_RESIDUES,
  KEY_SELE_RES_SEQ_NUMS,
  KEY_SELE_SERIALS,
  KEY_SELE_NOT_CHAIN_IDS,
  KEY_SELE_NOT_ELEMENTS,
  KEY_SELE_NOT_INDICES,
  KEY_SELE_NOT_NAMES,
  KEY_SELE_NOT_RESIDUES,
  KEY_SELE_NOT_RES_SEQ_NUMS,
  KEY_SELE_NOT_SERIALS,
  KEY_SELE_BONDED_CHAIN_IDS,
  KEY_SELE_BONDED_ELEMENTS,
  KEY_SELE_BONDED_INDICES,
  KEY_SELE_BONDED_NAMES,
  KEY_SELE_BONDED_RESIDUES,
  KEY_SELE_BONDED_RES_SEQ_NUMS,
  KEY_SELE_BONDED_SERIALS,
  KEY_SELE_NOT_BONDED_CHAIN_IDS,
  KEY_SELE_NOT_BONDED_ELEMENTS,
  KEY_SELE_NOT_BONDED_INDICES,
  KEY_SELE_NOT_BONDED_NAMES,
  KEY_SELE_NOT_BONDED_RESIDUES,
  KEY_SELE_NOT_BONDED_RES_SEQ_NUMS,
  KEY_SELE_NOT_BONDED_SERIALS,
  KEY_SELE_WITHIN_DISTANCE,
  KEY_SELE_NOT_WITHIN_DISTANCE,
  KEY_SELE_CELL_IDS,
  KEY_SELE_NOT_CELL_IDS,
  KEY_SELE_PRIMARY_CELL,
  KEY_CELL_TYPE,
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
/// `Schema` lists the kinds of tables in the input file.
#[derive(Debug,Clone,Copy,PartialEq)]
enum Schema{
  TopLevel,
  DetectedSpin,
  Orientations,
  Output,
  PairCutoffs,
  Pulse,
  Group,
  Selection,
  Isotope,
  Tensor,
  Axes,
  Axis,
}

const ALL_SCHEMAS: [Schema; 12] = [
  Schema::TopLevel,
  Schema::DetectedSpin,
  Schema::Orientations,
  Schema::Output,
  Schema::PairCutoffs,
  Schema::Pulse,
  Schema::Group,
  Schema::Selection,
  Schema::Isotope,
  Schema::Tensor,
  Schema::Axes,
  Schema::Axis,
];

impl Schema{
  //----------------------------------------------------------------------------
  // This function returns the fixed keys of the table.
  fn keys(&self) -> &'static [&'static str]{
    match self{
      Schema::TopLevel => &TOP_LEVEL_KEYS,
      Schema::DetectedSpin => &DETECTED_SPIN_KEYS,
      Schema::Orientations => &ORIENTATIONS_KEYS,
      Schema::Output => &OUTPUT_KEYS,
      Schema::PairCutoffs => &PAIR_CUTOFFS_KEYS,
      Schema::Pulse => &PULSE_KEYS,
      Schema::Group => &GROUP_KEYS,
      Schema::Selection => &SELECTION_KEYS,
      Schema::Isotope => &ISOTOPE_KEYS,
      Schema::Tensor => &TENSOR_KEYS,
      Schema::Axes => &AXES_KEYS,
      Schema::Axis => &AXIS_KEYS,
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
        "detected_spin" => Some((Schema::DetectedSpin,false)),
        "orientations" => Some((Schema::Orientations,false)),
        "output" => Some((Schema::Output,false)),
        "pair_cutoffs" => Some((Schema::PairCutoffs,false)),
        "groups" => Some((Schema::Group,true)),
        "pulses" => Some((Schema::Pulse,true)),
        _ => None,
      },
      Schema::DetectedSpin => match key{
        "g_matrix" | "electric_quadrupole" | "zerofield"
          => Some((Schema::Tensor,false)),
        _ => None,
      },
      Schema::Group => {
        if key == KEY_SELECTION{
          Some((Schema::Selection,false))
        }else if is_isotope_key(key){
          Some((Schema::Isotope,false))
        }else{
          None
        }
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
      Schema::Orientations | Schema::Output | Schema::PairCutoffs
        | Schema::Pulse | Schema::Selection | Schema::Axis => None,
    }
  }
  //----------------------------------------------------------------------------
  // This function names the table in messages.
  fn description(&self) -> &'static str{
    match self{
      Schema::TopLevel => "the top level",
      Schema::DetectedSpin => "[detected_spin]",
      Schema::Orientations => "[orientations]",
      Schema::Output => "[output]",
      Schema::PairCutoffs => "[pair_cutoffs]",
      Schema::Pulse => "[[pulses]]",
      Schema::Group => "[[groups]]",
      Schema::Selection => "the selection table of [[groups]]",
      Schema::Isotope => "an isotope table of [[groups]]",
      Schema::Tensor => "a tensor table (values, axes, \
spherical_coefficients, stevens_coefficients)",
      Schema::Axes => "the axes table of a tensor (x, y, z)",
      Schema::Axis => "an axis table of a tensor (from, from_bonded_to, \
from_same_molecule_as, to, to_bonded_to, to_same_molecule_as)",
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
/// path, for example `groups[2].selection.not_elemnts`;
/// `groups[i]` counts the [[groups]] entries from 0.
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
    let key_path = if path.is_empty(){
      key.to_string()
    }else{
      format!("{}.{}", path, key)
    };

    if !schema.allows(key){
      errors.push(not_allowed_message(&key_path, key, schema, group));
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
// This function writes the message for `key` at `key_path`, which is not
// allowed in a table of kind `schema`.
fn not_allowed_message(key_path: &str, key: &str, schema: Schema,
    group: Option<&str>) -> String
{
  let mut message = key_path.to_string();
  if let Some(name) = group{
    message = format!("{} (group \"{}\")", message, name);
  }

  let places = ALL_SCHEMAS.iter()
    .filter(|s| s.allows(key))
    .map(|s| s.description())
    .collect::<Vec::<&str>>();

  if places.is_empty(){
    format!("{}: unknown key", message)
  }else{
    format!("{}: not allowed in {}; this key belongs in {}", message,
        schema.description(), places.join(" or "))
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
  /// (`partitioning` with its `kmeans_size`, `orientations.file`, and a
  /// csv `detected_spin.position`) are checked by the parser.
  /// Errors that depend on the structure (group names, overlapping groups,
  /// empty selections) are found only by a run.
  pub fn check_completeness(&self) -> Result<(),CluEError>{
    let mut problems = Vec::<String>::new();

    // Settings without a default.
    let required = [
      (self.input_structure_file.is_some(), "input_structure_file"),
      (self.radius.is_some(), "radius"),
      (self.magnetic_field.is_some(), "magnetic_field"),
      (self.detected_spin_position.is_some(), "detected_spin.position"),
      (self.pulse_sequence.is_some(), "pulse_sequence"),
      (!self.number_timepoints.is_empty(), "number_timepoints"),
      (!self.tau_increments.is_empty(), "tau_increments"),
      (self.cluster_method.is_some(), "cluster_method"),
      (self.max_cluster_size.is_some(), "max_cluster_size"),
      // Set from unit_of_time or its default when the input is parsed;
      // it can be unset only in a `Config` built in code.
      (self.unit_of_time_to_seconds.is_some(), "unit_of_time"),
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
      problems.push(format!("temperature is required for populations = \
\"{}\"", KEY_DENSITY_MATRIX_THERMAL));
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
for cluster_method = \"gcce\"".to_string()),
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
for cluster_method = \"gcce\"".to_string()),
        Some(DetectedPopulation::Thermal) if !has_temperature =>
          problems.push(format!("temperature is required for \
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
        problems.push("temperature is required for \
thermalize_cluster_populations = true".to_string());
      }
    }

    // Files and directories.
    if let Some(filename) = &self.input_structure_file{
      if let Some(problem) = check_path(filename, false){
        problems.push(format!("input_structure_file: {}", problem));
      }
    }

    match &self.cluster_source{
      Some(ClusterSource::File(filename)) => {
        if let Some(problem) = check_path(filename, false){
          problems.push(format!("cluster_source: {}", problem));
        }
      },
      Some(ClusterSource::OrientationsDir(dirname)) => {
        if let Some(problem) = check_path(dirname, true){
          problems.push(format!("cluster_source: {}", problem));
        }
      },
      // The partition table is read only when clusters are built from
      // the structure.
      Some(ClusterSource::Structure) | None => {
        if let Some(PartitionTableConfig::CSV(filename))
          = &self.partition_table
        {
          if let Some(problem) = check_path(filename, false){
            problems.push(format!("partition_table: {}", problem));
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
        radius = 10
        [detected_spin]
        position = [0.0, 0.0, 0.0]
        g_matrix.values = [2.0, 2.0, 2.0]
        g_matrix.axes.x = {from = "a", to_bonded_to = "b"}
        g_matrix.axes.y = "random"
        zerofield.values = [1.0, 1.0, -2.0]
        zerofield.axes.x = [1.0, 0.0, 0.0]
        zerofield.axes.z.from_same_molecule_as = "a"
        zerofield.axes.z.to_same_molecule_as = "b"
        zerofield.spherical_coefficients = [[0.0, 0.0, 0.0, 0.0, 0.0]]
        [orientations]
        grid = "lebedev"
        number = 170
        [output]
        bath = true
        partition_table = true
        detected_spin = true
        [pair_cutoffs]
        distance = 10
        coupling = 1.0
        [[pulses]]
        name = "pi"
        matrix = [[0.0, 1.0], [1.0, 0.0]]
        [[groups]]
        name = "h"
        drop_probability = 0.1
        selection = {elements = ["H"], cell_type = "all_cells"}
        1H.abundance = 1.0
        1H.hyperfine.values = [1.0, 1.0, -2.0]
        1H.hyperfine.axes.x = {from_bonded_to = "a", to = "b"}
        1H.hyperfine.stevens_coefficients = [[0.0, 0.0, 0.0, 0.0, 0.0]]
        e.g_matrix = [2.0]
        [[groups]]
        name = "o"
        selection.not_bonded_residues = ["SOL"]
        selection.within_distance = 5.0
        "##);
    assert_eq!(errors, Vec::<String>::new());
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_unknown_keys_are_reported_together(){
    let errors = key_errors(r##"
        radius = 10
        radiuss = 10
        lebedev = 1
        [detected_spin]
        positon = [0.0, 0.0, 0.0]
        g_matrix.valeus = [2.0, 2.0, 2.0]
        g_matrix.axes.w = [1.0, 0.0, 0.0]
        g_matrix.axes.x = {form = "a"}
        [[groups]]
        name = "a"
        selection.elements = ["H"]
        [[groups]]
        name = "b"
        selection.not_elemnts = ["H"]
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
        "groups[1].selection.not_elemnts (group \"b\"): unknown key"
          .to_string(),
        "lebedev: unknown key".to_string(),
        "radiuss: unknown key".to_string(),
    ]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_misplaced_keys_name_their_table(){
    let errors = key_errors(r##"
        grid = "lebedev"
        1H.abundance = 1.0
        [detected_spin]
        radius = 10
        [output]
        number_runs = 2
        [[groups]]
        name = "a"
        elements = ["H"]
        selection.name = "a"
        "##);
    assert_eq!(errors, vec![
        "1H: not allowed in the top level; this key belongs in [[groups]]"
          .to_string(),
        "detected_spin.radius: not allowed in [detected_spin]; this key \
belongs in the top level".to_string(),
        "grid: not allowed in the top level; this key belongs in \
[orientations]".to_string(),
        "groups[0].elements (group \"a\"): not allowed in [[groups]]; this \
key belongs in the selection table of [[groups]]".to_string(),
        "groups[0].selection.name (group \"a\"): not allowed in the \
selection table of [[groups]]; this key belongs in [[pulses]] or \
[[groups]]".to_string(),
        "output.number_runs: not allowed in [output]; this key belongs in \
the top level".to_string(),
    ]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_key_check_runs_when_toml_is_parsed(){
    let result = Config::from_toml_string(r##"
        radius = 10
        [orientations]
        grid = "lebedev"
        numbr = 170
        "##);
    assert_eq!(result.err(), Some(CluEError::InvalidTOMLKeys(vec![
          "orientations.numbr: unknown key".to_string()])));

    // A TOML syntax error is reported before the key check.
    let result = Config::from_toml_string("radius = ");
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
        [[pulses]]
        name = "pi"
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
    for key in ["1h", "H", "H1", "", "selection", "lebedev"]{
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
        "input_structure_file is not set".to_string(),
        "radius is not set".to_string(),
        "magnetic_field is not set".to_string(),
        "detected_spin.position is not set".to_string(),
        "pulse_sequence is not set".to_string(),
        "number_timepoints is not set".to_string(),
        "tau_increments is not set".to_string(),
        "cluster_method is not set".to_string(),
        "max_cluster_size is not set".to_string(),
    ]);

    // unit_of_time is set when the input is parsed; 
    // a Config built in code can lack it.
    let mut config = Config::from_toml_file(
        "assets/1omp_K26R1_0.003228966616048703.toml").unwrap();
    config.unit_of_time_to_seconds = None;
    assert_eq!(problems_of(&config), vec![
        "unit_of_time is not set".to_string()]);

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
        "temperature is required for populations = \"thermal\"".to_string()]);

    // A field set after parsing: gCCE needs the detected spin operators.
    let mut config = Config::from_toml_file(cce_file).unwrap();
    config.cluster_method = Some(ClusterMethod::GCCE);
    assert_eq!(problems_of(&config), vec![
        "detected_spin.detection_operator is required for cluster_method = \
\"gcce\"".to_string(),
        "detected_spin.density_matrix is required for cluster_method = \
\"gcce\"".to_string(),
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
        "temperature is required for populations = \"thermal\"".to_string(),
        "temperature is required for detected_spin.density_matrix = \
\"thermal\"".to_string(),
        "thermalize_cluster_populations = true requires ensemble_cce = true"
          .to_string(),
        "temperature is required for thermalize_cluster_populations = true"
          .to_string(),
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
    // The asset uses ensemble_cce = true and mean_fields = false.
    let mut config = Config::from_toml_file(
        "assets/1omp_K26R1_0.003228966616048703.toml").unwrap();
    config.mean_fields = Some(true);
    assert_eq!(problems_of(&config), vec![
        "mean field averaging is not implemented for ensemble CCE"
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
        "input_structure_file: \"assets/no_such_file.pdb\" does not exist"
          .to_string(),
        "partition_table: \"assets\" is not a file".to_string(),
    ]);

    // The partition table is not read for clusters from a file.
    config.input_structure_file = Some(
        "assets/1omp_K26R1_0.003228966616048703.pdb".to_string());
    config.cluster_source = Some(ClusterSource::File(
          "assets/no_clusters.toml".to_string()));
    assert_eq!(problems_of(&config), vec![
        "cluster_source: \"assets/no_clusters.toml\" does not exist"
          .to_string()]);

    config.cluster_source = Some(ClusterSource::OrientationsDir(
          "assets/TEMPO.pdb".to_string()));
    assert_eq!(problems_of(&config), vec![
        "cluster_source: \"assets/TEMPO.pdb\" is not a directory".to_string()]);

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
            vec!["radius is not set".to_string()])));
    assert!(!Path::new("CluE-test_run_checks_completeness_first").exists());
  }
  //----------------------------------------------------------------------------
}
