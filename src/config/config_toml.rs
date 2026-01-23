use crate::clue_errors::CluEError;
use crate::io::FromTOMLString;

use std::fmt;
use std::collections::HashMap;

use serde::{Serialize,Deserialize};

const ALLOWED_KEYS: [&str;120] = [
  "abundance",
  "active",
  "auxiliary_signals",
  "axes",
  "bath",
  "bonded_elements",
  "bonded_indices",
  "bonded_names",
  "bonded_residues",
  "bonded_residue_sequence_numbers",
  "bonded_serials",
  "c3_tunnel_splitting",
  "cell_ids",
  "cell_type",
  "cosubstitute",
  "coupling",
  "coupling_xx_yy",
  "clash_distance",
  "clash_distance_pbc",
  "cluster_batch_size",
  "cluster_method",
  "cluster_source",
  "clusters",
  "config",
  "delta_hyperfine_zz",
  "detected_spin",
  "distance",
  "drop_probability",
  "electric_quadrupole",
  "elements",
  "exchange_coupling",
  "exchange_groups",
  "file",
  "from",
  "from_bonded_to",
  "from_same_molecule_as",
  "g_matrix",
  "grid",
  "groups",
  "hahn_mod_depth",
  "hahn_taylor_4",
  "hyperfine",
  "identity",
  "indices",
  "info",
  "input_structure_file",
  "isotpes",
  "lebedev",
  "magnetic_field",
  "max_cell_size",
  "max_cluster_size",
  "max_spins",
  "methyl_partitions",
  "min_cell_size",
  "name",
  "names",
  "not_bonded_elements",
  "not_bonded_indices",
  "not_bonded_names",
  "not_bonded_residues",
  "not_bonded_residue_sequence_numbers",
  "not_bonded_serials",
  "not_cell_ids",
  "not_elements",
  "not_indices",
  "not_isotopes",
  "not_names",
  "not_residues",
  "not_residue_sequence_numbers",
  "not_serials",
  "not_within_distance",
  "number",
  "number_runs",
  "number_timepoints",
  "orientation_grid",
  "orientations",
  "orientation_signals",
  "output",
  "output_directory",
  "pair_cutoffs",
  "partitioning",
  "partition_table",
  "pdb_model_index",
  "point_dipole_perpendicular",
  "populations",
  "position",
  "pulse_sequence",
  "radius",
  "random",
  "replicate_unit_cell",
  "residues",
  "residue_sequence_numbers",
  "rng_seed",
  "run_name",
  "sans_spin_signals",
  "save_dir",
  "selection",
  "serials",
  "singles",
  "spin_multiplicity",
  "structure_pdb",
  "tau_increments",
  "temperature",
  "tensors",
  "thermal",
  "to",
  "to_bonded_to",
  "to_same_molecule_as",
  "transition",
  "unit_of_distance",
  "unit_of_energy",
  "unit_of_magnetic_field",
  "unit_of_time",
  "values",
  "vector",
  "vector_grid",
  "within_distance",
  "x",
  "y",
  "z",
];

pub const DEFAULT_UNIT_ENERGY: &str = "MHz";
pub const DEFAULT_UNIT_DISTANCE: &str = "Å";
pub const DEFAULT_UNIT_MAGNETIC_FIELD: &str = "T";
pub const DEFAULT_UNIT_TIME: &str = "μs";

// General Keys
pub const KEY_CUTOFF_COUPLING: &str = "coupling_xx_yy";
pub const KEY_CUTOFF_DELTA_HF: &str = "delta_hyperfine_zz";
pub const KEY_CUTOFF_DIPOLE_PERP: &str = "point_dipole_perpendicular";
pub const KEY_CUTOFF_DISTANCE: &str = "distance";
pub const KEY_CUTOFF_HAHN_MOD_DEPTH: &str = "hahn_mod_depth";
pub const KEY_CUTOFF_HAHN_TAYLOR_4: &str = "hahn_taylor_4";

pub const KEY_OUT_AUX_SIGS: &str = "auxiliary_signals";
pub const KEY_OUT_BATH: &str = "bath";
pub const KEY_OUT_DET_SPIN: &str = "detected_spin";
pub const KEY_OUT_CLUSTERS: &str = "clusters";
pub const KEY_OUT_CONFIG: &str = "config";
pub const KEY_OUT_INFO: &str = "info";
pub const KEY_OUT_EXCHANGE_GROUPS: &str = "exchange_groups";
pub const KEY_OUT_METHYL_PARTITIONS: &str = "methyl_partitions";
pub const KEY_OUT_ORI_GRID: &str = "orientation_grid";
pub const KEY_OUT_ORI_SIGS: &str = "orientation_signals";
pub const KEY_OUT_PART_TAB: &str = "partition_table";
pub const KEY_OUT_SANS_SPIN_SIGS: &str = "sans_spin_signals";
pub const KEY_OUT_STRUC_PDB: &str = "structure_pdb";
pub const KEY_OUT_TENSORS: &str = "tensors";

pub const KEY_DENSITY_MATRIX_ID: &str  = "identity";
pub const KEY_DENSITY_MATRIX_THERMAL: &str  = "thermal";

pub const KEY_PARTITION_EX_GROUPS: &str = "exchange_groups";
pub const KEY_PARTITION_PARTICLE: &str = "singles";

pub const KEY_ORI_LEBEDEV: &str = "lebedev";
pub const KEY_ORI_RANDOM: &str = "random";
pub const KEY_ORI_FILE: &str = "file";
pub const KEY_ORI_VECTOR: &str = "vector";
pub const KEY_ORI_VECTORGRID: &str = "vector_grid";


pub const KEY_EIG_VALUES: &str = "values";                                           
pub const KEY_EIG_AXES: &str = "axes";                                               
pub const KEY_EIG_X_AXIS: &str = "x";                                                
pub const KEY_EIG_Y_AXIS: &str = "y";                                                
pub const KEY_EIG_Z_AXIS: &str = "z";

// Group Keys
pub const KEY_NAME: &str = "name";

pub const KEY_DROP_PROB: &str = "drop_probability";

pub const KEY_ISO_COSUBSTITUTE: &str = "cosubstitute";


// Filter Keys
pub const KEY_SELECTION: &str = "selection";
pub const KEY_CELL_TYPE: &str = "cell_type";


pub const KEY_SELE_INDICES: &str = "indices";
pub const KEY_SELE_NOT_INDICES: &str = "not_indices";

pub const KEY_SELE_CELL_IDS: &str = "cell_ids";
pub const KEY_SELE_NOT_CELL_IDS: &str = "not_cell_ids";

pub const KEY_SELE_ELEMENTS: &str = "elements";
pub const KEY_SELE_NOT_ELEMENTS: &str = "not_elements";

pub const KEY_SELE_SERIALS: &str = "serials";
pub const KEY_SELE_NOT_SERIALS: &str = "not_serials";

pub const KEY_SELE_NAMES: &str = "names";
pub const KEY_SELE_NOT_NAMES: &str = "not_names";

pub const KEY_SELE_RESIDUES: &str = "residues";
pub const KEY_SELE_NOT_RESIDUES: &str = "not_residues";

pub const KEY_SELE_RES_SEQ_NUMS: &str = "residue_sequence_numbers";
pub const KEY_SELE_NOT_RES_SEQ_NUMS: &str = "not_residue_sequence_numbers";

pub const KEY_SELE_ISOTOPES: &str = "isotpes";
pub const KEY_SELE_NOT_ISOTOPES: &str = "not_isotopes";

pub const KEY_SELE_BONDED_INDICES: &str = "bonded_indices";
pub const KEY_SELE_NOT_BONDED_INDICES: &str = "not_bonded_indices";

pub const KEY_SELE_WITHIN_DISTANCE: &str = "within_distance";
pub const KEY_SELE_NOT_WITHIN_DISTANCE: &str = "not_within_distance";

pub const KEY_SELE_BONDED_ELEMENTS: &str = "bonded_elements";
pub const KEY_SELE_NOT_BONDED_ELEMENTS: &str = "not_bonded_elements";

pub const KEY_SELE_BONDED_SERIALS: &str = "bonded_serials";
pub const KEY_SELE_NOT_BONDED_SERIALS: &str = "not_bonded_serials";

pub const KEY_SELE_BONDED_NAMES: &str = "bonded_names";
pub const KEY_SELE_NOT_BONDED_NAMES: &str = "not_bonded_names";

pub const KEY_SELE_BONDED_RESIDUES: &str = "bonded_residues";
pub const KEY_SELE_NOT_BONDED_RESIDUES: &str = "not_bonded_residues";

pub const KEY_SELE_BONDED_RES_SEQ_NUMS: &str = "bonded_residue_sequence_numbers";
pub const KEY_SELE_NOT_BONDED_RES_SEQ_NUMS: &str 
    = "not_bonded_residue_sequence_numbers";

// Isotope Key
pub const KEY_ISO_ABUNDACE: &str = "abundance";
pub const KEY_ISO_ACTIVE: &str = "active";
pub const KEY_ISO_G_MATRIX: &str = "g_matrix";
pub const KEY_ISO_HYPERFINE: &str = "hyperfine";
pub const KEY_ISO_ELEC_QUADRUPOLE: &str = "electric_quadrupole";
pub const KEY_ISO_EXCHANGE_COUPLING: &str = "exchange_coupling";
pub const KEY_ISO_C3_TUNNEL_SPLITTING: &str = "c3_tunnel_splitting";
pub const KEY_ISO_: &str = "";

pub const KEY_VEC_SPECIFIER_FROM: &str = "from";
pub const KEY_VEC_SPECIFIER_FROM_BONDED_TO: &str = "from_bonded_to";
pub const KEY_VEC_SPECIFIER_FROM_SAME_MOLECULE_AS: &str = "from_same_molecule_as";
pub const KEY_VEC_SPECIFIER_TO: &str = "to";
pub const KEY_VEC_SPECIFIER_TO_BONDED_TO: &str = "to_bonded_to";
pub const KEY_VEC_SPECIFIER_TO_SAME_MOLECULE_AS: &str = "to_same_molecule_as";
pub const KEY_VEC_SPECIFIER_RANDOM: &str = "random";
   
// moveed trait def to io.rs.
//pub trait FromTOMLString{
//  fn from_toml_string(s: &str) -> Result<Self,CluEError> where Self: Sized;
//}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct DetectedSpinTOML{
  pub multiplicity: Option<usize>,
  pub g_matrix: Option<toml::Value>,
  //pub g_values: Option<Vec::<f64>>,
  //pub gx: Option<toml::Value>,
  //pub gy: Option<toml::Value>,
  //pub gz: Option<toml::Value>,
  //pub position: Option<VectorSpecifierTOML>,
  pub position: Option<toml::Value>,
  pub transition: Option<[usize;2]>,
  pub identity: Option<String>,
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ParticleProperties` specifies custom particle properties.
/// 'cosubstitute` selects the set of particles that should always be the same
/// isotope when the isotopic distribution is randomized.
/// `isotopic_distribution` specifies how elements are assigned an isotope.
/// `isotope_properties` defines some physical properties of the spin.
#[derive(Debug,Clone,PartialEq,Default,Serialize,Deserialize)]
pub struct ParticlePropertiesTOML{
  //pub cosubstitute: Option<SecondaryParticleFilter>,
  //pub isotope_abundances: Option<Vec::<IsotopeAbundanceTOML>>,
  pub isotope_abundances: Option<HashMap::<String,f64>>,
  pub void_probability: Option<f64>,
}
//------------------------------------------------------------------------------
#[derive(Debug,Clone,PartialEq,Serialize,Deserialize)]
pub struct OrientationsTOML{
  pub grid: Option<String>,
  pub number: Option<usize>,
  pub file: Option<String>,
  pub vector: Option<Vec::<f64>>,
  pub vector_grid: Option<Vec::<Vec::<f64>>>,
}
//------------------------------------------------------------------------------
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct ConfigTOML{
  pub clash_distance: Option<f64>, 
  pub clash_distance_pbc: Option<f64>,
  pub cluster_batch_size: Option<usize>, 
  pub populations: Option<String>, 
  pub cluster_method: Option<String>,
  pub cluster_source: Option<String>,
  pub input_structure_file: Option<String>,
  pub magnetic_field: Option<f64>,
  pub max_cell_size: Option<usize>,
  pub max_cluster_size: Option<usize>,
  pub max_spins: Option<usize>,
  pub min_cell_size: Option<usize>,
  pub number_runs: Option<usize>, 
  pub number_timepoints: Option<Vec::<usize>>,
  pub replicate_unit_cell: Option<toml::Value>,
  pub partitioning: Option<String>, 
  pub partition_table: Option<toml::Value>,
  pub pdb_model_index: Option<usize>,
  pub pulse_sequence: Option<String>,  
  pub radius: Option<f64>,
  pub rng_seed: Option<u64>,
  pub output_directory: Option<String>,
  pub run_name: Option<String>,
  pub temperature: Option<f64>,  
  pub tau_increments: Option<Vec::<f64>>,
  pub unit_of_energy: Option<String>,
  pub unit_of_magnetic_field: Option<String>,
  pub unit_of_distance: Option<String>,
  pub unit_of_time: Option<String>,

  
  pub detected_spin: Option<DetectedSpinTOML>, 
  pub orientations: Option<OrientationsTOML>,
  pub output: Option<HashMap::<String,bool>>, 
  pub pair_cutoffs: Option<HashMap::<String,f64>>,
  pub groups: Option<Vec::<toml::Value>>, // TODO
}
impl ConfigTOML{
  fn set_default_units(&mut self){
    if self.unit_of_energy.is_none(){
      self.unit_of_energy = Some(DEFAULT_UNIT_ENERGY.to_string());
    }
    if self.unit_of_distance.is_none(){
      self.unit_of_distance = Some(DEFAULT_UNIT_DISTANCE.to_string());
    }
    if self.unit_of_magnetic_field.is_none(){
      self.unit_of_magnetic_field = Some(DEFAULT_UNIT_MAGNETIC_FIELD.to_string());
    }
    if self.unit_of_time.is_none(){
      self.unit_of_time = Some(DEFAULT_UNIT_TIME.to_string());
    }
  }
}

//------------------------------------------------------------------------------
fn check_toml_str(toml_str: &str) -> Result<(),CluEError>{
  let config: toml::Table = match toml::from_str(toml_str){
    Ok(cfg) => cfg,
    Err(err) => return Err(CluEError::CannotReadTOMLFile( format!("{}",err) )), 
  };

  check_toml_table(&config, 1)
}
//------------------------------------------------------------------------------
const MAX_DEPTH: usize = 10;

fn check_toml_table(config: &toml::Table, depth: usize) -> Result<(),CluEError>{  
  assert!(depth <= MAX_DEPTH);

  for (key, value) in config.iter(){
    if !ALLOWED_KEYS.contains( &&key[..] ){
      return Err(CluEError::CannotReadTOMLFile(format!("invalid key {}",key))); 
    }

    if let toml::Value::Table(table) = value{
      check_toml_table(table,depth +1)?;
    } 
  }

  Ok(())
}
//------------------------------------------------------------------------------


impl FromTOMLString for ConfigTOML{
  fn from_toml_string(toml_str: &str) -> Result<Self,CluEError>{
    check_toml_str(toml_str)?;

    let decoded: Result<ConfigTOML,_> = toml::from_str(toml_str);
    match decoded {
      Ok(mut config) => {
        config.set_default_units();
        Ok(config)
      },
      Err(err) => Err(CluEError::CannotReadTOMLFile( format!("{}",err) )), 
      //Err(err) => panic!("TODO: implement error: {}.",err)
    }
  }
}
impl fmt::Display for ConfigTOML{
   fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
     let toml = toml::to_string(self).unwrap();
     write!(f,"{}",toml)
   }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//==============================================================================
#[cfg(test)]
mod tests{
  use super::*;

  #[allow(non_snake_case)]
  #[test]
  fn test_ConfigTOML_from_string(){
    // Units:
    // Tims: μs,
    // Distance: Å,
    // Energy: Mhz,
    // Magnetic Field, T,
    let toml_str = r##"
        replicate_unit_cell = false
        clash_distance_pbc = 0.1
        cluster_batch_size = 20000

        populations = "thermal"
        temperature = 20

        cluster_method = "CCE"
        cluster_source = "clusters_file.toml"
        input_structure_file = "../../assets/TEMPO_wat_gly_70A.pdb"
        magnetic_field = 1.2
        max_cell_size = 2
        max_cluster_size = 4
        max_spins = 8
        min_cell_size = 1
        number_runs = 1
        partitioning = "exchange_groups"
        pdb_model_index = 0
        
        ##pulse_sequence = { CarrPurcell = 1 }
        pulse_sequence = "CP-1"

        radius = 80
        rng_seed = 0
        save_dir = "save_directory"

        number_timepoints = [40,60]
        tau_increments = [1, 500] # ns

        [orientations]
        grid = "lebedev" 
        number = 170

        #[orientations]
        #grid = "random" 
        #number = 170

        #[orientations]
        #grid = "file"
        #file = "xyzw.csv"

        #[orientations]
        #grid = "vector"
        #vector = [1,0,0]

        [detected_spin]
        spin_multiplicity = 2
        transition = [0,1]
        g_matrix.values = [2.0097, 2.0064, 2.0025]
        g_matrix.axes.x = [-1.1500, -0.4700, 0.7100]
        g_matrix.axes.y = { from = "tempo_c1", to = "tempo_c19"  }

        # single position (array of floats)
        #position = [0.0, 0.0, 0.0]
          
        # delocalized position (array of arrays of floats)
        position = [
          [0.0,0.0,0.0,0.5],
          [0.0,0.0,1.0,0.25],
          [0.0,0.0,-1.0,0.25],
        ]

        # delocalized position from file (string)
        #position = "xyzw.csv"

        # centroid over serials (array of ints)
        # position = [28, 29]

        # centroid over groups (array of string)
        # position = ["group1", "group2", "group3"]


        [pair_cutoffs]
        coupling = 1e+3
        delta_hyperfine_zz = 1e04
        point_dipole_perpendicular = 100
        distance = 10
        hahn_mod_depth = 1e-10
        hahn_taylor_4 = 1e-9

        [output]
        auxiliary_signals = true
        bath = true
        clusters = true
        exchange_groups = true
        info = true
        methyl_partitions = true
        orientation_signals = true
        sans_spin_signals = false
        structure_pdb = true
        tensors = false


        [[groups]]
        name = "hydrogens"
        
        drop_probability = 0.5
        1H.abundance = 0.9 
        2H.abundance = 0.1 

        1H.active = true

        2H.active = false

        [groups.selection] 
        elements = ["H"]
        not_elements = ["N",  "O"]
        within_distance = 25
        not_within_distance = 4


        [[groups]]

        name = "nitroxide_N"

        selection = {elements = ["N"],residues = ["R1M"]}

        14N.abundance = 0.9
        15N.abundance = 0.1


        14N.hyperfine_coupling.values = [14.7,14.7,101.4]
        14N.hyperfine_coupling.axes.x = { from = "self", to_bonded_to = "r1m_o" }
        14N.hyperfine_coupling.axes.y = { from_bonded_to = "r1m_c1", to_bonded_to = "r1m_c19" }

        14N.electric_quadrupole_coupling.values = [
          -0.6714,  0.4899, 1.4813, 
          0.4899,  1.1125, -0.1011, 
          1.4813, -0.1011, -0.4411
        ]
        14N.electric_quadrupole_coupling.axes.x = [1,0,0] 
        14N.electric_quadrupole_coupling.axes.y = [0,1,0] 


      "##;

    let config = ConfigTOML::from_toml_string(toml_str).unwrap();  

    
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_check_toml_str(){
    let toml_str = r##"
      this_is_not_a_key = True
      "##;
    assert!(check_toml_str(&toml_str).is_err());
  }
  //----------------------------------------------------------------------------

}
