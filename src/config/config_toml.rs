use crate::clue_errors::CluEError;
use crate::io::FromTOMLString;
use crate::config::validate::check_toml_keys;

use std::fmt;
use std::collections::HashMap;

use serde::{Serialize,Deserialize};
pub use crate::config::toml_keys::*;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct DetectedSpinTOML{
  pub multiplicity: Option<usize>,
  pub g_matrix: Option<toml::Value>,
  pub zerofield: Option<toml::Value>,
  pub density_matrix: Option<toml::Value>,
  pub detection_operator: Option<toml::Value>,
  pub frame: Option<String>,
  //pub g_values: Option<Vec::<f64>>,
  //pub gx: Option<toml::Value>,
  //pub gy: Option<toml::Value>,
  //pub gz: Option<toml::Value>,
  //pub position: Option<VectorSpecifierTOML>,
  pub position: Option<toml::Value>,
  pub transition: Option<[usize;2]>,
  pub weight_tolerance: Option<f64>,
  pub identity: Option<String>,
  pub max_spherical_tensor_rank: Option<usize>,
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
  pub weight_tolerance: Option<f64>,
}
//------------------------------------------------------------------------------
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `UnitsTOML` holds the [units] table: the units of the numbers in the input.
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct UnitsTOML{
  pub distance: Option<String>,
  pub energy: Option<String>,
  pub magnetic_field: Option<String>,
  pub time: Option<String>,
}
//------------------------------------------------------------------------------
/// `StructureTOML` holds the [structure] table: the structure file and how
/// the system is built from it.
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct StructureTOML{
  pub file: Option<String>,
  pub pdb_model_index: Option<usize>,
  pub radius: Option<f64>,
  pub replicate_unit_cell: Option<toml::Value>,
  pub clash_distance: Option<f64>, 
  pub clash_distance_pbc: Option<f64>,
}
//------------------------------------------------------------------------------
/// `ConditionsTOML` holds the [conditions] table: the sample conditions.
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct ConditionsTOML{
  pub magnetic_field: Option<f64>,
  pub temperature: Option<f64>,  
}
//------------------------------------------------------------------------------
/// `PulseSequenceTOML` holds the [pulse_sequence] table: the sequence and 
/// its time axes.
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct PulseSequenceTOML{
  pub sequence: Option<toml::Value>,  
  pub number_timepoints: Option<Vec::<usize>>,
  pub number_timepoints2: Option<Vec::<usize>>,
  pub tau_increments: Option<Vec::<f64>>,
  pub tau2_increments: Option<Vec::<f64>>,
}
//------------------------------------------------------------------------------
/// `ClusterExpansionTOML` holds the [cluster_expansion] table: the method,
/// the truncation order, and the starting state of the cluster calculations.
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct ClusterExpansionTOML{
  pub method: Option<String>,
  pub max_cluster_size: Option<usize>,
  pub ensemble: Option<bool>,
  pub populations: Option<String>, 
  pub mean_fields: Option<bool>,
  pub thermalize_cluster_populations: Option<bool>,
  pub clustering: Option<ClusteringTOML>,
}
//------------------------------------------------------------------------------
/// `ClusteringTOML` holds the [cluster_expansion.clustering] table: 
/// how the clusters are read or built.
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct ClusteringTOML{
  pub file: Option<String>,
  pub max_spins: Option<usize>,
  pub partitioning: Option<String>, 
  pub kmeans_size: Option::<usize>,
  pub kmeans_iter_cap: Option::<usize>,
  pub partition_table: Option<toml::Value>,
  pub connect_exchange_groups: Option<bool>,
  pub pair_cutoffs: Option<HashMap::<String,f64>>,
}
//------------------------------------------------------------------------------
/// `RunTOML` holds the [run] table: how the calculation is run.
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct RunTOML{
  pub rng_seed: Option<u64>,
  pub number_runs: Option<usize>, 
  pub parallel: Option<bool>,
  pub cluster_batch_size: Option<usize>, 
}
//------------------------------------------------------------------------------
/// `OutputTOML` holds the [output] table: where the output is written and
/// which files are written.
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct OutputTOML{
  pub directory: Option<String>,
  pub auxiliary_signals: Option<bool>,
  pub bath: Option<bool>,
  pub clusters: Option<bool>,
  pub detected_spin: Option<bool>,
  pub exchange_groups: Option<bool>,
  pub info: Option<bool>,
  pub methyl_partitions: Option<bool>,
  pub orientation_grid: Option<bool>,
  pub orientation_signals: Option<bool>,
  pub partition_table: Option<bool>,
  pub sans_spin_signals: Option<bool>,
  pub structure_pdb: Option<bool>,
  pub tensors: Option<bool>,
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ConfigTOML` holds the parsed input file. Its top level holds only tables.
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct ConfigTOML{
  pub units: Option<UnitsTOML>,
  pub structure: Option<StructureTOML>,
  pub detected_spin: Option<DetectedSpinTOML>, 
  pub groups: Option<Vec::<toml::Value>>, 
  pub conditions: Option<ConditionsTOML>,
  pub orientations: Option<OrientationsTOML>,
  pub pulse_sequence: Option<PulseSequenceTOML>,
  pub cluster_expansion: Option<ClusterExpansionTOML>,
  pub run: Option<RunTOML>,
  pub output: Option<OutputTOML>, 
}
impl ConfigTOML{
  fn set_default_units(&mut self){
    let units = self.units.get_or_insert_with(UnitsTOML::default);
    if units.energy.is_none(){
      units.energy = Some(DEFAULT_UNIT_ENERGY.to_string());
    }
    if units.distance.is_none(){
      units.distance = Some(DEFAULT_UNIT_DISTANCE.to_string());
    }
    if units.magnetic_field.is_none(){
      units.magnetic_field = Some(DEFAULT_UNIT_MAGNETIC_FIELD.to_string());
    }
    if units.time.is_none(){
      units.time = Some(DEFAULT_UNIT_TIME.to_string());
    }
  }
}

//------------------------------------------------------------------------------
// This function parses `toml_str` and checks that every key is allowed in
// the table that contains it (`check_toml_keys()`).
fn check_toml_str(toml_str: &str) -> Result<(),CluEError>{
  let config: toml::Table = match toml::from_str(toml_str){
    Ok(cfg) => cfg,
    Err(err) => return Err(CluEError::CannotReadTOMLFile( format!("{}",err) )), 
  };

  check_toml_keys(&config)
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
        [structure]
        file = "../../assets/TEMPO_wat_gly_70A.pdb"
        pdb_model_index = 0
        radius = 80
        replicate_unit_cell = false
        clash_distance_pbc = 0.1

        [conditions]
        magnetic_field = 1.2
        temperature = 20

        [pulse_sequence]
        #sequence = "CP-1"
        sequence = [
          ["pulse","pi/2"],
          ["delay", "tau"],
          ["pulse","pi"],
          ["delay", "tau"],
          ["detect"],
        ]
        number_timepoints = [40,60]
        tau_increments = [1, 500] # ns

        [cluster_expansion]
        method = "CCE"
        max_cluster_size = 4
        populations = "thermal"

        [cluster_expansion.clustering]
        file = "clusters_file.toml"
        connect_exchange_groups = true
        max_spins = 8
        partitioning = "exchange_groups"

        [cluster_expansion.clustering.pair_cutoffs]
        delta_hyperfine_zz = 1e04
        point_dipole_perpendicular = 100
        distance = 10
        hahn_mod_depth = 1e-10
        hahn_taylor_4 = 1e-9

        [run]
        cluster_batch_size = 20000
        number_runs = 1
        rng_seed = 0

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
        multiplicity = 2
        detection_operator = [[0,1],[0,0]]  
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


        [output]
        directory = "save_directory"
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

        [groups.include] 
        elements = ["H"]
        distance = [4, 25]

        [groups.exclude] 
        elements = ["N",  "O"]


        [[groups]]

        name = "nitroxide_N"

        include = {elements = ["N"],residues = ["R1M"]}

        14N.abundance = 0.9
        15N.abundance = 0.1


        14N.hyperfine.values = [14.7,14.7,101.4]
        14N.hyperfine.axes.x = { from = "self", to_bonded_to = "r1m_o" }
        14N.hyperfine.axes.y = { from_bonded_to = "r1m_c1", to_bonded_to = "r1m_c19" }

        14N.electric_quadrupole.values = [
          -0.6714,  0.4899, 1.4813, 
          0.4899,  1.1125, -0.1011, 
          1.4813, -0.1011, -0.4411
        ]
        14N.electric_quadrupole.axes.x = [1,0,0] 
        14N.electric_quadrupole.axes.y = [0,1,0] 


      "##;

    let config = ConfigTOML::from_toml_string(toml_str).unwrap();  
    let pulse_sequence = config.pulse_sequence.unwrap();
    assert!(pulse_sequence.sequence.is_some());
    assert_eq!(pulse_sequence.number_timepoints, Some(vec![40,60]));
    let clustering = config.cluster_expansion.unwrap().clustering.unwrap();
    assert_eq!(clustering.file, Some("clusters_file.toml".to_string()));
    assert_eq!(clustering.pair_cutoffs.unwrap().get("distance"), Some(&10.0));
    assert_eq!(config.output.unwrap().directory, 
        Some("save_directory".to_string()));
    // Missing units get their defaults.
    let units = config.units.unwrap();
    assert_eq!(units.distance, Some(DEFAULT_UNIT_DISTANCE.to_string()));
    assert_eq!(units.time, Some(DEFAULT_UNIT_TIME.to_string()));

    
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
