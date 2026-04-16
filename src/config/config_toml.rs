use crate::clue_errors::CluEError;
use crate::io::FromTOMLString;

use std::fmt;
use std::collections::HashMap;

use serde::{Serialize,Deserialize};
pub use crate::config::toml_keys::*;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone,Default,Serialize,Deserialize)]
pub struct DetectedSpinTOML{
  pub multiplicity: Option<usize>,
  pub g_matrix: Option<toml::Value>,
  pub electric_quadrupole: Option<toml::Value>,
  pub zerofield: Option<toml::Value>,
  pub density_matrix: Option<toml::Value>,
  pub detection_operator: Option<toml::Value>,
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
  pub connect_exchange_groups: Option<bool>,
  pub populations: Option<String>, 
  pub cluster_method: Option<String>,
  pub cluster_source: Option<String>,
  pub input_structure_file: Option<String>,
  pub kmeans_size: Option::<usize>,
  pub magnetic_field: Option<f64>,
  pub max_cell_size: Option<usize>,
  pub max_cluster_size: Option<usize>,
  pub max_spins: Option<usize>,
  pub min_cell_size: Option<usize>,
  pub number_runs: Option<usize>, 
  pub number_timepoints: Option<Vec::<usize>>,
  pub replicate_unit_cell: Option<toml::Value>,
  pub run_in_parallel: Option<bool>,
  pub partitioning: Option<String>, 
  pub partition_table: Option<toml::Value>,
  pub pdb_model_index: Option<usize>,
  pub pulse_sequence: Option<toml::Value>,  
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
    if let toml::Value::Array(array) = value{
      check_toml_array(array,depth +1)?;
    }
  }

  Ok(())
}
//------------------------------------------------------------------------------
fn check_toml_array(config: &[toml::Value], depth: usize) -> Result<(),CluEError>
{  
  assert!(depth <= MAX_DEPTH);

  for el in config.iter(){
    match el{
      toml::Value::Table(table) => check_toml_table(table,depth +1)?,
      toml::Value::Array(array) => check_toml_array(array,depth +1)?,
      _ => (),  
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
        connect_exchange_groups = true
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
        output_directory = "save_directory"

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

    let _config = ConfigTOML::from_toml_string(toml_str).unwrap();  

    
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
