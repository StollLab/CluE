
use crate::CluEError;
use crate::config::config_toml::*;
use crate::misc::{
  are_all_same_type,
  scaled_vec_cxmat_from_toml_array,
};
use crate::structure::particle_filter::{ParticleFilter,VectorSpecifier};
use crate::isotopes::Isotope;
use crate::physical_constants::{
  C3_TUNNEL_SPLITTING_TO_EXCHANGE_COUPLING,
  ONE,
};
use crate::space_3d::SymmetricTensor3D;

use std::collections::HashMap;
use strum::IntoEnumIterator;

use num_complex::Complex64;
use ndarray::Array2;
type CxMat = Array2::<Complex64>;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ParticleConfig` stores user setting for how to treat particles/
/// `label` is an identifier.
/// `filter` selects the set of particles to modify.
/// `properties` specify the custom properties.
/// `cell_type` (the group key `cells`) allows for treating the main cell and
/// PBC copies differently.
#[derive(Debug,Clone,PartialEq)]
pub struct ParticleConfig{
  pub label: String,
  pub filter: Option<ParticleFilter>,
  pub properties: Option<ParticleProperties>,
  pub cell_type: CellType,
}



pub fn set_particle_configs_from_toml_table(
    particle_configs: &mut Vec::<ParticleConfig>, 
    table: toml::Table,
    unit_of_distance: f64,
    unit_of_energy: f64)
  -> Result<(),CluEError>
{
  let label = if let Some(value) = table.get(KEY_NAME){
    let Some(name) = value.as_str() else{
      return Err(CluEError::InvalidTOMLValue(
          format!("groups.{}", KEY_NAME), "a string".to_string(),
          value.type_str().to_string()));
    };
    name
  } else{
    return Err(CluEError::MissingGroupName);
  };
 

  let mut found_particle_config: Option<usize> = None;
  for (idx, p_cfg) in particle_configs.iter().enumerate(){
    if *label == p_cfg.label{
      found_particle_config = Some(idx);
      break;
    }
  }

  if let Some(idx) = found_particle_config{
    particle_configs[idx].set_from_toml_table(table,
        unit_of_distance, unit_of_energy)?;
  }else{
    particle_configs.push(
        ParticleConfig::from_toml_table(table,
          unit_of_distance, unit_of_energy)?);
  }


  Ok(())
}

impl ParticleConfig{
  /// This function generates a default `ParticleConfig` with the input `label`.
  pub fn new(label: String) -> Self{
    ParticleConfig{
      label,
      filter: None,
      properties: None,
      cell_type: CellType::AllCells,
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_table(table: toml::Table, 
      unit_of_distance:f64, unit_of_energy: f64) 
      -> Result<Self,CluEError>
  {

    let label = if let Some(value) = table.get(KEY_NAME){
      let Some(name) = value.as_str() else{
        return Err(CluEError::InvalidTOMLValue(
          format!("groups.{}", KEY_NAME), "a string".to_string(),
          value.type_str().to_string()));
      };
      name
    } else {
      return Err(CluEError::MissingGroupName);
    };

    let mut particle_config = Self::new(label.to_string());  
    particle_config.set_from_toml_table(table,unit_of_distance, unit_of_energy)?;

    Ok(particle_config)
  }
  //----------------------------------------------------------------------------
  fn set_from_toml_table(&mut self, 
      table: toml::Table, 
      unit_of_distance: f64, unit_of_energy: f64) 
    -> Result<(),CluEError>
  {
    let label = if let Some(value) = table.get(KEY_NAME){
      let Some(name) = value.as_str() else{
        return Err(CluEError::InvalidTOMLValue(
          format!("groups.{}", KEY_NAME), "a string".to_string(),
          value.type_str().to_string()));
      };
      name.to_string()
    } else {
      return Err(CluEError::MissingGroupName);
    };

    if label != self.label{
      return Err(CluEError::MismatchedGroupNames(label,self.label.to_string()));
    }


    if let Some(value) = table.get(KEY_CELLS){
      let Some(cells) = value.as_str() else{
        return Err(CluEError::InvalidTOMLValue(
              format!("group \"{}\": {}", self.label, KEY_CELLS),
              "a string".to_string(), value.type_str().to_string()));
      };
      self.cell_type = CellType::from(cells)?;
    }

    // A group with neither include nor exclude has no filter and selects
    // no particles.
    let include = selection_table(&table, KEY_INCLUDE, &self.label)?;
    let exclude = selection_table(&table, KEY_EXCLUDE, &self.label)?;
    if include.is_some() || exclude.is_some(){

      // D and T are read as hydrogen.
      for (key, symbol) in hydrogen_isotope_symbols_in_selection(&table){
        println!("Group \"{}\": \"{}\" in {} is read as \
hydrogen, so it applies to all hydrogens.", self.label, symbol, key);
      }

      if let Some(filter) = &mut self.filter {
        filter.set_from_toml_tables(include, exclude, unit_of_distance)?;
      }else{
        self.filter = Some(ParticleFilter::from_toml_tables(
              include, exclude, unit_of_distance)?);
      }
    }
    if let Some(properties) = &mut self.properties {
      properties.set_from_toml_table(&table,unit_of_energy,&self.label)?;
    }else{
      self.properties = Some(ParticleProperties::from_toml_table(
            &table,unit_of_energy,&self.label)?);
    }

    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function looks at all spins that the `ParticleConfig` might apply to
  /// and returns the largest spin multiplicity, 2S+1.
  pub fn max_possible_spin_multiplicity(&self) -> Option<usize> {
    
    let properties = &self.properties.as_ref()?;

    if properties.isotopic_distribution.isotope_abundances.is_empty(){
      return None;
    }

    let mut max_spin_mult = 0;
    for isotope_abundance in properties.isotopic_distribution
      .isotope_abundances.iter()
    {
      max_spin_mult = usize::max(max_spin_mult,
          isotope_abundance.isotope.spin_multiplicity());
    }

    Some(max_spin_mult)
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This enum specifies the unit cells that a group is used for
/// (the group key `cells`): all cells, the primary cell, or the periodic 
/// copies of the primary cell (`Config::set_extracell_particles()`).
#[derive(Debug,Clone,PartialEq)]
pub enum CellType{
  AllCells,
  PrimaryCell,
  Extracells,
}
impl CellType{
  fn from(cell_type: &str) -> Result<Self,CluEError>{
    match cell_type{
      "all" => Ok(Self::AllCells),
      "primary" => Ok(Self::PrimaryCell),
      "extra" => Ok(Self::Extracells),
      _ => Err(CluEError::CannotParseCellType(cell_type.to_string())),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//------------------------------------------------------------------------------
// This function returns the table `key` (include or exclude) of the group 
// `table` with name `label`, if the group has one.
fn selection_table<'a>(table: &'a toml::Table, key: &str, label: &str)
  -> Result<Option<&'a toml::Table>,CluEError>
{
  let Some(value) = table.get(key) else{
    return Ok(None);
  };
  let Some(selection) = value.as_table() else{
    return Err(CluEError::InvalidTOMLValue(
          format!("group \"{}\": {}", label, key),
          "a table".to_string(), value.type_str().to_string()));
  };
  Ok(Some(selection))
}
//------------------------------------------------------------------------------
// This function returns the element selections of the group `table` 
// (include.elements, exclude.elements, and their bonded_to forms) that 
// contain the symbols "D" or "T", which are read as hydrogen, 
// as (key, symbol) pairs.
fn hydrogen_isotope_symbols_in_selection(table: &toml::Table)
  -> Vec::<(String, String)>
{
  let mut found = Vec::<(String, String)>::new();
  for key in [KEY_INCLUDE, KEY_EXCLUDE]{
    let Some(selection) = table.get(key).and_then(|v| v.as_table()) else{
      continue;
    };
    let bonded_to = selection.get(KEY_BONDED_TO).and_then(|v| v.as_table());
    for (path, sub) in [(key.to_string(), Some(selection)),
        (format!("{}.{}", key, KEY_BONDED_TO), bonded_to)]
    {
      let Some(value) = sub.and_then(|t| t.get(KEY_SELE_ELEMENTS)) else{
        continue;
      };
      // A single value is a list of one.
      let symbols = match value{
        toml::Value::Array(array) => array.iter().collect::<Vec::<_>>(),
        single => vec![single],
      };
      for symbol in symbols.iter().filter_map(|v| v.as_str()){
        if symbol == "D" || symbol == "T"{
          found.push((format!("{}.{}", path, KEY_SELE_ELEMENTS), 
                symbol.to_string()));
        }
      }
    }
  }
  found
}
//------------------------------------------------------------------------------

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `Cosubstitute` selects the particles that are cosubstituted with a 
/// particle, that is, that get the same random isotope and are kept or 
/// dropped together:
/// `Group` selects all particles of its group,
/// `Particle` selects only the particle itself,
/// and `SameMolecule` selects the particles of its group in the same molecule.
#[derive(Debug,Clone,PartialEq)]
pub enum Cosubstitute{
  Group,
  Particle,
  SameMolecule,
}
impl Cosubstitute{
  /// This function converts the value of "cosubstitute" in group `group`
  /// to a `Cosubstitute`.
  pub fn from(cosubstitute: &str, group: &str) -> Result<Self,CluEError>{
    match cosubstitute{
      "group" | "filter" => Ok(Self::Group),
      "particle" => Ok(Self::Particle),
      "same_molecule" => Ok(Self::SameMolecule),
      _ => Err(CluEError::CannotParseCosubstitute(group.to_string(),
            cosubstitute.to_string())),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ParticleProperties` specifies custom particle properties.
/// 'cosubstitute` selects the set of particles that should always be the same
/// isotope when the isotopic distribution is randomized.
/// `isotopic_distribution` specifies how elements are assigned an isotope.
/// `isotope_properties` defines some physical properties of the spin.
#[derive(Debug,Clone,PartialEq,Default)]
pub struct ParticleProperties{
  pub cosubstitute: Option<Cosubstitute>,
  // TODO: move the fields of isotopic_distribution to ParticleProperties,
  // to avoid the subtable.
  pub isotopic_distribution:  IsotopeDistribution,
  pub isotope_properties: HashMap::<String,IsotopeProperties>,
}

impl ParticleProperties{
  /// This function creates a new instance of `ParticleProperties`.
  pub fn new() -> Self{
    Default::default()
  }
  //----------------------------------------------------------------------------
  /// This function reads the properties of the group `group` 
  /// from a TOML table.
  pub fn from_toml_table(table: &toml::Table, unit_of_energy: f64,
      group: &str) 
      -> Result<Self,CluEError>
  {
    let mut properties = Self::new();  
    properties.set_from_toml_table(table,unit_of_energy,group)?;

    Ok(properties)
  }
  //----------------------------------------------------------------------------
  fn set_from_toml_table(&mut self, table: &toml::Table, unit_of_energy: f64,
      group: &str) 
    -> Result<(),CluEError>
  {
    if let Some(value) = table.get(KEY_ISO_COSUBSTITUTE){
      let Some(cosubstitute) = value.as_str() else{
        return Err(CluEError::InvalidTOMLValue(
              format!("group \"{}\": {}", group, KEY_ISO_COSUBSTITUTE),
              "a string".to_string(), value.type_str().to_string()));
      };
      self.cosubstitute = Some(Cosubstitute::from(cosubstitute,group)?);
    }

    if let Some(value) = table.get(KEY_DROP_PROB){
      let Some(p) = value.as_float() else{
        return Err(CluEError::InvalidTOMLValue(
              format!("group \"{}\": {}", group, KEY_DROP_PROB),
              "a float".to_string(), value.type_str().to_string()));
      };
      self.isotopic_distribution.void_probability = Some(p);
    }

    if let Some(value) = table.get(KEY_ABUNDANCE_TOL){
      let Some(tolerance) = value.as_float() else{
        return Err(CluEError::InvalidTOMLValue(
              format!("group \"{}\": {}", group, KEY_ABUNDANCE_TOL),
              "a float".to_string(), value.type_str().to_string()));
      };
      self.isotopic_distribution.abundance_tolerance = Some(tolerance);
    }

    for isotope in Isotope::iter(){
      let iso_str = isotope.to_string();
      let Some(iso_value) = table.get(&iso_str) else {
        continue;
      };
      let Some(iso_table) = iso_value.as_table() else {
        return Err(CluEError::InvalidTOMLValue(
              format!("group \"{}\": {}", group, iso_str),
              "a table".to_string(), iso_value.type_str().to_string()));
      };

      if let Some(value) = iso_table.get(KEY_ISO_ABUNDACE){
        let Some(abundance) = value.as_float() else{
          return Err(CluEError::InvalidTOMLValue(
                format!("group \"{}\": {}.{}", group, iso_str, 
                  KEY_ISO_ABUNDACE),
                "a float".to_string(), value.type_str().to_string()));
        };
        self.isotopic_distribution.isotope_abundances.push(
            IsotopeAbundance{isotope,abundance} );
      }

      let context = format!("group \"{}\": {}", group, iso_str);
      if let Some(iso_prop) = self.isotope_properties.get_mut(&iso_str){
        iso_prop.set_from_toml_table(iso_table,unit_of_energy,&context)?;
      } else{
        let iso_prop = IsotopeProperties::from_toml_table(iso_table,
            unit_of_energy,&context)?;
        self.isotope_properties.insert(iso_str.clone(), iso_prop);
      }

    }


    Ok(())
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `IsotopeProperties` defines some physical properties of a spin.
/// `active` specifies whether the spin should be included in simulations.
/// `electric_quadrupole_coupling` defines the coupling of the spin to the 
/// local electric field gradient.
/// `exchange_coupling` defines an effective coupling coming from
/// symmetry requirements of the wavefunction upon the exchange of 
/// identical particles.
/// `g_matrix` specifies the effective coupling between the spin and the applied
/// magnetic field.
/// `hyperfine_coupling` specifices the coupling to the detected electron.
#[derive(Debug,Clone,PartialEq,Default)]
pub struct IsotopeProperties{
  pub active: Option<bool>,
  pub electric_quadrupole_coupling: Option<TensorSpecifier>,
  pub zerofield_coupling: Option<TensorSpecifier>,
  pub exchange_coupling: Option<f64>,
  pub g_matrix: Option<TensorSpecifier>,
  pub hyperfine_coupling: Option<TensorSpecifier>,
}

impl IsotopeProperties{
  /// This function creates a new instance of `IsotopeProperties`.
  pub fn new() -> Self{ IsotopeProperties::default() }
  //----------------------------------------------------------------------------
  /// `context` names the group and isotope for error messages.
  pub fn from_toml_table(table: &toml::Table, unit_of_energy: f64,
      context: &str) 
      -> Result<Self,CluEError>
  {
    let mut properties = Self::new();  
    properties.set_from_toml_table(table,unit_of_energy,context)?;

    Ok(properties)
  }
  //----------------------------------------------------------------------------
  fn set_from_toml_table(&mut self, table: &toml::Table, unit_of_energy: f64,
      context: &str) 
    -> Result<(),CluEError>
  {
    if let Some(value) = table.get(KEY_ISO_ACTIVE){
      let Some(active) = value.as_bool() else{
        return Err(CluEError::InvalidTOMLValue(
              format!("{}.{}", context, KEY_ISO_ACTIVE),
              "a bool".to_string(), value.type_str().to_string()));
      };
      self.active = Some(active);
    }

    if let Some(value) = table.get(KEY_ISO_G_MATRIX){
      let tensor = TensorSpecifier::from_toml_value(value.clone(),1.0,
          &format!("{}.{}", context, KEY_ISO_G_MATRIX))?;
      self.g_matrix = Some(tensor);
    }

    if let Some(value) = table.get(KEY_ISO_ELEC_QUADRUPOLE){
      let tensor = TensorSpecifier::from_toml_value(
          value.clone(), unit_of_energy,
          &format!("{}.{}", context, KEY_ISO_ELEC_QUADRUPOLE))?;
      self.electric_quadrupole_coupling = Some(tensor);
    }

    if let Some(value) = table.get(KEY_ISO_ZEROFIELD){
      let tensor = TensorSpecifier::from_toml_value(
          value.clone(), unit_of_energy,
          &format!("{}.{}", context, KEY_ISO_ZEROFIELD))?;
      self.zerofield_coupling = Some(tensor);
    }

    if let Some(value) = table.get(KEY_ISO_HYPERFINE){
      let tensor = TensorSpecifier::from_toml_value(value.clone(),
          unit_of_energy, &format!("{}.{}", context, KEY_ISO_HYPERFINE))?;
      self.hyperfine_coupling = Some(tensor);
    }

    let n_ex: usize = [table.contains_key(KEY_ISO_EXCHANGE_COUPLING),
      table.contains_key(KEY_ISO_C3_TUNNEL_SPLITTING)].iter().map(|b| if *b{1}else{0} ).sum();

    if n_ex > 1{
      return Err(CluEError::MutuallyExclusiveKeys(
            KEY_ISO_EXCHANGE_COUPLING.to_string(),
            KEY_ISO_C3_TUNNEL_SPLITTING.to_string()));
    }

    if let Some(value) = table.get(KEY_ISO_EXCHANGE_COUPLING)
    {
      let Some(ex) = value.as_float() else{
        return Err(CluEError::InvalidTOMLValue(
              format!("{}.{}", context, KEY_ISO_EXCHANGE_COUPLING),
              "a float".to_string(), value.type_str().to_string()));
      };
      self.exchange_coupling = Some(ex*unit_of_energy);

    } else if let Some(value) = table.get(KEY_ISO_C3_TUNNEL_SPLITTING)
    {
      let Some(nut) = value.as_float() else{
        return Err(CluEError::InvalidTOMLValue(
              format!("{}.{}", context, KEY_ISO_C3_TUNNEL_SPLITTING),
              "a float".to_string(), value.type_str().to_string()));
      };
      self.exchange_coupling = Some(
        nut*C3_TUNNEL_SPLITTING_TO_EXCHANGE_COUPLING*unit_of_energy
      );
    }


    Ok(())
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `IsotopeDistribution` determine how isotopes are assigned to elements.
/// `isotope_abundances` lists the possible isotope and propabilities;
/// if any are given, they must sum to one within `abundance_tolerance`.
/// `abundance_tolerance` is the allowed difference between one and the sum
/// of the abundances (default `DEFAULT_ABUNDANCE_TOLERANCE`).
/// 'void_probability' is the probability that the particle will be removed.
#[derive(Debug,Clone,PartialEq,Default)]
pub struct IsotopeDistribution{
  // TODO: can isotope_abundances be a HashMap::<Isotope,f64> instead?
  pub isotope_abundances: Vec::<IsotopeAbundance>,
  pub abundance_tolerance: Option<f64>,
  pub void_probability: Option<f64>,
}
/// The default allowed difference between one and the sum of the isotope 
/// abundances of a group.
pub const DEFAULT_ABUNDANCE_TOLERANCE: f64 = 1e-9;
//------------------------------------------------------------------------------
/// `IsotopeAbundances` lists the possible isotope and abundances.
#[derive(Debug,Clone,PartialEq)]
pub struct IsotopeAbundance{
  pub isotope: Isotope,
  pub abundance: f64,         
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `TensorSpecifier` specifies a symmetric 3-by-3 coupling matrix.
#[derive(Debug,Clone,PartialEq,Default)]
pub enum TensorSpecifier{
  Eig(EigSpecifier),
  SymmetricTensor3D(SymmetricTensor3D),
  #[default]
  Unspecified,
}

impl TensorSpecifier{
  /// This function creates a new instance of `TensorSpecifier`.
  pub fn new() -> Self{
    Self::default()
  }
  /// This function reads a tensor specifier.
  /// `key` names the TOML key for error messages.
  pub fn from_toml_value(value: toml::Value, units: f64, key: &str) 
      -> Result<Self,CluEError>
  {
    match value{
      toml::Value::Array(array) => {
        let mut tensor = SymmetricTensor3D::from_toml_array(array)?;
        tensor.scale_mut(units);
        Ok(Self::SymmetricTensor3D(tensor))
      },
      toml::Value::Table(table) => Ok(Self::Eig(
            EigSpecifier::from_toml_table(table,units,key)?)),
      other => Err(CluEError::InvalidTOMLValue(key.to_string(),
            "an array or a table with values and axes".to_string(),
            other.type_str().to_string())),
    }
  
  }
  /*
  pub fn from_toml_value(value: toml::Value) -> Result<Self,CluEError> 
  {
    match value{
      Value::Array(array) => Self::from_toml_array(array),
      _ => Err(CluEError::TOMLArrayDoesNotSpecifyATensor),
    }
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_array(array: Vec::<toml::Value>) -> Result<Self,CluEError>
  {

    if array.is_empty(){
      return Err(CluEError::TOMLArrayIsEmpty);
    }
    if !are_all_same_type(&array){
      return Err(CluEError::TOMLArrayContainsMultipleTypes);
    }
    match array[0]{
      Value::Float(_) => {
        let vec: Vec::<f64> = array.iter()
          .filter_map(|v| v.as_float()).collect();

        match vec.len(){
          1 => (),
          2 => (),
          3 => (),
          6 => (),
          9 => ()
        }

        Ok(Self::XYZ(vec3d))

      },
      Value::Integer(_) => {
        let vec: Vec::<u32> = array.iter()
          .filter_map(|v| v.as_integer()).map(|n| n as u32).collect();

        Ok(Self::CentroidOverSerials(vec))
      },
      Value::String(_) => {
        let vec: Vec::<String> = array.iter()
          .filter_map(|v| v.as_str()).map(|s| s.to_string()).collect();

        Ok(Self::CentroidOverGroup(vec))
      },
      _ => Err(CluEError::TOMLArrayDoesNotSpecifyAVector),
    }
  }
  */
}

/// `EigSpecifier` specifies a symmetric 3-by-3 coupling matrix.
/// `values` contains the three eigenvalues.
/// `e_axis` for `e` in {'x','y','z'} specify eigenvectors.
/// Note that only two axes should be specified. 
#[derive(Debug,Clone,PartialEq,Default)]
pub struct EigSpecifier{
  pub values: Option<[f64; 3]>,
  pub x_axis: Option<VectorSpecifier>,
  pub y_axis: Option<VectorSpecifier>,
  pub z_axis: Option<VectorSpecifier>,
  pub spherical_tensor_coefficients: Vec::<CxMat>,
  pub stevens: bool,
}


impl EigSpecifier{
  /// This function creates a new instance of `EigSpecifier`.
  pub fn new() -> Self{
    Self::default()
  }
  //----------------------------------------------------------------------------
  /// `key` names the TOML key for error messages.
  pub fn from_toml_table(table: toml::Table,units: f64, key: &str) 
      -> Result<Self,CluEError>
  {
    let values_key = format!("{}.{}", key, KEY_EIG_VALUES);
    let values = match table.get(KEY_EIG_VALUES){
      Some(toml::Value::Array(array)) => {
        if array.len() != 3 || !are_all_same_type(array) 
          || !array[0].is_float(){
          return Err(CluEError::InvalidTOMLValue(values_key,
                "an array of three floats".to_string(),
                toml::Value::Array(array.clone()).to_string()));
        }
        array.iter().filter_map(|v| v.as_float()).collect::<Vec::<f64>>()
      },
      Some(value) => return Err(CluEError::InvalidTOMLValue(values_key,
            "an array of three floats".to_string(),
            value.type_str().to_string())),
      None => return Err(CluEError::MissingTOMLKey(values_key,
            format!("{} given as a table", key))),
    }; 
    let values = Some([values[0]*units, values[1]*units, values[2]*units]);

    let Some(axes) = table.get(KEY_EIG_AXES) else{
      return Err(CluEError::MissingTOMLKey(
            format!("{}.{}", key, KEY_EIG_AXES),
            format!("{} given as a table", key)));
    };

    let x_axis = if let Some(ax) = axes.get(KEY_EIG_X_AXIS){
      Some(VectorSpecifier::from_toml_value(ax.clone(),
            &format!("{}.{}.{}", key, KEY_EIG_AXES, KEY_EIG_X_AXIS))?)
    }else{
      None
    };
    let y_axis = if let Some(ax) = axes.get(KEY_EIG_Y_AXIS){
      Some(VectorSpecifier::from_toml_value(ax.clone(),
            &format!("{}.{}.{}", key, KEY_EIG_AXES, KEY_EIG_Y_AXIS))?)
    }else{
      None
    };
    let z_axis = if let Some(ax) = axes.get(KEY_EIG_Z_AXIS){
      Some(VectorSpecifier::from_toml_value(ax.clone(),
            &format!("{}.{}.{}", key, KEY_EIG_AXES, KEY_EIG_Z_AXIS))?)
    }else{
      None
    };

    let mut stevens = false;

    let spherical_tensor_coefficients = 
          match (table.get(KEY_IST_COEF), table.get(KEY_STEVENS_COEF)){
      (Some(toml::Value::Array(a)),None) => 
        scaled_vec_cxmat_from_toml_array(
          a.clone(), units*ONE)?,
      (None,Some(toml::Value::Array(a))) => {
        stevens = true;
        scaled_vec_cxmat_from_toml_array(a.clone(),units*ONE)?
      },
      (Some(v),None) => return Err(CluEError::InvalidIST(v.to_string())),   
      (None,Some(v)) => return Err(CluEError::InvalidStevens(v.to_string())),   
      (None,None) => Vec::<CxMat>::new(),
      (_,_) => return Err(CluEError::BothStevensAndIST),
    };

    Ok(Self{
        values,
        x_axis,
        y_axis,
        z_axis,
        spherical_tensor_coefficients,
        stevens,
    })
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  use crate::elements::Element;
  //----------------------------------------------------------------------------
  #[test]
  fn test_hydrogen_isotope_symbols_in_selection(){
    let group: toml::Table = toml::from_str(r##"
        include.elements = ["D", "C"]
        include.bonded_to.elements = "T"
        exclude.elements = ["O"]
        exclude.bonded_to.elements = ["N", "D"]
        "##).unwrap();
    let s = |x: &str| x.to_string();
    assert_eq!(hydrogen_isotope_symbols_in_selection(&group), vec![
        (s("include.elements"), s("D")),
        (s("include.bonded_to.elements"), s("T")),
        (s("exclude.bonded_to.elements"), s("D")),
    ]);

    // The selection still treats D and T as hydrogen.
    let include = group.get(KEY_INCLUDE).and_then(|v| v.as_table());
    let exclude = group.get(KEY_EXCLUDE).and_then(|v| v.as_table());
    let filter = ParticleFilter::from_toml_tables(include, exclude, 1.0)
      .unwrap();
    assert_eq!(filter.elements, vec![Element::Hydrogen, Element::Carbon]);
    assert_eq!(filter.bonded_elements, vec![Element::Hydrogen]);
    assert_eq!(filter.not_bonded_elements, 
        vec![Element::Nitrogen, Element::Hydrogen]);
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
