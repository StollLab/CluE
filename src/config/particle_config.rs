
use crate::CluEError;
use crate::config::config_toml::*;
use crate::misc::are_all_same_type;
use crate::structure::particle_filter::{ParticleFilter,VectorSpecifier,
  SecondaryParticleFilter};
use crate::isotopes::Isotope;
use crate::physical_constants::C3_TUNNEL_SPLITTING_TO_EXCHANGE_COUPLING;
use crate::space_3d::SymmetricTensor3D;

use std::collections::HashMap;
use strum::IntoEnumIterator;


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ParticleConfig` stores user setting for how to treat particles/
/// `label` is an identifier.
/// `filter` selects the set of particles to modify.
/// `properties` specify the custom properties.
/// `cell_type` allows for treating the main cell and PBC copies differently.
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
      return Err(CluEError::ExpectedTOMLString(value.type_str().to_string()));
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
        return Err(CluEError::ExpectedTOMLString(value.type_str().to_string()));
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
        return Err(CluEError::ExpectedTOMLString(value.type_str().to_string()));
      };
      name.to_string()
    } else {
      return Err(CluEError::MissingGroupName);
    };

    if label != self.label{
      return Err(CluEError::MismatchedGroupNames(label,self.label.to_string()));
    }


    if let Some(value) = table.get(KEY_SELECTION){
      let Some(selection) = value.as_table() else{
        return Err(CluEError::ExpectedTOMLTable(value.type_str().to_string()));
      };

      if let Some(value) = selection.get(KEY_CELL_TYPE){
        let Some(cell_type_str) = value.as_str() else{
          return Err(CluEError::ExpectedTOMLString(
                value.type_str().to_string()));
        };
        self.cell_type = CellType::from(cell_type_str)?;
      }

      if let Some(filter) = &mut self.filter {
        filter.set_from_toml_table(selection.clone(),unit_of_distance)?;
      }else{
        self.filter = Some(ParticleFilter::from_toml_table(
              selection.clone(), unit_of_distance)?);
      }
    }
    if let Some(properties) = &mut self.properties {
      properties.set_from_toml_table(&table,unit_of_energy)?;
    }else{
      self.properties = Some(ParticleProperties::from_toml_table(
            &table,unit_of_energy)?);
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
/// This enum specifies different sets of PBC cells.
#[derive(Debug,Clone,PartialEq)]
pub enum CellType{
  AllCells,
  PrimaryCell,
  Extracells,
}
impl CellType{
  fn from(cell_type: &str) -> Result<Self,CluEError>{
    match cell_type{
      "all_cells" => Ok(Self::AllCells),
      "primary_cell" => Ok(Self::PrimaryCell),
      "extracells" => Ok(Self::Extracells),
      _ => Err(CluEError::CannotParseCellType(cell_type.to_string())),
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
  pub cosubstitute: Option<SecondaryParticleFilter>,
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
  pub fn from_toml_table(table: &toml::Table, unit_of_energy: f64) 
      -> Result<Self,CluEError>
  {
    let mut properties = Self::new();  
    properties.set_from_toml_table(table,unit_of_energy)?;

    Ok(properties)
  }
  //----------------------------------------------------------------------------
  fn set_from_toml_table(&mut self, table: &toml::Table, unit_of_energy: f64) 
    -> Result<(),CluEError>
  {
    if let Some(value) = table.get(KEY_ISO_COSUBSTITUTE){
      let Some(cosubstitute) = value.as_str() else{
        return Err(CluEError::ExpectedTOMLString(value.type_str().to_string() ));
      };
      self.cosubstitute = Some(SecondaryParticleFilter::from(cosubstitute)?);
    }

    if let Some(value) = table.get(KEY_DROP_PROB){
      let Some(p) = value.as_float() else{
        return Err(CluEError::ExpectedTOMLFloat(value.type_str().to_string() ));
      };
      self.isotopic_distribution.void_probability = Some(p);
    }

    for isotope in Isotope::iter(){
      let iso_str = isotope.to_string();
      let Some(iso_value) = table.get(&iso_str) else {
        continue;
      };
      let Some(iso_table) = iso_value.as_table() else {
        return Err(CluEError::ExpectedTOMLTable(
              iso_value.type_str().to_string() ));
      };

      if let Some(value) = iso_table.get(KEY_ISO_ABUNDACE){
        let Some(abundance) = value.as_float() else{
          return Err(CluEError::ExpectedTOMLFloat(value.type_str().to_string() ));
        };
        self.isotopic_distribution.isotope_abundances.push(
            IsotopeAbundance{isotope,abundance} );
      }

      if let Some(iso_prop) = self.isotope_properties.get_mut(&iso_str){
        iso_prop.set_from_toml_table(iso_table,unit_of_energy)?;
      } else{
        let iso_prop 
            = IsotopeProperties::from_toml_table(iso_table,unit_of_energy)?;
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
  pub fn from_toml_table(table: &toml::Table, unit_of_energy: f64) 
      -> Result<Self,CluEError>
  {
    let mut properties = Self::new();  
    properties.set_from_toml_table(table,unit_of_energy)?;

    Ok(properties)
  }
  //----------------------------------------------------------------------------
  fn set_from_toml_table(&mut self, table: &toml::Table, unit_of_energy: f64) 
    -> Result<(),CluEError>
  {
    if let Some(value) = table.get(KEY_ISO_ACTIVE){
      let Some(active) = value.as_bool() else{
        return Err(CluEError::ExpectedTOMLBool(value.type_str().to_string() ));
      };
      self.active = Some(active);
    }

    if let Some(value) = table.get(KEY_ISO_G_MATRIX){
      let tensor = TensorSpecifier::from_toml_value(value.clone(),1.0)?;
      self.g_matrix = Some(tensor);
    }

    if let Some(value) = table.get(KEY_ISO_ELEC_QUADRUPOLE){
      let tensor = TensorSpecifier::from_toml_value(
          value.clone(), unit_of_energy)?;
      self.electric_quadrupole_coupling = Some(tensor);
    }

    if let Some(value) = table.get(KEY_ISO_ZEROFIELD){
      let tensor = TensorSpecifier::from_toml_value(
          value.clone(), unit_of_energy)?;
      self.zerofield_coupling = Some(tensor);
    }

    if let Some(value) = table.get(KEY_ISO_HYPERFINE){
      let tensor = TensorSpecifier::from_toml_value(value.clone(),
          unit_of_energy)?;
      self.hyperfine_coupling = Some(tensor);
    }

    let n_ex = [table.contains_key(KEY_ISO_EXCHANGE_COUPLING),
      table.contains_key(KEY_ISO_C3_TUNNEL_SPLITTING)].iter().map(|b| if *b{1}else{0} ).sum();

    if n_ex > 1{
      return Err(CluEError::TooManyExchangeCouplingsSpecified(n_ex));
    }

    if let Some(value) = table.get(KEY_ISO_EXCHANGE_COUPLING)
    {
      let Some(ex) = value.as_float() else{
        return Err(CluEError::ExpectedTOMLFloat(value.type_str().to_string() ));
      };
      self.exchange_coupling = Some(ex*unit_of_energy);

    } else if let Some(value) = table.get(KEY_ISO_C3_TUNNEL_SPLITTING)
    {
      let Some(nut) = value.as_float() else{
        return Err(CluEError::ExpectedTOMLFloat(value.type_str().to_string() ));
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
/// `isotope_abundances` lists the possible isotope and propabilities.
/// 'void_probability' is the probability that the particle will be included
/// at all.
#[derive(Debug,Clone,PartialEq,Default)]
pub struct IsotopeDistribution{
  // TODO: can isotope_abundances be a HashMap::<Isotope,f64> instead?
  pub isotope_abundances: Vec::<IsotopeAbundance>,
  pub void_probability: Option<f64>,
}
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
  pub fn from_toml_value(value: toml::Value, units: f64) 
      -> Result<Self,CluEError>
  {
    match value{
      toml::Value::Array(array) => {
        let mut tensor = SymmetricTensor3D::from_toml_array(array)?;
        tensor.scale_mut(units);
        Ok(Self::SymmetricTensor3D(tensor))
      },
      toml::Value::Table(table) => Ok(Self::Eig(
            EigSpecifier::from_toml_table(table,units)?)),
      _ => Err(CluEError::TOMLArrayDoesNotSpecifyATensor),
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
}


impl EigSpecifier{
  /// This function creates a new instance of `EigSpecifier`.
  pub fn new() -> Self{
    Self::default()
  }
  //----------------------------------------------------------------------------
  pub fn from_toml_table(table: toml::Table,units: f64) 
      -> Result<Self,CluEError>
  {
  
    let values = if let Some(toml::Value::Array(array)) 
        = table.get(KEY_EIG_VALUES){
      if !are_all_same_type(array){
        return Err(CluEError::TOMLArrayContainsMultipleTypes);
      }
      if array.len() != 3 || !array[0].is_float(){
        return Err(CluEError::TOMLArrayDoesNotSpecifyATensor);
      }
      array.iter().filter_map(|v| v.as_float()).collect::<Vec::<f64>>()

    }else{
      return Err(CluEError::TOMLArrayDoesNotSpecifyATensor);
    }; 
    let values = Some([values[0]*units, values[1]*units, values[2]*units]);

    let Some(axes) = table.get(KEY_EIG_AXES) else{
      return Err(CluEError::TOMLArrayDoesNotSpecifyATensor);
    };

    let x_axis = if let Some(ax) = axes.get(KEY_EIG_X_AXIS){
      Some(VectorSpecifier::from_toml_value(ax.clone())?)
    }else{
      None
    };
    let y_axis = if let Some(ax) = axes.get(KEY_EIG_Y_AXIS){
      Some(VectorSpecifier::from_toml_value(ax.clone())?)
    }else{
      None
    };
    let z_axis = if let Some(ax) = axes.get(KEY_EIG_Z_AXIS){
      Some(VectorSpecifier::from_toml_value(ax.clone())?)
    }else{
      None
    };

    Ok(Self{
        values,
        x_axis,
        y_axis,
        z_axis,
    })
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



