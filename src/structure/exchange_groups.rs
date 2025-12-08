use crate::CluEError;
use crate::physical_constants::ANGSTROM;
use crate::space_3d::Vector3D;
use crate::structure::Structure;

use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This trait allows the exchange group to be translated.
pub trait Translate{ fn translate(&mut self, r: &Vector3D); }

/// This trait gets the corrdinates to the exchange group centroid.
pub trait GetCentroid{fn centroid(&self) -> &Vector3D; }

/// This trait gets the exchange group normal vector.
pub trait GetNormal{fn normal(&self) -> &Vector3D; }

/// This trait gets the indices of the particles in the exchange group.
pub trait GetIndices{fn indices(&self) -> Vec::<usize>; }

/// This trait gets the exchange group's exchange coupling.
pub trait GetExchangeCoupling{fn exchange_coupling(&self) -> f64; }

//pub trait ToStringResult{fn to_string_result(&self, structure: &Structure) 
//  -> Result<String, CluEError>;}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// 'ExchangeGroupManager' contains the info needed to build the effective
/// spin Hamiltonian for identical spins with a non-neglegible exchange rate
/// via quantum tunneling.
/// `exchange_groups` contains a list of the exchange groups.
/// `exchange_group_ids` has an entry for every spin, `None` entries for 
/// particles not part of an exchange group and `Some(idx)` entries for
/// particles in the exchange group `exchange_groups[idx]`.
/// `exchange_couplings` contains the effective couplings for each group.
#[derive(Debug,Clone)]
pub struct ExchangeGroupManager{
  pub exchange_groups: Vec::<ExchangeGroup>,
  pub exchange_group_ids: Vec::<Option<usize>>,
  pub exchange_couplings: Vec::<f64>, // one entry per exchange_group
}

impl ExchangeGroupManager{
  /// This function write the exchange groups in a `Structure` to a csv file.
  pub fn to_csv(&self, filename: &str,structure: &Structure) 
    -> Result<(),CluEError>
  {
    let Ok(file) = File::create(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()) );
    };

    let n_char_f64 = 16;
    let bytes_per_char = 32;
    let n_bytes = 8*bytes_per_char*n_char_f64*(self.exchange_groups.len()+1);

    let mut stream = BufWriter::with_capacity(n_bytes,file);

    let line = "exchange_group,\
center_x,center_y,center_z,\
normal_x,normal_y,normal_z,\
exchange_coupling,res_seq_id,res\n".to_string();
    if stream.write(line.as_bytes()).is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }

    let pdb_origin = structure.pdb_origin.scale(1.0/ANGSTROM);
    for (ii,exchange_group) in self.exchange_groups.iter().enumerate(){

      let methyl_str = exchange_group.to_string();
      let r = exchange_group.centroid().scale(1.0/ANGSTROM);
      let n = exchange_group.normal();
      let j = self.exchange_couplings[ii];

      let indices = exchange_group.indices();

      let res_seq_id = if let Some(idx) = &structure.bath_particles[indices[0]]
        .residue_sequence_number
      {
        format!("{}",idx)
      }else{
        "".to_string()
      };
      let res = if let Some(r) = &structure.bath_particles[indices[0]]
        .residue
      {
        format!("{}",r)
      }else{
        "".to_string()
      };

      let line = format!("{},{},{},{},{},{},{},{},{},{}\n",
          methyl_str, 
          r.x() - pdb_origin.x(), 
          r.y() - pdb_origin.y(), 
          r.z() - pdb_origin.z(), 
          n.x(),n.y(), n.z(), 
          j, res_seq_id,res);
    
      if stream.write(line.as_bytes()).is_err(){
        return Err(CluEError::CannotWriteFile(filename.to_string()) );
      }
    }


    Ok(())
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `ExchangeGroup` lists the implemented exchange groups.
#[derive(Debug, Clone)]
pub enum ExchangeGroup{
  Methyl(C3Rotor),
  PrimaryAmonium(C3Rotor),
}


impl std::fmt::Display for ExchangeGroup {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    let string = match self{
      ExchangeGroup::Methyl(rotor) => format!("methyl_{}",rotor),
      ExchangeGroup::PrimaryAmonium(rotor) 
        => format!("primary_amonium_{}",rotor),
    };
    write!(f,"{}",string)
  }
}
//------------------------------------------------------------------------------
impl GetCentroid for ExchangeGroup{
  // This function implements `GetCentroid` for `ExchangeGroup`.
  fn centroid(&self) -> &Vector3D{
    match self{
      ExchangeGroup::Methyl(rotor) => rotor.centroid(),
      ExchangeGroup::PrimaryAmonium(rotor) => rotor.centroid(),
    }
  }
}

impl GetIndices for ExchangeGroup{
  // This function implements `GetIndices` for `ExchangeGroup`.
  fn indices(&self) -> Vec::<usize>{
    match self{
      ExchangeGroup::Methyl(rotor) => rotor.indices(),
      ExchangeGroup::PrimaryAmonium(rotor) => rotor.indices(),
    }
  }
}

impl GetNormal for ExchangeGroup{
  // This function implements `GetNormal` for `ExchangeGroup`.
  fn normal(&self) -> &Vector3D{
    match self{
      ExchangeGroup::Methyl(rotor) => rotor.normal(),
      ExchangeGroup::PrimaryAmonium(rotor) => rotor.normal(),
    }
  }
}

impl Translate for ExchangeGroup{
  // This function implements `Translate` for `ExchangeGroup`.
  fn translate(&mut self, r: &Vector3D){
    match self{
      ExchangeGroup::Methyl(rotor) => rotor.translate(r),
      ExchangeGroup::PrimaryAmonium(rotor) => rotor.translate(r),
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
/// `C3Rotor` contains the info of a C3 symmetric rotor exchange group.
/// 'center' is the coordinates of the centroid.
/// 'normal' is the normal vector.
/// 'indices' is a list of which particles, via index, make up the rotor.
#[derive(Debug, Clone)]
pub struct C3Rotor{
  pub center: Vector3D,
  pub normal: Vector3D,
  pub indices: [usize;3],
}

//------------------------------------------------------------------------------
impl std::fmt::Display for C3Rotor {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    // The output should be the reference indices, but the rotor stores the
    // structure indices, so a "+1" is used to convert from structure indices
    // to reference indices.
    write!(f,"{}_{}_{}", self.indices[0]+1,self.indices[1]+1,self.indices[2]+1)
  }
}
//------------------------------------------------------------------------------
impl GetCentroid for C3Rotor{
  // This function implements `GetCentroid` for `C3Rotor`.
  fn centroid(&self) -> &Vector3D {                                                
    &self.center                                                           
  }
}
impl C3Rotor{                                                                     
                                                                                 
//------------------------------------------------------------------------------
/// This function builds a `C3Rotor` from the positions and indices of the
/// hydrogens and the position of the carbon atom within a methyl group.
pub fn from(r_carbon: Vector3D,                                                   
    h0: Vector3D,                                                                 
    h1: Vector3D,                                                                 
    h2: Vector3D,
    indices: [usize;3],    
    ) -> Self{                                                                 
                                                                                 
  let center = (&( &h0 + &h1) + &h2).scale(1.0/3.0);                             
                                                                                 
  let x = (&h0 - &center).normalize();                                           
                                                                                 
  let y = (&h1 - &center).normalize();                                           
  let y = &y - &x.scale(x.dot(&y));                                              
                                                                                 
                                                                                 
  let mut normal = x.cross(&y).normalize();                                                  
                                                                                 
  let axis = &center -  &r_carbon;                                               
  if normal.dot(&axis) < 0.0 {                                                   
    normal = normal.scale(-1.0);                                                 
  }                                                                              
                                                                                 
  C3Rotor{ center, normal, indices}
}                                                                                

//------------------------------------------------------------------------------
/// This function return the unit vector normal to the plane of hydrogens in
/// a methyl group, pointing away from the carbon.
pub fn normal(&self) -> &Vector3D {                                                
  &self.normal                                                      
}
}

//------------------------------------------------------------------------------
impl GetIndices for C3Rotor{
  // This function implements `GetIndices` for `C3Rotor`.
  fn indices(&self) -> Vec::<usize>{
    
    let mut out = Vec::<usize>::with_capacity(3);

    for h in self.indices{
      out.push(h);
    }
    out
  }
}
//------------------------------------------------------------------------------
impl Translate for C3Rotor{
  // This function implements `Translate` for `C3Rotor`.
  fn translate(&mut self, r: &Vector3D){
   self.center = &self.center + r;
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


