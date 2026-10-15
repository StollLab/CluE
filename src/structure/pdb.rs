use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::elements::Element;
use crate::physical_constants::{ANGSTROM,PI};
use crate::isotopes::Isotope;
use crate::structure::{Structure, particle::Particle};
use crate::space_3d::Vector3D;

use std::io::BufRead;
use std::io::BufWriter;
use std::collections::HashMap;
use substring::Substring;
use std::fs::File;
use std::io::prelude::*;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This function reads a PDB file and tries to build a `Structure`.
/// If it cannot it will err.
pub fn parse_pdb(filename: &str,target_model: usize) 
  -> Result< Structure, CluEError>
{

   let (n_atoms, _n_models) = count_pdb_atoms(filename)?;
   
   let particle_map = parse_atoms(filename, n_atoms, target_model)?;
   let ParticleMap{ particles: bath_particles, map : serial_to_index,
     hydrogen_symbol_counts} = particle_map;

   // D and T are read as hydrogen.
   for (n, symbol, isotope, name) in [
       (hydrogen_symbol_counts[0], "D", "2H", "deuterium"),
       (hydrogen_symbol_counts[1], "T", "3H", "tritium")]
   {
     if n > 0{
       println!("{} atoms with element symbol \"{}\" were read as hydrogen, \
starting as 1H; to make them {}, put them in a group with {}.abundance = 1.0.",
           n, symbol, name, isotope);
     }
   }

   let connections = parse_connections(filename,serial_to_index)?;

   let (cell_offsets, space_group) = parse_cell(filename)?;
   
   let mut structure = Structure::new(
      bath_particles,
      connections,
      cell_offsets);
   structure.space_group = space_group;

   Ok(structure)
}
//------------------------------------------------------------------------------
// This function tries to parse the crystallographic data, and returns 
// the unit cell edges and the space group.
// If it cannot it will err.
fn parse_cell(filename: &str) 
  -> Result<(Vec::<Vector3D>,Option<String>),CluEError>
{

   let mut cell_offsets = Vec::<Vector3D>::new();
   let mut space_group: Option<String> = None;

   let Ok(file) = std::fs::File::open(filename) else{
     return Err(CluEError::CannotOpenFile(filename.to_string()));
   };
   let lines = std::io::BufReader::new(file).lines();


   for line_result in lines {
     let Ok(line) = line_result else{
       return Err(CluEError::CannotOpenFile(filename.to_string()));
     };
     if line.starts_with("CRYST1"){
       (cell_offsets, space_group) = parse_crystal_line(&line)?;
       break;
     }
   }
    
  Ok((cell_offsets, space_group))
}
//------------------------------------------------------------------------------
// This function tries to parse the connections data.
// If it cannot it will err.
fn parse_connections(filename: &str ,
    serial_to_index: HashMap::<u32,Option<usize>>) 
  -> Result<AdjacencyList,CluEError>
{
   let n_atoms: usize = serial_to_index.iter()
     .filter_map(|(_key,val)| if val.is_some(){Some(1)}else{None})
     .sum();

   let mut connections = AdjacencyList::with_capacity(n_atoms);

   let Ok(file) = std::fs::File::open(filename) else{
     return Err(CluEError::CannotOpenFile(filename.to_string()));
   };

   let lines = std::io::BufReader::new(file).lines();

   'line_loop: for line_result in lines {

     let Ok(line) = line_result else{
       return Err(CluEError::CannotOpenFile(filename.to_string()));
     };

     if line.starts_with("CONECT"){
       let bonds_serial = parse_connections_line(&line)?;

       let mut bonds_index =  Vec::<usize>::with_capacity(bonds_serial.len());
       
       for serial in bonds_serial.iter(){
         
         let Some(idx_opt) = serial_to_index.get(serial) else{
           return Err(CluEError::CannotConvertSerialToIndex(*serial));
         };

         if let Some(idx) = idx_opt{
           bonds_index.push(*idx);
         }else if bonds_index.is_empty(){
           // If the first index correspond to an invalid atom, 
           // all the connections that would be to the invalid atom can be 
           // dropped.
           continue 'line_loop;
         }
       }

       for idx in 1..bonds_index.len(){
         connections.connect(bonds_index[0],bonds_index[idx]);
       }
     }
   }
   Ok(connections)
}
//------------------------------------------------------------------------------
// `ParticleMap` contains a list of particles and a map from the PDB ID to 
// the list index.
// `hydrogen_symbol_counts` counts the atoms with element symbol "D" and "T",
// which are read as hydrogen.
struct ParticleMap{
  particles: Vec::<Particle>,
  map: HashMap::<u32,Option<usize>>,
  hydrogen_symbol_counts: [usize;2],
}
//------------------------------------------------------------------------------
// This function tries to parse the atom data.
// If it cannot it will err.
fn parse_atoms(filename: &str,n_atoms: usize, target_model:usize) 
  -> Result<ParticleMap,CluEError>
{
   let mut serial_to_index = HashMap::<u32,Option<usize>>::new(); 
   let mut unknown_elements = HashMap::<String,usize>::new(); 
   let mut n_elements_from_names = 0;
   let mut hydrogen_symbol_counts = [0,0];

   let mut bath_particles = Vec::<Particle>::with_capacity(n_atoms);

   let mut model_idx = 0;

   let Ok(file) = std::fs::File::open(filename) else{
     return Err(CluEError::CannotOpenFile(filename.to_string()));
   };
   
   let lines = std::io::BufReader::new(file).lines();

   for line_result in lines {

     let Ok(line) = line_result else{
       return Err(CluEError::CannotOpenFile(filename.to_string()));
     };

     if line.starts_with("ENDMDL"){
       model_idx += 1;
       if model_idx > target_model{
         break;
       }
       continue;
     }

     if model_idx != target_model
       || !(line.starts_with("ATOM") || line.starts_with("HETATM")) {
       continue;
     }

     let particle: Particle = match parse_atom_line(&line){
       Ok(p) => p,

       Err((CluEError::CannotParseElement(unk_el),serial)) => {
         
         serial_to_index.insert(serial,None);

         if let Some(value) = unknown_elements.get_mut(&unk_el){
           *value += 1;
         }else{
           unknown_elements.insert(unk_el,1);
         }

         continue;
       }

       Err((err,_)) => return Err(err),
     };

     // Count the elements read from atom names, and the symbols D and T.
     if let Ok((symbol, from_name)) = element_symbol_of_atom_line(&line){
       if from_name{
         n_elements_from_names += 1;
       }
       match symbol.as_deref(){
         Some("D") => hydrogen_symbol_counts[0] += 1,
         Some("T") => hydrogen_symbol_counts[1] += 1,
         _ => (),
       }
     }

     if let Some(serial) = particle.serial{
       if serial_to_index.contains_key(&serial) {
         return Err(CluEError::DuplicateSerial(serial));
       }
       serial_to_index.insert(serial, Some(bath_particles.len()));
     }

     bath_particles.push(particle);
   }

   for (unk_el, unk_count) in unknown_elements.iter(){
     println!("Could not parse {} instances of element \"{}\".",
         unk_count, unk_el);
   }
   if n_elements_from_names > 0{
     println!("Element symbols were taken from the atom names for {} atoms.",
         n_elements_from_names);
   }
   Ok( ParticleMap{
     particles: bath_particles,
     map: serial_to_index,
     hydrogen_symbol_counts})
}
//------------------------------------------------------------------------------
// This function tries to parse the crystallographic data line, and returns
// the unit cell edges and the space group.
// The unit cell is empty for the placeholder 1 A cubic cell that PDB files 
// use for structures without a unit cell, and for cells without volume.
// The space group is `None` if it is blank.
// If it cannot it will err.
// The expected format is as follows.
// https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html
/*
COLUMNS       DATA  TYPE    FIELD          DEFINITION
-------------------------------------------------------------
 1 -  6       Record name   "CRYST1"
 7 - 15       Real(9.3)     a              a (Angstroms).
16 - 24       Real(9.3)     b              b (Angstroms).
25 - 33       Real(9.3)     c              c (Angstroms).
34 - 40       Real(7.2)     alpha          alpha (degrees).
41 - 47       Real(7.2)     beta           beta (degrees).
48 - 54       Real(7.2)     gamma          gamma (degrees).
56 - 66       LString       sGroup         Space  group.
67 - 70       Integer       z              Z value.
*/
fn parse_crystal_line(line: &str) 
  -> Result<(Vec::<Vector3D>,Option<String>),CluEError>
{
  // The numbers are required and must be complete.
  let parse_columns = |first: usize, last: usize| -> Result<f64,CluEError>{
    if line.len() < last{
      return Err(CluEError::CannotParseLine(line.to_string()));
    }
    match pdb_columns(line, first, last)?.map(|s| s.trim().parse::<f64>()){
      Some(Ok(x)) => Ok(x),
      _ => Err(CluEError::CannotParseLine(line.to_string())),
    }
  };
  if let (Ok(a),Ok(b),Ok(c),Ok(mut alpha),Ok(mut beta),Ok(mut gamma)) = (
      parse_columns(7,15),
      parse_columns(16,24),
      parse_columns(25,33),
      parse_columns(34,40),
      parse_columns(41,47),
      parse_columns(48,54) ){

    let space_group = pdb_columns(line, 56, 66)?
      .map(|s| s.trim().to_string()).filter(|s| !s.is_empty());

    // Placeholder cell for structures without a unit cell.
    if [a,b,c] == [1.0,1.0,1.0] && [alpha,beta,gamma] == [90.0,90.0,90.0]{
      return Ok((Vec::<Vector3D>::new(), space_group));
    }

    alpha *= PI/180.0;
    beta *= PI/180.0;
    gamma *= PI/180.0;

    let a = a*ANGSTROM;
    let b = b*ANGSTROM;
    let c = c*ANGSTROM;

    let a_vec = Vector3D::from([a, 0.0, 0.0]);
    let b_vec = Vector3D::from([b*f64::cos(gamma), b*f64::sin(gamma), 0.0]);
    let cx = c * f64::cos(beta);
    let cy = c*(f64::cos(alpha) - f64::cos(beta)*f64::cos(gamma))/f64::sin(gamma);
    let cz = f64::sqrt( c*c -cx*cx - cy*cy);
    let c_vec = Vector3D::from([cx, cy, cz]);
  
    let cell_edges = vec![a_vec,b_vec,c_vec];
    if crate::structure::reciprocal_vectors(&cell_edges).is_none(){
      return Ok((Vec::<Vector3D>::new(), space_group));
    }

    return Ok((cell_edges, space_group));
  }
  Err(CluEError::CannotParseLine(line.to_string() ))
}

//------------------------------------------------------------------------------
// This function tries to parse an atom data line.
// If it cannot it will err.
// https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
/*
COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
*/
fn parse_atom_line(line: &str) -> Result<Particle,(CluEError,u32)>{

  let cannot_parse = || (CluEError::CannotParseLine(line.to_string()),0);
  let columns = |first: usize, last: usize| -> Result<Option<&str>,_>{
    pdb_columns(line, first, last).map_err(|_| cannot_parse())
  };
  // Required fields must be complete.
  let required_columns = |first: usize, last: usize| -> Result<&str,_>{
    if line.len() < last{
      return Err(cannot_parse());
    }
    columns(first,last)?.ok_or_else(cannot_parse)
  };

  // Required fields.
  let mut xyz = [0.0; 3];
  for (ii, first) in [31, 39, 47].into_iter().enumerate(){
    let Ok(x) = required_columns(first, first + 7)?.trim().parse::<f64>()
    else{
      return Err(cannot_parse());
    };
    xyz[ii] = x;
  }
  let coordinates = Vector3D::from(xyz).scale(ANGSTROM);

  let Ok(serial) = required_columns(7,11)?.trim().parse::<u32>() else{
    return Err(cannot_parse());
  };

  // Optional fields.
  let name_columns = columns(13,16)?;
  let name = name_columns.map(|s| s.trim().to_string());
  let residue = columns(18,20)?.map(|s| s.trim().to_string());
  let chain_id = columns(22,22)?.map(|s| s.to_string());
  let residue_sequence_number = columns(23,26)?
    .and_then(|s| s.trim().parse::<u32>().ok());

  // The element is in columns 77-78; if those are missing or blank,
  // it is read from the atom name.
  let (symbol, from_name) = element_symbol_of_atom_line(line)
    .map_err(|_| cannot_parse())?;
  let Some(Ok(element)) = symbol.as_deref().map(Element::from) else{
    // Report the symbol as written.
    let written = if from_name{
      name_columns.unwrap_or("").trim()
    }else{
      columns(77,78)?.unwrap_or("").trim()
    };
    return Err((CluEError::CannotParseElement(written.to_string()), serial));
  };

  Ok(Particle{
      element,
      isotope: Isotope::most_common_for(&element),
      coordinates,
      active: true,
      serial: Some(serial),
      name,
      residue,
      residue_sequence_number,
      chain_id,
      //state: Vec::<Complex64>::new(),
      })
}
//------------------------------------------------------------------------------
// This function returns the text in columns `first` to `last` (1-based, 
// inclusive, as in the PDB format specification) of `line`.
// It returns `Ok(None)` if the line does not reach column `first`, 
// and an error if the columns do not fall on character boundaries.
fn pdb_columns(line: &str, first: usize, last: usize) 
  -> Result<Option<&str>,CluEError>
{
  if line.len() < first{
    return Ok(None);
  }
  match line.get(first - 1..usize::min(last, line.len())){
    Some(s) => Ok(Some(s)),
    None => Err(CluEError::CannotParseLine(line.to_string())),
  }
}
//------------------------------------------------------------------------------
// This function writes an element symbol with an upper case first letter and
// lower case other letters, since PDB files write symbols in upper case,
// for example "FE" for "Fe".
fn capitalize_element_symbol(symbol: &str) -> String{
  let mut chars = symbol.chars();
  match chars.next(){
    Some(first) => first.to_uppercase().chain(chars.flat_map(|c| 
          c.to_lowercase())).collect(),
    None => String::new(),
  }
}
//------------------------------------------------------------------------------
// This function infers the element from the atom name (columns 13-16), 
// for PDB files without element columns.
// In the PDB format, the element symbol is right-justified in columns 13-14:
// if column 13 is blank or a digit, column 14 is the element (" CA " is 
// carbon, "1HB " is hydrogen); a name of at most two letters starting in 
// column 13 is a two-letter element ("CA  " is calcium, "FE  " is iron);
// otherwise column 13 is the element ("HG21" and "HW1 " are hydrogen).
// It returns the element symbol, or `None` if no element is found.
fn element_symbol_from_atom_name(name: &str) -> Option<String>{
  let chars = name.chars().collect::<Vec::<char>>();
  let col13 = *chars.first()?;
  let valid = |symbol: String| -> Option<String>{
    Element::from(&symbol).ok().map(|_| symbol)
  };

  if col13 == ' ' || col13.is_ascii_digit(){
    let col14 = *chars.get(1)?;
    return valid(capitalize_element_symbol(&col14.to_string()));
  }
  if !col13.is_ascii_alphabetic(){
    return None;
  }

  let rest_is_blank = chars.iter().skip(2).all(|c| *c == ' ');
  if let Some(&col14) = chars.get(1){
    if col14.is_ascii_alphabetic() && rest_is_blank{
      let symbol = capitalize_element_symbol(&format!("{}{}",col13,col14));
      if let Some(symbol) = valid(symbol){
        return Some(symbol);
      }
    }
  }
  valid(capitalize_element_symbol(&col13.to_string()))
}
//------------------------------------------------------------------------------
// This function returns the element symbol of an ATOM or HETATM line, 
// capitalized as in "Fe", and whether it was read from the atom name.
// The symbol is read from columns 77-78, or from the atom name (columns 
// 13-16) if those are missing or blank; it is `None` if the atom name gives 
// no element.  A symbol from columns 77-78 is returned even if it is not an 
// element.
fn element_symbol_of_atom_line(line: &str) 
  -> Result<(Option<String>, bool),CluEError>
{
  if let Some(symbol) = pdb_columns(line,77,78)?.map(|s| s.trim())
    .filter(|s| !s.is_empty())
  {
    return Ok((Some(capitalize_element_symbol(symbol)), false));
  }
  let name = pdb_columns(line,13,16)?.unwrap_or("");
  Ok((element_symbol_from_atom_name(name), true))
}
//------------------------------------------------------------------------------
// This function tries to parse a connections data line.
// If it cannot it will err.
// https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html
/*
COLUMNS       DATA  TYPE      FIELD        DEFINITION
-------------------------------------------------------------------------
 1 -  6        Record name    "CONECT"
 7 - 11        Integer        serial       Atom  serial number
12 - 16        Integer        serial       Serial number of bonded atom
17 - 21        Integer        serial       Serial  number of bonded atom
22 - 26        Integer        serial       Serial number of bonded atom
27 - 31        Integer        serial       Serial number of bonded atom
*/
fn parse_connections_line(conect_line: &str) -> Result<Vec::<u32>,CluEError>{
  let Some(line) = conect_line.get(6..) else{
    return Err(CluEError::CannotParseLine(conect_line.to_string() ));
  };
  let line = line.trim_end();

  let num = line.len();
  let n = 5;
  let mut serials = Vec::<u32>::with_capacity(num/n);
  let mut idx = 0;
  while idx < num {

    let idx_end = usize::min(idx+n,num);
    
    let connect: u32;
    
    if let Some(Ok(serial)) = line.get(idx .. idx_end)
      .map(|s| s.trim().parse::<u32>()){
      connect = serial;
    }else{
      return Err(CluEError::CannotParseLine(conect_line.to_string() ));
    }
    
    serials.push(connect);
    
    idx += n;
  }

  Ok(serials)

}
//------------------------------------------------------------------------------
// This function tries to count the number of atoms in the PDB.
// If it cannot it will err.
fn count_pdb_atoms(filename: &str) -> Result<(usize,usize),CluEError> {

   let mut count = 0;
   let mut model_count = 0;
   let mut increment = 1;
   let Ok(file) = File::open(filename) else{
     return Err(CluEError::CannotOpenFile(filename.to_string()));
   };
   let lines = std::io::BufReader::new(file).lines();
   for line_result in lines {
     let Ok(line) = line_result else{
       return Err(CluEError::CannotOpenFile(filename.to_string()));
     };
     if line.starts_with("MODEL"){
       model_count += 1;
     }else if line.starts_with("ENDMDL"){
       increment = 0;
     }else if line.starts_with("ATOM") || line.starts_with("HETATM") {
       count += increment;
     }
   }

   if model_count==0{
     model_count += 1;
   }
   Ok((count,model_count))
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl Structure{
  /// This function writes a `Structure` in PDB format.
  pub fn write_pdb(&self,filename: &str) -> Result<(),CluEError>{

    let n_active = self.number_active();

    let pdb_chars_per_line = 80;
    let bytes_per_char = 32;
    let n_bytes = (n_active +1)*pdb_chars_per_line*bytes_per_char;

    let Ok(file) = File::create(filename)
    else {
      return Err(CluEError::CannotOpenFile(filename.to_string()) );
    };
    let mut stream = BufWriter::with_capacity(n_bytes,file);

    let pdb_det_str = self.get_detected_particle_pdb_string()?;

    if stream.write(pdb_det_str.as_bytes()).is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }

    let mut serial_num = 1;
    let mut line  = String::new();
    for particle in self.bath_particles.iter(){
      if !particle.active {continue;}
      self.set_bath_particle_pdb_string(&mut line, particle, &mut serial_num)?;

      
      let stream_result = stream.write(line.as_bytes());
      if stream_result.is_err(){
        return Err(CluEError::CannotWriteFile(filename.to_string()) );
      }
      

    }

    let stream_result = stream.flush();
    if stream_result.is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }

    Ok(())
  }
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// This function write the detected particle in PDB format.
fn get_detected_particle_pdb_string(&self) 
  -> Result<String,CluEError>
{
  let pdb_str: String;

  // Write coordinate lines.
  if let Some(detected_particle) = &self.detected_particle{
    let det_isotope_str = &detected_particle.isotope.to_string();

    let serial = format!("{: >5}", 0.to_string());
    let name =  format!(" {: <4}",det_isotope_str);
    let alt_loc = " ".to_string();
    let res_name = String::from("DET");
    let chain_id = "  ".to_string();
    let res_seq = format!("{: >4}",0);
    let icode = "    ".to_string();
    let Some(r) = detected_particle.weighted_coordinates.weighted_mean() else{
      return Err(CluEError::InvalidSpinGridWeights(
            detected_particle.weighted_coordinates.weights.iter().fold(0.0, |sum, w| sum + w)));
    };
    let r0 = &self.pdb_origin;
    let x = format!("{: >8}",
      ((r[0] -r0.x())/ANGSTROM).to_string().substring(0,7));
    let y = format!("{: >8}",
      ((r[1]-r0.y())/ANGSTROM).to_string().substring(0,7));
    let z = format!("{: >8}",
      ((r[2]-r0.z())/ANGSTROM).to_string().substring(0,7));
    let occupancy = "  1.00";
    let temp_factor = "  0.00";
    let element =  format!("          {: >2}",det_isotope_str);
    pdb_str = format!("HETATM{}{}{}{}{}{}{}{}{}{}{}{}{}\n"
        ,serial, name, alt_loc, res_name, chain_id,
        res_seq, icode, x, y, z,
        occupancy, temp_factor, element);

    
  }else{
    pdb_str = String::new();
  }
  Ok(pdb_str)
}  
//------------------------------------------------------------------------------
// This function sets 'line' and a PDB formated line for a `Particle`.
fn set_bath_particle_pdb_string(&self, line: &mut String, particle: &Particle, 
    serial_num: &mut usize) 
  -> Result<(),CluEError>
{

  let serial = format!("{: >5}", serial_num.to_string());
  let name = format!(" {: <4}",particle.element.to_string());
  let alt_loc = " ".to_string();
  let res_name: String;
  if let Some(res) = &particle.residue{
    res_name = format!("{: >3}",res);
  }else{
    res_name = "UNK".to_string();
  }
  let chain_id = "  ".to_string();
  let res_seq: String;
  if let Some(res_seq_num) = particle.residue_sequence_number{
    res_seq = format!("{: >4}",res_seq_num);
  }else{
    res_seq = "   0".to_string();
  }
  let icode = "    ".to_string();
  let r = &particle.coordinates - &self.pdb_origin;
  let x = format!("{: >8}",
      (r.x()/ANGSTROM).to_string().substring(0,7));
  let y = format!("{: >8}",
      (r.y()/ANGSTROM).to_string().substring(0,7));
  let z = format!("{: >8}",
      (r.z()/ANGSTROM).to_string().substring(0,7));

  let occupancy = "  1.00";
  let temp_factor = "  0.00";
  let element = format!("          {: >2}",particle.element.to_string());

  *line = format!("HETATM{}{}{}{}{}{}{}{}{}{}{}{}{}\n"
      ,serial, name, alt_loc, res_name, chain_id,
      res_seq, icode, x, y, z,
      occupancy, temp_factor, element);

  *serial_num += 1;


  Ok(())
}
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests {
  use super::*;
  //----------------------------------------------------------------------------
  // This function writes a CRYST1 line in the PDB format.
  fn cryst1_line(edges: [f64;3], angles: [f64;3], space_group: &str) 
    -> String
  {
    format!("CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} {:<11}{:>4}",
        edges[0], edges[1], edges[2], angles[0], angles[1], angles[2],
        space_group, 1)
  }
  //----------------------------------------------------------------------------
  // This function writes an ATOM line in the PDB format, with the atom name 
  // `name` (columns 13-16) and the element `element` (columns 77-78).
  fn atom_line(name: &str, element: &str) -> String{
    assert_eq!(name.len(),4);
    assert_eq!(element.len(),2);
    format!("ATOM  {:>5} {}{}{:>3} {}{:>4}{}   {:>8.3}{:>8.3}{:>8.3}\
{:>6.2}{:>6.2}          {}",
        7, name, " ", "SOL", "A", 12, " ", 1.0, -2.5, 30.25, 1.0, 0.0, element)
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_atom_line_columns(){
    let line = atom_line(" OW ", " O");
    assert_eq!(line.len(),78);
    let particle = parse_atom_line(&line).unwrap();
    assert_eq!(particle.element, Element::Oxygen);
    assert_eq!(particle.serial, Some(7));
    assert_eq!(particle.name, Some("OW".to_string()));
    assert_eq!(particle.residue, Some("SOL".to_string()));
    assert_eq!(particle.chain_id, Some("A".to_string()));
    assert_eq!(particle.residue_sequence_number, Some(12));
    assert_eq!(particle.coordinates, 
        Vector3D::from([1.0, -2.5, 30.25]).scale(ANGSTROM));

    // The element is in column 77 and trailing spaces are removed.
    let line = atom_line(" CA ", "C ");
    let line = line.trim_end();
    assert_eq!(line.len(),77);
    assert_eq!(parse_atom_line(line).unwrap().element, Element::Carbon);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_lines_that_are_too_short(){
    // Lines that are too short are errors rather than panics.
    let line = atom_line(" OW ", " O");
    for n in [10, 30, 50]{
      assert!(matches!(parse_atom_line(&line[..n]),
            Err((CluEError::CannotParseLine(_),_))));
    }
    let line = cryst1_line([50.0, 60.0, 70.0], [90.0, 90.0, 90.0], "P 1");
    assert!(matches!(parse_crystal_line(&line[..40]),
          Err(CluEError::CannotParseLine(_))));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_lines_with_non_ascii_characters(){
    // Fixed columns that split characters are errors rather than panics.
    // The two-byte character Å replaces column 12, so that the atom name 
    // (columns 13-16) starts in the middle of it.
    let line = atom_line(" OW ", " O");
    let line = format!("{}Å{}",&line[..11],&line[12..]);
    assert!(parse_atom_line(&line).is_err());
    // Å replaces column 15, so that the edge a (columns 7-15) ends in the 
    // middle of it.
    let line = cryst1_line([50.0, 60.0, 70.0], [90.0, 90.0, 90.0], "P 1");
    let line = format!("{}Å{}",&line[..14],&line[15..]);
    assert!(parse_crystal_line(&line).is_err());
    // The second serial (columns 12-16) ends in the middle of Å.
    assert!(parse_connections_line("CONECT    1    Å").is_err());
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_atom_line_upper_case_elements(){
    // PDB files write element symbols in upper case.
    for (symbol, element) in [(" C", Element::Carbon), ("FE", Element::Iron),
        ("NA", Element::Sodium), ("CL", Element::Chlorine), 
        ("ZN", Element::Zinc), ("Fe", Element::Iron)]
    {
      let line = atom_line(" X  ", symbol);
      assert_eq!(parse_atom_line(&line).unwrap().element, element);
    }
    // Virtual sites are not elements.
    for symbol in ["Mw","MW"," M"]{
      let line = atom_line("MW  ", symbol);
      assert!(matches!(parse_atom_line(&line),
            Err((CluEError::CannotParseElement(_),7))));
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_atom_line_element_from_name(){
    // Without element columns, the element is read from the atom name, 
    // where it is right-justified in columns 13-14.
    for (name, element) in [(" CA ", Element::Carbon), ("CA  ", Element::Calcium),
        ("FE  ", Element::Iron), ("1HB ", Element::Hydrogen), 
        (" HW1", Element::Hydrogen), ("HW1 ", Element::Hydrogen),
        (" OW ", Element::Oxygen), ("HG21", Element::Hydrogen),
        ("HD11", Element::Hydrogen)]
    {
      let line = atom_line(name, "  ");
      assert_eq!(parse_atom_line(&line[..66]).unwrap().element, element,
          "{}", name);
      assert_eq!(parse_atom_line(&line).unwrap().element, element,
          "{}", name);
    }
    let line = atom_line("MW  ", "  ");
    assert!(matches!(parse_atom_line(&line[..66]),
          Err((CluEError::CannotParseElement(_),7))));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_detected_particle_pdb_string_at_weighted_mean(){
    // The detected spin is written at the weighted mean of its grid,
    // in the frame of the structure file.
    use crate::config::{Config, DetectedSpinCoordinates};
    use crate::integration_grid::IntegrationGrid;
    use rand::SeedableRng;
    let mut grid = IntegrationGrid::new(3);
    grid.push(vec![36.0e-10, 37.0e-10, 38.0e-10], 0.9).unwrap();
    grid.push(vec![46.0e-10, 37.0e-10, 38.0e-10], 0.1).unwrap();
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::ProbabilityDistribution(grid));
    config.set_defaults().unwrap();
    let mut structure = parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    let mut rng = rand_chacha::ChaCha20Rng::seed_from_u64(0);
    structure.build_primary_structure(&mut rng,&config).unwrap();

    let line = structure.get_detected_particle_pdb_string().unwrap();
    let xyz = [&line[30..38], &line[38..46], &line[46..54]].map(|s| 
        s.trim().parse::<f64>().unwrap());
    for (x, expected) in xyz.iter().zip([37.0, 37.0, 38.0]){
      assert!((x - expected).abs() < 1e-3, "{}", line);
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_hydrogen_isotope_symbols(){
    // D and T are read as hydrogen, starting as 1H, and are counted.
    let text = std::fs::read_to_string("./assets/TEMPO.pdb").unwrap();
    let mut n_h = 0;
    let lines = text.lines().map(|line|{
      if line.starts_with("ATOM") && line.len() >= 78 
        && line[76..78].trim() == "H"
      {
        n_h += 1;
        let symbol = if n_h == 1 {" T"} else {" D"};
        format!("{}{}{}", &line[..76], symbol, &line[78..])
      }else{
        line.to_string()
      }
    }).collect::<Vec::<String>>();
    assert_eq!(n_h, 18);
    let filename = std::env::temp_dir().join(format!(
          "clue_oxide_test_tempo_dt_{}.pdb", std::process::id()));
    std::fs::write(&filename, lines.join("\n")).unwrap();
    let filename = filename.to_str().unwrap().to_string();

    let particle_map = parse_atoms(&filename, 29, 0).unwrap();
    assert_eq!(particle_map.hydrogen_symbol_counts, [17,1]);
    let structure = parse_pdb(&filename,0).unwrap();
    std::fs::remove_file(&filename).unwrap();
    let hydrogens = structure.bath_particles.iter()
      .filter(|p| p.element == Element::Hydrogen).collect::<Vec::<_>>();
    assert_eq!(hydrogens.len(), 18);
    assert!(hydrogens.iter().all(|p| p.isotope == Isotope::Hydrogen1));

    // A deuteron named in the atom name, without element columns.
    let line = atom_line(" D1 ", "  ");
    assert_eq!(element_symbol_of_atom_line(&line[..66]),
        Ok((Some("D".to_string()), true)));
    assert_eq!(parse_atom_line(&line[..66]).unwrap().element, 
        Element::Hydrogen);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_filter_missing_residue_sequence_number(){
    // A residue sequence number that is not a number, as in hybrid-36, 
    // is missing, so it is not selected by residue_sequence_numbers.
    use crate::structure::particle_filter::ParticleFilter;
    let line = atom_line(" OW ", " O");
    let line = format!("{}A000{}", &line[..22], &line[26..]);
    let particle = parse_atom_line(&line).unwrap();
    assert_eq!(particle.residue_sequence_number, None);

    let mut structure = Structure::new(vec![particle], 
        AdjacencyList::with_capacity(1), Vec::<Vector3D>::new());
    structure.primary_cell_indices = vec![0];
    structure.cell_indices = vec![vec![Some(0)]];
    structure.cell_ids = vec![0];

    let mut filter = ParticleFilter::new();
    filter.residue_sequence_numbers = vec![1];
    assert!(filter.filter(&structure).is_empty());

    let mut filter = ParticleFilter::new();
    filter.not_residue_sequence_numbers = vec![1];
    assert_eq!(filter.filter(&structure), vec![0]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_pdb_fields_unchanged(){
    // Files with element columns are read as the columns specify.
    for filename in ["./assets/TEMPO.pdb", "./assets/water_gly_20pc_vol.pdb"]{
      let structure = parse_pdb(filename,0).unwrap();
      let text = std::fs::read_to_string(filename).unwrap();
      let lines = text.lines()
        .filter(|l| l.starts_with("ATOM") || l.starts_with("HETATM"))
        .filter(|l| Element::from(l[76..78].trim()).is_ok())
        .collect::<Vec::<&str>>();
      assert_eq!(lines.len(), structure.bath_particles.len());
      for (line, particle) in lines.iter().zip(structure.bath_particles.iter()){
        assert_eq!(particle.element, Element::from(line[76..78].trim()).unwrap());
        assert_eq!(particle.serial, line[6..11].trim().parse::<u32>().ok());
        assert_eq!(particle.name.as_deref(), Some(line[12..16].trim()));
        assert_eq!(particle.residue.as_deref(), Some(line[17..20].trim()));
        assert_eq!(particle.chain_id.as_deref(), Some(&line[21..22]));
        assert_eq!(particle.residue_sequence_number, 
            line[22..26].trim().parse::<u32>().ok());
        let r = Vector3D::from([
            line[30..38].trim().parse::<f64>().unwrap(),
            line[38..46].trim().parse::<f64>().unwrap(),
            line[46..54].trim().parse::<f64>().unwrap()]).scale(ANGSTROM);
        assert_eq!(particle.coordinates, r);
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_crystal_line_space_group(){
    let edges = [50.0, 60.0, 70.0];
    let angles = [90.0, 90.0, 90.0];
    for (space_group, expected) in [("P 21 21 21", Some("P 21 21 21")),
        ("P 1", Some("P 1")), ("", None)]
    {
      let line = cryst1_line(edges, angles, space_group);
      let (cell_edges, sg) = parse_crystal_line(&line).unwrap();
      assert_eq!(sg, expected.map(|s| s.to_string()));
      assert_eq!(cell_edges, vec![
          Vector3D::from([50.0, 0.0, 0.0]).scale(ANGSTROM),
          Vector3D::from([0.0, 60.0, 0.0]).scale(ANGSTROM),
          Vector3D::from([0.0, 0.0, 70.0]).scale(ANGSTROM),
      ]);
    }

    let structure = parse_pdb("./assets/TEMPO.pdb",0).unwrap();
    assert_eq!(structure.space_group, Some("P 1".to_string()));
    let structure = parse_pdb("./assets/water_gly_20pc_vol.pdb",0).unwrap();
    assert_eq!(structure.space_group, None);
    assert_eq!(structure.cell_offsets.len(), 3);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_crystal_line_without_unit_cell(){
    // The placeholder cell of structures without a unit cell, 
    // and cells without volume, are ignored.
    for (edges, angles) in [
      ([1.0, 1.0, 1.0], [90.0, 90.0, 90.0]),
      ([0.0, 0.0, 0.0], [90.0, 90.0, 90.0]),
      ([10.0, 10.0, 10.0], [90.0, 90.0, 180.0]),
    ]{
      let line = cryst1_line(edges, angles, "P 1");
      let (cell_edges, space_group) = parse_crystal_line(&line).unwrap();
      assert!(cell_edges.is_empty(), "{}", line);
      assert_eq!(space_group, Some("P 1".to_string()));
    }
  }
  #[test]
  fn test_parse_pdb_water(){
    let filename = "./assets/water.pdb";
    let structure = parse_pdb(&filename,0).unwrap();
    assert_eq!(structure.bath_particles.len(),3);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_pdb_tempo_wat_gly_7nm(){
    // n    : Molecules    : Hydrons 
    // 1    : TEMPO        :    18
    // 1500 : glycerols    : 12000
    // 7469 : waters       : 14938
    //
    // total hydrons =  26956.
    //
    let filename = "./assets/TEMPO_wat_gly_70A.pdb";
    let structure = parse_pdb(&filename,0).unwrap();
    assert_eq!(structure.bath_particles.len(),43436);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_connections(){

    let filename = "./assets/water.pdb";
    let mut serial_to_index = HashMap::<u32,Option<usize>>::new();
    serial_to_index.insert(21030,Some(0));
    serial_to_index.insert(21031,Some(1));
    serial_to_index.insert(21032,Some(2));
    serial_to_index.insert(21033,None);
    let adjacency = parse_connections(&filename,serial_to_index).unwrap();

    assert!(adjacency.are_connected(0,1));
    assert!(adjacency.are_connected(0,2));
    assert!(!adjacency.are_connected(0,3));
    assert!(!adjacency.are_connected(1,3));
    assert!(!adjacency.are_connected(1,2));
    assert!(!adjacency.are_connected(2,3));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_parse_pdb_tempo(){
    let filename = "./assets/TEMPO.pdb";
    let structure = parse_pdb(&filename,0).unwrap();

    assert_eq!(structure.bath_particles.len(),29);
    assert_eq!(structure.bath_particles[27].element,Element::Nitrogen);
    assert_eq!(structure.bath_particles[28].element,Element::Oxygen);
    assert!(structure.connections.are_connected(27,28));

    assert_eq!(structure.bath_particles[27].serial, Some(28));
    assert_eq!(structure.bath_particles[27].residue,
        Some("TEM".to_string()));
    assert_eq!(structure.bath_particles[27].residue_sequence_number,
        Some(1));
    assert_eq!(structure.bath_particles[27].coordinates.x(),
        36.440*ANGSTROM);
    assert_eq!(structure.bath_particles[27].coordinates.y(),
        36.900*ANGSTROM);
    assert_eq!(structure.bath_particles[27].coordinates.z(),
        37.100*ANGSTROM);

    assert_eq!(structure.cell_offsets[0].x(),72.5676*ANGSTROM);
    assert_eq!(structure.cell_offsets[0].y(),0.0);
    assert_eq!(structure.cell_offsets[0].z(),0.0);

    assert!((structure.cell_offsets[1].x()-0.0).abs() < 1e12);
    assert_eq!(structure.cell_offsets[1].y(),72.5676*ANGSTROM);
    assert!((structure.cell_offsets[1].z()-0.0).abs() < 1e12);

    assert!((structure.cell_offsets[2].x()-0.0).abs() < 1e12);
    assert!((structure.cell_offsets[2].y()-0.0).abs() < 1e12);
    assert_eq!(structure.cell_offsets[2].z(),72.5676*ANGSTROM);


    //         O
    //   3HC   N  CH3
    //  3HC-C    C-CH3
    //    2HC    CH2
    //        CH2

    assert!(structure.connections.are_connected(0,1));
    assert!(structure.connections.are_connected(0,5));
    assert!(structure.connections.are_connected(0,9));
    assert!(structure.connections.are_connected(0,27));
    assert!(structure.connections.are_connected(1,2));
    assert!(structure.connections.are_connected(1,3));
    assert!(structure.connections.are_connected(1,4));
    assert!(structure.connections.are_connected(5,6));
    assert!(structure.connections.are_connected(5,7));
    assert!(structure.connections.are_connected(5,8));
    assert!(structure.connections.are_connected(9,10));
    assert!(structure.connections.are_connected(9,11));
    assert!(structure.connections.are_connected(9,12));
    assert!(structure.connections.are_connected(12,13));
    assert!(structure.connections.are_connected(12,14));
    assert!(structure.connections.are_connected(12,15));
    assert!(structure.connections.are_connected(15,16));
    assert!(structure.connections.are_connected(15,17));
    assert!(structure.connections.are_connected(15,18));
    assert!(structure.connections.are_connected(18,19));
    assert!(structure.connections.are_connected(18,23));
    assert!(structure.connections.are_connected(18,27));
    assert!(structure.connections.are_connected(19,20));
    assert!(structure.connections.are_connected(19,21));
    assert!(structure.connections.are_connected(19,22));
    assert!(structure.connections.are_connected(23,24));
    assert!(structure.connections.are_connected(23,25));
    assert!(structure.connections.are_connected(23,26));
    assert!(structure.connections.are_connected(27,28));

    assert!(!structure.connections.are_connected(0,28));
  }
  //----------------------------------------------------------------------------
}
