use crate::clue_errors::CluEError;
use crate::cluster::adjacency::AdjacencyList;
use crate::elements::Element;
use crate::physical_constants::{ANGSTROM,PI};
use crate::isotopes::Isotope;
use crate::structure::{Structure, particle::Particle};
use crate::space_3d::Vector3D;

use num_complex::Complex64;

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
   let ParticleMap{ particles: bath_particles, map : serial_to_index} 
     = particle_map;

   let connections = parse_connections(filename,serial_to_index)?;

   let cell_offsets = parse_cell(filename)?;
   
   Ok( Structure::new(
      bath_particles,
      connections,
      cell_offsets)
    )
}
//------------------------------------------------------------------------------
// This function tries to parse the crystallographic data.
// If it cannot it will err.
fn parse_cell(filename: &str) -> Result<Vec::<Vector3D>,CluEError>
{

   let mut cell_offsets = Vec::<Vector3D>::new();

   let Ok(file) = std::fs::File::open(filename) else{
     return Err(CluEError::CannotOpenFile(filename.to_string()));
   };
   let lines = std::io::BufReader::new(file).lines();


   for line_result in lines {
     let Ok(line) = line_result else{
       return Err(CluEError::CannotOpenFile(filename.to_string()));
     };
     if line.contains("CRYST1"){
       cell_offsets = parse_crystal_line(&line)?;
       break;
     }
   }
    
  Ok(cell_offsets)
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

     if line.contains("CONECT"){
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
struct ParticleMap{
  particles: Vec::<Particle>,
  map: HashMap::<u32,Option<usize>>
}
//------------------------------------------------------------------------------
// This function tries to parse the atom data.
// If it cannot it will err.
fn parse_atoms(filename: &str,n_atoms: usize, target_model:usize) 
  -> Result<ParticleMap,CluEError>
{
   let mut serial_to_index = HashMap::<u32,Option<usize>>::new(); 
   let mut unknown_elements = HashMap::<String,usize>::new(); 

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

     if line.contains("ENDMDL"){
       model_idx += 1;
       if model_idx > target_model{
         break;
       }
       continue;
     }

     if model_idx != target_model
       || !(line.contains("ATOM") || line.contains("HETATM")) {
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

     if let Some(serial) = particle.serial{
       serial_to_index.entry(serial)
         .or_insert_with(|| Some(bath_particles.len()));
     }

     bath_particles.push(particle);
   }

   for (unk_el, unk_count) in unknown_elements.iter(){
     println!("Could not parse {} instances of element \"{}\".",
         unk_count, unk_el);
   }

   Ok( ParticleMap{
     particles: bath_particles,
     map: serial_to_index})
}
//------------------------------------------------------------------------------
// This function tries to parse the crystallographic data line.
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
  -> Result<Vec::<Vector3D>,CluEError>
{
  if let (Ok(a),Ok(b),Ok(c),Ok(mut alpha),Ok(mut beta),Ok(mut gamma)) = (
      line[6..=14].trim().parse::<f64>(),
      line[15..=23].trim().parse::<f64>(),
      line[24..=32].trim().parse::<f64>(),
      line[33..=39].trim().parse::<f64>(),
      line[40..=46].trim().parse::<f64>(),
      line[47..=53].trim().parse::<f64>() ){

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
  

    return Ok(vec![a_vec,b_vec,c_vec]);
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

  let serial_result = line[6..=10].trim().parse::<u32>();
  let name = Some(line[12..=15].trim().to_string());
  let x_coor = line[30..=37].trim().parse::<f64>();
  let y_coor = line[38..=45].trim().parse::<f64>();
  let z_coor = line[46..=53].trim().parse::<f64>();
  let residue = Some(line[17..=19].trim().to_string());
  let residue_sequence_number = line[22..=25].trim().parse::<u32>();
  let chain_id = Some(line[21..=21].to_string());



  let coordinates: Vector3D;
  if let (Ok(x),Ok(y), Ok(z)) = (x_coor, y_coor, z_coor) {
    coordinates = Vector3D::from([x*ANGSTROM,y*ANGSTROM,z*ANGSTROM]);
  }else{
    return Err((CluEError::CannotParseLine(line.to_string()),0));
  }

  let serial: u32;
  if let Ok(ser) = serial_result{
    serial = ser;
  }else{
    return Err((CluEError::CannotParseLine(line.to_string()),0));
  }

  let element: Element = match Element::from(line[76..=77].trim()) {
    Ok(elmt) => elmt,
    Err(err) => return Err((err,serial))
  };

  Ok(Particle{
      element,
      isotope: Isotope::most_common_for(&element),
      coordinates,
      active: true,
      serial: Some(serial),
      name,
      residue,
      residue_sequence_number: residue_sequence_number.ok(),
      chain_id,
      //state: Vec::<Complex64>::new(),
      })
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
  let line = conect_line[6..].trim_end();

  let num = line.len();
  let n = 5;
  let mut serials = Vec::<u32>::with_capacity(num/n);
  let mut idx = 0;
  while idx < num {

    let idx_end = usize::min(idx+n,num);
    
    let connect: u32;
    
    if let Ok(serial) = line[idx .. idx_end].trim().parse::<u32>(){
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
     if line.contains("MODEL"){
       model_count += 1;
     }else if line.contains("ENDMDL"){
       increment = 0;
     }else if line.contains("ATOM") || line.contains("HETATM") {
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
    let r = detected_particle.weighted_coordinates.mean();
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
