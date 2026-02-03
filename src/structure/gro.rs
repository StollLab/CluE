use crate::clue_errors::CluEError;
use crate::physical_constants::NANOMETER;
use crate::structure::{Structure, particle::Particle};

use std::io::BufWriter;
use substring::Substring;
use std::fs::File;
use std::io::prelude::*;


const GRO_PBC_INDENT: &str  = "    ";

impl Structure{
  //----------------------------------------------------------------------------
  // TODO: pub fn read_gro(file_name: &str) -> Result<Structure,CluEError>{}
  //----------------------------------------------------------------------------
  /// This function writes a `Structure` in GRO format.
  pub fn write_gro(&self,filename: &str) -> Result<(),CluEError>{

    let n_active = self.number_active();

    let chars_per_line = 68;
    let bytes_per_char = 32;
    let n_bytes = (n_active +1)*chars_per_line*bytes_per_char;

    let Ok(file) = File::create(filename)
    else {
      return Err(CluEError::CannotOpenFile(filename.to_string()) );
    };
    let mut stream = BufWriter::with_capacity(n_bytes,file);

    let mut atom_num = 1;
    let mut line: String;

    // Write title.
    line = "CluE Oxide System\n".to_string();
    let stream_result = stream.write(line.as_bytes());
    if stream_result.is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }

    // Write number of atoms.
    line = format!("{: >5}\n",self.primary_cell_indices.len());
    let stream_result = stream.write(line.as_bytes());
    if stream_result.is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }

    for (idx,particle) in self.bath_particles.iter().enumerate(){
      if self.cell_id(idx) != Ok(0){
        break;
      }
      self.set_bath_particle_gro_string(&mut line, particle, &mut atom_num)?;

      
      let stream_result = stream.write(line.as_bytes());
      if stream_result.is_err(){
        return Err(CluEError::CannotWriteFile(filename.to_string()) );
      }
      

    }

    // Write PBC.
    if self.cell_offsets.len() ==3{
      line = format!("{}{}{}{}{}{}{}{}{}{}",GRO_PBC_INDENT,
          self.cell_offsets[0].x(),
          self.cell_offsets[1].y(),
          self.cell_offsets[2].z(),
          self.cell_offsets[0].y(),
          self.cell_offsets[0].z(),
          self.cell_offsets[1].x(),
          self.cell_offsets[2].z(),
          self.cell_offsets[2].x(),
          self.cell_offsets[2].y(),
      );
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
/*
residue number (5 positions, integer)

residue name (5 characters)

atom name (5 characters)

atom number (5 positions, integer)

position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)

velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions 
with 4 decimal places   
*/   
fn set_bath_particle_gro_string(&self, line: &mut String, particle: &Particle, 
    atom_num: &mut usize) 
  -> Result<(),CluEError>
{

  // atom number (5 positions, integer)
  let atom_num_str = format!("{: >5}", atom_num.to_string());
  
  // atom name (5 characters)
  let name = if let Some(part_name) = &particle.name{
    format!("{: >5}",part_name)
  }else{
    format!("{: >5}",particle.element.to_string())
  };
  
  // residue name (5 characters)
  let res_name = if let Some(res) = &particle.residue{
    format!("{: >5}",res)
  }else{
    "  UNK".to_string()
  };

  // residue number (5 positions, integer)
  let res_seq = if let Some(res_seq_num) = particle.residue_sequence_number{
    format!("{: >5}",res_seq_num)
  }else{
    "    0".to_string()
  };
 
  //position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
  let r = &particle.coordinates - &self.pdb_origin;
 
  let x = format!("{: >8}",
      ((1000.0*r.x()/NANOMETER).round()/1000.0).to_string().substring(0,7));
  let y = format!("{: >8}",
      ((1000.0*r.y()/NANOMETER).round()/1000.0).to_string().substring(0,7));
  let z = format!("{: >8}",
      ((1000.0*r.z()/NANOMETER).round()/1000.0).to_string().substring(0,7));

  // velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions 
  // with 4 decimal places   
  let v = "  0.0000".to_string();

  *line = format!("{}{}{}{}{}{}{}{}{}{}\n"
      ,//GRO_ATOM_INDENT, 
      res_seq,res_name, 
      name, atom_num_str, 
      x, y, z,
      v, v, v,
      );

  *atom_num += 1;


  Ok(())
}  
//------------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod test{

  /*
  fn get_gro_string() -> String{
    // https://manual.gromacs.org/2025.1/reference-manual/file-formats.html#gro
    "\
MD of 2 waters, t= 0.0\n\
    6\n\
    1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434\n\
    1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791\n\
    1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180\n\
    2WATER  OW1    4   1.275   0.053   0.622  0.2519  0.3140 -0.1734\n\
    2WATER  HW2    5   1.337   0.002   0.680 -1.0641 -1.1349  0.0257\n\
    2WATER  HW3    6   1.326   0.120   0.568  1.9427 -0.8216 -0.0244\n\
   1.82060   1.82060   1.82060".to_string()
  }
  */

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
