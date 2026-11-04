use crate::clue_errors::CluEError;
use crate::physical_constants::NANOMETER;
use crate::space_3d::Vector3D;
use crate::structure::Structure;

use std::fs::File;
use std::io::prelude::*;


// GRO atom and residue numbers have five digits and wrap around.
const GRO_NUMBER_MODULUS: usize = 100_000;
// Characters in an atom line without velocities, including the newline.
const GRO_CHARS_PER_LINE: usize = 45;

impl Structure{
  //----------------------------------------------------------------------------
  // TODO: pub fn read_gro(file_name: &str) -> Result<Structure,CluEError>{}
  //----------------------------------------------------------------------------
  /// This function writes the primary cell of a `Structure` in GRO format.
  /// The detected spin, if there is one, is written first as atom 0 of
  /// residue "DET" at the weighted mean of its coordinates.
  /// All bath particles of the primary cell follow, whether active or not.
  /// Periodic copies are not written, since the box line describes the
  /// primary cell. Velocities are not written.
  pub fn write_gro(&self,filename: &str) -> Result<(),CluEError>{

    let gro = self.gro_string()?;

    let mut file = match File::create(filename){
      Ok(file) => file,
      Err(err) => return Err(CluEError::CannotCreateFile(
          filename.to_string(), err.to_string())),
    };

    if let Err(err) = file.write_all(gro.as_bytes()){
      return Err(CluEError::CannotWriteFile(filename.to_string(),
            err.to_string()));
    }

    Ok(())
  }
  //----------------------------------------------------------------------------
  // This function returns the contents of the GRO file written by
  // `write_gro()`.
  fn gro_string(&self) -> Result<String,CluEError>{

    // A particle without a cell ID is treated as being in the primary cell.
    let primary_indices = (0..self.number())
      .filter(|&idx| self.cell_ids.get(idx).map_or(true, |&id| id == 0))
      .collect::<Vec::<usize>>();

    let detected_line = self.get_detected_particle_gro_line()?;

    let n_atoms = primary_indices.len() + detected_line.iter().count();

    let mut gro = String::with_capacity((n_atoms + 3)*GRO_CHARS_PER_LINE);

    // Write title.
    gro.push_str("CluE Oxide System\n");

    // Write number of atoms.
    gro.push_str(&format!("{: >5}\n",n_atoms));

    if let Some(line) = detected_line{
      gro.push_str(&line);
    }

    for (ii,&idx) in primary_indices.iter().enumerate(){
      let particle = &self.bath_particles[idx];

      let res_seq = particle.residue_sequence_number.unwrap_or(0) as usize;

      let res_name = particle.residue.clone()
        .unwrap_or_else(|| "UNK".to_string());

      let name = particle.name.clone()
        .unwrap_or_else(|| particle.element.to_string());

      let r = &particle.coordinates - &self.pdb_origin;

      gro.push_str(&gro_atom_line(res_seq, &res_name, &name, ii + 1, &r));
    }

    // Write PBC.
    gro.push_str(&self.get_gro_box_line()?);

    Ok(gro)
  }
  //----------------------------------------------------------------------------
  // This function returns the GRO line of the detected spin,
  // or `None` if there is no detected spin.
  fn get_detected_particle_gro_line(&self)
    -> Result<Option<String>,CluEError>
  {
    let Some(detected_particle) = &self.detected_particle else{
      return Ok(None);
    };

    let Some(r) = detected_particle.weighted_coordinates.weighted_mean() else{
      return Err(CluEError::InvalidSpinGridWeights(
            detected_particle.weighted_coordinates.weights.iter()
            .fold(0.0, |sum, w| sum + w)));
    };
    let r = &Vector3D::from([r[0],r[1],r[2]]) - &self.pdb_origin;

    Ok(Some(gro_atom_line(0, "DET", &detected_particle.isotope.to_string(),
            0, &r)))
  }
  //----------------------------------------------------------------------------
  // This function returns the GRO box line of the unit cell, in nm:
  // v1(x) v2(y) v3(z) for a rectangular cell, otherwise
  // v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y).
  // A structure without a unit cell gets a zero box.
  fn get_gro_box_line(&self) -> Result<String,CluEError>{

    let (a,b,c) = match self.cell_edges.as_slice(){
      [] => return Ok(format!("{:10.5}{:10.5}{:10.5}\n",0.0,0.0,0.0)),
      [a,b,c] => (a,b,c),
      edges => return Err(CluEError::WrongNumberOfCellEdges(edges.len())),
    };

    // Values that are written as zero.
    let is_zero = |x: f64| (x/NANOMETER).abs() < 5e-6;

    // GRO requires v1(y) = v1(z) = v2(z) = 0.
    if !(is_zero(a.y()) && is_zero(a.z()) && is_zero(b.z())){
      return Err(CluEError::InvalidGroBox(format!("{}, {}, {}",
              a.scale(1.0/NANOMETER), b.scale(1.0/NANOMETER),
              c.scale(1.0/NANOMETER))));
    }

    let mut line = format!("{:10.5}{:10.5}{:10.5}",
        a.x()/NANOMETER, b.y()/NANOMETER, c.z()/NANOMETER);

    if !(is_zero(b.x()) && is_zero(c.x()) && is_zero(c.y())){
      for x in [0.0, 0.0, b.x(), 0.0, c.x(), c.y()]{
        line.push_str(&format!("{:10.5}",x/NANOMETER));
      }
    }
    line.push('\n');

    Ok(line)
  }
//------------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/*
residue number (5 positions, integer)

residue name (5 characters)

atom name (5 characters)

atom number (5 positions, integer)

position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)

velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions
with 4 decimal places
*/
// This function returns a GRO atom line without velocities.
// Longer names are cut to five characters, and numbers wrap at 100000.
fn gro_atom_line(res_seq: usize, res_name: &str, name: &str, atom_num: usize,
    r: &Vector3D) -> String
{
  format!("{: >5}{: <5.5}{: >5.5}{: >5}{:8.3}{:8.3}{:8.3}\n",
      res_seq % GRO_NUMBER_MODULUS, res_name,
      name, atom_num % GRO_NUMBER_MODULUS,
      r.x()/NANOMETER, r.y()/NANOMETER, r.z()/NANOMETER)
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod test{
  use super::*;
  use crate::cluster::adjacency::AdjacencyList;
  use crate::config::Config;
  use crate::elements::Element;
  use crate::io::FromTOMLString;
  use crate::physical_constants::ANGSTROM;
  use crate::structure::particle::Particle;
  use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

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
  //----------------------------------------------------------------------------
  // This function checks the layout of GRO file contents and returns
  // the atom lines and the box values.
  fn check_gro_layout(gro: &str) -> (Vec::<String>, Vec::<f64>){
    let lines = gro.lines().map(|s| s.to_string()).collect::<Vec::<String>>();
    assert!(gro.ends_with('\n'));
    assert_eq!(lines[0],"CluE Oxide System");

    let n_atoms = lines[1].trim().parse::<usize>().unwrap();
    assert_eq!(lines.len(), n_atoms + 3);

    let atom_lines = lines[2..2+n_atoms].to_vec();
    for line in atom_lines.iter(){
      // Fixed columns with three decimal places.
      assert_eq!(line.len(), 44, "{}", line);
      for &col in [24, 32, 40].iter(){
        assert_eq!(&line[col..col+1], ".", "{}", line);
      }
    }

    let box_line = &lines[n_atoms + 2];
    assert_eq!(box_line.len() % 10, 0);
    let box_values = (0..box_line.len()/10)
      .map(|ii| box_line[10*ii..10*(ii+1)].trim().parse::<f64>().unwrap())
      .collect::<Vec::<f64>>();

    (atom_lines, box_values)
  }
  //----------------------------------------------------------------------------
  // This function reads the coordinates in meters from a GRO atom line.
  fn gro_coordinates(line: &str) -> [f64;3]{
    let x = |ii: usize| line[20 + 8*ii..28 + 8*ii].trim()
      .parse::<f64>().unwrap()*NANOMETER;
    [x(0), x(1), x(2)]
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_write_gro_primary_cell(){
    for replicate in ["false","true"]{
      let mut config = Config::from_toml_string(&format!(r##"
        input_structure_file = "assets/TEMPO.pdb"
        radius = 80
        detected_spin.position = [28,29]
        number_timepoints = [101]
        tau_increments = [1e-1]
        cluster_method = "cce"
        max_cluster_size = 2
        magnetic_field = 1.2
        replicate_unit_cell = {}

        [[groups]]
        name = "tempo_h"
        selection.elements = ["H"]
        "##, replicate)).unwrap();
      config.set_defaults().unwrap();

      let mut rng = ChaCha20Rng::seed_from_u64(0);
      let structure = Structure::build_structure(&mut rng,&config).unwrap();
      if replicate == "true"{
        assert!(structure.number() > 29);
      }

      let filename = std::env::temp_dir().join(format!(
            "clue_oxide_test_write_gro_{}_{}.gro", replicate,
            std::process::id()));
      let filename = filename.to_str().unwrap().to_string();
      structure.write_gro(&filename).unwrap();
      let gro = std::fs::read_to_string(&filename).unwrap();
      std::fs::remove_file(&filename).unwrap();

      let (atom_lines, box_values) = check_gro_layout(&gro);

      // The detected spin and the 29 atoms of TEMPO.pdb.
      assert_eq!(atom_lines.len(), 30);
      assert_eq!(&atom_lines[0][0..20],"    0DET      e    0");

      let r_det = structure.detected_particle.as_ref().unwrap()
        .weighted_coordinates.weighted_mean().unwrap();
      let r_det = &Vector3D::from([r_det[0],r_det[1],r_det[2]])
        - &structure.pdb_origin;
      let r_det = [r_det.x(), r_det.y(), r_det.z()];
      let r = gro_coordinates(&atom_lines[0]);
      for ii in 0..3{
        assert!((r[ii] - r_det[ii]).abs() <= 5e-4*NANOMETER);
      }

      for (ii,line) in atom_lines[1..].iter().enumerate(){
        let particle = &structure.bath_particles[ii];
        assert_eq!(line[15..20].trim().parse::<usize>().unwrap(), ii + 1);
        assert_eq!(line[10..15].trim(), particle.name.as_ref().unwrap());

        let r0 = &particle.coordinates - &structure.pdb_origin;
        let r0 = [r0.x(), r0.y(), r0.z()];
        let r = gro_coordinates(line);
        for jj in 0..3{
          assert!((r[jj] - r0[jj]).abs() <= 5e-4*NANOMETER);
        }
      }

      // CRYST1  72.5676 72.5676 72.5676  90.00  90.00  90.00
      assert_eq!(box_values.len(), 3);
      for x in box_values{
        assert!((x*NANOMETER - 72.5676*ANGSTROM).abs() <= 5e-6*NANOMETER);
      }
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_gro_string_from_new_structure(){
    let particles = vec![
      Particle::new(Element::Hydrogen, 0.0, 0.0, 0.0),
      Particle::new(Element::Hydrogen, 1.5*NANOMETER, -2.0*NANOMETER,
          1e-4*NANOMETER),
    ];
    let structure = Structure::new(particles, AdjacencyList::with_capacity(2),
        Vec::<Vector3D>::new());

    let gro = structure.gro_string().unwrap();
    let (atom_lines, box_values) = check_gro_layout(&gro);

    assert_eq!(atom_lines, vec![
        "    0UNK      H    1   0.000   0.000   0.000".to_string(),
        "    0UNK      H    2   1.500  -2.000   0.000".to_string(),
    ]);
    assert_eq!(box_values, vec![0.0, 0.0, 0.0]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_gro_box_line(){
    let nm = NANOMETER;
    let triclinic = vec![
      Vector3D::from([1.0*nm, 0.0, 0.0]),
      Vector3D::from([0.5*nm, 0.8*nm, 0.0]),
      Vector3D::from([0.1*nm, 0.2*nm, 0.9*nm]),
    ];
    let structure = Structure::new(Vec::<Particle>::new(),
        AdjacencyList::with_capacity(0), triclinic);
    assert_eq!(structure.get_gro_box_line().unwrap(),
        "   1.00000   0.80000   0.90000   0.00000   0.00000   0.50000   \
0.00000   0.10000   0.20000\n");

    // The first edge must be along x and the second in the xy plane.
    for (a,b) in [
      ([1.0*nm, 0.1*nm, 0.0], [0.0, 1.0*nm, 0.0]),
      ([1.0*nm, 0.0, 0.0], [0.0, 1.0*nm, 0.1*nm]),
    ]{
      let edges = vec![Vector3D::from(a), Vector3D::from(b),
        Vector3D::from([0.0, 0.0, 1.0*nm])];
      let structure = Structure::new(Vec::<Particle>::new(),
          AdjacencyList::with_capacity(0), edges);
      assert!(matches!(structure.get_gro_box_line(),
            Err(CluEError::InvalidGroBox(_))));
    }

    let structure = Structure::new(Vec::<Particle>::new(),
        AdjacencyList::with_capacity(0),
        vec![Vector3D::from([1.0*nm, 0.0, 0.0])]);
    assert!(matches!(structure.get_gro_box_line(),
          Err(CluEError::WrongNumberOfCellEdges(1))));
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
