
use crate::clue_errors::CluEError;
use crate::config::Config;
use crate::config::particle_config::{EigSpecifier,TensorSpecifier};
use crate::isotopes::Isotope;
use crate::physical_constants::{HBAR, JOULES_TO_HERTZ,
  MU0,MUB,MUN,PI};
use crate::space_3d::{SymmetricTensor3D,UnitSpherePoint,Vector3D};
use crate::structure::{DetectedSpin,Structure};
use crate::structure::particle::Particle;
use crate::structure::exchange_groups::ExchangeGroup;
use crate::symmetric_list_2d::SymList2D;


use rand_chacha::ChaCha20Rng;

//use std::fs;
//use std::fs::File;
//use std::io::{Error, Write};
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// `HamiltonianTensors` contians the information needed to write the spin
/// Hamiltonian a the sum of spin operators times coupling tensors.
/// `spin_multiplicities` lists the spin multiplicities, (2S+1), for each spin.
/// 'spin1_tensors' the specifies the terms that contain only a
/// single spin operator, such as the Zeeman coupling.
/// 'spin2_tensors' the specifies the terms that contain two spin operators
/// either between two different spins or from the same spin.
/// Higher order terms are not supported.
/// `detected_gamma_matrix` is the gamma matrix matrix for the detected spin.
/// `magnetic_field` is the applied magnetic field, assumed to be uniform.
#[derive(Debug,Clone)] 
pub struct HamiltonianTensors{
  pub spin_multiplicities: Vec::<usize>,
  pub spin1_tensors: Spin1Tensors, // O(S)
  pub spin2_tensors: Spin2Tensors, // O(S^2)
  pub detected_gamma_matrix: SymmetricTensor3D,
  pub magnetic_field: Vector3D,
}
impl HamiltonianTensors{
  //----------------------------------------------------------------------------
  /// This function returns the number of spins covered by the Hamiltonian.
  pub fn len(&self) -> usize{
    self.spin1_tensors.len()
  }
  //----------------------------------------------------------------------------
  /// This function return `true` iff there are no spins covered by the 
  /// Hamiltonian.
  pub fn is_empty(&self) -> bool{
    self.spin1_tensors.is_empty()
  }
  //----------------------------------------------------------------------------
  /// This function performs a passive rotation with respect to the system.
  pub fn rotate_pasive(&mut self, dir: &UnitSpherePoint){
    let gamma_matrix = self.detected_gamma_matrix.rotate_pasive(dir);
    self.spin1_tensors.set(0,
        construct_zeeman_tensor(&gamma_matrix,&self.magnetic_field));

    self.spin2_tensors.rotate_pasive(dir);
  }
  //----------------------------------------------------------------------------
  /// This function builds the spin Hamiltonian from the input structure
  /// and user specifications.
  pub fn generate(rng: &mut ChaCha20Rng,
      structure: &Structure, config: &Config) 
    -> Result<Self,CluEError>
  {

    let n_spins = structure.number_active() + 1;
    let mut spin_multiplicities = Vec::<usize>::with_capacity(n_spins);
    let mut spin1_tensors = Spin1Tensors::new(n_spins);
    let mut spin2_tensors = Spin2Tensors::new(n_spins);

    let Some(magnetic_field) = &config.magnetic_field else{
      return Err(CluEError::NoMagneticField);
    };

    let Some(detected_particle) = &structure.detected_particle else {
      return Err(CluEError::NoCentralSpin);
    };

    let gamma_matrix = detected_particle.gyromagnetic_ratio_matrix();

    spin_multiplicities.push(detected_particle.spin_multiplicity());

    let det_zeeman = construct_zeeman_tensor(gamma_matrix,magnetic_field);
    if det_zeeman.any_nan(){
      return Err(CluEError::NANTensorDetectedZeeman);
    }
    spin1_tensors.set(0, det_zeeman);

    let eye = SymmetricTensor3D::eye();

    let mut tensor_indices = Vec::<Option<usize>>::with_capacity(
        structure.bath_particles.len());

    let mut idx0 = 0;
    for (particle_idx0, particle0) in structure.bath_particles.iter()
      .enumerate(){

      if !particle0.active {
        tensor_indices.push(None);
        continue;
      }
      idx0 += 1;
      tensor_indices.push(Some(idx0));

      let gamma0 = particle0.isotope.gyromagnetic_ratio();

      spin_multiplicities.push(particle0.isotope.spin_multiplicity());


      // nuclear Zeeman
      let bath_zeeman = construct_bath_zeeman_tensor(rng,particle_idx0,
          magnetic_field, structure, config)?;
      //let bath_zeeman = construct_zeeman_tensor(&(gamma0*&zz),magnetic_field);

      if bath_zeeman.any_nan(){
        return Err(CluEError::NANTensorBathZeeman(
              particle_idx0,particle0.isotope.to_string()));
      }
      spin1_tensors.set(idx0, bath_zeeman);


      // hyperfine
      let hf_ten = construct_hyperfine_tensor(rng,detected_particle, particle0, 
          particle_idx0, structure, config)?; 

      spin2_tensors.set(0,idx0, hf_ten);

      
      // electric quadrupole coupling
      let quadrupole_opt = construct_electric_quadrupole_tensor(rng,particle0, 
          particle_idx0, structure, config)?;

      if let Some(quadrupole_ten) = quadrupole_opt{
        spin2_tensors.set(idx0,idx0, quadrupole_ten);
      }


      // nucleus-nucleus dipole-dipole
      let mut idx1 = 0;
      for (particle_idx1,particle1) in structure.bath_particles.iter()
        .enumerate()
      {
        if !particle1.active {continue;}
         idx1 += 1;
         if idx1 == idx0 {break;}

         let gamma1 = particle1.isotope.gyromagnetic_ratio();

         let delta_r01 = &particle0.coordinates - &particle1.coordinates;

         let dd_ten = construct_point_dipole_dipole_tensor(gamma1,gamma0,
              &delta_r01);

         if dd_ten.any_nan(){
           return Err(CluEError::NANTensorBathDipoleDipole(
                 particle_idx0,particle0.isotope.to_string(),
                 particle_idx1,particle1.isotope.to_string(),
               ));
         }

         spin2_tensors.set(idx1,idx0, dd_ten);
      }
    }


    if let Some(exchange_group_manager) = &structure.exchange_groups{
      'exchange_loop: for (ex_id, exchange_group) in exchange_group_manager
        .exchange_groups.iter().enumerate()
      {
        // Add tensors bassed on the identity of the exchange group.
        match exchange_group{
          ExchangeGroup::Methyl(c3rotor) | 
            ExchangeGroup::PrimaryAmonium(c3rotor) =>{
              let j = exchange_group_manager.exchange_couplings[ex_id];

              let j_tensor = eye.scale(j);

              // Loop through all hydrogens.
              for &h0 in c3rotor.indices.iter(){
                if !structure.bath_particles[h0].active{
                  continue 'exchange_loop;
                }
              }

              for &h0 in c3rotor.indices.iter(){
                //let h0_indices = &structure.cell_indices[h0];
                
                // Loop through all hydrogens with indices less than the first.
                for &h1 in c3rotor.indices.iter(){
                  if h1 == h0 { break; }

                  let Some(h0_idx) = tensor_indices[h0] else{
                    return Err(CluEError::TensorNotSet(h0));
                  };
                  let Some(h1_idx) = tensor_indices[h1] else{
                    return Err(CluEError::TensorNotSet(h1));
                  };

                  if j_tensor.any_nan(){
                    return Err(CluEError::NANTensorExchangeCoupling(
                      h0, structure.bath_particles[h0].isotope.to_string(),
                      h1, structure.bath_particles[h1].isotope.to_string(),
                      ));
                  }

                  spin2_tensors.add(h0_idx,h1_idx, j_tensor.clone());

                }
              }

            },
        }

       
      }
    }


    Ok(HamiltonianTensors{
      spin_multiplicities,
      spin1_tensors,
      spin2_tensors,
      detected_gamma_matrix: gamma_matrix.clone(),
      magnetic_field: magnetic_field.clone(),
      })

  }
  //----------------------------------------------------------------------------
  /// This function writes the coupling tensors to a text file.
  pub fn save(&self, filename: &str, structure: &Structure) 
    -> Result<(),CluEError>
  {
    let Ok(file) = File::create(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()) );
    };

    let n_char_spin_mult = 30;
    let n_char_f64 = 16;
    let n_char_s1 = 30 + 3*n_char_f64;
    let n_char_s2 = 35 + 6*n_char_f64;


    let n_spins = self.spin1_tensors.len();

    let bytes_per_char = 32;

    let n_bytes = 3200 
      + bytes_per_char*(n_spins+1)*(n_char_spin_mult + n_char_s1 + n_char_s2);

    let mut stream = BufWriter::with_capacity(n_bytes,file);
    
    let line = format!("#[tensors(number = {})]\n",n_spins);
    if stream.write(line.as_bytes()).is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }

    // Write spin multiplicities.
    let line = "\n// Spin Multiplicities\n".to_string();
    if stream.write(line.as_bytes()).is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }

    for (ii, spin_mult) in self.spin_multiplicities.iter().enumerate(){

      let idx = structure.get_reference_index_of_nth_active(ii)?;

      let line = format!("spin_mutiplicity[{}] = {};\n", idx, spin_mult);

      if stream.write(line.as_bytes()).is_err(){
        return Err(CluEError::CannotWriteFile(filename.to_string()) );
      }
    }

    // Write O(S^1) tensors.
    let line = "\n// One Spin Tensors\n".to_string();
    if stream.write(line.as_bytes()).is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }

    for ii in 0..n_spins{
      if let Some(tensor) = self.spin1_tensors.get(ii){
        let line = format!("tensor[{}] = {}; // Hz.\n",
            structure.get_reference_index_of_nth_active(ii)?,
            tensor);
       
        if stream.write(line.as_bytes()).is_err(){
          return Err(CluEError::CannotWriteFile(filename.to_string()) );
        } 
      }
    }

    // Write O(S^2) tensors.
    let line = "\n// Two Spin Tensors\n".to_string();
    if stream.write(line.as_bytes()).is_err(){
      return Err(CluEError::CannotWriteFile(filename.to_string()) );
    }

    for ii in 0..n_spins{
      for jj in ii..n_spins{

        if let Some(tensor) = self.spin2_tensors.get(ii,jj){
          let line = format!("tensor[{},{}] = {}; // Hz.\n",
              structure.get_reference_index_of_nth_active(ii)?,
              structure.get_reference_index_of_nth_active(jj)?, 
              tensor);
         
          if stream.write(line.as_bytes()).is_err(){
            return Err(CluEError::CannotWriteFile(filename.to_string()) );
          } 
        }
      }
    }

    Ok(())
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)] 
pub struct Spin1Tensors{
  tensors: Vec::< Option<Vector3D> >,
}

impl<'a> Spin1Tensors{

  pub fn new(number: usize) -> Spin1Tensors {
    
    let mut tensors 
     = Vec::<Option<Vector3D>>::with_capacity(number);
    for _ii in 0..number {
      tensors.push(None);
    }

    Spin1Tensors{
      tensors,
    }
  }
  //----------------------------------------------------------------------------
  pub fn set(&mut self, n: usize, ten: Vector3D) {
    self.tensors[n] = Some(ten);
  }
  //----------------------------------------------------------------------------
  pub fn get(&'a self, n: usize) -> Option< &'a Vector3D> {
    match &self.tensors[n] {
      Some(ten) => Some(ten),
      None => None,  
    }
  }
  //----------------------------------------------------------------------------
  pub fn remove(&mut self, n: usize){
    self.tensors[n] = None;
  }
  //----------------------------------------------------------------------------
  pub fn len(&self) -> usize{
    self.tensors.len()
  }
  //----------------------------------------------------------------------------
  pub fn is_empty(&self) -> bool{
    self.tensors.is_empty()
  }
  //----------------------------------------------------------------------------
  /*
  pub fn find_nans(&self) -> Vec::<usize>{

    let mut nan_indices = Vec::<usize>::new();
    for (idx, ten_opt) in self.tensors.iter().enumerate(){
      if let Some(ten) = &ten_opt{
        if ten.any_nan(){
          nan_indices.push(idx);
        }
      }
    }

    nan_indices
  }
  */
  //----------------------------------------------------------------------------
  /*
  pub fn rotate_active(&mut self, dir: &UnitSpherePoint){
    for tensor in self.tensors.iter_mut().flatten(){
      *tensor = tensor.rotate_active(dir);
    }
  }
  */
  //----------------------------------------------------------------------------

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug,Clone)] 
pub struct Spin2Tensors{
  tensors: SymList2D::<SymmetricTensor3D>,
}

impl<'a> Spin2Tensors{

  pub fn new(number: usize) -> Self {
    Spin2Tensors{
      tensors: SymList2D::new(number),
    }
  }
  //----------------------------------------------------------------------------
  pub fn get(&'a self, m: usize, n: usize) -> Option<&'a SymmetricTensor3D >{
    self.tensors.get(m,n)
  }
  //----------------------------------------------------------------------------
  pub fn set(&mut self, m: usize, n: usize, ten: SymmetricTensor3D){
    self.tensors.set(m,n,ten);
  }
  //----------------------------------------------------------------------------
  pub fn add(&mut self, m: usize, n: usize, mut ten: SymmetricTensor3D){
    if let Some(ten0) = self.tensors.get(m,n){
      ten = &ten + ten0;
    }
    self.tensors.set(m,n,ten);
  }
  //----------------------------------------------------------------------------
  pub fn remove(&mut self, m: usize, n: usize){
    self.tensors.remove(m,n);
  }
  //----------------------------------------------------------------------------
  /*
  pub fn rotate_active(&mut self, dir: &UnitSpherePoint){

    let dim = self.tensors.dim();
    for m in 0..dim{
      for n in 0..=m{
        if let Some(ten) = self.tensors.get(m,n){
          let rot_ten = ten.rotate_active(dir);
          self.tensors.set(m,n, rot_ten );
        } 
      }
    }
  }
  */
  //----------------------------------------------------------------------------
  pub fn rotate_pasive(&mut self, dir: &UnitSpherePoint){

    let dim = self.tensors.dim();
    for m in 0..dim{
      for n in 0..=m{
        if let Some(ten) = self.tensors.get(m,n){
          let rot_ten = ten.rotate_pasive(dir);
          self.tensors.set(m,n, rot_ten );
        } 
      }
    }
  }
  //----------------------------------------------------------------------------
  /*
  pub fn find_nans(&self) -> Vec::<(usize,usize)>{

    let mut nan_indices = Vec::<(usize,usize)>::new();
    let dim = self.tensors.dim();
    for m in 0..dim{
      for n in 0..=m{
        if let Some(ten) = self.tensors.get(m,n){
          if ten.any_nan(){
            nan_indices.push((m,n));
          }
        } 
      }
    }

    nan_indices
  }
  */
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
pub fn construct_zeeman_tensor(
    gyromagnetic_ratio: &SymmetricTensor3D, 
    magnetic_field: &Vector3D) -> Vector3D
{

  (-0.5/PI)*&(magnetic_field * gyromagnetic_ratio)
}
//------------------------------------------------------------------------------
fn construct_bath_zeeman_tensor(rng: &mut ChaCha20Rng,
    particle_index: usize, magnetic_field: &Vector3D,
    structure: &Structure, config: &Config)
  -> Result<Vector3D, CluEError>
{
  let particle = &structure.bath_particles[particle_index];

  let gamma_matrix = if let Some(tensor_specifier) = structure
    .extract_g_matrix_specifier(particle_index,config)
  {
    let g_matrix = construct_symmetric_tensor_from_tensor_specifier(rng,
        tensor_specifier, Some(particle_index),structure, config)?;
     
    let mu: f64 = if particle.isotope == Isotope::Electron{
      -MUB
    }else{
      MUN
    };

    g_matrix.scale(mu/HBAR)

  }else{
    let gamma0 = particle.isotope.gyromagnetic_ratio();

    let zz = SymmetricTensor3D::from([0.0,0.0,0.0,
                                          0.0,0.0,
                                              1.0]);
    
    gamma0*&zz
  };

  Ok(construct_zeeman_tensor(&gamma_matrix,magnetic_field))
}
//------------------------------------------------------------------------------
fn construct_electric_quadrupole_tensor(rng: &mut ChaCha20Rng,
    particle0: &Particle,particle_index: usize,
    structure: &Structure, config: &Config)
  -> Result<Option<SymmetricTensor3D>, CluEError>
{
  if particle0.isotope.spin_multiplicity() < 3 {
    return Ok(None);
  }

  match structure.extract_electric_quadrupole_specifier(particle_index,config)
  {
    Some(tensor_specifier) => {
      let tensor = construct_symmetric_tensor_from_tensor_specifier(rng,
        tensor_specifier, Some(particle_index),structure, config)?;

      if tensor.any_nan(){
        return Err(CluEError::NANTensorQuadrupole(
          particle_index,particle0.isotope.to_string()));
      }
      Ok(Some(tensor))
    },
    None => Ok(None),
  }
}
//------------------------------------------------------------------------------
fn construct_hyperfine_tensor(rng: &mut ChaCha20Rng,
    detected_particle: &DetectedSpin, 
    particle0: &Particle, particle_index: usize, 
    structure: &Structure, config: &Config) 
  -> Result<SymmetricTensor3D, CluEError>
{

  let tensor: SymmetricTensor3D;
  if let Some(tensor_specifier) = structure
    .extract_hyperfine_specifier(particle_index,config)
  {
    tensor = construct_symmetric_tensor_from_tensor_specifier(rng,
        tensor_specifier, Some(particle_index),structure, config)?;

  }else{

    let gamma_e = detected_particle.isotope.gyromagnetic_ratio();

    let gamma0 = particle0.isotope.gyromagnetic_ratio();
    
    let mut hf_ten = SymmetricTensor3D::zeros();
    for ii in 0..detected_particle.weighted_coordinates.len(){
      let Ok(r) = detected_particle.weighted_coordinates.xyz(ii) else{
        return Err(CluEError::NoCentralSpinCoor);
      }; 
      let delta_r = &r - &particle0.coordinates;
    
      let w = detected_particle.weighted_coordinates.weight(ii);
    
      let ten = construct_point_dipole_dipole_tensor(gamma_e,gamma0,&delta_r);
      hf_ten = &hf_ten + &(w*&ten);
    }
    tensor = hf_ten;
        
  }

  if tensor.any_nan(){
    return Err(CluEError::NANTensorHyperfine(
          particle_index,particle0.isotope.to_string()));
  }
  Ok(tensor)
}
//------------------------------------------------------------------------------
pub fn get_perpendicular_dipole_dipole_frequency(
    gyromagnetic_ratio_1: f64,
    gyromagnetic_ratio_2: f64,
    r: f64
    ) -> f64 {
  JOULES_TO_HERTZ*
  MU0/(4.0*PI)*gyromagnetic_ratio_1*gyromagnetic_ratio_2   
  *HBAR*HBAR/r/r/r
}
//------------------------------------------------------------------------------
pub fn construct_point_dipole_dipole_tensor(
    gyromagnetic_ratio_1: f64,
    gyromagnetic_ratio_2: f64,
    delta_r: &Vector3D
    ) -> SymmetricTensor3D {

  let r = delta_r.norm();
  
  let n = (1.0/r) * delta_r;
 
  let n3nt = 3.0 * &n.self_outer(); 

  let h_perp = get_perpendicular_dipole_dipole_frequency(gyromagnetic_ratio_1,
      gyromagnetic_ratio_2,r);

  h_perp*&(&SymmetricTensor3D::eye() - &n3nt)

}
//------------------------------------------------------------------------------
/// This function takes three values and three vectors and constructs
/// ___T___ = ___XE___trans(___X___), where ___E___ is the digonal matrix of
/// the values and ___X___ is a 3×3 matrix with each vector as a column.
pub fn construct_symmetric_tensor_from_values_and_vectors(
    evals: &[f64;3], evecs: &[Vector3D;3]) -> SymmetricTensor3D
{

  let mut ten = SymmetricTensor3D::zeros();
  for (ii,v) in evecs.iter().enumerate(){
    ten = &ten + &(v.self_outer()).scale(evals[ii]);
  }

  ten

}
//------------------------------------------------------------------------------
/// This function takes `&TensorSpecifier` and generates a `SymmetricTensor3D`.
/// Depending on the `TensorSpecifier`, additional information may be 
/// required, which is supplied by the other arguments.
pub fn construct_symmetric_tensor_from_tensor_specifier(rng: &mut ChaCha20Rng,
    tensor_specifier: &TensorSpecifier, particle_index_opt: Option<usize>,
    structure: &Structure, config: &Config) 
  -> Result<SymmetricTensor3D, CluEError>
{
  match tensor_specifier{
    TensorSpecifier::Unspecified => Err(CluEError::NoTensorSpecifier),
    TensorSpecifier::Eig(eig_specifier) 
        => construct_symmetric_tensor_from_eig_specifier(rng,eig_specifier,
            particle_index_opt, structure, config),
    TensorSpecifier::SymmetricTensor3D(tensor) => Ok(tensor.clone()),
  }
}
//------------------------------------------------------------------------------
fn construct_symmetric_tensor_from_eig_specifier(rng: &mut ChaCha20Rng,
    tensor_specifier: &EigSpecifier, particle_index_opt: Option<usize>,
    structure: &Structure, config: &Config) 
  -> Result<SymmetricTensor3D, CluEError>
{
  let Some(values) = tensor_specifier.values else{
    return Err(CluEError::NoTensorValues);
  }; 

  let mut axes = Vec::<Vector3D>::with_capacity(2);
  let mut axis_dims = Vec::<usize>::with_capacity(2);

  const X: usize = 0;
  const Y: usize = 1;
  const Z: usize = 2;
  if let Some(axis_specifier) = &tensor_specifier.z_axis{
    let mut axis = axis_specifier.to_vector3d(rng,particle_index_opt,structure,
        config)?;
    axis = axis.normalize();
    axes.push(axis);
    axis_dims.push(Z);
  }

  if let Some(axis_specifier) = &tensor_specifier.x_axis{
    let mut axis = axis_specifier.to_vector3d(rng,particle_index_opt,structure,
        config)?;

    if axes.len() == 1{
      axis = &axis - &axes[0].scale(axes[0].dot(&axis));
    }
    axis = axis.normalize();
    axes.push(axis);
    axis_dims.push(X);
  }

  if axes.len() < 2 {
    if let Some(axis_specifier) = &tensor_specifier.y_axis{
      let mut axis = axis_specifier.to_vector3d(rng,
          particle_index_opt, structure, 
          config)?;
      if axes.len() == 1{
        axis = &axis - &axes[0].scale(axes[0].dot(&axis));
      }
      axis = axis.normalize();
      axes.push(axis);
      axis_dims.push(Y);
    }
  }


  if axis_dims.len() != 2{
    return Err(CluEError::WrongNumberOfAxes(axes.len(),2))
  }
  let axis_dims = [axis_dims[0],axis_dims[1]];

  let all_3_axes: [Vector3D; 3]; 
  match axis_dims{
    [Z,X] => {
       let mut axis = axes[0].cross(&axes[1]);
       axis = axis.scale(1.0/axis.norm());
       all_3_axes = [axes[1].clone(),axis,axes[0].clone()];
    },
    [Z,Y] => {
       let mut axis = axes[1].cross(&axes[0]);
       axis = axis.scale(1.0/axis.norm());
       all_3_axes = [axis,axes[1].clone(),axes[0].clone()];
    },
    [X,Y] => {
       let mut axis = axes[0].cross(&axes[1]);
       axis = axis.scale(1.0/axis.norm());
       all_3_axes = [axes[0].clone(),axes[1].clone(),axis];
    },
    _ => return Err(CluEError::InvalidAxes),
  }

  let ten = construct_symmetric_tensor_from_values_and_vectors(
      &values, &all_3_axes);

  Ok(ten)

}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;
  use crate::physical_constants::*;
  use crate::elements::Element;
  use crate::space_3d::{SymmetricTensor3D, Vector3D};
  use crate::config::particle_config::ParticleConfig;
  use crate::structure::pdb;
  use crate::Config;
  use crate::config::DetectedSpinCoordinates;
  use crate::structure::particle_filter::{ParticleFilter, 
    SecondaryParticleFilter,VectorSpecifier};
  //use crate::config::lexer::get_tokens_from_line;
  use crate::io::FromTOMLString;

  use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

  //----------------------------------------------------------------------------
  #[allow(non_snake_case)]
  #[test]
  fn test_HamiltonianTensors_generate(){

    let mut config = Config::from_toml_string(r##"
        input_structure_file = "assets/TEMPO.pdb"
        radius = 18 
        number_timepoints = [101]
        tau_increments = [1e-1] 
        cluster_method = "cce"
        max_cluster_size = 2
        magnetic_field = 1.2
        replicate_unit_cell = true

        [detected_spin]
        position = [28,29]

        [[groups]]
        name = "tempo_h"
        selection.residues = ["TEM"]
        selection.elements = ["H"]

        1H.c3_tunnel_splitting = 80e-3

        [[groups]]
        name = "tempo_o"
        selection.residues = ["TEM"]
        selection.elements = ["O"]

        [[groups]]
        name = "tempo_c"
        selection.residues = ["TEM"]
        selection.elements = ["C"]

        [[groups]]
        name = "tempo_n"
        selection.residues = ["TEM"]
        selection.elements = ["N"]

        14N.electric_quadrupole.values = [-1.12, -5.88, 7.0]
        14N.electric_quadrupole.axes.x = {to_same_molecule_as = "tempo_o"}
        14N.electric_quadrupole.axes.y.from_bonded_to = "tempo_c"
        14N.electric_quadrupole.axes.y.to_bonded_to = "tempo_c"
        

        14N.hyperfine.values = [16.086, 16.086, 103.356]
        14N.hyperfine.axes.x = [-1.15, -0.47,  0.71]
        14N.hyperfine.axes.y = [-1.09,  2.23, -0.49]

        # unrealistic, but useful for testing
        14N.g_matrix.values = [1.0, 1.0, 1.0]
        14N.g_matrix.axes.x = [-1.15, -0.47,  0.71]
        14N.g_matrix.axes.y = [-1.09,  2.23, -0.49]


        "##).unwrap();
    config.set_defaults().unwrap();

    let mut rng = ChaCha20Rng::from_entropy();

    let structure = Structure::build_structure(&mut rng,&config).unwrap();

    assert_eq!(structure.number_active(),19);

    let tensors = HamiltonianTensors::generate(&mut rng, &structure, &config)
        .unwrap();

    assert_eq!(tensors.spin_multiplicities.len(),20);


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Test nitrogen Zeeman

    let tensor = tensors.spin1_tensors.get(19).unwrap().clone();

    let ref_value = -3693247.037734574/0.403761;
    assert!( (tensor.z() - ref_value ).abs() < 1e-12);
    assert!( (tensor.x()/ref_value).abs() < 1e-12);
    assert!( (tensor.y()/ref_value).abs() < 1e-12);
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Test nitrogen electric quadrupole coupling.
    let tensor = tensors.spin2_tensors.get(19,19).unwrap().clone();
    let values = [-1120000.0, -5880000.0, 7000000.0];

    let r_n = Vector3D::from([36.440e-10, 36.900e-10,  37.100e-10]);
    let r_o = Vector3D::from([35.290e-10,  36.430e-10,  37.810e-10]);
    let r_c1 = Vector3D::from([37.700e-10, 36.150e-10, 37.340e-10]);
    let r_c19 = Vector3D::from([36.610e-10, 38.380e-10, 36.850e-10]);

    let r_no = &r_o - &r_n;
    let r_cc = &r_c19 - &r_c1;

    let x = r_no.normalize();
    let y = (&r_cc - &x.scale(x.dot(&r_cc))).normalize();

    assert_eq!(&tensor*&r_no,r_no.scale(values[0]));
    assert!((&(&tensor*&x) 
          - &x.scale(values[0])).norm()/values[0].abs() < 1e-9 );

    assert!((&(&tensor*&y) 
          - &y.scale(values[1])).norm()/values[1].abs() < 1e-9 );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Test nitrogen hyperfine coupling.
    let values = [16086000.0, 16086000.0, 103356000.0];
    let tensor = tensors.spin2_tensors.get(0,19).unwrap().clone();

    assert_eq!(&tensor*&r_no,r_no.scale(values[0]));
    assert!((&(&tensor*&x) 
          - &x.scale(values[0])).norm()/values[0].abs() < 1e-9 );

    assert!((&(&tensor*&y) 
          - &y.scale(values[1])).norm()/values[1].abs() < 1e-9 );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Test methyl exchange coupling.
    let methyl_tensor_indices = [[1,2,3], [4,5,6], [13,14,15], [16,17,18]];
    let j = -2.0*80e3/3.0;
    let tol = 1e-12;
    for methyl in methyl_tensor_indices.iter(){
      for &h in methyl.iter(){
        for idx in 0..20{
          if h == idx {
            assert!(tensors.spin2_tensors.get(h,idx).is_none());
            continue;
          }
          let ten = tensors.spin2_tensors.get(h,idx).unwrap();

          if methyl.contains(&idx){
            assert!( (ten.trace()/(3.0*j) - 1.0).abs() < tol)
          }else{
            assert!((ten.trace()/ten.xx()).abs() < tol)
          }
        }
      }
    }


  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_construct_symmetric_tensor_from_tensor_specifier(){
    let mut rng = ChaCha20Rng::from_entropy();
    let filename = "./assets/TEMPO.pdb";
    let mut structure = pdb::parse_pdb(&filename,0).unwrap();
    let mut config = Config::new();
    config.detected_spin_position = Some(
        DetectedSpinCoordinates::CentroidOverSerials(vec![28,29]) );
    config.set_defaults().unwrap();
    structure.build_primary_structure(&mut rng, &config).unwrap();

    config.particles.push( ParticleConfig::new("nitrogen".to_string()) );
    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Nitrogen];
    config.particles[0].filter = Some(filter);

    let label = String::from("oxygen");
    config.particles.push( ParticleConfig::new(label.clone()) );
    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Oxygen];
    config.particles[1].filter = Some(filter);

    let label = String::from("carbon");
    config.particles.push( ParticleConfig::new(label.clone()) );
    let mut filter = ParticleFilter::new();
    filter.elements = vec![Element::Carbon];
    config.particles[2].filter = Some(filter);

    let vector_specifier_no = VectorSpecifier::Diff(
        SecondaryParticleFilter::Particle,"nitrogen".to_string(),
        SecondaryParticleFilter::Bonded,"oxygen".to_string());

    let vector_specifier_cc = VectorSpecifier::Diff(
        SecondaryParticleFilter::Bonded,"carbon".to_string(),
        SecondaryParticleFilter::Bonded,"carbon".to_string());

    let r_n = Vector3D::from([36.440e-10, 36.900e-10,  37.100e-10]);
    let r_o = Vector3D::from([35.290e-10,  36.430e-10,  37.810e-10]);
    let r_c1 = Vector3D::from([37.700e-10, 36.150e-10, 37.340e-10]);
    let r_c19 = Vector3D::from([36.610e-10, 38.380e-10, 36.850e-10]);

    let delta_r_no = &r_o - &r_n;
    let delta_r_cc = &r_c19 - &r_c1;

    let particle_index = 27;
    assert_eq!(structure.bath_particles[particle_index].element,
        Element::Nitrogen);
    let r_no = vector_specifier_no.to_vector3d(&mut rng,
        Some(particle_index),&structure,
        &config).unwrap();
    assert_eq!(r_no,delta_r_no);

    let r_cc = vector_specifier_cc.to_vector3d(&mut rng,
        Some(particle_index),&structure,
        &config).unwrap();
    assert_eq!(r_cc,delta_r_cc);


    // de Oliveira, M.; Knitsch, R.; Sajid, M.; Stute, A.; Elmer, L.-M.; 
    // Kehr, G.; Erker, G.; Magon, C. J.; Jeschke, G.; Eckert, H.
    // Aminoxyl Radicals of B/P Frustrated Lewis Pairs:
    // Refinement of the Spin-Hamiltonian Parameters by Field- and
    // Temperature-Dependent Pulsed EPR Spectroscopy.
    // PLoS ONE 2016, 11 (6), e0157944.
    // https://doi.org/10.1371/journal.pone.0157944.
    let e2qqh = 3.5*1e6; // Hz
    let eta = 0.68;
    let values = [e2qqh*(eta-1.0), e2qqh*(-eta-1.0), 2.0*e2qqh];
    assert_eq!(values,[-1119999.9999999998, -5880000.000000001, 7000000.0]);

    let tensor_specifier = TensorSpecifier::Eig(EigSpecifier{
      values: Some(values.clone()),
      x_axis: Some(vector_specifier_no),
      y_axis: Some(vector_specifier_cc),
      z_axis: None,
    });

    let tensor = construct_symmetric_tensor_from_tensor_specifier(&mut rng,
        &tensor_specifier, Some(particle_index),&structure, &config).unwrap();

    let x = r_no.normalize();
    let y = (&r_cc - &x.scale(x.dot(&r_cc))).normalize();

    assert_eq!(&tensor*&r_no,r_no.scale(values[0]));
    assert!((&(&tensor*&x) 
          - &x.scale(values[0])).norm()/values[0].abs() < 1e-9 );

    assert!((&(&tensor*&y) 
          - &y.scale(values[1])).norm()/values[1].abs() < 1e-9 );
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_construct_symmetric_tensor_from_values_and_vectors(){
    let x = Vector3D::from([1.0, 0.0, 0.0]);
    let y = Vector3D::from([0.0, 1.0, 0.0]);
    let z = Vector3D::from([0.0, 0.0, 1.0]);

    let xyz = [x,y,z];
    let evals = [-1.0,2.0,10.0];

    let ten = construct_symmetric_tensor_from_values_and_vectors(
        &evals,&xyz);

    let xyz_ten = SymmetricTensor3D::from([evals[0],0.0,0.0,
                                                    evals[1], 0.0,
                                                              evals[2]]);

    assert_eq!(ten.xx(),xyz_ten.xx());
    assert_eq!(ten.xy(),xyz_ten.xy());
    assert_eq!(ten.xz(),xyz_ten.xz());
    assert_eq!(ten.yy(),xyz_ten.yy());
    assert_eq!(ten.yz(),xyz_ten.yz());
    assert_eq!(ten.zz(),xyz_ten.zz());

    let v1 = Vector3D::from([1.0,1.0,1.0]).scale(1.0/3.0f64.sqrt());
    let v2 = Vector3D::from([-1.0,-1.0,2.0]).scale(1.0/6.0f64.sqrt());
    let v3 = Vector3D::from([1.0,-1.0,0.0]).scale(1.0/2.0f64.sqrt());

    assert!( (v1.norm() - 1.0) < 1e-12);
    assert!( (v2.norm() - 1.0) < 1e-12);
    assert!( (v3.norm() - 1.0) < 1e-12);
    assert!( v1.dot(&v2) < 1e-12);
    assert!( v1.dot(&v3) < 1e-12);
    assert!( v2.dot(&v3) < 1e-12);
    
    let evecs = [v1,v2,v3];


    let ten = construct_symmetric_tensor_from_values_and_vectors(
        &evals,&evecs);

    for (ii,v) in evecs.iter().enumerate(){
      assert_eq!(&ten*v , v.scale(evals[ii]));
    }
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_construct_zeeman_tensor() {
    let gmr = -ELECTRON_G*MUB/HBAR;
    assert!((gmr+1.76085963023e11).abs() < 0.00000000053e11);
    let gyromagnetic_ratio = SymmetricTensor3D::from([gmr,0.0, 0.0,
                                                          gmr, 0.0,
                                                               gmr]); 
    let magnetic_field = Vector3D::from([0.0, 0.0, 1.2]);

    let zeeman = construct_zeeman_tensor(&gyromagnetic_ratio,&magnetic_field);
    assert_eq!(zeeman.x(),0.0);
    assert_eq!(zeeman.y(),0.0);
    assert!( (zeeman.z() - 33629941709.0).abs() < 0.1);

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_get_perpendicular_dipole_dipole_frequency(){
    let gmr = Isotope::Hydrogen1.gyromagnetic_ratio();
    let gmr0 = 2.6752218744e8;
    assert!( (gmr-gmr0).abs() < 0.000_000_0011e8);
    let r = 1.5e-10;
    let freq = get_perpendicular_dipole_dipole_frequency(gmr,gmr,r);
    
    let ref_freq = 35591.15890;
    assert!((freq-ref_freq).abs()/ref_freq<1e-9);
  
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_construct_point_dipole_dipole_tensor(){
  
    let gmr = Isotope::Hydrogen1.gyromagnetic_ratio();
    let r = 1.5e-10;
    let freq = get_perpendicular_dipole_dipole_frequency(gmr,gmr,r);

    let phi_list = (0..=10).map(|x| (x as f64)/10.0*PI*2.0)
      .collect::<Vec::<f64>>();

    let theta_list = (0..=10).map(|x| (2.0*(x as f64)/10.0 - 1.0f64).acos())
      .collect::<Vec::<f64>>();

    for &theta in theta_list.iter(){
      let ct = theta.cos();
      let st = theta.sin();

      for &phi in phi_list.iter(){
        let cp = phi.cos();
        let sp = phi.sin();


        let delta_r = r*&Vector3D::from([st*cp,st*sp,ct]);

        assert!((delta_r.norm()-r).abs()/r<1e-12);

        let ten = construct_point_dipole_dipole_tensor(gmr,gmr,&delta_r);
        
        let nnt = SymmetricTensor3D::from(
            [cp*cp*st*st, cp*sp*st*st, cp*st*ct,
                          sp*sp*st*st, sp*st*ct,
                                          ct*ct]);
        let ref_ten = freq*&(&SymmetricTensor3D::eye() - &(3.0*&nnt));

        let err_ten = &ten - &ref_ten; 

        let tol = freq*1e-9;
        assert!(err_ten.xx().abs() < tol);
        assert!(err_ten.xy().abs() < tol);
        assert!(err_ten.xz().abs() < tol);
        assert!(err_ten.yy().abs() < tol);
        assert!(err_ten.yz().abs() < tol);
        assert!(err_ten.zz().abs() < tol);


      }
    }
  }

  //----------------------------------------------------------------------------
  #[test]
  #[allow(non_snake_case)]
  fn test_Spin1Tensors() {
    let number = 3;
    let mut sf_tens = Spin1Tensors::new(number);

    for ii in 0..number {

      let mut ten = Vector3D::zeros();
      let val = 1.0 + (ii as f64);
      ten.set_x(val);

      sf_tens.set(ii,ten);

    }

    let expected_values = vec![ 1.0, 2.0, 3.0];

    for ii in 0..(number as usize) {

      let ten = sf_tens.get(ii).unwrap();
      let val = ten.x();
      assert!( (val-expected_values[ii]).abs() < 1e-12 ); 

    }


  }
  //----------------------------------------------------------------------------
  #[test]
  #[allow(non_snake_case)]
  fn test_Spin2Tensors() {
    let number = 3;
    let mut ss_tens = Spin2Tensors::new(number);

    let mut val = 0.0;
    for ii in 0..(number as usize) {
      for jj in 0..(number as usize) {

        let ten = SymmetricTensor3D::eye();
        val += 1.0;
        let ten = ten.scale(val);

        ss_tens.set(ii,jj,ten);
      }
    }
    let expected_values = vec![ 
      1.0, 4.0, 7.0, 
      4.0, 5.0, 8.0,
      7.0, 8.0, 9.0];

    let mut idx = 0;
    for ii in 0..number {
      for jj in 0..number {
        let ten = ss_tens.get(ii,jj).unwrap();
        let val = ten.xx();
        assert!( (val-expected_values[idx]).abs() < 1e-12 ); 
        idx += 1;
      }
    }
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
