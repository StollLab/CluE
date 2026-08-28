use crate::config::{
  ClusterPopulations,
  Config,
  DetectedPopulation,  
  pulse_sequence::PulseSequence,
};
use crate::clue_errors::CluEError;
use crate::signal::Signal;
use crate::HamiltonianTensors;
use crate::quantum::spin_hamiltonian::*;
use crate::quantum::gcce_hamiltonian::{
  build_spin_hamiltonian,
  get_electron_cluster_thermal_density_matrix,
  propagate_custom_pulse_sequence,
  propagate_pulse_sequence,
};
use crate::quantum::pulse_sequences::{
  get_standard_pulses,
  generate_pulse_sequence,  
};
use crate::quantum::cluster_operators::ClusterSpinOperators;
use crate::quantum::spin_states::SpinStates;

use ndarray::linalg::kron;
use ndarray_linalg::Trace;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/// This function calculates the cluster correlation expansion (CCE)
/// approximation of the echo decay.
/// The CCE method is developed  in
/// W. Yang and R.-B. Liu, “Quantum many-body theory of qubit decoherence
/// in a finite-size spin bath baths,” Phys. Rev. B 78, 085315 (2008).
/// and 
/// W. Yang and R.-B. Liu, “Quantum many-body theory of qubit decoherence
/// in a finite-size spin bath. II. Ensemble dynamics,” Phys. Rev. B 79, 115320
/// (2009).
//------------------------------------------------------------------------------
// This function calculate the cluster signal for the cluster specified by
// tensor_indices.
pub fn cce(tensor_indices: &[usize],states: &SpinStates, 
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors, 
    config: &Config) 
  -> Result<Option<Signal>,CluEError>
{

  if tensor_indices.is_empty(){
    let number_timepoints = &config.number_timepoints;
    if number_timepoints.is_empty(){
      return Err(CluEError::NoTimepoints);
    }
    let n_tot = config.get_total_number_timesteps();
    return Ok(Some(Signal::ones(n_tot)));
  }

  let hamiltonian = build_block_diag_hamiltonian(
      tensor_indices,spin_ops, tensors,config)?;

  let Some(cluster_populations) = &config.cluster_populations else{
    return Err(CluEError::NoClusterDensityMatrixMethod);
  };
  let density_matrix = match cluster_populations{
    ClusterPopulations::Thermal  => 
        get_cluster_thermal_density_matrix(&hamiltonian, config)?,
    _ => states.density_matrix_for(tensor_indices)?,
  };

  let signal = propagate_pulse_sequence_block_diag(
      &density_matrix, &hamiltonian, config)?;
  Ok(Some(signal))
}
//------------------------------------------------------------------------------
// This function calculate the cluster signal for the cluster specified by
// tensor_indices.
pub fn gcce(tensor_indices: &[usize], states: &SpinStates, 
    spin_ops: &ClusterSpinOperators, tensors: &HamiltonianTensors, 
    config: &Config) 
  -> Result<Option<Signal>,CluEError>
{

  let mut spin_indices = Vec::<usize>::with_capacity(1 + tensor_indices.len());
  spin_indices.push(0);
  for idx in tensor_indices.iter(){
    spin_indices.push(*idx);
  }

  let spin_multiplicities: Vec::<usize> =
      spin_indices.iter().map(|idx| tensors.spin_multiplicities[*idx])
      .collect();

  let spin_multiplicity = match spin_multiplicities.len(){
    0 => return Err(CluEError::NoDetectedSpinMultiplicity),
    1 => spin_multiplicities[0], // TODO: check if this is right.
    _ => spin_multiplicities[1],
  };

  let cluster_size = spin_multiplicities.len();

  let (h_eigvals, h_eigvecs) = build_spin_hamiltonian(
      &spin_indices,spin_ops,tensors)?;


  /*
  let mut density_matrix = match states.density_matrix_for(tensor_indices)?{
    None => {
      get_electron_cluster_thermal_density_matrix(
          &h_eigvals, &h_eigvecs, config)?},
    Some(rho) => {
      let detected_spin_density_matrix = spin_ops.get_density_matrix(
          spin_multiplicity, 1)?;  
      kron(&detected_spin_density_matrix,&rho)
    },     
  };
  */
  let Some(detected_population) = &config.detected_population else{
    return Err(CluEError::NoDetectedSpinDensityMatrix);
  };
  let Some(cluster_populations) = &config.cluster_populations else{
    return Err(CluEError::NoClusterDensityMatrixMethod);
  };
  let mut density_matrix = match (detected_population,cluster_populations){
    (DetectedPopulation::Thermal,ClusterPopulations::Thermal) => {
        get_electron_cluster_thermal_density_matrix(
            &h_eigvals, &h_eigvecs, config)?
    },
    (_, _) => {
          let rho = states.density_matrix_for(tensor_indices)?;
          let detected_spin_density_matrix = spin_ops.get_density_matrix(
              spin_multiplicity, 1)?;  
          kron(&detected_spin_density_matrix,&rho)
    },
  };

  /*
  let mut density_matrix = match states.density_matrix_for(tensor_indices)?{

    None => get_electron_cluster_thermal_density_matrix(
          &h_eigvals, &h_eigvecs, config)?,

    Some(rho) => {
      let detected_spin_density_matrix = match config.detected_population{
        Some(DetectedPopulation::Thermal) 
            => states.incoherent_density_matrix_for(&[0])?,
        Some(_) => spin_ops.get_density_matrix(spin_multiplicity, 1)?,
        None => return Err(CluEError::DetectedSpinDensityMatrix),
      };  
      kron(&detected_spin_density_matrix,&rho)
    },     

  };
  */

  /* Normalization is not needed here and will fail for ρ = Sz.
  let Ok(z) = density_matrix.trace() else{ 
    return Err(CluEError::CannotTakeTrace(format!("{}",density_matrix)));
  };
  density_matrix /= z;
  */

  let Some(pulse_sequence) = &config.pulse_sequence else{
    return Err(CluEError::NoPulseSequence);
  }; 

  let signal = match pulse_sequence{
    PulseSequence::Custom(pulse_sequence_specifier) =>{
      let pulse_sequence = generate_pulse_sequence(
          pulse_sequence_specifier, spin_ops,
          spin_multiplicity,cluster_size)?;

      propagate_custom_pulse_sequence(
          &pulse_sequence, &density_matrix,&h_eigvals, &h_eigvecs,config
      )?
    }  
    PulseSequence::CarrPurcell(_) | PulseSequence::RefocusedHahnEcho =>{ 
      let pulses = get_standard_pulses(
          spin_ops,spin_multiplicity,cluster_size)?;

      propagate_pulse_sequence(&pulses,&density_matrix,
          &h_eigvals, &h_eigvecs, config,
          )?
    },  
    PulseSequence::FreeEvolution => {
      let pulses = get_standard_pulses(
          spin_ops,spin_multiplicity,cluster_size)?;

      propagate_pulse_sequence(&pulses,&density_matrix,
          &h_eigvals, &h_eigvecs, config,
          )?
    },  
  };
  Ok(Some(signal))
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests{
  use super::*;

  use crate::quantum::tensors::*;
  use crate::space_3d::{SymmetricTensor3D,Vector3D};
  use crate::cluster_methods::appa;

  use rand_chacha::ChaCha20Rng;
  use rand::SeedableRng;
  //----------------------------------------------------------------------------
  #[test]
  fn test_cce(){


    let z0 = 33.0e9;
    let z1 = 80.0e6;
    let a1 = 10.0e6;
    let a2 = -10.0e6;
    let b = 10.0e3;
    let tensors = build_restricted_three_spin_tensors(z0, z1, a1, a2, b);
  
    let spin_indices = vec![1,2];

    let mut config = Config::new();
    config.set_defaults().unwrap();
    
    let spin_ops = ClusterSpinOperators::new(1,&vec![2],2,&config).unwrap();

    let mut config = Config::new();
    config.number_timepoints = vec![21];
    let delta_hf = a1 - a2;
    let freq = appa::appa_hahn_frequency(delta_hf,b);
    config.tau_increments = vec![0.05/freq];
    config.pulse_sequence = Some(PulseSequence::CarrPurcell(1));

    config.set_defaults().unwrap();
    config.set_tau_axis().unwrap();

    let mut rng = ChaCha20Rng::from_entropy();
    let states = SpinStates::generate(&mut rng, &tensors,&config).unwrap();
  
    let signal_opt  = cce(&vec![1,2], &states,&spin_ops, &tensors, 
        &config).unwrap();

    let Some(signal) = signal_opt else{
      panic!("Could not calculate signal.");
    }; 
    assert_eq!(signal.data.len(),21);


    let ref_signal_opt = appa::appa_hahn(
        &spin_indices,&tensors,&config).unwrap();
    let Some(ref_signal) = ref_signal_opt else{
      panic!("Could not calculate reference signal.");
    };

    for (ii,v) in signal.data.iter().enumerate(){
      let v0 = ref_signal.data[ii];
      let err: f64;
      if (v+v0).norm() < 1e12{
        err = (v-v0).norm();
      }else{
        err = (2.0+(v-v0)/(v+v0)).norm();
      }
      assert!(err < 1e-9);
    }
  }
  //----------------------------------------------------------------------------
  fn build_restricted_three_spin_tensors(z0: f64, z1: f64, a1: f64, a2: f64, 
      b: f64) -> HamiltonianTensors{
    let spin_multiplicities = vec![2,2,2,2,2];

    let mut spin1_tensors = Spin1Tensors::new(5);
    let zeeman0 = Vector3D::from([0.0, 0.0, z0]);
    let zeeman1 = Vector3D::from([0.0, 0.0, z1]);
    spin1_tensors.set(0,zeeman0);
    spin1_tensors.set(1,zeeman1.clone());
    spin1_tensors.set(2,zeeman1.clone());
    spin1_tensors.set(3,zeeman1.clone());
    spin1_tensors.set(4,zeeman1);

    let mut spin2_tensors = Spin2Tensors::new(5);
    let hf1 = SymmetricTensor3D::from([ 0.0, 0.0, 0.0,
                                             0.0, 0.0,
                                                    a1]);
    let hf2 = SymmetricTensor3D::from([ 0.0, 0.0, 0.0,
                                             0.0, 0.0,
                                                   a2]);

    let dip = SymmetricTensor3D::from([ -b/2.0,    0.0, 0.0,
                                                -b/2.0, 0.0,
                                                         b]);

    spin2_tensors.set(0,1,hf1.clone());
    spin2_tensors.set(0,2,hf2.clone());
    spin2_tensors.set(1,2,dip.clone());

    spin2_tensors.set(0,3,hf1);
    spin2_tensors.set(0,4,hf2);
    spin2_tensors.set(3,4,dip);

    let ge = -1.7609e11;
    HamiltonianTensors{
      spin_multiplicities,
      spin1_tensors,
      spin2_tensors,
      detected_gamma_matrix: SymmetricTensor3D::from([ge, 0.0, 0.0,
                                                           ge, 0.0,
                                                                ge]),
      magnetic_field: Vector3D::from([0.0,0.0,1.2]),
      mean_field_couplings: None,
      }
  }  
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
