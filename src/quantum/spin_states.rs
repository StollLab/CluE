use crate::clue_errors::CluEError;
use crate::config::{ClusterPopulations,Config};

use crate::physical_constants::PI;
use crate::HamiltonianTensors;
use crate::physical_constants::{ONE,I,ZERO};

use num_complex::Complex64;
use ndarray::{Array1,Array2};
use ndarray::linalg::kron;
use ndarray_linalg::Trace;

use rand_chacha::ChaCha20Rng;
use rand_distr::uniform::SampleRange;
use rand::Rng;
use rand::distributions::Uniform;
use rand_distr::Distribution;
use ndarray_linalg::Norm;

type CxVec = Array1::<Complex64>;
type CxMat = Array2::<Complex64>;

pub struct SpinStates{
  pub states: Vec::<CxMat>,
  pub pure_states: bool,          
}

impl SpinStates{

  //----------------------------------------------------------------------------
  pub fn new() -> Self{
    Self{
      states: Vec::<CxMat>::new(),
      pure_states: true,
    }
  }
  //----------------------------------------------------------------------------
  pub fn generate(
      rng: &mut ChaCha20Rng,
      tensors: &HamiltonianTensors,
      config: &Config) 
      -> Result<Self,CluEError>
  {
    let Some(populations) =  &config.cluster_populations else{
      return Err(CluEError::NoClusterDensityMatrixMethod);
    };

    let Some(ensemble_cce) = config.ensemble_cce else{
      return Err(CluEError::NoEnsembleCCE);
    };

    let states = match populations{
      ClusterPopulations::Uniform => {
        Self::uniform(&tensors.spin_multiplicities,!ensemble_cce)
      },
      ClusterPopulations::Random => {
        let mut s = Self::zeros(&tensors.spin_multiplicities,!ensemble_cce);
        if let Some(temperature) = config.temperature{
          let weights = tensors.get_zeeman_boltzmann_weights(temperature);
          s.randomize_state_weighted(rng,&weights);
        }else{
          s.randomize_state(rng);
        }
        s
      },
      ClusterPopulations::Thermal => {
        let Some(temperature) = config.temperature else{
          return Err(CluEError::NoTemperature);
        };
        let weights = tensors.get_zeeman_boltzmann_weights(temperature);
        Self::from_weights(weights,!ensemble_cce, Some(rng)) 
      },
      ClusterPopulations::Zeeman => {
        if let Some(temperature) = config.temperature{
          let weights = tensors.get_zeeman_boltzmann_weights(temperature);
          Self::from_weights(weights,!ensemble_cce, None) 
        } else{
          Self::uniform(&tensors.spin_multiplicities,!ensemble_cce)
        }
      }
    };

    Ok(states)
  }
  //----------------------------------------------------------------------------
  pub fn density_matrix_for(&self, indices: &[usize]) 
      -> Result<CxMat,CluEError>
  {
    if self.pure_states{
      self.pure_state_density_matrix_for(indices)
    }else{
      self.incoherent_density_matrix_for(indices)
    }
  }
  //----------------------------------------------------------------------------
  /// The function returns
  /// |ψ_idx1> ... |ψ_idxn><ψ_idxn| ... <ψ_idx|
  /// for ind in `indices`.
  pub fn pure_state_density_matrix_for(&self, indices: &[usize]) 
      -> Result<CxMat,CluEError>
  {
    let mut psi = CxMat::ones([1,1]);
    for &idx in indices.iter(){
      psi = kron(&psi, &self.states[idx])
    }
 
    let psi_t = psi.t().map(|psi_i| psi_i.conj() );
    let mut rho = psi.dot(&psi_t);
    let Ok(z) = rho.trace() else{
      return Err(CluEError::CannotTakeTrace(format!("{}",rho)));
    };
    rho = (ONE/z)*rho;

    Ok(rho)
  }
  //----------------------------------------------------------------------------
  /// The function returns
  /// diag(diag(|ψ_idx1> ... |ψ_idxn><ψ_idxn| ... <ψ_idx|))
  /// for ind in `indices`.
  pub fn incoherent_density_matrix_for(&self, indices: &[usize]) 
      -> Result<CxMat,CluEError>
  {
    let mut psi = CxMat::ones([1,1]);
    for &idx in indices.iter(){
      psi = kron(&psi, &self.states[idx])
    }
 
    let mut rho = CxMat::from_diag(&CxVec::from_vec(
        psi.map(|psi_i| psi_i.norm_sqr()*ONE)
        .into_raw_vec()));

    let Ok(z) = rho.trace() else{
      return Err(CluEError::CannotTakeTrace(format!("{}",rho)));
    };
    rho = (ONE/z)*rho;

    Ok(rho)
  }
  //----------------------------------------------------------------------------
  pub fn zeros(spin_multiplicities: &[usize],pure_states:bool) -> Self{
    let n = spin_multiplicities.len();
    let mut states = Vec::<CxMat>::with_capacity(n);

    for &s in spin_multiplicities.iter(){
      states.push(CxMat::zeros([s,1]))
    }

    Self{ states,pure_states}
  }
  //----------------------------------------------------------------------------
  pub fn uniform(spin_multiplicities: &[usize], pure_states: bool) -> Self{
    let n = spin_multiplicities.len();
    let mut states = Vec::<CxMat>::with_capacity(n);

    for &s in spin_multiplicities.iter(){
      let u0 = ONE*( 1.0/(s as f64) ).sqrt();
      states.push(u0*CxMat::ones([s,1]))
    }

    Self{ states,pure_states}
  }
  //----------------------------------------------------------------------------
  pub fn from_weights(weights_list: Vec::<Vec::<f64>>, 
      pure_states: bool, mut randomize_phase: Option<&mut ChaCha20Rng>) 
    -> Self
  {
    let range = Uniform::new(0.0f64, 1.0);
    let n = weights_list.len();
    let mut states = Vec::<CxMat>::with_capacity(n);

    for (ii,weights) in weights_list.iter().enumerate(){

      let mult = weights.len();

      states.push(CxMat::zeros([mult,1]));
      for (jj,w) in weights.iter().enumerate(){
        states[ii][[jj,0]] = ONE*w;

        if let Some(rng) = &mut randomize_phase{
          let phi = 2.0*PI*range.sample(rng);
          states[ii][[jj,0]] *= Complex64{re: phi.cos(), im: phi.sin()};
        }

        let sum_sqr: f64 = states[ii]
            .iter().map(|z| z.norm_sqr()).sum();
        states[ii][[jj,0]] *= Complex64{re: 1.0/sum_sqr.sqrt(), im: 0.0};
      }
    }

    Self{ states,pure_states}
  }
  //----------------------------------------------------------------------------
  pub fn randomize_state_weighted(&mut self,rng: &mut ChaCha20Rng,
      weights: &[Vec::<f64>]){

    let range = Uniform::new(0.0f64, 1.0);
    for (ii,psi) in self.states.iter_mut().enumerate(){

      let mult = psi.dim().0;
      'rng_lp: loop{
        let mut psi2 = 0.0;

        for idx in 0..mult{ 
          let phi = 2.0*PI*range.sample(rng);
          let r = range.sample(rng);
          let z = Complex64{re: r*phi.cos(), im: r*phi.sin()};
          psi2 += z.norm_sqr();
          if psi2 > 1.0{ continue 'rng_lp };  
      
          let w = weights[ii][idx];
          psi[[idx,0]] = w*z; 
        }
        let sum_sqr: f64 = psi.iter().map(|z| z.norm_sqr()).sum();
        *psi *= Complex64{re: 1.0/sum_sqr.sqrt(), im: 0.0};

        break;
      }
    }

  }
  //----------------------------------------------------------------------------
  pub fn randomize_state(&mut self,rng: &mut ChaCha20Rng){
    let range = Uniform::new(0.0f64, 1.0);
    for psi in self.states.iter_mut(){
      let mult = psi.dim().0;
      'rng_lp: loop{
        let mut psi2 = 0.0;
        for idx in 0..mult{ 
          let phi = 2.0*PI*range.sample(rng);
          let r = range.sample(rng);
          let z = Complex64{re: r*phi.cos(), im: r*phi.sin()};
          psi2 += z.norm_sqr();
          if psi2 > 1.0{ continue 'rng_lp };  
      
          psi[[idx,0]] = z; 
        }

        *psi *= Complex64{re: 1.0/psi2.sqrt(), im: 0.0};

        break;
      }
    }
  }
  //----------------------------------------------------------------------------
}


#[cfg(test)]
mod tests{
  use super::*;
  use ndarray::array;
  use rand::SeedableRng;

  //----------------------------------------------------------------------------
  #[test]
  fn test_pure_state_density_matrix_for(){
    let states = SpinStates{
      states: vec![
          array![[ONE],[ZERO]],
          array![[ONE],[I]],
          array![[ZERO],[I]],
      ],
      pure_states: true,
    };

    let density_matrix = states.pure_state_density_matrix_for(&[0,2])
        .unwrap();
    let answer = array![
        [ZERO,ZERO,ZERO,ZERO],
        [ZERO,ONE,ZERO,ZERO],
        [ZERO,ZERO,ZERO,ZERO],
        [ZERO,ZERO,ZERO,ZERO],
    ];
    assert_eq!(density_matrix, answer);

    let states = SpinStates{
      states: vec![
          array![[ONE],[I]],
          array![[ONE],[I]],
      ],
      pure_states: true,
    };

    let density_matrix = states.pure_state_density_matrix_for(&[0])
        .unwrap();
    let answer = array![
        [ONE,-I],
        [I,ONE],
    ]/2.0;
    assert_eq!(density_matrix, answer);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_incoherent_density_matrix_for(){
    let states = SpinStates{
      states: vec![
          array![[ONE],[ZERO]],
          array![[ONE],[I]],
          array![[ZERO],[I]],
      ],
      pure_states: false,
    };

    let density_matrix = states.density_matrix_for(&[0,2]).unwrap();
    let answer = array![
        [ZERO,ZERO,ZERO,ZERO],
        [ZERO,ONE,ZERO,ZERO],
        [ZERO,ZERO,ZERO,ZERO],
        [ZERO,ZERO,ZERO,ZERO],
    ];
    assert_eq!(density_matrix, answer);

    let states = SpinStates{
      states: vec![
          array![[ONE],[I]],
          array![[ONE],[I]],
      ],
      pure_states: false,
    };

    let density_matrix = states.density_matrix_for(&[0]).unwrap();
    let answer = array![
        [ONE,ZERO],
        [ZERO,ONE],
    ]/2.0;
    assert_eq!(density_matrix, answer);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_SpinStates_zeros(){
    let s = SpinStates::zeros(&[2,3],true);
    assert_eq!(s.states[0].dim().0,2);
    assert_eq!(s.states[0].dim().1,1);
    assert_eq!(s.states[0], array![[ZERO],[ZERO]]);
    assert_eq!(s.states[1], array![[ZERO],[ZERO], [ZERO]]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_randomize_state_weighted(){
    let mut rng = ChaCha20Rng::from_entropy();
    let mut s = SpinStates::zeros(&[2,3], true);

    let mut x = 0.0;
    let mut y = 0.0;
    let n = 1000;
    let weights = vec![vec![0.8,0.2], vec![0.45, 0.55, 0.0]];
    for _ in 0..n{
      s.randomize_state_weighted(&mut rng, &weights);
      for psi in s.states.iter(){
        assert!( (psi.norm()-1.0).abs() < 1e-12);
      }
      x += s.states[0][[0,0]].norm_sqr();
      y += s.states[1][[0,0]].norm_sqr();
    }
    x /= n as f64;
    y /= n as f64;

    assert!( x > 0.75);
    assert!( x < 0.85);
    assert!( y > 0.4);
    assert!( y < 0.5);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_randomize_state(){
    let mut rng = ChaCha20Rng::from_entropy();
    let mut s = SpinStates::zeros(&[2,3],true);

    let mut x = 0.0;
    let mut y = 0.0;
    let n = 1000;
    for _ in 0..n{
      s.randomize_state(&mut rng);
      x += s.states[0][[0,0]].norm_sqr();
      y += s.states[1][[0,0]].norm_sqr();
      for psi in s.states.iter(){
        assert!( (psi.norm()-1.0).abs() < 1e-12);
      }

    }
    x /= n as f64;
    y /= n as f64;

  }
}
