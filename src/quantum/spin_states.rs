use crate::clue_errors::CluEError;
use crate::config::{ClusterMethod,Config};

use crate::HamiltonianTensors;
use crate::ndarray_linalg_eigh_workaround::eigh_hermitian;
use crate::physical_constants::{BOLTZMANN,JOULES_TO_HERTZ,ONE};
use crate::quantum::general_spin_hamiltonian::build_single_spin_hamiltonian;

use num_complex::Complex64;
use ndarray::{Array1,Array2};
use ndarray::linalg::kron;
use ndarray_linalg::{Trace,UPLO};

use rand_chacha::ChaCha20Rng;
use rand::distributions::Uniform;
use rand_distr::Distribution;

type CxMat = Array2::<Complex64>;

/// `SpinStates` holds one pure state per spin, as a column vector in the
/// Zeeman basis (index i is m = S - i). `states[0]` belongs to the detected
/// spin; it is the first basis vector and is not used (P4 in
/// plan_populations.md).
pub struct SpinStates{
  pub states: Vec::<CxMat>,
}

impl SpinStates{

  //----------------------------------------------------------------------------
  pub fn new() -> Self{
    Self{
      states: Vec::<CxMat>::new(),
    }
  }
  //----------------------------------------------------------------------------
  /// This function draws the initial state of every bath spin n ≥ 1: an
  /// eigenstate |k⟩ of its interaction-free Hamiltonian
  /// H_n = Σ_k E_k |k⟩⟨k| (`build_single_spin_hamiltonian`), chosen with
  /// probability p_k = exp(-βE_k)/Z_n, with β = h/(k_B T), or with
  /// p_k = 1/d_n (T = ∞) if `conditions.temperature` is not set.
  /// H_n has no mean field and no coupling to other spins. It includes the
  /// spherical-tensor terms only for gCCE, the method that propagates them.
  ///
  /// How the cluster initial state is built from these states (`cce`, `gcce`):
  /// | method | ensemble | cluster initial state                          |
  /// | CCE    | false    | ⊗ |k_n⟩⟨k_n| (this function)                   |
  /// | CCE    | true     | exp(-βH_C⁰)/Z, or 1/d without T                |
  /// | gCCE   | false    | ρ_S ⊗ (⊗ |k_n⟩⟨k_n|)                           |
  /// | gCCE   | true     | ρ_S ⊗ exp(-βH_C⁰)/Z (or 1/d without T) if
  ///                       detected_spin.density_matrix is set,
  ///                       else exp(-βH_{S+C})/Z                          |
  ///
  /// One uniform random number is drawn per bath spin, whatever the method
  /// and the ensemble setting, so that the random number stream does not
  /// depend on them.
  /// The average of ⊗|k_n⟩⟨k_n| over the draws is ⊗ exp(-βH_n)/Z_n, so
  /// averaging over samples gives the signal of the product thermal state.
  pub fn generate(
      rng: &mut ChaCha20Rng,
      tensors: &HamiltonianTensors,
      config: &Config)
      -> Result<Self,CluEError>
  {
    let include_spherical = config.cluster_method == Some(ClusterMethod::GCCE);

    // β in 1/Hz; the energies are in Hz.
    let beta_opt = config.temperature
      .map(|temperature| 1.0/(JOULES_TO_HERTZ*BOLTZMANN*temperature));

    let range = Uniform::new(0.0f64, 1.0);

    let n_spins = tensors.spin_multiplicities.len();
    let mut states = Vec::<CxMat>::with_capacity(n_spins);

    for (index, &multiplicity) in tensors.spin_multiplicities.iter()
      .enumerate()
    {
      if index == 0{
        // The detected spin: not sampled (P4).
        let mut psi = CxMat::zeros((multiplicity,1));
        if multiplicity > 0{ psi[[0,0]] = ONE; }
        states.push(psi);
        continue;
      }

      let ham = build_single_spin_hamiltonian(tensors, index,
          include_spherical)?;

      let Ok((energies, eigvecs)) = eigh_hermitian(&ham, UPLO::Lower) else{
        return Err(CluEError::CannotDiagonalizeHamiltonian(
              crate::math::cx_matrix_summary(&ham)));
      };

      let probabilities = level_probabilities(&energies, beta_opt);

      let u = range.sample(rng);
      let k = choose_level(&probabilities, u);

      let mut psi = CxMat::zeros((multiplicity,1));
      psi.column_mut(0).assign(&eigvecs.column(k));
      states.push(psi);
    }

    Ok(Self{states})
  }
  //----------------------------------------------------------------------------
  /// The function returns the product state density matrix
  /// |ψ_idx1> ... |ψ_idxn><ψ_idxn| ... <ψ_idx1|
  /// for idx in `indices`, normalized to unit trace.
  pub fn density_matrix_for(&self, indices: &[usize])
      -> Result<CxMat,CluEError>
  {
    self.pure_state_density_matrix_for(indices)
  }
  //----------------------------------------------------------------------------
  /// The function returns
  /// |ψ_idx1> ... |ψ_idxn><ψ_idxn| ... <ψ_idx1|
  /// for idx in `indices`, normalized to unit trace.
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
      return Err(CluEError::CannotTakeTrace(format!("a {}x{} matrix",
              rho.dim().0, rho.dim().1)));
    };
    rho = (ONE/z)*rho;

    Ok(rho)
  }
  //----------------------------------------------------------------------------
}

impl Default for SpinStates {
    fn default() -> Self {
        Self::new()
    }
}
//------------------------------------------------------------------------------
// This function returns p_k = exp(-β(E_k - E_min))/Σ_j exp(-β(E_j - E_min)),
// or 1/d if `beta_opt` is `None` (T = ∞). The shift by E_min avoids overflow
// and does not change p_k.
fn level_probabilities(energies: &Array1::<f64>, beta_opt: Option<f64>)
  -> Vec::<f64>
{
  let d = energies.len();
  let Some(beta) = beta_opt else{
    return vec![1.0/(d as f64); d];
  };
  let e_min = energies.iter().cloned().fold(f64::INFINITY, f64::min);
  let weights: Vec::<f64> = energies.iter()
    .map(|e| (-beta*(e - e_min)).exp()).collect();
  let z: f64 = weights.iter().sum();
  weights.iter().map(|w| w/z).collect()
}
//------------------------------------------------------------------------------
// This function returns the first k with Σ_{j≤k} p_j > u, or the last index
// if rounding leaves the cumulative sum at or below u.
fn choose_level(probabilities: &[f64], u: f64) -> usize{
  let mut cumulative = 0.0;
  for (k, p) in probabilities.iter().enumerate(){
    cumulative += p;
    if cumulative > u{
      return k;
    }
  }
  probabilities.len().saturating_sub(1)
}
//------------------------------------------------------------------------------


#[cfg(test)]
mod tests{
  use super::*;
  use crate::physical_constants::{I,PLANCK,ZERO};
  use crate::quantum::tensors::*;
  use crate::space_3d::{SymmetricTensor3D,Vector3D};
  use ndarray::array;
  use ndarray_linalg::Norm;
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
    assert_eq!(states.density_matrix_for(&[0,2]).unwrap(), answer);

    let states = SpinStates{
      states: vec![
          array![[ONE],[I]],
          array![[ONE],[I]],
      ],
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
  fn test_choose_level(){
    let p = [0.2, 0.5, 0.3];
    assert_eq!(choose_level(&p, 0.0), 0);
    assert_eq!(choose_level(&p, 0.19), 0);
    assert_eq!(choose_level(&p, 0.2), 1);
    assert_eq!(choose_level(&p, 0.69), 1);
    assert_eq!(choose_level(&p, 0.71), 2);
    assert_eq!(choose_level(&p, 0.999_999), 2);
    // Rounding: the probabilities sum to less than u.
    assert_eq!(choose_level(&[0.5, 0.4999], 0.99995), 1);

    let p = level_probabilities(&array![0.0, 1.0, 2.0], None);
    assert_eq!(p, vec![1.0/3.0; 3]);
    // Without the shift by E_min, exp(-βE) would underflow to 0 here.
    // E = 2^20 + (0, 1, 2)·2^-10 and β = 2^10 are exact in floating point.
    let e0 = 2f64.powi(20);
    let de = 2f64.powi(-10);
    let p = level_probabilities(&array![e0, e0 + de, e0 + 2.0*de], 
        Some(2f64.powi(10)));
    let z = 1.0 + (-1.0f64).exp() + (-2.0f64).exp();
    assert!((p[0] - 1.0/z).abs() < 1e-9);
    assert!((p[2] - (-2.0f64).exp()/z).abs() < 1e-9);
    let p = level_probabilities(&array![5.0, 6.0, 7.0], Some(1.0));
    let z = 1.0 + (-1.0f64).exp() + (-2.0f64).exp();
    assert!((p[0] - 1.0/z).abs() < 1e-15);
    assert!((p[1] - (-1.0f64).exp()/z).abs() < 1e-15);
    assert!((p[2] - (-2.0f64).exp()/z).abs() < 1e-15);
  }
  //----------------------------------------------------------------------------
  // Tensors for the detected electron (index 0) and one bath spin (index 1)
  // with Zeeman vector `nu` and, if multiplicity 3, the axial tensor
  // Q = d(n n^T - 1/3), n at angle theta from z in the xz-plane.
  fn one_bath_spin_tensors(multiplicity: usize, nu: Vector3D, d: f64,
      theta: f64) -> HamiltonianTensors
  {
    let mut spin1_tensors = Spin1Tensors::new(2);
    spin1_tensors.set(0, Vector3D::from([0.0, 0.0, 33.0e9]));
    spin1_tensors.set(1, nu);
    let mut spin2_tensors = Spin2Tensors::new(2);
    if multiplicity == 3{
      let (st, ct) = (theta.sin(), theta.cos());
      spin2_tensors.set(1, 1, SymmetricTensor3D::from([
          d*(st*st - 1.0/3.0), 0.0, d*st*ct,
                          -d/3.0, 0.0,
                                  d*(ct*ct - 1.0/3.0)]));
    }
    HamiltonianTensors{
      spin_multiplicities: vec![2, multiplicity],
      // No bath gamma matrices: every bath Zeeman vector is fixed.
      bath_gamma_matrices: vec![None; spin1_tensors.len()],
      spin1_tensors,
      spin2_tensors,
      spin_spherical_tensors: SpinSphericalTensors::new(),
      max_spherical_tensor_rank: 0,
      detected_gamma_matrix: SymmetricTensor3D::eye(),
      magnetic_field: Vector3D::from([0.0,0.0,1.2]),
      mean_field_couplings: None,
    }
  }
  //----------------------------------------------------------------------------
  fn config_with(temperature: Option<f64>) -> Config{
    let mut config = Config::new();
    config.cluster_method = Some(ClusterMethod::CCE);
    config.temperature = temperature;
    config
  }
  //----------------------------------------------------------------------------
  // exp(a) for a matrix a, from its Taylor series (scaled and squared);
  // independent of eigh.
  fn expm(a: &CxMat) -> CxMat{
    let n = a.dim().0;
    let squarings = 10;
    let a = a/(2f64.powi(squarings));
    let mut term = CxMat::eye(n);
    let mut sum = CxMat::eye(n);
    for k in 1..30{
      term = term.dot(&a)/(k as f64);
      sum += &term;
    }
    for _ in 0..squarings{ sum = sum.dot(&sum); }
    sum
  }
  //----------------------------------------------------------------------------
  // (a) 1H at 1.2 T: ν = -γB/2π = -51.09 MHz (positive gyromagnetic ratio),
  // so m = +1/2 (index 0) is the lower level. At T = h|ν|/k_B the population
  // ratio is e, and the fraction in the lower level is e/(1+e).
  #[test]
  fn test_generate_spin_half(){
    let nu = -51.09e6;
    let tensors = one_bath_spin_tensors(2, Vector3D::from([0.0,0.0,nu]),
        0.0, 0.0);
    let temperature = PLANCK*nu.abs()/BOLTZMANN;
    let config = config_with(Some(temperature));
    let mut rng = ChaCha20Rng::seed_from_u64(7);

    let n = 20_000;
    let mut n_lower = 0;
    for _ in 0..n{
      let states = SpinStates::generate(&mut rng, &tensors, &config).unwrap();
      assert_eq!(states.states.len(), 2);
      assert_eq!(states.states[0], array![[ONE],[ZERO]]);
      let psi = &states.states[1];
      assert!((psi.norm() - 1.0).abs() < 1e-12);
      // ±|m⟩ up to a phase.
      let up = psi[[0,0]].norm();
      let down = psi[[1,0]].norm();
      assert!((up - 1.0).abs() < 1e-12 && down < 1e-12
          || (down - 1.0).abs() < 1e-12 && up < 1e-12);
      if up > 0.5 { n_lower += 1; }
    }
    let e = 1f64.exp();
    let p = e/(1.0 + e);
    let sigma = (p*(1.0 - p)/(n as f64)).sqrt();
    let fraction = n_lower as f64/(n as f64);
    assert!((fraction - p).abs() < 4.0*sigma,
        "fraction {} expected {} ± {}", fraction, p, sigma);
  }
  //----------------------------------------------------------------------------
  // (b) Spin 1 with a tilted axial Q: each sample is an eigenstate, and the
  // sample mean of |ψ⟩⟨ψ| matches exp(-βH)/Z (computed with a Taylor series,
  // not with eigh) within 4σ per element.
  // (c) Without a temperature, the three levels are equally likely.
  #[test]
  fn test_generate_spin_one(){
    let tensors = one_bath_spin_tensors(3, Vector3D::from([0.0,0.0,3.0e6]),
        2.0e6, 0.7);
    let ham = build_single_spin_hamiltonian(&tensors, 1, false).unwrap();
    let ham_norm = ham.norm();
    let temperature = PLANCK*3.0e6/BOLTZMANN;
    let beta = 1.0/(JOULES_TO_HERTZ*BOLTZMANN*temperature);

    let mut rho_ref = expm(&(&ham*(ONE*(-beta))));
    let z = rho_ref.trace().unwrap();
    rho_ref /= z;

    for temperature_opt in [Some(temperature), None]{
      let config = config_with(temperature_opt);
      let mut rng = ChaCha20Rng::seed_from_u64(11);
      let n = 20_000;
      let mut sum = CxMat::zeros((3,3));
      let mut sum_sq_re = Array2::<f64>::zeros((3,3));
      let mut sum_sq_im = Array2::<f64>::zeros((3,3));
      let mut energies = Vec::<f64>::with_capacity(n);
      for _ in 0..n{
        let states = SpinStates::generate(&mut rng, &tensors, &config)
          .unwrap();
        let psi = &states.states[1];
        let h_psi = ham.dot(psi);
        let energy = psi.t().map(|x| x.conj()).dot(&h_psi)[[0,0]];
        assert!(energy.im.abs() < 1e-12*ham_norm);
        assert!((&h_psi - &(psi*energy)).norm() <= 1e-12*ham_norm);
        energies.push(energy.re);

        let rho = psi.dot(&psi.t().map(|x| x.conj()));
        sum += &rho;
        sum_sq_re += &rho.map(|x| x.re*x.re);
        sum_sq_im += &rho.map(|x| x.im*x.im);
      }
      let nf = n as f64;
      let mean = &sum/nf;

      match temperature_opt{
        Some(_) => {
          for i in 0..3{
            for j in 0..3{
              let m = mean[[i,j]];
              let var_re = sum_sq_re[[i,j]]/nf - m.re*m.re;
              let var_im = sum_sq_im[[i,j]]/nf - m.im*m.im;
              let s_re = (var_re.max(0.0)/nf).sqrt();
              let s_im = (var_im.max(0.0)/nf).sqrt();
              let r = rho_ref[[i,j]];
              assert!((m.re - r.re).abs() <= 4.0*s_re + 1e-12,
                  "({},{}) re {} vs {} ± {}", i, j, m.re, r.re, s_re);
              assert!((m.im - r.im).abs() <= 4.0*s_im + 1e-12,
                  "({},{}) im {} vs {} ± {}", i, j, m.im, r.im, s_im);
            }
          }
        },
        None => {
          // Count the levels by energy.
          let (levels, _) = eigh_hermitian(&ham, UPLO::Lower).unwrap();
          let mut counts = [0usize; 3];
          for e in energies.iter(){
            let k = (0..3).min_by(|&a,&b| (levels[a] - e).abs()
                .partial_cmp(&(levels[b] - e).abs()).unwrap()).unwrap();
            assert!((levels[k] - e).abs() <= 1e-9*ham_norm);
            counts[k] += 1;
          }
          let p = 1.0/3.0;
          let sigma = (p*(1.0 - p)/nf).sqrt();
          for c in counts.iter(){
            assert!((*c as f64/nf - p).abs() < 4.0*sigma);
          }
          // The mean is then 1/3.
          assert!((mean - CxMat::eye(3)/3.0).norm() < 0.05);
        },
      }
    }
  }
  //----------------------------------------------------------------------------
  // The random number stream does not depend on the ensemble setting, and
  // the same seed gives the same states.
  #[test]
  fn test_generate_is_reproducible(){
    let tensors = one_bath_spin_tensors(3, Vector3D::from([0.0,0.0,3.0e6]),
        2.0e6, 0.7);
    let mut config = config_with(Some(1e-4));
    config.ensemble_cce = Some(false);
    let a = SpinStates::generate(&mut ChaCha20Rng::seed_from_u64(3),
        &tensors, &config).unwrap();
    config.ensemble_cce = Some(true);
    let b = SpinStates::generate(&mut ChaCha20Rng::seed_from_u64(3),
        &tensors, &config).unwrap();
    assert_eq!(a.states, b.states);
  }
  //----------------------------------------------------------------------------
}
