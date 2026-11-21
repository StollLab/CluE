//! Workaround for wrong eigenvectors from `ndarray_linalg::Eigh::eigh` for
//! complex Hermitian matrices in C (row-major) layout.
//!
//! The bug, in ndarray-linalg 0.17.0 (used by CluE) and in 0.18.1 (the latest
//! release, checked 2026-11-20):
//! `EighInplace::eigh_inplace` (ndarray-linalg `src/eigh.rs`) handles a
//! C-layout matrix by transposing it (`self.swap_axes(0, 1)`, comment
//! "XXX Force layout to be Fortran (see #146)") and passes the buffer to
//! LAPACK `zheev` without conjugating it. For a Hermitian A, the transpose is
//! the complex conjugate, A^T = A*, so the returned eigenvectors are those
//! of A*, not of A. The eigenvalues are correct, because A and A* have the
//! same spectrum. Every `CxMat` in CluE is in C layout. Real symmetric
//! matrices are not affected.
//! The issue is reported upstream as
//! https://github.com/rust-ndarray/ndarray-linalg/issues/307
//! ("wrong eigh eigenvectors for complex valued matrix", open, with no fix).
//!
//! The workaround calls `eigh` on a Fortran-layout copy of the matrix. That
//! path is not transposed, so LAPACK sees A itself. It also stays correct if
//! a future ndarray-linalg conjugates C-layout input properly.
//!
//! The test `test_upstream_eigh_still_conjugates_c_layout` checks that the
//! bug is still present. If it fails after an ndarray-linalg update, the bug
//! is probably fixed upstream, and this module may be removed (replace
//! `eigh_hermitian(&a, uplo)` by `a.eigh(uplo)` again).

use ndarray::{Array1, Array2, ShapeBuilder};
use ndarray_linalg::{Eigh, UPLO};
use ndarray_linalg::error::Result as LinalgResult;
use num_complex::Complex64;

type CxMat = Array2::<Complex64>;

//------------------------------------------------------------------------------
/// This function returns the eigenvalues (ascending) and the eigenvectors
/// (columns) of the Hermitian matrix `mat`, such that
/// `mat · V = V · diag(E)`.
/// It is `mat.eigh(uplo)` applied to a Fortran-layout copy of `mat`; see the
/// module documentation for why the copy is needed.
pub fn eigh_hermitian(mat: &CxMat, uplo: UPLO)
  -> LinalgResult<(Array1::<f64>, CxMat)>
{
  let mut mat_f = CxMat::zeros(mat.raw_dim().f());
  mat_f.assign(mat);
  mat_f.eigh(uplo)
}
//------------------------------------------------------------------------------

//==============================================================================
#[cfg(test)]
mod tests{
  use super::*;
  use ndarray::array;
  use ndarray_linalg::Norm;

  //----------------------------------------------------------------------------
  fn c(re: f64, im: f64) -> Complex64{ Complex64{re, im} }
  //----------------------------------------------------------------------------
  // Complex Hermitian test matrices in C layout: the matrix of upstream
  // issue #307, and a 4×4 matrix with complex entries in every off-diagonal
  // element.
  fn test_matrices() -> Vec::<CxMat>{
    let a = array![
      [c(1.0, 0.0), c(0.0,-2.0)],
      [c(0.0, 2.0), c(5.0, 0.0)],
    ];
    let b = array![
      [c( 2.0, 0.0), c( 0.3,-1.1), c(-0.7, 0.4), c( 0.2, 0.9)],
      [c( 0.3, 1.1), c(-1.0, 0.0), c( 0.5, 0.5), c(-1.3,-0.2)],
      [c(-0.7,-0.4), c( 0.5,-0.5), c( 0.5, 0.0), c( 0.8,-0.6)],
      [c( 0.2,-0.9), c(-1.3, 0.2), c( 0.8, 0.6), c( 3.0, 0.0)],
    ];
    vec![a, b]
  }
  //----------------------------------------------------------------------------
  // max_k ||A v_k - λ_k v_k||, relative to ||A||.
  fn residual(mat: &CxMat, eigvals: &Array1::<f64>, eigvecs: &CxMat) -> f64{
    let mut res = 0.0f64;
    for k in 0..eigvals.len(){
      let v = eigvecs.column(k).to_owned();
      let r = mat.dot(&v) - &v*c(eigvals[k], 0.0);
      res = res.max(r.norm());
    }
    res/mat.norm()
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_eigh_hermitian(){
    for mat in test_matrices(){
      assert!(mat.is_standard_layout());
      let hermiticity = (&mat - &mat.t().map(|x| x.conj())).norm();
      assert_eq!(hermiticity, 0.0);

      for uplo in [UPLO::Lower, UPLO::Upper]{
        let (eigvals, eigvecs) = eigh_hermitian(&mat, uplo).unwrap();
        let n = eigvals.len();

        // A V = V diag(E).
        assert!(residual(&mat, &eigvals, &eigvecs) < 1e-14);

        // V is unitary.
        let vdv = eigvecs.t().map(|x| x.conj()).dot(&eigvecs);
        assert!((vdv - CxMat::eye(n)).norm() < 1e-14);

        // Ascending eigenvalues.
        for k in 1..n{ assert!(eigvals[k-1] <= eigvals[k]); }
      }
    }

    // The 2×2 matrix of issue #307 has the eigenvalues 3 ∓ √8.
    let (eigvals, _) = eigh_hermitian(&test_matrices()[0], UPLO::Lower)
      .unwrap();
    assert!((eigvals[0] - (3.0 - 8f64.sqrt())).abs() < 1e-14);
    assert!((eigvals[1] - (3.0 + 8f64.sqrt())).abs() < 1e-14);

    // The input is not changed.
    let mat = test_matrices()[1].clone();
    let _ = eigh_hermitian(&mat, UPLO::Lower).unwrap();
    assert_eq!(mat, test_matrices()[1]);
  }
  //----------------------------------------------------------------------------
  // This test checks that the upstream bug is still present. If it fails,
  // ndarray-linalg's eigh now returns the eigenvectors of A for complex
  // C-layout input, and this workaround module may be removed.
  #[test]
  fn test_upstream_eigh_still_conjugates_c_layout(){
    for mat in test_matrices(){
      assert!(mat.is_standard_layout());
      let (eigvals, eigvecs) = mat.eigh(UPLO::Lower).unwrap();
      let res = residual(&mat, &eigvals, &eigvecs);
      assert!(res > 1e-3,
        "ndarray-linalg eigh now returns correct eigenvectors for a complex \
Hermitian matrix in C layout (max ||Av - λv||/||A|| = {:e}). The upstream \
bug (rust-ndarray/ndarray-linalg issue #307) appears to be fixed: the \
workaround in src/ndarray_linalg_eigh_workaround.rs may be removed.", res);

      // The returned vectors are the eigenvectors of A*.
      let mat_conj = mat.map(|x| x.conj());
      assert!(residual(&mat_conj, &eigvals, &eigvecs) < 1e-14);
    }
  }
  //----------------------------------------------------------------------------
}
