/// @file MatrixNormal.h
///
/// Density evaluation and random number generation for the Matrix Normal distribution.

#ifndef MatrixNormal_h
#define MatrixNormal_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
#include "TriUtils.h"

/// The Matrix Normal distribution.
class MatrixNormal {
 private:
  // storage
  int p_;
  int q_;
  MatrixXd Z_;
  LLT<MatrixXd> cholRowV_;
  // double ldRowV;
  LLT<MatrixXd> cholColV_;
  // double ldColV;
 public:
  /// Constructor.
  MatrixNormal(int p, int q);
  /// Log-density evaluation.
  double LogDens(const Ref<const MatrixXd>& X, const Ref<const MatrixXd>& Mu,
		 const Ref<const MatrixXd>& RowV, const Ref<const MatrixXd>& ColV);
  /// Log-density evaluation with precomputations.
  double LogDens(const Ref<const MatrixXd>& X,
		 const Ref<const MatrixXd>& Mu,
		 LLT<MatrixXd>& cholRowV, double ldRowV, 
		 LLT<MatrixXd>& cholColV, double ldColV);
  /// Random draw with RowV/ColV on the variance/variance scale.
  void GenerateRowSColS(Ref<MatrixXd> X,
			const Ref<const MatrixXd>& Lambda,
			const Ref<const MatrixXd>& RowSigmaL,
			const Ref<const MatrixXd>& ColSigmaU);
  /// Random number generation with RowV/ColV on the variance/precision scale.
  void GenerateRowSColO(Ref<MatrixXd> X,
			const Ref<const MatrixXd>& Lambda,
			const Ref<const MatrixXd>& RowSigmaL,
			const Ref<const MatrixXd>& ColOmegaL);
  /// Random number generation with RowV/ColV on the precision/precision scale.
  void GenerateRowOColO(Ref<MatrixXd> X,
			const Ref<const MatrixXd>& Lambda,
			const Ref<const MatrixXd>& RowOmegaU,
			const Ref<const MatrixXd>& ColOmegaL);
};

/// @param [in] p Number of rows of Matrix Normal distribution.
/// @param [in] q Number of columns of Matrix Normal distribution.
inline MatrixNormal::MatrixNormal(int p, int q) {
  p_ = p;
  q_ = q;
  Z_ = MatrixXd::Zero(p,q);
  cholRowV_.compute(MatrixXd::Identity(q_,q_));
  cholColV_.compute(MatrixXd::Identity(p_,p_));
}

// --- log-density -------------------------------------------------------------

/// @param [in] X Observation matrix of size `p x q`.
/// @param [in] Mu Mean matrix of size `p x q`.
/// @param [in] RowV Row-variance matrix of size `p x p`.
/// @param [in] ColV Column-variance matrix of size `q x q`.
///
/// @return The log-density evaluated at `X`.
inline double MatrixNormal::LogDens(const Ref<const MatrixXd>& X,
				    const Ref<const MatrixXd>& Mu,
				    const Ref<const MatrixXd>& RowV,
				    const Ref<const MatrixXd>& ColV) {
  cholRowV_.compute(RowV);
  cholColV_.compute(ColV);
  return LogDens(X, Mu,
		 cholRowV_, logDetCholV(cholRowV_),
		 cholColV_, logDetCholV(cholColV_));
}

/// Identical to the shorter `LogDens`, except with pre-computed Cholesky factors and log-determinants for `RowV` and `ColV` (faster in a for-loop where one of these is held fixed).
///
/// @param [in] X Observation matrix of size `p x q`.
/// @param [in] Mu Mean matrix of size `p x q`.
/// @param [in] cholRowV Lower Cholesky factor of row-variance matrix, as an `Eigen::LLT` object.
/// @param [in] ldRowV Log-determinant of `cholRowV`.
/// @param [in] cholColV Lower Cholesky factor of column-variance matrix, as an `Eigen::LLT` object.
///
/// @param [in] ldColV Log-determinant of `cholColV`.
///
/// @return The log-density evaluated as `X`.
inline double MatrixNormal::LogDens(const Ref<const MatrixXd>& X,
				    const Ref<const MatrixXd>& Mu,
				    LLT<MatrixXd>& cholRowV, double ldRowV, 
				    LLT<MatrixXd>& cholColV, double ldColV) {
  //double ldens;
  // double p = X.rows();
  // double q = X.cols();
  Z_ = X-Mu;
  // if(CalcRowV) {
  //   cholRowV.compute(RowV);
  //   ldRowV = 0.0;
  //   for(int ii=0; ii<p; ii++) {
  //     ldRowV += log(cholRowV.matrixL()(ii,ii));
  //   }
  // }
  // if(CalcColV) {
  //   cholColV.compute(ColV);
  //   ldColV = 0.0;
  //   for(int ii=0; ii<q; ii++) {
  //     ldColV += log(cholColV.matrixL()(ii,ii));
  //   }
  // }
  cholRowV.matrixL().solveInPlace(Z_);
  cholColV.matrixU().solveInPlace<OnTheRight>(Z_);
  return -(.5 * Z_.squaredNorm() + q_ * ldRowV + p_ * ldColV + p_*q_ * M_LN_SQRT_2PI);
}

// --- rng ---------------------------------------------------------------------

/// @param [out] X Matrix of size `p x q` in which to store the random draw.
/// @param [in] Lambda Mean matrix of size `p x q`.
/// @param [in] RowSigmaL Lower Cholesky factor of the row-variance matrix `RowV` (a matrix of size `p x p`).
/// @param [in] ColSigmaU Upper Cholesky factor of the column-variance matrix `ColV` (a matrix of size `q x q`).
inline void MatrixNormal::GenerateRowSColS(Ref<MatrixXd> X,
					   const Ref<const MatrixXd>& Lambda,
					   const Ref<const MatrixXd>& RowSigmaL,
					   const Ref<const MatrixXd>& ColSigmaU) {
  int ii, jj;
  // int p = Lambda.rows();
  // int q = Lambda.cols();
  // populate Z_ with iid normals
  for(ii=0; ii<p_; ii++) {
    for(jj=0; jj<q_; jj++) {
      Z_(ii,jj) = norm_rand();
    }
  }
  // scaling
  triMultXU(X, Z_, ColSigmaU);
  triMultLX(Z_, RowSigmaL, X);
  X = Z_ + Lambda;
  return;
}

/// @param [out] X Matrix of size `p x q` in which to store the random draw.
/// @param [in] Lambda Mean matrix of size `p x q`.
/// @param [in] RowSigmaL Lower Cholesky factor of the row-variance matrix `RowV` (a matrix of size `p x p`).
/// @param [in] ColOmegaL Lower Cholesky factor of the column-precision matrix `ColV^{-1}` (a matrix of size `q x q`).
inline void MatrixNormal::GenerateRowSColO(Ref<MatrixXd> X,
					   const Ref<const MatrixXd>& Lambda,
					   const Ref<const MatrixXd>& RowSigmaL,
					   const Ref<const MatrixXd>& ColOmegaL) {
  int ii, jj;
  // int p = Lambda.rows();
  // int q = Lambda.cols();
  // populate Z_ with iid normals
  for(ii=0; ii<p_; ii++) {
    for(jj=0; jj<q_; jj++) {
      Z_(ii,jj) = norm_rand();
    }
  }
  // scaling
  triMultXLi(Z_, ColOmegaL);
  triMultLX(X, RowSigmaL, Z_);
  X += Lambda;
  return;
}

/// @param [out] X Matrix of size `p x q` in which to store the random draw.
/// @param [in] Lambda Mean matrix of size `p x q`.
/// @param [in] RowOmegaU Upper Cholesky factor of the row-precision matrix `RowV^{-1}` (a matrix of size `p x p`).
/// @param [in] ColOmegaL Lower Cholesky factor of the column-precision matrix `ColV^{-1}` (a matrix of size `q x q`).
inline void MatrixNormal::GenerateRowOColO(Ref<MatrixXd> X,
					   const Ref<const MatrixXd>& Lambda,
					   const Ref<const MatrixXd>& RowOmegaU,
					   const Ref<const MatrixXd>& ColOmegaL) {
  int ii, jj;
  // int p = Lambda.rows();
  // int q = Lambda.cols();
  // populate Z with iid normals
  for(ii=0; ii<p_; ii++) {
    for(jj=0; jj<q_; jj++) {
      Z_(ii,jj) = norm_rand();
    }
  }
  // scaling
  triMultXLi(Z_, ColOmegaL);
  triMultUiX(RowOmegaU, Z_);
  X = Z_ + Lambda;
  return;
}

#endif
