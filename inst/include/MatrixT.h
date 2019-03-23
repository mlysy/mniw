/// @file MatrixT.h
///
/// Density evaluation and random number generation for the Matrix-t distribution.

#ifndef MatrixT_h
#define MatrixT_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
#include "TriUtils.h"
#include "Wishart.h"
#include "MatrixNormal.h"

/// The Matrix-t distribution.
class MatrixT {
 private:
  // storage
  int p_;
  int q_;
  int pq_; // p * q
  bool pMin_; // isTRUE(p < q)
  MatrixXd Z_;
  LLT<MatrixXd> cholRowV_;
  LLT<MatrixXd> cholColV_;
  LLT<MatrixXd> lltd_;
  MatrixXd A_;
  MatrixXd B_;
  MatrixXd CL_;
  Wishart *wish_;
  MatrixNormal *matnorm_;
 public:
  /// Constructor.
  MatrixT(int p, int q);
  /// Destructor.
  ~MatrixT();
  /// Log-density evaluation.
  double LogDens(const Ref<const MatrixXd>&  X,
		 const Ref<const MatrixXd>& Mu,
		 const Ref<const MatrixXd>& RowV,
		 const Ref<const MatrixXd>& ColV,
		 double nu);
  /// Log-density evaluation with precomputations.
  double LogDens(const Ref<const MatrixXd>& X,
		 const Ref<const MatrixXd>& Mu,
		 const Ref<const MatrixXd>& RowV,
		 LLT<MatrixXd>& cholRowV, double ldRowV,
		 const Ref<const MatrixXd>& ColV,
		 LLT<MatrixXd>& cholColV, double ldColV, double nu);
  /// Random draw with RowV/ColV on the variance/precision scale.
  void GenerateRowSColO(Ref<MatrixXd> X,
			const Ref<const MatrixXd>& Mu,
			const Ref<const MatrixXd>& RowVL,
			const Ref<const MatrixXd>& iColVL, double nu);
  /// Random draw with RowV/ColV on the precision/precision scale.
  void GenerateRowOColO(Ref<MatrixXd> X,
			const Ref<const MatrixXd>& Mu,
			const Ref<const MatrixXd>& iRowVU,
			const Ref<const MatrixXd>& iColVL, double nu);
};

/// @param [in] p Number of rows of Matrix-t distribution.
/// @param [in] q Number of columns of Matrix-t distribution.
inline MatrixT::MatrixT(int p, int q) {
  // problem dimensions
  p_ = p;
  q_ = q;
  // precomputations
  pMin_ = p_ < q_;
  pq_ = p_*q_; 
  // storage
  Z_ = MatrixXd::Zero(p_,q_);
  cholRowV_.compute(MatrixXd::Identity(p_,p_));
  cholColV_.compute(MatrixXd::Identity(q_,q_));
  if(pMin_) {
    A_ = MatrixXd::Zero(q_,p_);
    B_ = MatrixXd::Zero(p_,p_);
    lltd_.compute(MatrixXd::Identity(p_,p_));
  }
  else {
    A_ = MatrixXd::Zero(p_,q_);
    B_ = MatrixXd::Zero(q_,q_);
    lltd_.compute(MatrixXd::Identity(q_,q_));
  }
  CL_ = MatrixXd::Zero(q_,q_);
  wish_ = new Wishart(q_);
  matnorm_ = new MatrixNormal(p_,q_);
}

inline MatrixT::~MatrixT() {
  delete wish_;
  delete matnorm_;
}

// --- log-density --------------------------------------------------------

/// @param [in] X Observation matrix of size `p x q`.
/// @param [in] Mu Mean matrix of size `p x q`.
/// @param [in] RowV Row-variance matrix of size `p x p`.
/// @param [in] ColV Column-variance matrix of size `q x q`.
/// @param [in] nu Shape parameter.
///
/// @return The log-density evaluated at `X`.
inline double MatrixT::LogDens(const Ref<const MatrixXd>&  X,
			       const Ref<const MatrixXd>& Mu,
			       const Ref<const MatrixXd>& RowV,
			       const Ref<const MatrixXd>& ColV,
			       double nu) {
  cholRowV_.compute(RowV);
  cholColV_.compute(ColV);
  return LogDens(X, Mu,
		 RowV, cholRowV_, logDetCholV(cholRowV_),
		 ColV, cholColV_, logDetCholV(cholColV_), nu);
}

/// Identical to the shorter `LogDens`, except with pre-computed Cholesky factor and log-determinants for `RowV` and `ColV` (faster in a for-loop where one of these is held fixed).
///
/// @param [in] X Observation matrix of size `p x q`.
/// @param [in] Mu Mean matrix of size `p x q`.
/// @param [in] RowV Row-variance matrix of size `p x p`.
/// @param [in] cholRowV Lower Cholesky factor of row-variance matrix, as an `Eigen::LLT` object.
/// @param [in] ldRowV Log-determinant of `cholRowV`.
/// @param [in] ColV Column-variance matrix of size `q x q`.
/// @param [in] cholColV Lower Cholesky factor of column-variance matrix, as an `Eigen::LLT` object.
/// @param [in] ldColV Log-determinant of `cholColV`.
/// @param [in] nu Shape parameter.
///
/// @return The log-density evaluated at `X`.
inline double MatrixT::LogDens(const Ref<const MatrixXd>& X,
			       const Ref<const MatrixXd>& Mu,
			       const Ref<const MatrixXd>& RowV,
			       LLT<MatrixXd>& cholRowV, double ldRowV,
			       const Ref<const MatrixXd>& ColV,
			       LLT<MatrixXd>& cholColV, double ldColV,
			       double nu) {
  double ldens;
  double nuq = nu + q_ - 1.0;
  double nupq = nuq + p_;
  // calculate log-determinant by Sylvester formula
  Z_.noalias() = X - Mu;
  if(pMin_) {
    A_ = cholColV.solve(Z_.adjoint());
    B_.noalias() = RowV + Z_ * A_;
    lltd_.compute(B_);
    ldens = logDetCholV(lltd_) - ldRowV;
  }
  else {
    A_ = cholRowV.solve(Z_);
    B_.noalias() = ColV + Z_.adjoint() * A_;
    lltd_.compute(B_);
    ldens = logDetCholV(lltd_) - ldColV;
  }
  // add components
  ldens = nupq * ldens + q_ * ldRowV + p_ * ldColV + pq_ * M_LN_SQRT_PI;
  return -ldens + logMultiGamma(.5 * nupq, q_) - logMultiGamma(.5 * nuq, q_);
}

// --- rng ----------------------------------------------------------------

/// @param [out] X Matrix of size `p x q` in which to store the random draw.
/// @param [in] Mu Mean matrix of size `p x q`.
/// @param [in] RowVL Lower Cholesky factor of the row-variance matrix `RowV` (a matrix of size `p x p`).
/// @param [in] iColVL Lower Cholesky factor of the column-precision matrix: `ColV^{-1} = iColVL * t(iColVL)` (a matrix of size `q x q`).
/// @param [in] nu Shape parameter.
inline void MatrixT::GenerateRowSColO(Ref<MatrixXd> X,
				      const Ref<const MatrixXd>& Mu,
				      const Ref<const MatrixXd>& RowVL,
				      const Ref<const MatrixXd>& iColVL,
				      double nu) {
  wish_->GenerateLowerTriXi(CL_, iColVL, nu);
  matnorm_->GenerateRowSColO(X, Mu, RowVL, CL_);
  return;
}

/// @param [out] X Matrix of size `p x q` in which to store the random draw.
/// @param [in] Mu Mean matrix of size `p x q`.
/// @param [in] iRowVU Upper Cholesky factor of the row-precision matrix: `RowV^{-1} = t(iRowVU) * iRowVU` (a matrix of size `p x p`).
/// @param [in] iColVL Lower Cholesky factor of the column-precision matrix: `ColV^{-1} = iColVL * t(iColVL)` (a matrix of size `q x q`).
/// @param [in] nu Shape parameter.
inline void MatrixT::GenerateRowOColO(Ref<MatrixXd> X,
				      const Ref<const MatrixXd>& Mu,
				      const Ref<const MatrixXd>& iRowVU,
				      const Ref<const MatrixXd>& iColVL,
				      double nu) {
  wish_->GenerateLowerTriXi(CL_, iColVL, nu);
  matnorm_->GenerateRowOColO(X, Mu, iRowVU, CL_);
  return;
}

#endif
