/// @file MatrixT.h
/// @author Martin Lysy (mlysy@uwaterloo.ca)
///
/// Density evaluation and random number generation for the Matrix-t distribution

#ifndef MatrixT_h
#define MatrixT_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
#include "TriUtils.h"
#include "Wishart.h"

/// The Matrix-t distribution
class MatrixT {
 private:
  // storage
  int p_;
  int q_;
  int pq_; // p * q
  // int d_; // min(p, q)
  bool pMin_; // isTRUE(p < q)
  MatrixXd Z_;
  LLT<MatrixXd> cholRowV_;
  LLT<MatrixXd> cholColV_;
  LLT<MatrixXd> lltd_;
  MatrixXd A_;
  MatrixXd B_;
  // MatrixXd Mqp_;
  // MatrixXd Mpq_;
  // MatrixXd Mp_;
  // VectorXd Ip_;
 public:
  /// Constructor
  MatrixT(int p, int q);
  /// Log-density evaluation
  double LogDens(const Ref<const MatrixXd>&  X,
		 const Ref<const MatrixXd>& Mu,
		 const Ref<const MatrixXd>& RowV,
		 const Ref<const MatrixXd>& ColV,
		 double nu);
  /// Log-density evaluation with precomputations
  double LogDens(const Ref<const MatrixXd>& X,
		 const Ref<const MatrixXd>& Mu,
		 const Ref<const MatrixXd>& RowV,
		 LLT<MatrixXd>& cholRowV, double ldRowV,
		 const Ref<const MatrixXd>& ColV,
		 LLT<MatrixXd>& cholColV, double ldColV, double nu);
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
  // lltd_.compute(MatrixXd::Identity(d_,d_));
  // Mpq_ = MatrixXd::Zero(p_,q_);
  // Mqp_ = MatrixXd::Zero(q_,p_);
  // Mp_ = MatrixXd::Zero(p_,p_);
  // Ip_ = VectorXd::Ones(p_);
}

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

/// Identical to MatrixT::LogDens(const Ref<const MatrixXd>&,const Ref<const MatrixXd>&,const Ref<const MatrixXd>&,const Ref<const MatrixXd>&,double), except with pre-computed Cholesky factor and log-determinants for `RowV` and `ColV` (faster in a for-loop where one of these is held fixed).
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
  // return ldens;
  // ldens += pq_ * M_LN_SQRT_PI;
  return -ldens + logMultiGamma(.5 * nupq, q_) - logMultiGamma(.5 * nuq, q_);
}

#endif
