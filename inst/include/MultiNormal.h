/// @file MultiNormal.h
/// @author Martin Lysy (mlysy@uwaterloo.ca)
///
/// Density evaluation and random number generation for the Multivariate Normal distribution

#ifndef MultiNormal_h
#define MultiNormal_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
#include "TriUtils.h"

/// The Multivariate Normal distribution
class MultiNormal {
 private:
  int q_;
  LLT<MatrixXd> cholV_;
  MatrixXd VL_;
  VectorXd z_;
  MatrixXd Z_;
 public:
  MultiNormal(int q);
  double LogDens(const Ref<const VectorXd>& x,
		 const Ref<const VectorXd>& mu,
		 LLT<MatrixXd>& cholV, double ldV);
  double LogDens(const Ref<const VectorXd>& x,
		 const Ref<const VectorXd>& mu,
		 const Ref<const MatrixXd>& V);
  void Generate(Ref<VectorXd> x, const Ref<const VectorXd>& mu,
		LLT<MatrixXd>& cholV, Ref<MatrixXd> VL);
  void Generate(Ref<VectorXd> x, const Ref<const VectorXd>& mu,
		const Ref<const MatrixXd>& V);
};

/// Constructor
///
/// @param [in] q Number of dimensions of the Multivariate Normal.
inline MultiNormal::MultiNormal(int q) {
  q_ = q;
  z_ = VectorXd::Zero(q_);
  VL_ = MatrixXd::Zero(q_,q_);
  cholV_.compute(MatrixXd::Identity(q_,q_));
}

/// Log-density with precomputations
///
/// Identical to the shorter version of `LogDens`, but with the Cholesky decomposition and its log-determinant pre-computed.
///
/// @param [in] x Vector of observations.
/// @param [in] mu Mean vector.
/// @param [in] cholV Cholesky decomposition of the variance matrix, supplied as an `Eigen::LLT` object.
/// @param [in] ldV Log-determinant of the Cholesky factor of the variance.
///
/// @return The log-density evaluated at `x`.
inline double MultiNormal::LogDens(const Ref<const VectorXd>& x,
				   const Ref<const VectorXd>& mu,
				   LLT<MatrixXd>& cholV, double ldV) {
  //double ldens;
  // double q = x.size();
  z_ = x-mu;
  // if(CalcV) {
  //   cholV.compute(V);
  //   ldV = 0.0;
  //   for(int ii=0; ii<q; ii++) {
  //     ldV += log(cholV.matrixL()(ii,ii));
  //   }
  // }
  cholV.matrixL().solveInPlace(z_);
  return -(.5 * z_.squaredNorm() + ldV + q_ * M_LN_SQRT_2PI);
}

/// Log-density
///
/// @param [in] x Vector of observations.
/// @param [in] mu Mean vector.
/// @param [in] V Variance matrix.
///
/// @return The log-density evaluated at `x`.
inline double MultiNormal::LogDens(const Ref<const VectorXd>& x,
				   const Ref<const VectorXd>& mu,
				   const Ref<const MatrixXd>& V) {
  cholV_.compute(V);
  return LogDens(x, mu, cholV_, logDetCholV(cholV_));
}

/// Random number generation with pre-computations.
///
/// @param [out] x Vector to which random draw is assigned.
/// @param [in] mu Mean vector.
/// @param [in] cholV Cholesky decomposition of the variance, supplied as an `Eigen::LLT` object.
/// @param [in] VL Dense representation of this object.  This is because `TriMultLX` does not accept `Eigen::TriangularView` objects and hopefully can be removed in the future...
inline void MultiNormal::Generate(Ref<VectorXd> x,
				  const Ref<const VectorXd>& mu,
				  LLT<MatrixXd>& cholV,
				  Ref<MatrixXd> VL) {
  for(int ii=0; ii<q_; ii++) {
    z_(ii) = norm_rand(); // generate iid standard normals
  }
  triMultLX(x, VL, z_); // multiply by Sigma^{1/2}
  x += mu; // add mu
  return;
}

/// Random number generation
///
/// @param [out] x Vector to which random draw is assigned.
/// @param [in] mu Mean vector.
/// @param [in] V Variance matrix.
inline void MultiNormal::Generate(Ref<VectorXd> x,
				  const Ref<const VectorXd>& mu,
				  const Ref<const MatrixXd>& V) {
  cholV_.compute(V);
  VL_ = cholV_.matrixL();
  Generate(x, mu, cholV_, VL_);
}

#endif
