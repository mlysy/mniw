/// @file Hierarchical.h
/// @author Martin Lysy (mlysy@uwaterloo.ca)
///
/// Gibbs sampler for posterior distribution of Hierarchical Normal model with unequal variances

#ifndef Hierarchical_h
#define Hierarchical_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
#include "TriUtils.h"
#include "Wishart.h"
#include "MatrixNormal.h"
#include "RandomEffects.h"

/// Gibbs sampler for posterior distribution of Hierarchical Normal model with unequal variances
///
/// @note TODO:
///
/// * Consistent transpose notation for `Y` and `Mu`.
/// * Rename `OmegaL` to `iSigmaL`, or clarify notation.
/// * Name parameters more carefully.
///
class HierNormal {
 private:
  // dimensions of problem
  int N_; // number of observations
  int p_; // number of row dimensions
  int q_; // number of column dimensions
  // data needed for updates
  MatrixXd X_;
  MatrixXd Psi_;
  // random variables in problem
  bool flatBeta_; // whether Beta has a flat prior
  MatrixXd Mu_; // random effect coefficients
  MatrixXd Beta_; // regression coefficients
  MatrixXd OmegaL_; // lower cholesky factor of inverse error-variance
  // storage variables
  double nuHat_;
  MatrixXd Yt_;
  MatrixXd iXtXO_;
  MatrixXd RowiOmegaHatL_;
  MatrixXd IH_;
  MatrixXd LtO_;
  MatrixXd LtOL_;
  MatrixXd PsiHat_;
  MatrixXd XiHatL_;
  MatrixXd MuXLtO_;
  MatrixXd IHMut_;
  MatrixXd LambdaHat_;
  MatrixXd x_;
  MatrixXd iV_;
  MatrixXd Lq_;
  MatrixXd Uq_;
  MatrixXd Iq_;
  LLT<MatrixXd> lltq_;
  // distributions
  Wishart *wish_;
  MatrixNormal *matnorm_;
  RanfxNormal *renorm_;
 public:
  /// Constructor
  HierNormal(const Ref<const MatrixXd>& Y,
	     const Ref<const MatrixXd>& X,
	     const Ref<const MatrixXd>& V,
	     const Ref<const MatrixXd>& Lambda,
	     const Ref<const MatrixXd>& Omega,
	     const Ref<const MatrixXd>& Psi, double nu);
  /// Destructor
  ~HierNormal();
  /// Set random effects
  void setMu(const Ref<const MatrixXd>& Mu);
  /// Set lower Cholesky factor of random effects precision matrix
  void setOmegaL(const Ref<const MatrixXd>& OmegaL);
  /// Set random effects precision matrix
  void setISigma(const Ref<const MatrixXd>& iSigma);
  /// Set random effects regression coefficients
  void setBeta(const Ref<const MatrixXd>& Beta);
  /// Get random effects
  void getMu(Ref<MatrixXd> Mu);
  /// Get lower Cholesky factor of random effects precision matrix
  void getOmegaL(Ref<MatrixXd> OmegaL);
  /// Get random effects variance matrix
  void getSigma(Ref<MatrixXd> Sigma);
  /// Get random effects regression coefficients
  void getBeta(Ref<MatrixXd> Beta);
  /// Update `Beta` and `Sigma`
  void UpdateBetaSigma();
  /// Update `Mu`
  void UpdateMu();
};

/// @param [in] Y Response matrix of size `N x q`, each row is an observation.
/// @param [in] X Covariate matrix of size `N x q`, each row is an observation.
/// @param [in] V Matrix of size `q x (Nq)` containing the variance of each response about its mean.
/// @param [in] Lambda Prior mean matrix of size `p x q`.
/// @param [in] Omega Prior row-precision matrix of size `p x p`.
/// @param [in] Psi Prior column-variance matrix of size `q x q`.
/// @param [in] nu Scale parameter.
inline HierNormal::HierNormal(const Ref<const MatrixXd>& Y,
			      const Ref<const MatrixXd>& X,
			      const Ref<const MatrixXd>& V,
			      const Ref<const MatrixXd>& Lambda,
			      const Ref<const MatrixXd>& Omega,
			      const Ref<const MatrixXd>& Psi, double nu) {
  // dimensions of the problem
  N_ = Y.rows();
  p_ = X.cols();
  q_ = Y.cols();
  flatBeta_ = Omega.squaredNorm() == 0.0;
  // data needed for updates
  X_ = X;
  Psi_ = Psi;
  // precomputed values
  nuHat_ = nu + N_ - (flatBeta_ ? p_ : 0);
  Yt_ = Y.adjoint();
  // initialize random variables
  Mu_ = MatrixXd::Zero(q_, N_*q_);
  OmegaL_ = MatrixXd::Identity(q_, q_);
  Beta_ = MatrixXd::Zero(p_, q_);
  // internal values
  Lq_ = MatrixXd::Zero(q_,q_);
  Uq_ = MatrixXd::Zero(q_,q_);
  Iq_ = MatrixXd::Identity(q_,q_);
  lltq_.compute(Iq_);
  iXtXO_ = (X.adjoint() * X + Omega).eval().llt().solve(MatrixXd::Identity(p_,p_));
  RowiOmegaHatL_ = iXtXO_.llt().matrixL();
  IH_ = MatrixXd::Identity(N_,N_) - (X * (X.adjoint() * X).llt().solve(X.adjoint())).eval();
  LtO_ = Lambda.adjoint() * Omega;
  LtOL_ = LtO_ * Lambda;
  PsiHat_ = MatrixXd::Zero(q_,q_);
  XiHatL_ = MatrixXd::Zero(q_,q_);
  MuXLtO_ = MatrixXd::Zero(q_,p_);
  IHMut_ = MatrixXd::Zero(N_,q_);
  LambdaHat_ = MatrixXd::Zero(p_,q_);
  x_ = MatrixXd::Zero(q_,N_);
  // inverses of the variances
  iV_ = MatrixXd::Zero(q_,N_*q_);
  for(int jj=0; jj<N_; jj++) {
    iV_.block(0,jj*q_,q_,q_) = V.block(0,jj*q_,q_,q_).llt().solve(Iq_);
  }
  // distributions
  wish_ = new Wishart(q_);
  matnorm_ = new MatrixNormal(p_, q_);
  renorm_ = new RanfxNormal(q_);
}

inline HierNormal::~HierNormal() {
  delete wish_;
  delete matnorm_;
  delete renorm_;
}

/// @param [in] Beta Random effects coefficient matrix of size `p x q`.
inline void HierNormal::setBeta(const Ref<const MatrixXd>& Beta) {
  Beta_ = Beta;
  return;
}

/// @param [in] OmegaL Lower Cholesky factor of the random effects precision matrix of size `q x q`.
inline void HierNormal::setOmegaL(const Ref<const MatrixXd>& OmegaL) {
  OmegaL_ = OmegaL.triangularView<Eigen::Lower>();
  return;
}

/// @param [in] iSigma Random effects precision matrix of size `q x q`.
inline void HierNormal::setISigma(const Ref<const MatrixXd>& iSigma) {
  OmegaL_ = iSigma.llt().matrixL();
  return;
}

/// @param [in] Mu Random effects matrix of size `N x q`.
void HierNormal::setMu(const Ref<const MatrixXd>& Mu) {
  Mu_ = Mu.adjoint();
  return;
}

/// @param Mu [out] Matrix of size `q x N` in which to output current value of random effects.
inline void HierNormal::getMu(Ref<MatrixXd> Mu) {
  Mu = Mu_;
  return;
}

/// @param OmegaL [out] Matrix of size `q x q` in which to output current value of the lower Cholesky factor of the random effects precision matrix.
inline void HierNormal::getOmegaL(Ref<MatrixXd> OmegaL) {
  OmegaL = OmegaL_;
  return;
}

/// @param Sigma [out] Matrix of size `q x q` in which to output current value of the random effects variance matrix.
inline void HierNormal::getSigma(Ref<MatrixXd> Sigma) {
  InverseLLt(Sigma, OmegaL_, Lq_, Uq_, Iq_);
  return;
}

/// @param Beta [out] Matrix of size `p x q` in which to output current value of the random effects regression coefficients.
inline void HierNormal::getBeta(Ref<MatrixXd> Beta) {
  Beta = Beta_;
}


/// Performs one draw of \f$(\boldsymbol(\beta), \boldsymbol(\Sigma)) \sim p(\boldsymbol(\beta), \boldsymbol(\Sigma) \mid \boldsymbol(\mu), \boldsymbol(Y), \boldsymbol(X))\f$ and assigns it to the internal values of `Beta` and `Sigma`.
inline void HierNormal::UpdateBetaSigma() {
  // new MNIW parameters
  MuXLtO_.noalias() = Mu_ * X_ + LtO_;
  LambdaHat_.noalias() = iXtXO_ * MuXLtO_.adjoint();
  if(!flatBeta_) {
    PsiHat_.noalias() = Mu_ * Mu_.adjoint();
    PsiHat_.noalias() -= MuXLtO_ * LambdaHat_;
    PsiHat_ += LtOL_ + Psi_;
  }
  else {
    IHMut_.noalias() = IH_ * Mu_.adjoint();
    PsiHat_.noalias() = Mu_ * IHMut_;
    PsiHat_ += Psi_;
  }
  // reverse cholesky decomposition PsiHat = XiHatU * XiHatL
  ReverseCholesky(XiHatL_, PsiHat_, lltq_);
  // generate Beta and OmegaL, where Sigma^{-1} = OmegaL * OmegaL'
  wish_->GenerateLowerTriXi(OmegaL_, XiHatL_, nuHat_);
  matnorm_->GenerateRowSColO(Beta_, LambdaHat_, RowiOmegaHatL_, OmegaL_);
  return;
}

inline void HierNormal::UpdateMu() {
  x_ = (X_ * Beta_).adjoint();
  x_ -= Yt_;
  for(int jj=0; jj<N_; jj++) {
    renorm_->GenerateOL(Mu_.col(jj), x_.col(jj), iV_.block(0,jj*q_,q_,q_), OmegaL_);
  }
  Mu_ += Yt_;
  return;
}

#endif
