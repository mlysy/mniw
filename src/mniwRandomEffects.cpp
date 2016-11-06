///////////////////////////////////////////////////////////////////

// Random Effects Normal Models

#include "mniwSetLib.h"

#ifdef R_FUNCTION_LIBRARY
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#endif

#ifdef MATLAB_FUNCTION_LIBRARY
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <Eigen/Dense>
#include "RmathUtils.h"
#endif

//#include <iostream>
using namespace Eigen;
#include "mniwRandomEffects.h"

///////////////////////////////////////////////////////////////////

// Generates a draw from p(mu | x) from the random effects model:
// x | mu ~ N(mu, V)
// mu ~ N(0, Sigma)

// Generate mu = [(I-B)*V]^{1/2} * z + B * x,
// where z ~iid N(0,1) and B = V*(V+Sigma)^{-1}.
// To be compatible with MNIW simulation, the parameters are:
// C = V^{-1}
// OmegaL, where OmegaL * OmegaU = Sigma^{-1}.
void GenerateRandomEffectsOL(Ref<VectorXd> mu, const Ref<VectorXd>& x,
			     const Ref<MatrixXd>& C, const Ref<MatrixXd>& OmegaL,
			     TempPQ *tmp) {
  int q = mu.size();
  // mu = Q * x
  CrossProdLLt(tmp->Mq, OmegaL, tmp->Uq);
  mu.noalias() = tmp->Mq * x;
  // mu = G_L^{-1} z
  tmp->Mq += C;
  tmp->lltq.compute(tmp->Mq);
  tmp->lltq.matrixL().solveInPlace(mu);
  // mu = mu + N(0,I)
  for(int ii=0; ii<q; ii++) {
    mu[ii] += norm_rand();
  }
  // mu = G_U^{-1} z + y
  tmp->lltq.matrixU().solveInPlace(mu);
  return;
}

// Generate mu = [(I-B)*V]^{1/2} * z + B * x,
// where z ~iid N(0,1) and B = V*(V+Sigma)^{-1}.
// the parameters are:
// C = V^{-1}
// Omega = Sigma^{-1}.
void GenerateRandomEffectsO(Ref<VectorXd> mu, const Ref<VectorXd>& x,
			    const Ref<MatrixXd>& C, const Ref<MatrixXd>& Omega,
			    TempPQ *tmp) {
  int q = mu.size();
  // mu = Q * x
  //CrossProdLLt(tmp->Mq, OmegaL, tmp->Uq);
  mu.noalias() = Omega * x;
  // mu = G_L^{-1} z
  tmp->Mq = Omega + C;
  tmp->lltq.compute(tmp->Mq);
  tmp->lltq.matrixL().solveInPlace(mu);
  // mu = mu + N(0,I)
  for(int ii=0; ii<q; ii++) {
    mu[ii] += norm_rand();
  }
  // mu = G_U^{-1} z + y
  tmp->lltq.matrixU().solveInPlace(mu);
  return;
}

