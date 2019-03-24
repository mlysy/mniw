/// @file RandomEffectsExports.cpp
///
/// Exported Rcpp functions for the Multivariate Random-Effects Normal distribution.

// #include <Rcpp.h>
// using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
//#include <iostream>
// #include "mniwWishart.h"
// #include "mniwMatNorm.h"
// #include "mniwRandomEffects.h"
#include "TriUtils.h"
#include "RandomEffects.h"
using namespace mniw;

//////////////////////////////////////////////////////////////////

/// Generate a random sample from the Random-Effects Normal distribution.
///
/// Generate `N` independent draws from \f$p(\boldsymbol{mu} \mid \boldsymbol{y})\f$, where
/// \f[
/// \begin{aligned}
/// \boldsymbol{\mu} & \sim \mathcal N(\boldsymbol{\lambda}, \boldsymbol{A}) \\
/// \boldsymbol{y} \mid \boldsymbol{\mu} & \sim \mathcal N(\boldsymbol{\mu}, \boldsymbol{V}).
/// \end{aligned}
/// \f]
/// Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] N Integer number of random draws.
/// @param [in] lambda Matrix of size `q x n` of prior means.
/// @param [in] y Matrix of size `q x n` of observations.
/// @param [in] V Matrix of size `q x nq` of observation variances.
/// @param [in] A Matrix of size `q x nq` of prior variances.
///
/// @return Matrix of size `q x Nq` of random draws of observation means `mu`.
//[[Rcpp::export]]
Eigen::MatrixXd GenerateRandomEffectsNormal(int N,
					    Eigen::MatrixXd lambda,
					    Eigen::MatrixXd y,
					    Eigen::MatrixXd V,
					    Eigen::MatrixXd A) {
  int q = lambda.rows();
  bool singleLambda = (lambda.cols() == 1);
  bool singleY = (y.cols() == 1);
  bool singleV = (V.cols() == q);
  bool singleA = (A.cols() == q);
  // output variables
  MatrixXd mu(q,N);
  // internal variables
  // TempPQ *tmp = new TempPQ(1,q);
  LLT<MatrixXd> lltq(q);
  MatrixXd Iq = MatrixXd::Identity(q,q);
  MatrixXd C = MatrixXd::Zero(q,q);
  MatrixXd Omega = MatrixXd::Zero(q,q);
  VectorXd x = VectorXd::Zero(q);
  RanfxNormal renorm(q);
  if(singleV) {
    // tmp->lltq.compute(V);
    // C = tmp->lltq.solve(tmp->Iq);
    lltq.compute(V);
    C = lltq.solve(Iq);
  }
  if(singleA) {
    // tmp->lltq.compute(A);
    // Omega = tmp->lltq.solve(tmp->Iq);
    lltq.compute(A);
    Omega = lltq.solve(Iq);
  }
  if(singleLambda && singleY) {
    x = lambda - y;
  }
  for(int ii=0; ii<N; ii++) {
    if(!singleV) {
      // tmp->lltq.compute(V.block(0,ii*q,q,q));
      // C = tmp->lltq.solve(tmp->Iq);
      lltq.compute(V.block(0,ii*q,q,q));
      C = lltq.solve(Iq);
    }
    if(!singleA) {
      // tmp->lltq.compute(A.block(0,ii*q,q,q));
      // Omega = tmp->lltq.solve(tmp->Iq);
      lltq.compute(A.block(0,ii*q,q,q));
      Omega = lltq.solve(Iq);
    }
    if(!(singleLambda && singleY)) {
      x = lambda.col(ii*(!singleLambda)) - y.col(ii*(!singleY));
    }
    // GenerateRandomEffectsO(mu.col(ii), x, C, Omega, tmp);
    renorm.GenerateO(mu.col(ii), x, C, Omega);
    mu.col(ii) += y.col(ii*(!singleY));
  }
  // delete tmp;
  return(mu);
}


