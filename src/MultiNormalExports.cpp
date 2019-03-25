/// @file MultiNormalExports.cpp
///
/// Exported Rcpp functions for the Multivariate Normal distribution.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
//#include <iostream>
#include "MultiNormal.h"
// #include "mniwMatNorm.h"
using namespace mniw;


//////////////////////////////////////////////////////////////////

/// Log-density of the Multivariate Normal distribution.
///
/// Evaluate the log-density of `N` observations of a `q` dimensional Multivariate Normal distribution.  Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] X Matrix of `q x n` random variable observations.
/// @param [in] Mu Matrix of `q x n` mean vectors.
/// @param [in] V Matrix of `q x nq` variance matrices.
///
/// @return Vector of `N` log-density evaluations.
//[[Rcpp::export]]
Eigen::VectorXd LogDensityMultivariateNormal(Eigen::MatrixXd X,
					     Eigen::MatrixXd Mu,
					     Eigen::MatrixXd V) {
  int q = X.rows();
  int N = X.cols();
  N = std::max<int>(N, Mu.cols());
  N = std::max<int>(N, V.cols()/q);
  // output variables
  VectorXd logDens(N);
  // internal variables
  // VectorXd z(q);
  LLT<MatrixXd> cholV(q);
  double ldV;
  MultiNormal mnorm(q);
  bool singleX = X.cols() == 1;
  bool singleV = V.cols() == q;
  bool singleMu = Mu.cols() == 1;
  if(singleV) {
    cholV.compute(V);
    ldV = logDetCholV(cholV);
    // ldV = 0.0;
    // for(int ii=0; ii<q; ii++) {
    //   ldV += log(cholV.matrixL()(ii,ii));
    // }
  }
  for(int ii=0; ii<N; ii++) {
    if(!singleV) {
      cholV.compute(V.block(0,ii*q,q,q));
      ldV = logDetCholV(cholV);
    }
    logDens(ii) = mnorm.LogDens(X.col(ii*(!singleX)),
				Mu.col(ii*(!singleMu)), cholV, ldV);
    // logDens(ii) = LogDensNormal(X.col(ii*(!singleX)),
    // 				Mu.col(ii*(!singleMu)), V.block(0,ii*(!singleV)*q,q,q),
    // 				!singleV,
    // 				z, ldV, cholV);
  }
  return logDens;
}

/// Generate a random sample from the Multivariate Normal distribution.
///
/// Generate `N` independent draws from a `q` dimensional Multivariate Normal distribution.  Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] N Integer number of random draws.
/// @param [in] Lambda Matrix of size `q x n` of mean vectors.
/// @param [in] Sigma Matrix of size `q x nq` of variance matrices.
///
/// @return Matrix of size `q x N` of random draws.
//[[Rcpp::export]]
Eigen::MatrixXd GenerateMultivariateNormal(int N, Eigen::MatrixXd Lambda,
					   Eigen::MatrixXd Sigma) {
  int q = Lambda.rows();
  bool singleLambda = (Lambda.cols() == 1);
  bool singleSigma = (Sigma.cols() == q);
  int ii,jj;
  // output variables
  MatrixXd X(q,N);
  // internal variables
  LLT<MatrixXd> lltq(q);
  MatrixXd SigmaL(q,q);
  VectorXd lambda(q);
  MatrixXd Z(q,N);
  // generate iid normals
  for(ii=0; ii<N; ii++) {
    for(jj=0; jj<q; jj++) {
      Z(jj,ii) = norm_rand();
    }
  }
  // multiply by Sigma^{1/2}
  if(singleSigma) {
    lltq.compute(Sigma);
    SigmaL = lltq.matrixL();
    triMultLX(X, SigmaL, Z);
  } else {
    for(ii=0; ii<N; ii++) {
      lltq.compute(Sigma.block(0,ii*q,q,q));
      SigmaL = lltq.matrixL();
      triMultLX(X.col(ii), SigmaL, Z.col(ii));
    }
  }
  // add Lambda
  if(singleLambda) {
    //lambda = Lambda;
    X.colwise() += Lambda.col(0);
  } else {
    X += Lambda;
  }
  return X;
}
