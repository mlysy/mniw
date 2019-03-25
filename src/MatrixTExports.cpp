/// @file MatrixTExports.cpp
///
/// Exported Rcpp functions for the Matrix-T distribution.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
#include "MatrixT.h"
using namespace mniw;

/// Log-density of the Matrix-T distribution.
///
/// Evaluate the log-density of `N` observations of a `p x q` dimensional Matrix-T distribution.  Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] X Matrix of `p x nq` random matrix observations.
/// @param [in] Mu Matrix of `p x nq` mean matrices.
/// @param [in] RowV Matrix of `p x np` row-wise variance matrices.
/// @param [in] ColV Matrix of `q x nq` column-wise variance matrices.
/// @param [in] nu Vector of `n` degrees-of-freedom parameters.
///
/// @return Vector of `N` log-density evaluations.
//[[Rcpp::export]]
Eigen::VectorXd LogDensityMatrixT(Eigen::MatrixXd X,
				  Eigen::MatrixXd Mu,
				  Eigen::MatrixXd RowV,
				  Eigen::MatrixXd ColV,
				  Eigen::VectorXd nu) {
  // dimensions of the problem
  int p = RowV.rows();
  int q = ColV.rows();
  int N = X.cols()/q;
  N = std::max<int>(N, Mu.cols()/q);
  N = std::max<int>(N, RowV.cols()/p);
  N = std::max<int>(N, ColV.cols()/q);
  N = std::max<int>(N, nu.size());
  // output variables
  VectorXd logDens(N);
  // internal variables
  LLT<MatrixXd> cholRowV(p);
  LLT<MatrixXd> cholColV(q);
  double ldRowV, ldColV;
  bool singleX = (X.cols() == q);
  bool singleMu = (Mu.cols() == q);
  bool singleRowV = (RowV.cols() == p);
  bool singleColV = (ColV.cols() == q);
  bool singleNu = (nu.size() == 1);
  MatrixT matT(p,q);
  // precompute Cholesky factors and log-determinants
  if(singleRowV) {
    cholRowV.compute(RowV);
    ldRowV = logDetCholV(cholRowV);
  }
  if(singleColV) {
    cholColV.compute(ColV);
    ldColV = logDetCholV(cholColV);
  }
  // main loop
  for(int ii=0; ii<N; ii++) {
    // compute Chol/log-det if not precomputed
    if(!singleRowV) {
      cholRowV.compute(RowV.block(0,p*ii,p,p));
      ldRowV = logDetCholV(cholRowV);
    }
    if(!singleColV) {
      cholColV.compute(ColV.block(0,q*ii,q,q));
      ldColV = logDetCholV(cholColV);
    }
    // density evaluation
    logDens(ii) = matT.LogDens(X.block(0,q*ii*(!singleX),p,q),
			       Mu.block(0,q*ii*(!singleMu),p,q),
			       RowV.block(0,p*ii*(!singleRowV),p,p),
			       cholRowV, ldRowV,
			       ColV.block(0,q*ii*(!singleColV),q,q),
			       cholColV, ldColV,
			       nu(ii*(!singleNu)));
  }
  return logDens;
}

/// Generate a random sample from the Matrix-T distribution.
///
/// Generate `N` independent draws from a `p x q` dimensional Matrix-T distribution.  Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] N Integer number of random draws
/// @param [in] Lambda Matrix of `p x nq` mean matrices.
/// @param [in] RowSigma Matrix of `p x np` row-wise variance or precision matrices.
/// @param [in] ColSigma Matrix of `q x nq` column-wise variance matrices.
/// @param [in] nu Vector of `n` degrees-of-freedom parameters.
/// @param [in] inverse Boolean; `true/false` indicates that `RowSigma` is on the precision/variance scale.
/// @return Matrix of `p x Nq` random draws.
// [[Rcpp::export]]
Eigen::MatrixXd GenerateMatrixT(int N, Eigen::MatrixXd Lambda,
				Eigen::MatrixXd RowSigma,
				Eigen::MatrixXd ColSigma,
				Eigen::VectorXd nu,
				bool inverse = false) {
  // problem dimensions
  int p = RowSigma.rows();
  int q = ColSigma.rows();
  bool singleLambda = (Lambda.cols() == q);
  bool singleRowSigma = (RowSigma.cols() == p);
  bool singleColSigma = (ColSigma.cols() == q);
  bool singleNu = (nu.size() == 1);
  // internal variables
  LLT<MatrixXd> lltq(q);
  LLT<MatrixXd> lltp(p);
  MatrixXd RowSigmaL = MatrixXd::Zero(p,p);
  MatrixXd RowISigmaU = MatrixXd::Zero(p,p);
  // MatrixXd CL = MatrixXd::Zero(q,q);
  MatrixXd ColISigmaL = MatrixXd::Zero(q,q);
  MatrixT matT(p,q);
  // output variables
  MatrixXd X(p,N*q);
  // precomputations
  if(singleColSigma) {
    ReverseCholesky(ColISigmaL, ColSigma, lltq);
  }
  if(singleRowSigma) {
    lltp.compute(RowSigma);
    if(!inverse) {
      RowSigmaL = lltp.matrixL();
    }
    else {
      RowISigmaU = lltp.matrixU();
    }
  }
  // main loop
  for(int ii=0; ii<N; ii++) {
    if(!singleColSigma) {
      ReverseCholesky(ColISigmaL, ColSigma.block(0,ii*q,q,q), lltq);
    }
    // wish.GenerateLowerTriXi(CL, XiL, nu(ii*(!singleNu)));
    if(!singleRowSigma) {
      lltp.compute(RowSigma.block(0,ii*p,p,p));
    }
    if(!inverse) {
      if(!singleRowSigma) {
	RowSigmaL = lltp.matrixL();
      }
      matT.GenerateRowSColO(X.block(0,ii*q,p,q),
			       Lambda.block(0,ii*q*(!singleLambda),p,q),
			    RowSigmaL, ColISigmaL, nu(ii*(!singleNu)));
      // matnorm.GenerateRowSColO(X.block(0,ii*q,p,q),
      // 			       Lambda.block(0,ii*q*(!singleLambda),p,q),
      // 			       RowSigmaL, CL);
    }
    else {
      if(!singleRowSigma) {
	RowISigmaU = lltp.matrixU();
      }
      matT.GenerateRowOColO(X.block(0,ii*q,p,q),
			    Lambda.block(0,ii*q*(!singleLambda),p,q),
			    RowISigmaU, ColISigmaL, nu(ii*(!singleNu)));
      // matnorm.GenerateRowOColO(X.block(0,ii*q,p,q),
      // 			       Lambda.block(0,ii*q*(!singleLambda),p,q),
      // 			       RowISigmaU, CL);
    }
  }
  return X;
}
