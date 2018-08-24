///////////////////////////////////////////////////////////////////

// Exported Matrix-T functions

//////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "MatrixT.h"
using namespace Eigen;

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
