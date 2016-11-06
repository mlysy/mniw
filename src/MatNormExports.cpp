///////////////////////////////////////////////////////////////////

// Exported Matrix Normal functions

//////////////////////////////////////////////////////////////////

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

//#include "CoreLibs.h"
//using namespace Rcpp;
//using namespace Eigen;
//#include <iostream>
#include "mniwMatNorm.h"

//////////////////////////////////////////////////////////////////

// log-density of multivariate normal
// X ~ N(Mu, V)
// note that Mu and V can be different for each X
// an R version of this wrapper will do transpositions automatically
//[[Rcpp::export]]
Eigen::VectorXd LogDensityMultivariateNormal(Eigen::MatrixXd X,
					     Eigen::MatrixXd Mu, Eigen::MatrixXd V) {
  int q = X.rows();
  int N = X.cols();
  N = std::max<int>(N, Mu.cols());
  N = std::max<int>(N, V.cols()/q);
  // output variables
  VectorXd logDens(N);
  // internal variables
  VectorXd z(q);
  LLT<MatrixXd> cholV(q);
  double ldV;
  bool singleX = X.cols() == 1;
  bool singleV = V.cols() == q;
  bool singleMu = Mu.cols() == 1;
  if(singleV) {
    cholV.compute(V);
    ldV = 0.0;
    for(int ii=0; ii<q; ii++) {
      ldV += log(cholV.matrixL()(ii,ii));
    }
  }
  for(int ii=0; ii<N; ii++) {
    logDens(ii) = LogDensNormal(X.col(ii*(!singleX)),
				Mu.col(ii*(!singleMu)), V.block(0,ii*(!singleV)*q,q,q),
				!singleV,
				z, ldV, cholV);
  }
  return logDens;
}

// log-density of matrix normal
// X ~ N(Mu, RowV, ColV)
// note that Mu, RowV, ColV can be different for each X
// an R version of this wrapper will do transpositions automatically
//[[Rcpp::export]]
Eigen::VectorXd LogDensityMatrixNormal(Eigen::MatrixXd X, Eigen::MatrixXd Mu,
				       Eigen::MatrixXd RowV, Eigen::MatrixXd ColV) {
  // dimensions of the problem
  int q = ColV.rows();
  int p = RowV.rows();
  int N = X.cols()/q;
  N = std::max<int>(N, Mu.cols()/q);
  N = std::max<int>(N, RowV.cols()/p);
  N = std::max<int>(N, ColV.cols()/q);
  // output variables
  VectorXd logDens(N);
  // internal variables
  int ii;
  MatrixXd Z(p,q);
  LLT<MatrixXd> cholRowV(p);
  LLT<MatrixXd> cholColV(q);
  double ldRowV, ldColV;
  bool singleX = X.cols() == q;
  bool singleMu = Mu.cols() == q;
  bool singleRowV = RowV.cols() == p;
  bool singleColV = ColV.cols() == q;
  if(singleRowV) {
    cholRowV.compute(RowV);
    ldRowV = 0.0;
    for(ii=0; ii<p; ii++) {
      ldRowV += log(cholRowV.matrixL()(ii,ii));
    }
  }
  if(singleColV) {
    cholColV.compute(ColV);
    ldColV = 0.0;
    for(ii=0; ii<q; ii++) {
      ldColV += log(cholColV.matrixL()(ii,ii));
    }
  }
  for(ii=0; ii<N; ii++) {
    logDens(ii) = LogDensMatrixNormal(X.block(0,q*ii*(!singleX),p,q),
				      Mu.block(0,q*ii*(!singleMu),p,q),
				      RowV.block(0,p*ii*(!singleRowV),p,p), !singleRowV,
				      ColV.block(0,q*ii*(!singleColV),q,q), !singleColV,
				      Z, ldRowV, cholRowV, ldColV, cholColV);
  }
  return logDens;
}

// Generate Multivariate Normal
// supports multple instance sof each parameter
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

// Generate Matrix Normal
// X ~ MNorm(Lambda, RowSigma, ColSigma)
// supports multiple instances of each parameter
//[[Rcpp::export]]
Eigen::MatrixXd GenerateMatrixNormal(int N, Eigen::MatrixXd Lambda,
				     Eigen::MatrixXd RowSigma,
				     Eigen::MatrixXd ColSigma) {
  int p = RowSigma.rows();
  int q = ColSigma.rows();
  bool singleLambda = (Lambda.cols() == q);
  bool singleRowSigma = (RowSigma.cols() == p);
  bool singleColSigma = (ColSigma.cols() == q);
  int ii;
  // output variables
  MatrixXd X(p,N*q);
  // internal variables
  TempPQ *tmp = new TempPQ(p,q);
  MatrixXd RowSigmaL = MatrixXd::Zero(p,p);
  MatrixXd ColSigmaU = MatrixXd::Zero(q,q);
  LLT<MatrixXd> lltp(p);
  if(singleRowSigma) {
    lltp.compute(RowSigma);
    RowSigmaL = lltp.matrixL();
  }
  if(singleColSigma) {
    tmp->lltq.compute(ColSigma);
    ColSigmaU = tmp->lltq.matrixU();
  }
  for(ii=0; ii<N; ii++) {
    if(!singleRowSigma) {
      lltp.compute(RowSigma.block(0,ii*p,p,p));
      RowSigmaL = lltp.matrixL();
    }
    if(!singleColSigma) {
      tmp->lltq.compute(ColSigma.block(0,ii*q,q,q));
      ColSigmaU = tmp->lltq.matrixU();
    }
    GenerateMatrixNormalRowSColS(X.block(0,ii*q,p,q),
				 Lambda.block(0,ii*q*(!singleLambda),p,q),
				 RowSigmaL, ColSigmaU, tmp->Xpq);
  }
  delete tmp;
  return X;
}
