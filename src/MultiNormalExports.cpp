///////////////////////////////////////////////////////////////////

// Exported Multivariate Normal functions

//////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
//#include <iostream>
#include "MultiNormal.h"
// #include "mniwMatNorm.h"


//////////////////////////////////////////////////////////////////

// log-density of multivariate normal
// X ~ N(Mu, V)
// note that Mu and V can be different for each X
// an R version of this wrapper will do transpositions automatically
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
