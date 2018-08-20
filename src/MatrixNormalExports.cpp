///////////////////////////////////////////////////////////////////

// Exported Matrix Normal functions

//////////////////////////////////////////////////////////////////

// #include <Rcpp.h>
// using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
//#include <iostream>
#include "MatrixNormal.h"
// #include "mniwMatNorm.h"

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
  MatrixNormal matnorm(p,q);
  if(singleRowV) {
    cholRowV.compute(RowV);
    ldRowV = logDetCholV(cholRowV);
    // ldRowV = 0.0;
    // for(ii=0; ii<p; ii++) {
    //   ldRowV += log(cholRowV.matrixL()(ii,ii));
    // }
  }
  if(singleColV) {
    cholColV.compute(ColV);
    ldColV = logDetCholV(cholColV);
    // ldColV = 0.0;
    // for(ii=0; ii<q; ii++) {
    //   ldColV += log(cholColV.matrixL()(ii,ii));
    // }
  }
  for(ii=0; ii<N; ii++) {
    if(!singleRowV) {
      cholRowV.compute(RowV.block(0,p*ii,p,p));
      ldRowV = logDetCholV(cholRowV);
    }
    if(!singleColV) {
      cholColV.compute(ColV.block(0,q*ii,q,q));
      ldColV = logDetCholV(cholColV);
    }
    logDens(ii) = matnorm.LogDens(X.block(0,q*ii*(!singleX),p,q),
				  Mu.block(0,q*ii*(!singleMu),p,q),
				  cholRowV, ldRowV, cholColV, ldColV);
    // logDens(ii) = LogDensMatrixNormal(X.block(0,q*ii*(!singleX),p,q),
    // 				      Mu.block(0,q*ii*(!singleMu),p,q),
    // 				      RowV.block(0,p*ii*(!singleRowV),p,p), !singleRowV,
    // 				      ColV.block(0,q*ii*(!singleColV),q,q), !singleColV,
    // 				      Z, ldRowV, cholRowV, ldColV, cholColV);
  }
  return logDens;
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
  // TempPQ *tmp = new TempPQ(p,q);
  MatrixXd RowSigmaL = MatrixXd::Zero(p,p);
  MatrixXd ColSigmaU = MatrixXd::Zero(q,q);
  LLT<MatrixXd> lltq(q);
  LLT<MatrixXd> lltp(p);
  MatrixNormal matnorm(p,q);
  if(singleRowSigma) {
    lltp.compute(RowSigma);
    RowSigmaL = lltp.matrixL();
  }
  if(singleColSigma) {
    // tmp->lltq.compute(ColSigma);
    // ColSigmaU = tmp->lltq.matrixU();
    lltq.compute(ColSigma);
    ColSigmaU = lltq.matrixU();
  }
  for(ii=0; ii<N; ii++) {
    if(!singleRowSigma) {
      lltp.compute(RowSigma.block(0,ii*p,p,p));
      RowSigmaL = lltp.matrixL();
    }
    if(!singleColSigma) {
      // tmp->lltq.compute(ColSigma.block(0,ii*q,q,q));
      // ColSigmaU = tmp->lltq.matrixU();
      lltq.compute(ColSigma.block(0,ii*q,q,q));
      ColSigmaU = lltq.matrixU();
    }
    matnorm.GenerateRowSColS(X.block(0,ii*q,p,q),
			     Lambda.block(0,ii*q*(!singleLambda),p,q),
			     RowSigmaL, ColSigmaU);
    // GenerateMatrixNormalRowSColS(X.block(0,ii*q,p,q),
    // 				 Lambda.block(0,ii*q*(!singleLambda),p,q),
    // 				 RowSigmaL, ColSigmaU, tmp->Xpq);
  }
  // delete tmp;
  return X;
}
