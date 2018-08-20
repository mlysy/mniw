///////////////////////////////////////////////////////////////////

// Exported Random Effects functions

//////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////

// Generate Random-Effects Normal distribution
// p(mu | y), where
// y | mu ~ N(mu, V) and  mu ~ N(lambda, A)
// supports multiple instances of each parameter (including y)
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


