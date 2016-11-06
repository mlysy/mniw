///////////////////////////////////////////////////////////////////

// Exported Wishart functions

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

//#include <iostream>
using namespace Eigen;
#include "mniwWishart.h"

//////////////////////////////////////////////////////////////////

// log-density of Wishart or inverse-Wishart:
// X ~ Wish(Psi, nu) or X ~ iWish(Psi, nu)
// can specify multiple Psi and nu
// in an R version the transposition will be done automatically
//[[Rcpp::export]]
Eigen::VectorXd LogDensityWishart(Eigen::MatrixXd X, Eigen::MatrixXd Psi,
				  Eigen::VectorXd nu, bool inverse = false) {
  int q = X.rows();
  int N = X.cols()/q;
  N = std::max<int>(N, Psi.cols()/q);
  N = std::max<int>(N, nu.size());
  // output variables
  VectorXd logDens(N);
  // internal variables
  MatrixXd Z(q,q);
  LLT<MatrixXd> cholX(q);
  LLT<MatrixXd> cholPsi(q);
  double ldPsi;
  double ldX;
  bool singleX = X.cols() == q;
  bool singlePsi = Psi.cols() == q;
  bool singleNu = nu.size() == 1;
  if(singleX) {
    cholX.compute(X);
    ldX = 0.0;
    for(int ii=0; ii<q; ii++) {
      ldX += log(cholX.matrixL()(ii,ii));
    }
  }
  if(singlePsi) {
    cholPsi.compute(Psi);
    ldPsi = 0.0;
    for(int ii=0; ii<q; ii++) {
      ldPsi += log(cholPsi.matrixL()(ii,ii));
    }
  }
  for(int ii=0; ii<N; ii++) {
    logDens(ii) = LogDensWishart(X.block(0,ii*(!singleX)*q,q,q),
				 Psi.block(0,ii*(!singlePsi)*q,q,q),
				 nu(ii*(!singleNu)), inverse, !singleX, !singlePsi,
				 Z, ldX, cholX, ldPsi, cholPsi);
  }
  return logDens;
}



// generate wishart or inverse wishart
// V ~ Wish(Psi, nu) or V ~ iWish(Psi, nu).
// supports multiple parameter arguments
//[[Rcpp::export]]
Eigen::MatrixXd GenerateWishart(int N, Eigen::MatrixXd Psi, Eigen::VectorXd nu,
				bool inverse = false) {
  int q = Psi.rows();
  bool singlePsi = (Psi.cols() == q);
  bool singleNu = (nu.size() == 1);
  // output variables
  MatrixXd X(q,N*q);
  // internal variables
  TempPQ *tmp = new TempPQ(1,q);
  MatrixXd PsiL = MatrixXd::Zero(q,q);
  MatrixXd VL = MatrixXd::Zero(q,q);
  // wishart lower triangular factor
  if(singlePsi) {
    if(!inverse) {
      tmp->lltq.compute(Psi);
      PsiL = tmp->lltq.matrixL();
    }
    else {
      ReverseCholesky(PsiL, Psi, tmp);
    }
  }
  // for-loop
  for(int ii=0; ii<N; ii++) {
    if(!inverse) {
      if(!singlePsi) {
	tmp->lltq.compute(Psi.block(0,ii*q,q,q));
	PsiL = tmp->lltq.matrixL();
      }
      GenerateWishartLowerTri(VL, PsiL, nu(ii*(!singleNu)), tmp->Lq);
      CrossProdLLt(X.block(0,ii*q,q,q), VL, tmp->Uq);
    }
    else {
      if(!singlePsi) {
	ReverseCholesky(PsiL, Psi.block(0,ii*q,q,q), tmp);
      }
      GenerateWishartLowerTriXi(VL, PsiL, nu(ii*(!singleNu)));
      InverseLLt(X.block(0,ii*q,q,q), VL, tmp->Lq, tmp->Uq, tmp->Iq);
    }
  }
  delete tmp;
  return X;
}
