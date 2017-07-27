///////////////////////////////////////////////////////////////////

// Cross-Products

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
#include "mniwUtils.h"

//////////////////////////////////////////////////////////////////

// inner products t(X) %*% V %*% X for possibly multiple X and V
// can specify inverse = TRUE which takes V <- V^{-1}
// each X is (p x q), such that each V is (p x p)
//[[Rcpp::export]]
Eigen::MatrixXd CrossProdVXX(Eigen::MatrixXd X, Eigen::MatrixXd V,
			     int p, int q, bool inverse = false) {
  bool singleX = (X.cols() == q);
  bool singleV = (V.cols() == p);
  int N = singleX ? V.cols()/p : X.cols()/q;
  // output variables
  MatrixXd Y(q,q*N);
  // internal variables
  LLT<MatrixXd> lltp(p);
  MatrixXd VX(p,q);
  MatrixXd Xt(q,p);
  if(singleX) {
    Xt = X.adjoint();
  }
  if(inverse & singleV) {
    lltp.compute(V);
  }
  for(int ii=0; ii<N; ii++) {
    if(!inverse) {
      VX.noalias() = V.block(0,ii*p*(!singleV),p,p) * X.block(0,ii*q*(!singleX),p,q);
    }
    else {
      if(!singleV) {
	lltp.compute(V.block(0,ii*p,p,p));
      }
      VX = lltp.solve(X.block(0,ii*q*(!singleX),p,q));
    }
    if(singleX) {
      Y.block(0,ii*q,q,q).noalias() = Xt * VX;
    }
    else {
      Y.block(0,ii*q,q,q).noalias() = X.block(0,ii*q,p,q).adjoint() * VX;
    }
  }
  return(Y);
}

// inner products t(X) %*% V %*% Y for possibly multiple X, Y, and V
// can specify inverse = TRUE which takes V <- V^{-1}
// each X is (p x q) and Y is (p x r), such that each V is (p x p)
// each output block is (q x r)
//[[Rcpp::export]]
Eigen::MatrixXd CrossProdVXY(Eigen::MatrixXd X, Eigen::MatrixXd Y,
			     Eigen::MatrixXd V,
			     int p, int q, int r, bool inverse = false) {
  bool singleX = (X.cols() == q);
  bool singleY = (Y.cols() == r);
  bool singleV = (V.cols() == p);
  int N = singleX ? V.cols()/p : X.cols()/q;
  N = singleY ? N : Y.cols()/r;
  // output variables
  MatrixXd W(q,r*N);
  // internal variables
  LLT<MatrixXd> lltp(p);
  MatrixXd VY(p,r);
  MatrixXd Xt(q,p);
  // if both X,V or both Y,V are single, compute that product once
  if(inverse & singleV) {
    lltp.compute(V);
  }
  if(singleX) {
    if(singleV) {
      // X'V or X'V^{-1}
      if(!inverse) {
	Xt.noalias() = X.adjoint() * V;
      }
      else {
	Xt = lltp.solve(X).adjoint();
      }
    } else {
      // just transpose X
      Xt = X.adjoint();
    }
  } else if (singleY & singleV) {
    // VY or V^{-1}Y
    if(!inverse) {
      VY.noalias() = V * Y;
    }
    else {
      VY = lltp.solve(Y);
    }
  }
  //std::cout << "t(X) = " << Xt << std::endl;
  //std::cout << "V = " << V << std::endl;
  //std::cout << "Y = " << Y << std::endl;
  for(int ii=0; ii<N; ii++) {
    if(singleV & singleX) {
      // single X'V
      W.block(0,ii*r,q,r).noalias() = Xt * Y.block(0,ii*r*(!singleY),p,r);
    } else if(singleV & singleY) {
      // single VY
      W.block(0,ii*r,q,r).noalias() = X.block(0,ii*q*(!singleX),p,q).adjoint() * VY;
    } else {
      // V is not single, or neither X nor Y are single
      if(!inverse) {
	VY.noalias() = V.block(0,ii*p*(!singleV),p,p) * Y.block(0,ii*r*(!singleY),p,r);
      } else {
	if(!singleV) {
	  lltp.compute(V.block(0,ii*p,p,p));
	}
	VY = lltp.solve(Y.block(0,ii*r*(!singleY),p,r));
      }
      if(singleX) {
	W.block(0,ii*r,q,r).noalias() = Xt * VY;
      } else {
	W.block(0,ii*r,q,r).noalias() = X.block(0,ii*q,p,q).adjoint() * VY;
      }
    }
  }
  return(W);
}
