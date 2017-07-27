#include <iostream>
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

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

/***R

crossprodV <- function(X, Y = NULL, V, inverse = FALSE) {
  # dimensions of X and V
  if(is.vector(X)) X <- as.matrix(X)
  p <- dim(X)[1]
  q <- dim(X)[2]
  Xnames <- dimnames(X)[2]
  X <- matrix(X,p)
  if(!all(dim(V)[1:2] == p)) stop("X and V have incompatible dimensions.")
  V <- matrix(V,p)
  if(is.null(Y)) {
    # no Y
    r <- q
    Ynames <- Xnames
    # check lengths
    n <- unique(sort(c(ncol(X)/q, ncol(V)/p)))
    n <- c(1, n[n>1])
    if(length(n) > 2) stop("X and V have different lengths.")
    n <- max(n)
    # C code
    W <- CrossProdVXX(X, V, p, q, inverse)
  } else {
    # dimensions of Y
    if(is.vector(Y)) Y <- as.matrix(Y)
    if(dim(Y)[1] != p) stop("Y and V have incompatible dimensions.")
    Ynames <- dimnames(Y)[2]
    r <- dim(Y)[2]
    Y <- matrix(Y,p)
    # check lengths
    n <- unique(sort(c(ncol(X)/q, ncol(V)/p, ncol(Y)/r)))
    n <- c(1, n[n>1])
    if(length(n) > 2) stop("X, Y, and V have different lengths.")
    n <- max(n)
    # C code
    W <- CrossProdVXY(X, Y, V, p, q, r, inverse)
  }
  W <- array(W, dim = c(q,r,n))
  if(!is.null(Xnames) || !is.null(Ynames)) {
    dimnames(W) <- c(Xnames, Ynames, NULL)
  }
  W
}

*/
