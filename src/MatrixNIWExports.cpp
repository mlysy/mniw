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
#include "TriUtils.h"
#include "mniwWishart.h"
#include "mniwMatNorm.h"

//////////////////////////////////////////////////////////////////

// Generate Matrix Normal-Inverse-Wishart
// V ~ iWish(Psi, nu)
// X | V ~ MNorm(Lambda, Sigma, V)
// allow for multiple instances of each parameter
// also, inverse = TRUE means that the inverses of Sigma are provided instead
//[[Rcpp::export]]
List GenerateMatrixNIW(int N,
		       Eigen::MatrixXd Lambda, Eigen::MatrixXd Sigma,
		       Eigen::MatrixXd Psi, Eigen::VectorXd nu,
		       bool inverse = false) {
  int p = Lambda.rows();
  int q = Psi.rows();
  bool singleLambda = (Lambda.cols() == q);
  bool singleSigma = (Sigma.cols() == p);
  bool singlePsi = (Psi.cols() == q);
  bool singleNu = (nu.size() == 1);
  int ii;
  // internal variables
  TempPQ *tmp = new TempPQ(p,q);
  MatrixXd CL = MatrixXd::Zero(q,q);
  MatrixXd SigmaL = MatrixXd::Zero(p,p);
  MatrixXd OmegaU = MatrixXd::Zero(p,p);
  MatrixXd XiL = MatrixXd::Zero(q,q);
  LLT<MatrixXd> lltp(p);
  if(singlePsi) {
    ReverseCholesky(XiL, Psi, tmp);
  }
  if(singleSigma) {
    lltp.compute(Sigma);
    if(!inverse) {
      SigmaL = lltp.matrixL();
    }
    else {
      OmegaU = lltp.matrixU();
    }
  }
  // output variables
  MatrixXd X(p,N*q);
  MatrixXd V(q,N*q);
  for(ii=0; ii<N; ii++) {
    if(!singlePsi) {
      ReverseCholesky(XiL, Psi.block(0,ii*q,q,q), tmp);
    }
    GenerateWishartLowerTriXi(CL, XiL, nu(ii*(!singleNu)));
    if(!singleSigma) {
      lltp.compute(Sigma.block(0,ii*p,p,p));
    }
    if(!inverse) {
      if(!singleSigma) {
	SigmaL = lltp.matrixL();
      }
      GenerateMatrixNormalRowSColO(X.block(0,ii*q,p,q),
				   Lambda.block(0,ii*q*(!singleLambda),p,q),
				   SigmaL, CL, tmp->Xpq);
    }
    else {
      if(!singleSigma) {
	OmegaU = lltp.matrixU();
      }
      GenerateMatrixNormalRowOColO(X.block(0,ii*q,p,q),
				   Lambda.block(0,ii*q*(!singleLambda),p,q),
				   OmegaU, CL, tmp->Xpq);
    }
    InverseLLt(V.block(0,ii*q,q,q), CL, tmp->Lq, tmp->Uq, tmp->Iq);
  }
  delete tmp;
  return List::create(_["X"] = wrap(X),
		      _["V"] = wrap(V));
}
