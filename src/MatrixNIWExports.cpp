///////////////////////////////////////////////////////////////////

// Exported Matrix Normal Inverse Wishart functions

//////////////////////////////////////////////////////////////////

// #include <Rcpp.h>
// using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
//#include <iostream>
#include "TriUtils.h"
#include "Wishart.h"
#include "MatrixNormal.h"
// #include "mniwWishart.h"
// #include "mniwMatNorm.h"

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
  // TempPQ *tmp = new TempPQ(p,q);
  LLT<MatrixXd> lltq(q);
  MatrixXd Lq = MatrixXd::Zero(q,q);
  MatrixXd Uq = MatrixXd::Zero(q,q);
  MatrixXd Iq = MatrixXd::Identity(q,q);
  MatrixXd CL = MatrixXd::Zero(q,q);
  MatrixXd SigmaL = MatrixXd::Zero(p,p);
  MatrixXd OmegaU = MatrixXd::Zero(p,p);
  MatrixXd XiL = MatrixXd::Zero(q,q);
  LLT<MatrixXd> lltp(p);
  Wishart wish(q);
  MatrixNormal matnorm(p,q);
  if(singlePsi) {
    ReverseCholesky(XiL, Psi, lltq);
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
      ReverseCholesky(XiL, Psi.block(0,ii*q,q,q), lltq);
    }
    wish.GenerateLowerTriXi(CL, XiL, nu(ii*(!singleNu)));
    if(!singleSigma) {
      lltp.compute(Sigma.block(0,ii*p,p,p));
    }
    if(!inverse) {
      if(!singleSigma) {
	SigmaL = lltp.matrixL();
      }
      matnorm.GenerateRowSColO(X.block(0,ii*q,p,q),
			       Lambda.block(0,ii*q*(!singleLambda),p,q),
			       SigmaL, CL);
      // GenerateMatrixNormalRowSColO(X.block(0,ii*q,p,q),
      // 				   Lambda.block(0,ii*q*(!singleLambda),p,q),
      // 				   SigmaL, CL, tmp->Xpq);
    }
    else {
      if(!singleSigma) {
	OmegaU = lltp.matrixU();
      }
      matnorm.GenerateRowOColO(X.block(0,ii*q,p,q),
			       Lambda.block(0,ii*q*(!singleLambda),p,q),
			       OmegaU, CL);
      // GenerateMatrixNormalRowOColO(X.block(0,ii*q,p,q),
      // 				   Lambda.block(0,ii*q*(!singleLambda),p,q),
      // 				   OmegaU, CL, tmp->Xpq);
    }
    InverseLLt(V.block(0,ii*q,q,q), CL, Lq, Uq, Iq);
  }
  // delete tmp;
  return List::create(_["X"] = wrap(X),
		      _["V"] = wrap(V));
}
