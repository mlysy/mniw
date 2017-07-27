///////////////////////////////////////////////////////////////////

// Density and Simulation of the Wishart and Inverse Wishart distributions

#include "mniwSetLib.h"

#ifdef R_FUNCTION_LIBRARY
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#define gammaln R::lgammafn
#define chisq_rand Rf_rchisq
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

///////////////////////////////////////////////////////////////////

// multivariate log-gamma function
double logMultiGamma(double alpha, int q) {
  double lmg = 0.0;
  for(int ii=0; ii<q; ii++) {
    lmg += gammaln(alpha - .5*ii);
  }
  lmg += .5*q*(q-1) * M_LN_SQRT_PI;
  return lmg;
}

///////////////////////////////////////////////////////////////////

// Density of Wishart or Inverse Wishart
// CalcPsi determines whether to calculate cholPsi and ldPsi
// similarly for CalcX
double LogDensWishart(const Ref<MatrixXd>& X, const Ref<MatrixXd>& Psi,
		      double nu, bool inv, bool CalcX, bool CalcPsi,
		      Ref<MatrixXd> Z, double &ldX, LLT<MatrixXd> cholX,
		      double &ldPsi, LLT<MatrixXd> cholPsi) {
  double ldens;
  int q = X.cols();
  if(CalcX) {
    cholX.compute(X);
  }
  if(CalcPsi) {
    cholPsi.compute(Psi);
  }
  // trace
  if(!inv) {
    Z.noalias() = cholPsi.solve(X);
  }
  else {
    Z.noalias() = cholX.solve(Psi);
  }
  // determinants
  if(CalcX) {
    ldX = 0.0;
  }
  if(CalcPsi) {
    ldPsi = 0.0;
  }
  for(int ii=0; ii<q; ii++) {
    if(CalcX) {
      ldX += log(cholX.matrixL()(ii,ii));
    }
    if(CalcPsi) {
      ldPsi += log(cholPsi.matrixL()(ii,ii));
    }
  }
  if(!inv) {
    ldens = (nu-q-1) * ldX - nu * ldPsi;
  }
  else {
    ldens = -(nu+q+1) * ldX + nu * ldPsi;
  }
  ldens -= .5*(Z.trace() + q*nu * M_LN2) +  + logMultiGamma(.5*nu, q);
  return ldens;
}

// simulation of lower triangular Bartlett decomposition of Wisharts
// VL: output
// PsiL: precision matrix
// nu: degrees of freedom
// XL: scratch space
// testing: no typedef
// How is there no tri * tri multiplication???
void GenerateWishartLowerTri(Ref<MatrixXd> VL,
			     const Ref<MatrixXd>& PsiL,
			     double nu,
			     Ref<MatrixXd> XL) {
  int q = PsiL.cols();
  // fill XL with normals and chi-squared's
  int ii, jj;
  for(ii=0; ii<q; ii++) {
    // diagonal
    XL(ii,ii) = sqrt(chisq_rand(nu-ii));
    // off-diagonals
    for(jj=0; jj<ii; jj++) {
      XL(ii,jj) = norm_rand();
    }
  }
  // multiply by precision matrix VL = PsiL * XL
  triMultLL(VL, PsiL, XL);
  return;
}

// same thing but with XiL = PsiL^{-1}
void GenerateWishartLowerTriXi(Ref<MatrixXd> VL,
			       const Ref<const MatrixXd>& XiL,
			       double nu) {
  int q = VL.cols();
  // fill VL with normals and chi-squared's
  int ii, jj;
  for(ii=0; ii<q; ii++) {
    // diagonal
    VL(ii,ii) = sqrt(chisq_rand(nu-ii));
    // off-diagonals
    for(jj=0; jj<ii; jj++) {
      VL(ii,jj) = norm_rand();
    }
  }
  // multiply by precision matrix
  triMultLiX(VL, XiL);
  return;
}
