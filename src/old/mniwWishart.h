///////////////////////////////////////////////////////////////////

// Density and Simulation of the Wishart and Inverse Wishart distributions

#ifndef mniwWishart_h
#define mniwWishart_h

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

///////////////////////////////////////////////////////////////////

// multivariate log-gamma function
double logMultiGamma(double alpha, int q);

// Density of Wishart or Inverse Wishart
double LogDensWishart(const Ref<MatrixXd>& X, const Ref<MatrixXd>& Psi,
		      double nu, bool inv, bool CalcX, bool CalcPsi,
		      Ref<MatrixXd> Z, double &ldX, LLT<MatrixXd> cholX,
		      double &ldPsi, LLT<MatrixXd> cholPsi);

// simulation of lower triangular Bartlett decomposition of Wisharts
void GenerateWishartLowerTri(Ref<MatrixXd> VL,
			     const Ref<MatrixXd>& PsiL,
			     double nu,
			     Ref<MatrixXd> XL);

// same thing but with XiL = PsiL^{-1}
void GenerateWishartLowerTriXi(Ref<MatrixXd> VL,
			       const Ref<const MatrixXd>& XiL,
			       double nu);

#endif
