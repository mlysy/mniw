///////////////////////////////////////////////////////////////////

// Density and Simulation of the Matrix Normal Distribution

#ifndef mniwMatNorm_h
#define mniwMatNorm_h 1

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

// log-densities

// multivariate normal
double LogDensNormal(const Ref<VectorXd>& x, const Ref<VectorXd>& mu,
		     const Ref<MatrixXd>& V, bool CalcV,
		     Ref<VectorXd> z,
		     double &ldV, LLT<MatrixXd> cholV);

// matrix normal
double LogDensMatrixNormal(const Ref<const MatrixXd>& X,
			   const Ref<const MatrixXd>& Mu,
			   const Ref<const MatrixXd>& RowV, bool CalcRowV,
			   const Ref<const MatrixXd>& ColV, bool CalcColV,
			   Ref<MatrixXd> Z,
			   double& ldRowV, LLT<MatrixXd> cholRowV,
			   double& ldColV, LLT<MatrixXd> cholColV);

///////////////////////////////////////////////////////////////////

// simulation of matrix normal distribution

// both variances on the covariance scale.
void GenerateMatrixNormalRowSColS(Ref<MatrixXd> X,
				  const Ref<const MatrixXd>& Lambda,
				  const Ref<const MatrixXd>& RowSigmaL,
				  const Ref<const MatrixXd>& ColSigmaU,
				  Ref<MatrixXd> Z);

// between-row on covariance scale, between-col on precision scale.
void GenerateMatrixNormalRowSColO(Ref<MatrixXd> X,
				  const Ref<const MatrixXd>& Lambda,
				  const Ref<const MatrixXd>& RowSigmaL,
				  const Ref<const MatrixXd>& ColOmegaL,
				  Ref<MatrixXd> Z);

// both variances are on the precision scale.
void GenerateMatrixNormalRowOColO(Ref<MatrixXd> X,
				  const Ref<MatrixXd>& Lambda,
				  const Ref<const MatrixXd>& RowOmegaU,
				  const Ref<const MatrixXd>& ColOmegaL,
				  Ref<MatrixXd> Z);

#endif
