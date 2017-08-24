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
inline double LogDensNormal(const Ref<VectorXd>& x, const Ref<VectorXd>& mu,
		     const Ref<MatrixXd>& V, bool CalcV,
		     Ref<VectorXd> z,
		     double &ldV, LLT<MatrixXd> cholV);

// matrix normal
inline double LogDensMatrixNormal(const Ref<const MatrixXd>& X,
			   const Ref<const MatrixXd>& Mu,
			   const Ref<const MatrixXd>& RowV, bool CalcRowV,
			   const Ref<const MatrixXd>& ColV, bool CalcColV,
			   Ref<MatrixXd> Z,
			   double& ldRowV, LLT<MatrixXd> cholRowV,
			   double& ldColV, LLT<MatrixXd> cholColV);

///////////////////////////////////////////////////////////////////

// simulation of matrix normal distribution

// both variances on the covariance scale.
inline void GenerateMatrixNormalRowSColS(Ref<MatrixXd> X,
				  const Ref<const MatrixXd>& Lambda,
				  const Ref<const MatrixXd>& RowSigmaL,
				  const Ref<const MatrixXd>& ColSigmaU,
				  Ref<MatrixXd> Z);

// between-row on covariance scale, between-col on precision scale.
inline void GenerateMatrixNormalRowSColO(Ref<MatrixXd> X,
				  const Ref<const MatrixXd>& Lambda,
				  const Ref<const MatrixXd>& RowSigmaL,
				  const Ref<const MatrixXd>& ColOmegaL,
				  Ref<MatrixXd> Z);

// both variances are on the precision scale.
inline void GenerateMatrixNormalRowOColO(Ref<MatrixXd> X,
				  const Ref<MatrixXd>& Lambda,
				  const Ref<const MatrixXd>& RowOmegaU,
				  const Ref<const MatrixXd>& ColOmegaL,
				  Ref<MatrixXd> Z);

///////////////////////////////////////////////////////////////////

// log-densities

// multivariate normal x ~ N(0, V)
// CalcV indicates whether to calculate the Cholesky factor and the log-determinant.
double LogDensNormal(const Ref<VectorXd>& x, const Ref<VectorXd>& mu,
                     const Ref<MatrixXd>& V, bool CalcV,
                     Ref<VectorXd> z,
                     double &ldV, LLT<MatrixXd> cholV) {
  //double ldens;
  double q = x.size();
  z = x-mu;
  if(CalcV) {
    cholV.compute(V);
    ldV = 0.0;
    for(int ii=0; ii<q; ii++) {
      ldV += log(cholV.matrixL()(ii,ii));
    }
  }
  cholV.matrixL().solveInPlace(z);
  return -(.5 * z.squaredNorm() + ldV + q * M_LN_SQRT_2PI);
}

// matrix normal X ~ N(0, RowV, ColV)
// calcRowV and calcColV indicate whether to calculate the row and column cholesky
// factors and determinants
double LogDensMatrixNormal(const Ref<const MatrixXd>& X,
                           const Ref<const MatrixXd>& Mu,
                           const Ref<const MatrixXd>& RowV, bool CalcRowV,
                           const Ref<const MatrixXd>& ColV, bool CalcColV,
                           Ref<MatrixXd> Z,
                           double& ldRowV, LLT<MatrixXd> cholRowV,
                           double& ldColV, LLT<MatrixXd> cholColV) {
  //double ldens;
  double p = X.rows();
  double q = X.cols();
  Z = X-Mu;
  if(CalcRowV) {
    cholRowV.compute(RowV);
    ldRowV = 0.0;
    for(int ii=0; ii<p; ii++) {
      ldRowV += log(cholRowV.matrixL()(ii,ii));
    }
  }
  if(CalcColV) {
    cholColV.compute(ColV);
    ldColV = 0.0;
    for(int ii=0; ii<q; ii++) {
      ldColV += log(cholColV.matrixL()(ii,ii));
    }
  }
  cholRowV.matrixL().solveInPlace(Z);
  cholColV.matrixU().solveInPlace<OnTheRight>(Z);
  return -(.5 * Z.squaredNorm() + q * ldRowV + p * ldColV + p*q * M_LN_SQRT_2PI);
}

///////////////////////////////////////////////////////////////////

// simulation of matrix normal distribution
// there are several variants of this function.
// the ending Row{S/O}Col{S/O} determines whether the between-row or between-column
// variance matrix should be provided on the variance (S) or precision (O) scale
// to be compatible with the MNIW distribution,

// both variances on the covariance scale.
// Lambda: mean matrix
// RowSigmaL: _lower_ triangular cholesky of the between-row _covariance_ matrix
// ColSigmaU: _upper_ triangular cholesky of the between-column _covariance_ matrix
void GenerateMatrixNormalRowSColS(Ref<MatrixXd> X,
                                  const Ref<const MatrixXd>& Lambda,
                                  const Ref<const MatrixXd>& RowSigmaL,
                                  const Ref<const MatrixXd>& ColSigmaU,
                                  Ref<MatrixXd> Z) {
  int ii, jj;
  int p = Lambda.rows();
  int q = Lambda.cols();
  // populate Z with iid normals
  for(ii=0; ii<p; ii++) {
    for(jj=0; jj<q; jj++) {
      Z(ii,jj) = norm_rand();
    }
  }
  // scaling
  triMultXU(X, Z, ColSigmaU);
  triMultLX(Z, RowSigmaL, X);
  X = Z + Lambda;
  return;
}

// between-row on covariance scale, between-col on precision scale.
// Lambda: mean matrix
// RowSigmaL: _lower_ triangular cholesky of the between-row _covariance_ matrix
// ColOmegaL: _lower_ triangular cholesky of the between-column _precision_ matrix
void GenerateMatrixNormalRowSColO(Ref<MatrixXd> X,
                                  const Ref<const MatrixXd>& Lambda,
                                  const Ref<const MatrixXd>& RowSigmaL,
                                  const Ref<const MatrixXd>& ColOmegaL,
                                  Ref<MatrixXd> Z) {
  int ii, jj;
  int p = Lambda.rows();
  int q = Lambda.cols();
  // populate Z with iid normals
  for(ii=0; ii<p; ii++) {
    for(jj=0; jj<q; jj++) {
      Z(ii,jj) = norm_rand();
    }
  }
  // scaling
  triMultXLi(Z, ColOmegaL);
  triMultLX(X, RowSigmaL, Z);
  X += Lambda;
  return;
}

// both variances are on the precision scale.
// Lambda: mean matrix
// RowOmegaU: _upper_ triangular cholesky of the between-row _precision_ matrix
// ColOmegaU: _lower_ triangular cholesky of the between-column _precision_ matrix
void GenerateMatrixNormalRowOColO(Ref<MatrixXd> X,
                                  const Ref<MatrixXd>& Lambda,
                                  const Ref<const MatrixXd>& RowOmegaU,
                                  const Ref<const MatrixXd>& ColOmegaL,
                                  Ref<MatrixXd> Z) {
  int ii, jj;
  int p = Lambda.rows();
  int q = Lambda.cols();
  // populate Z with iid normals
  for(ii=0; ii<p; ii++) {
    for(jj=0; jj<q; jj++) {
      Z(ii,jj) = norm_rand();
    }
  }
  // scaling
  triMultXLi(Z, ColOmegaL);
  triMultUiX(RowOmegaU, Z);
  X = Z + Lambda;
  return;
}

#endif
