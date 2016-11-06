///////////////////////////////////////////////////////////////////

// Random Effects Normal Models

#ifndef mniwRandomEffects_h
#define mniwRandomEffects_h 1

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

// Generate mu = [(I-B)*V]^{1/2} * z + B * x,
// OmegaL * OmegaU = Sigma^{-1}.
void GenerateRandomEffectsOL(Ref<VectorXd> mu, const Ref<VectorXd>& x,
			     const Ref<MatrixXd>& C, const Ref<MatrixXd>& OmegaL,
			     TempPQ *tmp);

// Generate mu = [(I-B)*V]^{1/2} * z + B * x,
// Omega = Sigma^{-1}.
void GenerateRandomEffectsO(Ref<VectorXd> mu, const Ref<VectorXd>& x,
			    const Ref<MatrixXd>& C, const Ref<MatrixXd>& Omega,
			    TempPQ *tmp);

#endif
