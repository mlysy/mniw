///////////////////////////////////////////

// Default libraries to include
// this is for multi-software portability (R/MATLAB).
//

// which software is being used: R or MATLAB
#define R_FUNCTION_LIBRARY
//#define MATLAB_FUNCTION_LIBRARY

// R includes (only includes)
#ifdef R_FUNCTION_LIBRARY
#include <Rcpp.h>
#include <RcppEigen.h>
#endif

// MATLAB includes
#ifdef MATLAB_FUNCTION_LIBRARY
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <Eigen/Dense>
#include "RmathUtils.h"
#endif
