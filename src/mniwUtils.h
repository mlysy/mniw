///////////////////////////////////////////////////////////////////

// Utilities for the other MNIW functions:
// * Class for memory allocation
// * various cross products
// * reverse cholesky decomposition

//////////////////////////////////////////////////////////////////

#ifndef mniwUtils_h
#define mniwUtils_h 1

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

////////////////////////////////////////////////////////////

// a scratch space class containing elements that shouldn't be created all the time
class TempPQ {
public:
  TempPQ(int,int);
  MatrixXd Mq; // (q x q) matrix
  MatrixXd Iq; // (q x q) identity matrix
  MatrixXd Lq; // (q x q) lower triangular matrix
  MatrixXd Uq; // (q x q) upper triangular matrix
  MatrixXd Xpq; // (p x q) matrix
  //VectorXd vq; // q-vector
  LLT<MatrixXd> lltq; // (q x q) cholesky solver
};

// Reverse Cholesky decomposition
void AntiTransposeLowerTri(Ref<MatrixXd> Y, const Ref<MatrixXd>& X);
void ReverseCholesky(Ref<MatrixXd> L, const Ref<MatrixXd>& V, TempPQ *tmp);

// log|V|
double logDetV(MatrixXd V, LLT<MatrixXd> cholV);


#endif
