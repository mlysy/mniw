/// @file mniwUtils.h
/// @author Martin Lysy (mlysy@uwaterloo.ca)
///
/// Utilities for MNIW functions
///
/// - Class for memory allocation
/// - log-determinant of a positive-definite matrix

#ifndef mniwUtils_h
#define mniwUtils_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
//#include <iostream>
using namespace Eigen;
#include "TriUtils.h"

//------------------------------------------------------------------------------

/// Temporary storage class
///
/// Allocates space for various matrices used for intermediate calculations.
class TempPQ {
public:
  TempPQ(int p,int q); ///< Constructor
  MatrixXd Mq; ///< Storage matrix of size `q x q`.
  MatrixXd Iq; ///< Identity matrix of size `q x q`.
  MatrixXd Lq; ///< Storage lower triangular matrix of size `q x q`.
  MatrixXd Uq; ///< Storage upper triangular matrix of size `q x q`.
  MatrixXd Xpq; ///< Storage matrix of size `p x q`.
  //VectorXd vq; // q-vector
  LLT<MatrixXd> lltq; ///< Storage Cholesky solver of size `q x q`.
};

/// @param p First dimension of storage variables.
/// @param q Second dimension of storage variables.
inline TempPQ::TempPQ(int p, int q) {
  Mq = MatrixXd::Zero(q,q);
  Iq = MatrixXd::Identity(q,q);
  Lq = MatrixXd::Zero(q,q);
  Uq = MatrixXd::Zero(q,q);
  Xpq = MatrixXd::Zero(p,q);
  //vq = VectorXd::Zero(q);
  lltq.compute(Iq);
}

/// Logarithm of the determinant of a positive-definite matrix
///
/// Calculates `log|V|`, where `V` is a positive-definite matrix.
///
/// @param [in] V Positive-definite matrix of size `n x n`.
/// @param [in] cholV Cholesky solver of size `n x n` required for intermediate calculations.
/// @return The logarithm of the determinant of `V`.
inline double logDetV(MatrixXd V, LLT<MatrixXd> cholV) {
  double ldV = 0.0;
  cholV.compute(V);
  for(int ii=0; ii<V.cols(); ii++) {
    ldV += log(cholV.matrixL()(ii,ii));
  }
  return 2.0*ldV;
}

#endif
