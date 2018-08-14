/// @file mniwUtils.h
/// @author Martin Lysy (mlysy@uwaterloo.ca)
///
/// Utilities for MNIW functions
///
/// - Class for memory allocation
/// - Reverse Cholesky decomposition

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

/// Anti-transpose of a matrix
///
/// Transpose the lower triangular elements of a square matrix across the anti-diagonal (the remaining elements of the output are left untouched).  So for example, we would have
///
/// \f[
/// \begin{bmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{bmatrix}
/// \qquad \longrightarrow \qquad
/// \begin{bmatrix} 9 & & \\ 8 & 5 & \\ 7 & 4 & 1 \end{bmatrix}.
/// \f]
///
/// @param [out] Y Matrix of size `n x n` to which anti-transpose is output.
/// @param [in] X Matrix of size `n x n` of which the lower triangular elements will be anti-transposed.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
inline void AntiTransposeLowerTri(Ref<MatrixXd> Y, const Ref<const MatrixXd>& X) {
  int q = X.cols();
  int ii, jj;
  // anti-transpose lower part of X into Y
  for(ii=0; ii<q; ii++) {
    for(jj=0; jj<=ii; jj++) {
      Y(q-1-jj,q-1-ii) = X(ii,jj);
    }
  }
  return;
}

/// Reverse-Cholesky decomposition of a positive-definite matrix
///
/// Calculates the lower triangular matrix `L` satisfying `V = L'L`, where `V` is a positive-definite matrix.
///
/// @param [out] L Lower triangular matrix of size `n x n`.
/// @param [in] V Positive-definite matrix of size `n x n`.
/// @param [in] tmp Object of class `TempPQ` used to store intermediate calculations.
inline void ReverseCholesky(Ref<MatrixXd> L, const Ref<const MatrixXd>& V, TempPQ *tmp) {
  AntiTransposeLowerTri(tmp->Mq, V);
  tmp->lltq.compute(tmp->Mq);
  tmp->Lq = tmp->lltq.matrixL();
  AntiTransposeLowerTri(L, tmp->Lq);
  return;
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
