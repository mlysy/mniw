///////////////////////////////////////////////////////////////////

// Utilities for the other MNIW functions:
// * Class for memory allocation
// * various cross products
// * reverse cholesky decomposition

//////////////////////////////////////////////////////////////////

#ifndef mniwUtils_h
#define mniwUtils_h 1

/* #include <Rcpp.h> */
/* using namespace Rcpp; */
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
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
inline void AntiTransposeLowerTri(Ref<MatrixXd> Y, const Ref<MatrixXd>& X);
inline void ReverseCholesky(Ref<MatrixXd> L, const Ref<MatrixXd>& V, TempPQ *tmp);

// log|V|
inline double logDetV(MatrixXd V, LLT<MatrixXd> cholV);

////////////////////////////////////////////////////////////

// a scratch space class containing elements that shouldn't be created all the time

/*
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
*/
inline TempPQ::TempPQ(int p, int q) {
  Mq = MatrixXd::Zero(q,q);
  Iq = MatrixXd::Identity(q,q);
  Lq = MatrixXd::Zero(q,q);
  Uq = MatrixXd::Zero(q,q);
  Xpq = MatrixXd::Zero(p,q);
  //vq = VectorXd::Zero(q);
  lltq.compute(Iq);
}

////////////////////////////////////////////////////////////////////

// reverse cholesky decomposition

// transpose the lower triangular elements of square matrix X across the anti-diagonal.
// does not operate in place.
void AntiTransposeLowerTri(Ref<MatrixXd> Y, const Ref<MatrixXd>& X) {
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

// reverse cholesky decomposition V = L'L.
void ReverseCholesky(Ref<MatrixXd> L, const Ref<MatrixXd>& V, TempPQ *tmp) {
  AntiTransposeLowerTri(tmp->Mq, V);
  tmp->lltq.compute(tmp->Mq);
  tmp->Lq = tmp->lltq.matrixL();
  AntiTransposeLowerTri(L, tmp->Lq);
  return;
}

// log|V|
double logDetV(MatrixXd V, LLT<MatrixXd> cholV) {
  double ldV = 0.0;
  cholV.compute(V);
  for(int ii=0; ii<V.cols(); ii++) {
    ldV += log(cholV.matrixL()(ii,ii));
  }
  return 2.0*ldV;
}


#endif
