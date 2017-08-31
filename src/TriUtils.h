//////////////////////////////////////////////////////////////////

// Utility functions for Triangular Matrices
// triangular views in Eigen need to be templated inside functions, so can't separate
// declaration and implementation in the usual way.
// instead include them both in this small header and proceed as usual in other headers.

//////////////////////////////////////////////////////////////////

#ifndef TriUtils_h
#define TriUtils_h 1

/* #include <Rcpp.h> */
/* using namespace Rcpp; */
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
//#include <iostream>
using namespace Eigen;

//////////////////////////////////////////////////////////////////

// various types of triangular multiplications
// NOTE: none of these work properly if source and destination refer to same space in memory
// NOTE: For safest results set all "off-triangular" arguments to zero.

// X = L1 * L2
template <typename T1, typename T2>
  void triMultLL(const Eigen::MatrixBase<T1>& X,
		 const Ref<MatrixXd>& L1,
		 const Eigen::MatrixBase<T2>& L2) {
  // hack to use triangularView: see Eigen documentation
  //"Writing Functions Taking Eigen Types as Parameters"
  #define _VL const_cast<Eigen::MatrixBase<T1>& >(VL)
  #define _XL const_cast<Eigen::MatrixBase<T2>& >(XL)
  #define _X const_cast<Eigen::MatrixBase<T1>& >(X)
  #define _L2 const_cast<Eigen::MatrixBase<T2>& >(L2)
  _X.template triangularView<Eigen::Lower>() = L1 * _L2.template triangularView<Eigen::Lower>();
  #undef _X
  #undef _L2
  return;
}

// X = L^{-1} * X
template <typename T1>
void triMultLiX(Ref<MatrixXd> X, const Eigen::MatrixBase<T1>& L) {
  #ifdef _MSC_VER
  L.template triangularView<Eigen::Lower>().solveInPlace(X);
  #else
  L.template triangularView<Eigen::Lower>().template solveInPlace(X);
  #endif
  return;
}

// Y = X * U;
template <typename T1>
void triMultXU(Ref<MatrixXd> Y, const Ref<const MatrixXd>& X,
	       const Eigen::MatrixBase<T1>& U) {
  Y.noalias() = X * U.template triangularView<Eigen::Upper>();
  return;
}

// Y = L * X
template <typename T1>
void triMultLX(Ref<MatrixXd> Y, const Eigen::MatrixBase<T1>& L,
		const Ref<const MatrixXd>& X) {
  Y.noalias() = L.template triangularView<Eigen::Lower>() * X;
  return;
}

// X = X * L^{-1}
template <typename T1>
void triMultXLi(Ref<MatrixXd> X, const Eigen::MatrixBase<T1>& L) {
  #ifdef _MSC_VER
  L.template triangularView<Eigen::Lower>().solveInPlace<OnTheRight>(X);
  #else
  L.template triangularView<Eigen::Lower>().template solveInPlace<OnTheRight>(X);
  #endif
  return;
}

// X = U^{-1} * X
template <typename T1>
void triMultUiX(const Eigen::MatrixBase<T1>& U, Ref<MatrixXd> X) {
  #ifdef _MSC_VER
  U.template triangularView<Eigen::Upper>().solveInPlace(X);
  #else
  U.template triangularView<Eigen::Upper>().template solveInPlace(X);
  #endif
  return;
}


//////////////////////////////////////////////////////////////////

// various types of "crossproducts"
// NOTE: none of these work properly if source and destination refer to same space in memory
// NOTE: For safest results set all "off-triangular" arguments to zero.

// X = U' * U for upper triangular U
// L is scratch space
template <typename T1, typename T2>
void crossProdUtU(const Eigen::MatrixBase<T1>& X,
		  const Ref<MatrixXd>& U,
		  const Eigen::MatrixBase<T2>& L) {
  // hack
  #define _X const_cast<Eigen::MatrixBase<T1>& >(X)
  #define _L const_cast<Eigen::MatrixBase<T2>& >(L)
  _L.template triangularView<Eigen::Lower>() = U.adjoint();
  _X.template triangularView<Eigen::Upper>() = L.template triangularView<Eigen::Lower>() * U;
  _X.template triangularView<Eigen::Lower>() = X.adjoint();
  #undef _X
  #undef _L
  return;
}

// X = L * L' for lower triangular L
// U is scratch space
template <typename T1, typename T2>
void CrossProdLLt(const Eigen::MatrixBase<T1>& X,
		  const Ref<MatrixXd>& L,
		  const Eigen::MatrixBase<T2>& U) {
  // hack
  #define _X const_cast<Eigen::MatrixBase<T1>& >(X)
  #define _U const_cast<Eigen::MatrixBase<T2>& >(U)
  _U.template triangularView<Eigen::Upper>() = L.adjoint();
  _X.template triangularView<Eigen::Upper>() = L * U.template triangularView<Eigen::Upper>();
  _X.template triangularView<Eigen::Lower>() = X.adjoint();
  #undef _X
  #undef _U
  return;
}

// X = L' * L for lower triangular L
// U is scratch space
template <typename T1, typename T2>
void CrossProdLtL(const Eigen::MatrixBase<T1>& X,
		  const Ref<MatrixXd>& L,
		  const Eigen::MatrixBase<T2>& U) {
  // hack
  #define _X const_cast<Eigen::MatrixBase<T1>& >(X)
  #define _U const_cast<Eigen::MatrixBase<T2>& >(U)
  _U.template triangularView<Eigen::Upper>() = L.adjoint();
  _X.template triangularView<Eigen::Upper>() = U.template triangularView<Eigen::Upper>() * L;
  _X.template triangularView<Eigen::Lower>() = X.adjoint();
  #undef _X
  #undef _U
  return;
}

// X = (LL')^{-1} for lower triangular L
template <typename T1>
void InverseLLt(Ref<MatrixXd> X,
		const Eigen::MatrixBase<T1>& L,
		Ref<MatrixXd> L2,
		Ref<MatrixXd> U,
		const Ref<MatrixXd>& I) {
  // invert L
  L2 = L.template triangularView<Eigen::Lower>().solve(I);
  //calculate L2' * L2
  CrossProdLtL(X, L2, U);
  return;
}


#endif
