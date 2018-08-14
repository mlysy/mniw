/// @file TriUtils.h
/// @author Martin Lysy (mlysy@uwaterloo.ca)
///
/// Utility functions for triangular matrices
///
/// Various types of products and solvers for linear systems involving lower/upper triangular matrices.  All arguments to these functions can be subsets of larger matrices; see "Writing Functions Taking Eigen Types as Parameters" in the Eigen documentation.
///
/// @note TODO:
///
/// - Put in namespace
/// - Remove `using namespace Eigen`
/// - Remove macro documentation (
///

#ifndef TriUtils_h
#define TriUtils_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
//#include <iostream>
using namespace Eigen;

//------------------------------------------------------------------------------

/// Multiplication of two lower triangular matrices
///
/// Performs the matrix multiplication `X = L1 * L2`, where `L1` and `L2` are lower triangular matrices.
///
/// @param [in] L1 First `n x n` lower triangular matrix.
/// @param [in] L2 Second `n x n` lower triangular matrix.
/// @param [out] X Matrix of size `n x n` containing the product of `L1` and `L2`.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1, typename T2>
  void triMultLL(const Eigen::MatrixBase<T1>& X,
		 const Ref<const MatrixXd>& L1,
		 const Eigen::MatrixBase<T2>& L2) {
  // hack to use triangularView: see Eigen documentation
  //"Writing Functions Taking Eigen Types as Parameters"
  // #define _VL const_cast<Eigen::MatrixBase<T1>& >(VL)
  // #define _XL const_cast<Eigen::MatrixBase<T2>& >(XL)
  #define _X const_cast<Eigen::MatrixBase<T1>& >(X)
  #define _L2 const_cast<Eigen::MatrixBase<T2>& >(L2)
  _X.template triangularView<Eigen::Lower>() = L1 * _L2.template triangularView<Eigen::Lower>();
  #undef _X
  #undef _L2
  return;
}

/// In-place solution of a lower triangular system
///
/// Performs the matrix multiplication `X = L^{-1} * X`, where `L` is a lower triangular matrix.
///
/// @param [in,out] X Matrix of size `n x p` on RHS of linear system.  This is also where solution will be returned, i.e., in-place.
/// @param [in] L Lower triangular matrix of size `n x n` on LHS of linear system.
///
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1>
void triMultLiX(Ref<MatrixXd> X, const Eigen::MatrixBase<T1>& L) {
  #ifdef _MSC_VER
  L.template triangularView<Eigen::Lower>().solveInPlace(X);
  #else
  L.template triangularView<Eigen::Lower>().template solveInPlace(X);
  #endif
  return;
}

/// Right-multiplication by upper triangular matrix
///
/// Performs the matrix multiplication `Y = X * U`, where `U` is an upper triangular matrix.
///
/// @param [in] X Matrix of size `n x p` to be multiplied.
/// @param [in] U Lower triangular matrix of size `n x n` right-multiplying `X`.
/// @param [out] Y Matrix of size `n x p` containing the product of `X` and `U`.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1>
void triMultXU(Ref<MatrixXd> Y, const Ref<const MatrixXd>& X,
	       const Eigen::MatrixBase<T1>& U) {
  Y.noalias() = X * U.template triangularView<Eigen::Upper>();
  return;
}

/// Left-multiplication by a lower triangular matrix
///
/// Performs the matrix multiplication `Y = L * X`, where `L` is a lower triangular matrix.
///
/// @param [in] X Matrix of size `n x p` to be multiplied.
/// @param [in] L Lower triangular matrix of size `n x n` left-multiplying `X`.
/// @param [out] Y Matrix of size `n x p` containing the product of `L` and `X`.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1>
void triMultLX(Ref<MatrixXd> Y, const Eigen::MatrixBase<T1>& L,
		const Ref<const MatrixXd>& X) {
  Y.noalias() = L.template triangularView<Eigen::Lower>() * X;
  return;
}

/// In-place solution of a reverse lower triangular system
///
/// Performs the multiplication `X = X * L^{-1}`, where `L` is a lower triangular matrix.
///
/// @param [in,out] X Matrix of size `n x p` on RHS of linear system.  This is also where solution will be returned, i.e., in-place.
/// @param [in] L Lower triangular matrix of size `n x n` on LHS of linear system.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1>
void triMultXLi(Ref<MatrixXd> X, const Eigen::MatrixBase<T1>& L) {
  #ifdef _MSC_VER
  L.template triangularView<Eigen::Lower>().solveInPlace<OnTheRight>(X);
  #else
  L.template triangularView<Eigen::Lower>().template solveInPlace<OnTheRight>(X);
  #endif
  return;
}

/// In-place solution of an upper triangular system
/// Performs the multiplication `X = U^{-1} * X`, where `U` is an upper triangular matrix.
///
/// @param [in,out] X Matrix of size `n x p` on RHS of linear system.  This is also where solution will be returned, i.e., in-place.
/// @param [in] U Upper triangular matrix of size `n x n` on LHS of linear system.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1>
void triMultUiX(const Eigen::MatrixBase<T1>& U, Ref<MatrixXd> X) {
  #ifdef _MSC_VER
  U.template triangularView<Eigen::Upper>().solveInPlace(X);
  #else
  U.template triangularView<Eigen::Upper>().template solveInPlace(X);
  #endif
  return;
}


//------------------------------------------------------------------------------

/// Transpose-product of upper triangular matrices
///
/// Performs the multiplication `X = U' * U`, where `U` is an upper triangular matrix.
///
/// @param [out] X Matrix of size `n x n` containing the transpose-product of `U`.
/// @param [in] U Upper triangular matrix of size `n x n`.
/// @param [in] L Lower triangular matrix of size `n x n` used for intermediate calculations.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1, typename T2>
void crossProdUtU(const Eigen::MatrixBase<T1>& X,
		  const Ref<const MatrixXd>& U,
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

/// Reverse transpose-product of lower triangular matrices
///
/// Performs the multiplication `X = L * L'`, where `L` is a lower triangular matrix.
///
/// @param [out] X Matrix of size `n x n` containing the reverse transpose-product of `L`.
/// @param [in] L Lower triangular matrix of size `n x n`.
/// @param [in] U Upper triangular matrix of size `n x n` used for intermediate calculations.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1, typename T2>
void CrossProdLLt(const Eigen::MatrixBase<T1>& X,
		  const Ref<const MatrixXd>& L,
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

/// Transpose-product of lower triangular matrices
///
/// Performs the multiplication `X = L' * L`, where `L` is a lower triangular matrix.
///
/// @param [out] X Matrix of size `n x n` containing the transpose-product of `L`.
/// @param [in] L Lower triangular matrix of size `n x n`.
/// @param [in] U Upper triangular matrix of size `n x n` used for intermediate calculations.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1, typename T2>
void CrossProdLtL(const Eigen::MatrixBase<T1>& X,
		  const Ref<const MatrixXd>& L,
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

/// Inverse of reverse transpose-product
///
/// Performs the calculation `X = (L * L')^{-1}`, where `L` is a lower triangular matrix.
///
/// @param [out] X Matrix of size `n x n` containing the inverse of the reverse transpose-product of `L`.
/// @param [in] L Lower triangular matrix of size `n x n`.
/// @param [in] U Upper triangular matrix of size `n x n` used for intermediate calculations.
///
/// @param [in] L2 Lower triangular matrix of size `n x n` used for intermediate calculations.
/// @param [in] I Identity matrix of size `n x n` used for intermediate calculations.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1>
void InverseLLt(Ref<MatrixXd> X,
		const Eigen::MatrixBase<T1>& L,
		Ref<MatrixXd> L2,
		Ref<MatrixXd> U,
		const Ref<const MatrixXd>& I) {
  // invert L
  L2 = L.template triangularView<Eigen::Lower>().solve(I);
  //calculate L2' * L2
  CrossProdLtL(X, L2, U);
  return;
}

#endif
