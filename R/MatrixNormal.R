#--- matrix normal distribution --------------------------------------------

#' @name MatrixNormal
#' @title Matrix-Normal Distribution.
#' @description Density and Random sampling for the Matrix-Normal Distribution.
#' @aliases dMNorm rMNorm
#' @param X Argument to the density function.  Either a \code{p x q} matrix or a \code{p x q x n} array.
#' @param n Number of random matrices to generate.
#' @param Mu Mean matrix or matrices.  Either a \code{p x q} matrix or a \code{p x q x n} array.
#' @param RowV Between-Row covariance matrix or matrices.  Either a \code{p x p} matrix or a \code{p x p x n} array.
#' @param ColV Between-Column covariance matrix or matrices.  Either a \code{q x q} matrix or a \code{q x q x n} array.
#' @param log Whether or not to compute the log-density.
#' @details  \code{dMNorm} and \code{rMNorm} both accept single or multiple values for each argument.
#'
#' The distribution of a \code{p x q} Matrix-Normal random matrix is
#' \deqn{f(X) = const * \exp{-0.5 Tr( ColV^{-1}(X-Mu)''RowV^{-1}(X-Mu) )}.}
#' Linear combinations of Matrix-Normals are normal.  That is, for fixed \code{a} and \code{b} and \code{X ~ MNorm(Mu, RowV, ColV)}, then
#' \deqn{a' X b \sim N( a' Mu b, (a' RowV a) * (b' ColV b) ).}
#' @return A vector for densities, or a \code{p x q x n} array for random sampling.

#--- lower-level functions -------------------------------------------------

#' @rdname MatrixNormal
#' @export
dMNorm <- function(X, Mu, RowV, ColV, log = FALSE, debug = FALSE) {
  # get dimensions
  PQ <- .getPQ(X = X, Lambda = Mu, Sigma = RowV, Psi = ColV)
  p <- PQ[1]
  q <- PQ[2]
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  # format arguments
  X <- .setDims(X, p = p, q = q)
  if(anyNA(X)) stop("Something went wrong.  Please report bug.")
  Mu <- .setDims(Mu, p = p, q = q)
  if(anyNA(Mu)) stop("Mu and X have incompatible dimensions.")
  RowV <- .setDims(RowV, p = p)
  if(anyNA(RowV)) stop("RowV and X have incompatible dimensions.")
  ColV <- .setDims(ColV, q = q)
  if(anyNA(ColV)) stop("ColV and X have incompatible dimensions.")
  # check lengths
  N <- .getN(p = p, q = q, X = X, Lambda = Mu, Sigma = RowV,
             Psi = ColV)
  if(length(N) > 2) stop("Arguments have different lengths.")
  if(debug) browser()
  ans <- LogDensityMatrixNormal(X, Mu, RowV, ColV)
  if(!log) ans <- exp(ans)
  ans
}

#' @rdname MatrixNormal
#' @export
rMNorm <- function(n, Mu, RowV, ColV, debug = FALSE) {
  # get dimensions
  PQ <- .getPQ(Lambda = Mu, Sigma = RowV, Psi = ColV)
  p <- PQ[1]
  q <- PQ[2]
  if(any(anyNA(PQ))) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  # format arguments
  Mu <- .setDims(Mu, p = p, q = q)
  if(anyNA(Mu)) stop("Something went wrong.  Please report bug.")
  RowV <- .setDims(RowV, p = p)
  if(anyNA(RowV)) stop("RowV and Mu have incompatible dimensions.")
  ColV <- .setDims(ColV, q = q)
  if(anyNA(ColV)) stop("ColV and Mu have incompatible dimensions.")
  # check lengths
  N <- .getN(p = p, q = q, Lambda = Mu, Sigma = RowV,
             Psi = ColV)
  if(length(N) > 2 || (length(N) == 2 && N[2] != n))
    stop("Arguments don't all have length n.")
  if(debug) browser()
  X <- GenerateMatrixNormal(n, Mu, RowV, ColV)
  if(n > 1) X <- array(X, dim = c(p,q,n))
  X
}
