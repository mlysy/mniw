#--- matrix-t distribution ------------------------------------------------

#' The Matrix-t distribution.
#'
#' Density and sampling for the Matrix-t distribution.
#'
#' @name MatrixT
#' @aliases dMT
#' @param X Argument to the density function.  Either a \code{p x q} matrix or a \code{p x q x n} array.
#' @param n Number of random matrices to generate.
#' @param Mu Mean parameter  Either a \code{p x q} matrix or a \code{p x q x n} array.
#' @param RowV Between-row covariance matrix.  Either a \code{p x p} matrix or a \code{p x p x n} array.
#' @param ColV Between-column covariance matrix  Either a \code{q x q} matrix or a \code{q x q x n} array.
#' @param nu Degrees-of-freedom parameter.  A scalar or vector.
#' @param log Logical; whether or not to compute the log-density.
#'
#' @return A vector of length \code{n} for density evaluation, or a \code{p x q x n} array for random sampling.
#'
#' @template details-matrixt

#--- lower-level functions ------------------------------------------------

#' @rdname MatrixT
#' @export
dMT <- function(X, Mu, RowV, ColV, nu, log = FALSE) {
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
  nu <- c(nu)
  # check lengths
  N <- .getN(p = p, q = q, X = X, Lambda = Mu, Sigma = RowV,
             Psi = ColV, nu = nu)
  if(length(N) > 2) stop("Arguments have different lengths.")
  ans <- LogDensityMatrixT(X, Mu, RowV, ColV, nu)
  if(!log) ans <- exp(ans)
  ans
}

#' @rdname MatrixT
#' @export
rMT <- function(n, Mu, RowV, ColV, nu) {
  # get dimensions
  PQ <- .getPQ(Lambda = Mu, Sigma = RowV, Psi = ColV)
  p <- PQ[1]
  q <- PQ[2]
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  Mu <- .setDims(Mu, p = p, q = q)
  if(anyNA(Mu)) stop("Something went wrong.  Please report bug.")
  RowV <- .setDims(RowV, p = p)
  if(anyNA(RowV)) stop("RowV and Mu have incompatible dimensions.")
  ColV <- .setDims(ColV, q = q)
  if(anyNA(ColV)) stop("ColV and Mu have incompatible dimensions.")
  nu <- c(nu)
  # check lengths
  N <- .getN(p = p, q = q, Lambda = Mu, Sigma = RowV,
             Psi = ColV, nu = nu)
  if(length(N) > 2 || (length(N) == 2 && N[2] != n)) {
    stop("Arguments don't all have length n.")
  }
  X <- GenerateMatrixT(n, Mu, RowV, ColV, nu)
  if(n > 1) X <- array(X, dim = c(p,q,n))
  X
}
