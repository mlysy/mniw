#--- matrix-t distribution ------------------------------------------------

#' The Matrix-t distribution.
#'
#' Density and sampling for the Matrix-t distribution.
#'
#' @name MatrixT
#' @aliases dMT
#' @param X argument to the density function.  Either a \code{p x q} matrix or a \code{p x q x n} array.
#' @param Mu mean matrix or matrices.  Either a \code{p x q} matrix or a \code{p x q x n} array.
#' @param RowV between-row covariance matrix or matrices.  Either a \code{p x p} matrix or a \code{p x p x n} array.
#' @param ColV between-column covariance matrix or matrices.  Either a \code{q x q} matrix or a \code{q x q x n} array.
#' @param nu degrees-of-freedom parameter.  A scalar or vector.
#' @param log logical. Whether or not to compute the log-density.
#' @return A vector for densities, or a \code{p x q x n} array for random sampling.

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
