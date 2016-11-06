#--- multivariate normal distribution -------------------------------------------

#' @name mNorm
#' @title Multivariate Normal Distribution.
#' @description Density and Random sampling for the Multivariate Normal Distribution.
#' @param x Argument to the density function.  A vector of length \code{q} or an \code{n x q} matrix.
#' @param n Number of random vectors to generate.
#' @param mu Mean vector(s).  Either a vector of length \code{q} or an \code{n x q} matrix.
#' @param V Covariance matrix or matrices.  Either a \code{q x q} matrix or a \code{q x q x n} array.
#' @param log Whether or not to compute the log-density.
#' @details  \code{dmNorm} and \code{rmNorm} both accept single or multiple values for each argument.
#'
#' @return A vector for densities, or a \code{n x q} array for random sampling.

#' @rdname mNorm
#' @export
dmNorm <- function(x, mu, V, log = FALSE) {
  debug <- FALSE
  if(FALSE) {
    .getPQ <- mniw:::.getPQ
    .getN <- mniw:::.getN
    .setDims <- mniw:::.setDims
  }
  # get dimensions
  # first convert to appropriate MN format
  if(is.vector(x)) {
    x <- matrix(x, nrow = 1)
  } else {
    x <- t(x)
    x <- array(x, dim = c(1, dim(x)))
  }
  if(!missing(mu)) {
    if(is.vector(mu)) {
      mu <- matrix(mu, nrow = 1)
    } else {
      mu <- t(mu)
      mu <- array(mu, dim = c(1,dim(mu)))
    }
  }
  PQ <- .getPQ(X = x, Lambda = mu, Psi = V)
  p <- PQ[1]
  q <- PQ[2]
  if(p != 1) stop("p has wrong dimension.")
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  # format arguments
  x <- .setDims(x, p = p, q = q)
  if(anyNA(x)) stop("Something went wrong.  Please report bug.")
  mu <- .setDims(mu, p = p, q = q)
  if(anyNA(mu)) stop("mu and x have incompatible dimensions.")
  V <- .setDims(V, q = q)
  if(anyNA(V)) stop("V and x have incompatible dimensions.")
  # check lengths
  N <- .getN(p = p, q = q, X = x, Lambda = mu, Psi = V)
  if(length(N) > 2) stop("Arguments have different lengths.")
  x <- matrix(x, nrow = q) # format for mN
  mu <- matrix(mu, nrow = q) # format for mN
  if(debug) browser()
  ans <- .Call('mniw_LogDensityMultivariateNormal', PACKAGE = 'mniw',
               x, mu, V)
  if(!log) ans <- exp(ans)
  ans
}

#' @rdname mNorm
#' @export
rmNorm <- function(n, mu, V, debug = FALSE) {
  # get dimensions
  # first convert to appropriate MN format
  if(!missing(mu)) {
    if(is.vector(mu)) {
      mu <- matrix(mu, nrow = 1)
    } else {
      mu <- t(mu)
      mu <- array(mu, dim = c(1,dim(mu)))
    }
  }
  PQ <- .getPQ(Lambda = mu, Psi = V)
  if(is.na(PQ[1])) PQ[1] <- 1
  p <- PQ[1]
  q <- PQ[2]
  if(p != 1) stop("p has wrong dimension.")
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  # format arguments
  mu <- .setDims(mu, p = p, q = q)
  if(anyNA(mu)) stop("mu and V have incompatible dimensions.")
  V <- .setDims(V, q = q)
  if(anyNA(V)) stop("V and mu have incompatible dimensions.")
  # check lengths
  N <- .getN(p = p, q = q, Lambda = mu, Psi = V)
  if(length(N) > 2 || (length(N) == 2 && N[2] != n))
    stop("Arguments don't all have length n.")
  if(debug) browser()
  mu <- matrix(mu, nrow = q) # format for mN
  X <- .Call('mniw_GenerateMultivariateNormal', PACKAGE = 'mniw',
             n, mu, V)
  if(n > 1) {
    X <- t(X)
  } else {
    X <- c(X)
  }
  X
}
