#--- multivariate normal distribution --------------------------------------

#' @name MultiNormal
#' @title The Multivariate Normal distribution.
#' @description Density and random sampling for the Multivariate Normal distribution.
#' @aliases dmNorm rmNorm
#' @param x argument to the density function.  A vector of length \code{q} or an \code{n x q} matrix.
#' @param n number of random vectors to generate.
#' @param mu mean vector(s).  Either a vector of length \code{q} or an \code{n x q} matrix.
#' @param V covariance matrix or matrices.  Either a \code{q x q} matrix or a \code{q x q x n} array.
#' @param log logical. Whether or not to compute the log-density.
#' @details  \code{dmNorm} and \code{rmNorm} both accept single or multiple values for each argument.
#' @return A vector for densities, or a \code{n x q} array for random sampling.
#' @examples
#' ## Multivariate Normal random sample and subsequent density calculation
#' n = 100
#' q = 2
#' mu = c(1,2)
#' V = toeplitz(exp(-seq(1:q)))
#'
#' # Random sample from Multivariate Normal distribution
#' mnormData = rmNorm(n, mu, V)
#'
#' # Calculate log density for each sampled vector
#' dmNorm(mnormData, mu, V, log=TRUE)
#' @rdname MultiNormal
#' @export
dmNorm <- function(x, mu, V, log = FALSE) {
  ## debug <- FALSE
  ## if(FALSE) {
  ##   .getPQ <- mniw:::.getPQ
  ##   .getN <- mniw:::.getN
  ##   .setDims <- mniw:::.setDims
  ## }
  # get dimensions
  # first convert to appropriate MN format
  x <- .vec2mn(x)
  if(!missing(mu)) mu <- .vec2mn(mu)
  ## if(is.vector(x)) {
  ##   x <- matrix(x, nrow = 1)
  ## } else {
  ##   x <- t(x)
  ##   x <- array(x, dim = c(1, dim(x)))
  ## }
  ## if(!missing(mu)) {
  ##   if(is.vector(mu)) {
  ##     mu <- matrix(mu, nrow = 1)
  ##   } else {
  ##     mu <- t(mu)
  ##     mu <- array(mu, dim = c(1,dim(mu)))
  ##   }
  ## }
  PQ <- .getPQ(X = x, Lambda = mu, Sigma = V)
  p <- PQ[1]
  q <- PQ[2]
  if(q != 1) stop("x and mu must be vectors or matrices.")
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  # format arguments
  x <- .setDims(x, p = p, q = q)
  if(anyNA(x)) stop("Something went wrong.  Please report bug.")
  mu <- .setDims(mu, p = p, q = q)
  if(anyNA(mu)) stop("mu and x have incompatible dimensions.")
  V <- .setDims(V, p = p)
  if(anyNA(V)) stop("V and x have incompatible dimensions.")
  # check lengths
  N <- .getN(p = p, q = q, X = x, Lambda = mu, Sigma = V)
  if(length(N) > 2) stop("Arguments have different lengths.")
  x <- matrix(x, nrow = q) # format for mN
  mu <- matrix(mu, nrow = q) # format for mN
  if(debug) browser()
  ans <- LogDensityMultivariateNormal(x, mu, V)
  if(!log) ans <- exp(ans)
  ans
}

#' @rdname MultiNormal
#' @export
rmNorm <- function(n, mu, V, debug = FALSE) {
  # get dimensions
  # first convert to appropriate MN format
  if(!missing(mu)) mu <- .vec2mn(mu)
  ## if(!missing(mu)) {
  ##   if(is.vector(mu)) {
  ##     mu <- matrix(mu, nrow = 1)
  ##   } else {
  ##     mu <- t(mu)
  ##     mu <- array(mu, dim = c(1,dim(mu)))
  ##   }
  ## }
  PQ <- .getPQ(Lambda = mu, Sigma = V)
  #if(is.na(PQ[1])) PQ[1] <- 1 # what is this for???
  p <- PQ[1]
  q <- PQ[2]
  if(q != 1) stop("mu must be a vector or matrix.")
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined: must provide mu or V.")
  }
  # format arguments
  mu <- .setDims(mu, p = p, q = q)
  if(anyNA(mu)) stop("mu and V have incompatible dimensions.")
  V <- .setDims(V, p = q)
  if(anyNA(V)) stop("V and mu have incompatible dimensions.")
  # check lengths
  N <- .getN(p = p, q = q, Lambda = mu, Sigma = V)
  if(length(N) > 2 || (length(N) == 2 && N[2] != n)) {
    stop("Arguments don't all have length n.")
  }
  if(debug) browser()
  mu <- matrix(mu, nrow = q) # format for mN
  X <- GenerateMultivariateNormal(n, mu, V)
  if(n > 1) {
    X <- t(X)
  } else {
    X <- c(X)
  }
  X
}
