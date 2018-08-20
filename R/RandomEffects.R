#--- random effects models -------------------------------------------------

#' Conditional sampling for Multivariate Normal Random-Effects model.
#'
#' Sample from the conditional parameter distribution given the data and hyperparameters of the Multivariate Normal Random Effects (mNormRE) model.
#'
#' @param n the number of random samples to generate.
#' @param y vector of length \code{q} or \code{n x q} matrix of observations.  In the latter case each row is a different vector (see Details).
#' @param V matrix of size \code{q x q} or \code{q x q x n} array of observation variances.
#' @param lambda vector of length \code{q} or \code{n x q} matrix of prior means.  In the latter case each row is a different mean.  Defaults to zeros.
#' @param A matrix of size \code{q x q} or \code{q x q x n} array of prior variances.  Defaults to identity matrix.
#' @details The normal-normal random effects model is:
#' \deqn{y | \mu \sim N(\mu, V),}
#' \deqn{\mu \sim N(\lambda, A).}
#' The posterior distribution for the random effects \eqn{\mu} is
#' \deqn{\mu | y \sim N(B(\lambda - y) + y,(I - B)V),}
#' where \eqn{B = V(V+A)^{-1}}.
#' @examples
#' # Conditional sampling with default prior
#' n = 100
#' q = 10
#' y = rnorm(q)
#' V = rwish(1, diag(q), q+1)
#' lambda = rep(0,q)
#' A = diag(q)
#' rmNormRE(n, y, V, lambda, A)
#'
#' @export
rmNormRE <- function(n, y, V, lambda, A) {
  # convert to MN format
  y <- .vec2mn(y)
  if(!missing(lambda)) lambda <- .vec2mn(lambda)
  # problem dimensions
  PQ1 <- .getPQ(X = y, Lambda = lambda, Sigma = V)
  PQ2 <- .getPQ(X = y, Sigma = A)
  p <- PQ1[1]
  q <- PQ1[2]
  if(q != 1) stop("y and lambda must be vectors or matrices.")
  # format arguments
  y <- .setDims(y, p = p, q = q)
  lambda <- .setDims(lambda, p = p, q = q)
  if(anyNA(lambda)) stop("lambda and y have incompatible dimensions.")
  V <- .setDims(V, p = p)
  if(anyNA(V)) stop("V and y have incompatible dimensions.")
  A <- .setDims(A, p = p)
  if(anyNA(A)) stop("A and y have incompatible dimensions.")
  # check lengths
  N1 <- .getN(p = p, q = q, X = y, Lambda = lambda, Sigma = V)
  N2 <- .getN(p = p, q = q, X = y, Sigma = A)
  N <- unique(sort(c(N1, N2)))
  N <- c(1, N[N>1])
  if(length(N) > 2 || (length(N) == 2 && N[2] != n)) {
    stop("Arguments don't all have length n.")
  }
  ## if(is.vector(y)) {
  ##   y <- as.matrix(y)
  ## } else {
  ##   x <- t(x)
  ##   x <- array(x, dim = c(1, dim(x)))
  ## }
  ## if(is.vector(lambda)) {
  ## }
  ## # problem dimensions
  ## if(is.vector(y)) y <- as.matrix(y)
  ## q <- dim(y)[1]
  ## y <- matrix(y, q)
  ## if(!all(dim(V)[1:2] == q)) {
  ##   stop("y and V have incompatible dimensions.")
  ## }
  ## V <- matrix(V, q)
  ## if(is.vector(lambda)) lambda <- as.matrix(lambda)
  ## if(dim(lambda)[1] != q) {
  ##   stop("y and lambda have incompatible dimensions.")
  ## }
  ## lambda <- matrix(lambda, q)
  ## if(!all(dim(A)[1:2] == q)) {
  ##   stop("y and A have incompatible dimensions.")
  ## }
  ## V <- matrix(A, q)
  # check lengths
  ## N <- unique(sort(c(ncol(y), ncol(V)/q, ncol(lambda), ncol(A)/q)))
  ## N <- c(1, N[N > 1])
  ## if(length(N) > 2 || (length(N) == 2 && N[2] != n))
  ##   stop("Arguments don't all have length n.")
  Mu <- GenerateRandomEffectsNormal(n, lambda, y, V, A)
  if(n > 1) {
    Mu <- t(Mu)
  } else {
    Mu <- c(Mu)
  }
  Mu
}
