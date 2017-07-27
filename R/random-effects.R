#--- random effects models -------------------------------------------------

#' @title Random Effects Sampling for Multivariate Normal-Normal Model.
#' @description Sample from the multivariate normal conditional distribution of the random effects given the data.
#' @param n The number of random samples to generate.
#' @param y Vector of length \code{q} or \code{q x n} matrix of observations.  In the latter case each column is a different vector (see details).
#' @param V Matrix of size \code{q x q} or \code{q x q x n} array of observation variances.
#' @param lambda Vector of length \code{q} or \code{q x n} matrix of prior means.
#' @param A Matrix of size \code{q x q} or \code{q x q x n} array of prior variances.
#' @details The normal-normal random effects model is:
#' \deqn{y | mu \sim N(mu, V), \qquad mu \sim \N(lambda, A).}
#' The posterior distribution for the random effects \code{mu} is
#' \deqn{mu | y \sim N(,).}
#' @export
rmNormRE <- function(n, y, V, lambda, A) {
  # problem dimensions
  if(is.vector(y)) y <- as.matrix(y)
  q <- dim(y)[1]
  y <- matrix(y, q)
  if(!all(dim(V)[1:2] == q)) {
    stop("y and V have incompatible dimensions.")
  }
  V <- matrix(V, q)
  if(is.vector(lambda)) lambda <- as.matrix(lambda)
  if(dim(lambda)[1] != q) {
    stop("y and lambda have incompatible dimensions.")
  }
  lambda <- matrix(lambda, q)
  if(!all(dim(A)[1:2] == q)) {
    stop("y and A have incompatible dimensions.")
  }
  V <- matrix(A, q)
  # check lengths
  N <- unique(sort(c(ncol(y), ncol(V)/q, ncol(lambda), ncol(A)/q)))
  N <- c(1, N[N > 1])
  if(length(N) > 2 || (length(N) == 2 && N[2] != n))
    stop("Arguments don't all have length n.")
  Mu <- GenerateRandomEffectsNormal(N, lambda, y, V, A)
  Mu
}
