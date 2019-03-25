#--- random effects models -------------------------------------------------

#' Conditional sampling for Multivariate-Normal Random-Effects model.
#'
#' Sample from the conditional parameter distribution given the data and hyperparameters of the Multivariate-Normal Random-Effects (mNormRE) model (see \strong{Details}).
#'
#' @template param-n
#' @param y Data observations.  Either a vector of length \code{q} or a \code{n x q} matrix.  In the latter case each row is a different vector.
#' @param V Observation variances.  Either a matrix of size \code{q x q} or a \code{q x q x n} array.
#' @param lambda Prior means.  Either a vector of length \code{q} or an \code{n x q} matrix.  In the latter case each row is a different mean.  Defaults to zeros.
#' @param A Prior variances.  Either a matrix of size \code{q x q} or a \code{q x q x n} array.  Defaults to identity matrix.
#' @details Consider the hierarchical multivariate normal model
#' \deqn{
#' \arraycolsep=1.4pt
#' \begin{array}{rcl}
#' \boldsymbol{\mu} & \sim & \mathcal{N}(\boldsymbol{\lambda}, \boldsymbol{A}) \\
#' \boldsymbol{y} \mid \boldsymbol{\mu} & \sim & \mathcal{N}(\boldsymbol{\mu}, \boldsymbol{V}).
#' \end{array}
#' }{
#' \mu ~ N(\lambda, A)
#' }
#' \deqn{
#' \vspace{-2em}
#' }{
#' y | \mu ~ N(\mu, V).
#' }
#' The Multivariate-Normal Random-Effects model \eqn{\boldsymbol{\mu} \sim \textrm{NormRX}(\boldsymbol{y}, \boldsymbol{V}, \boldsymbol{\lambda}, \boldsymbol{A})}{\mu ~ NormRX(y, V, \lambda, A)} on the random vector \eqn{\boldsymbol{\mu}_q}{\mu_q} is defined as the posterior distribution \eqn{p(\boldsymbol{\mu} \mid \boldsymbol{y}, \boldsymbol{\lambda}, \boldsymbol{A})}{p(\mu | y, \lambda, A)}.  This distribution is multivariate normal; for the mathematical specification of its parameters please see \code{vignette("mniw-distributions", package = "mniw")}.
#'
#' @example examples/RandomEffects.R
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
  Mu <- GenerateRandomEffectsNormal(n, lambda, y, V, A)
  if(n > 1) {
    Mu <- t(Mu)
  } else {
    Mu <- c(Mu)
  }
  Mu
}
