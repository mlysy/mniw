#--- random effects models -------------------------------------------------

#' Conditional Sampling for Multivariate Normal Random-Effects Model.
#'
#' Sample from the conditional parameter distribution given the data and hyperparameters of the Multivariate Normal Random Effects (mNormRE) model.
#'
#' @param n The number of random samples to generate.
#' @param y Vector of length \code{q} or \code{n x q} matrix of observations.  In the latter case each column is a different vector (see details).
#' @param V Matrix of size \code{q x q} or \code{q x q x n} array of observation variances.
#' @param lambda Vector of length \code{q} or \code{n x q} matrix of prior means.  In the latter case each column is a different mean.  Defaults to zeros.
#' @param A Matrix of size \code{q x q} or \code{q x q x n} array of prior variances.  Defaults to identity matrix.
#' @details The normal-normal random effects model is:
#' \deqn{y | mu \sim N(mu, V), \qquad mu \sim \N(lambda, A).}
#' The posterior distribution for the random effects \code{mu} is
#' \deqn{mu | y \sim N(,).}
#' @export
rmNormRE <- function(n, y, V, lambda, A) {
  # convert to MN format
  y <- .vec2mn(y)
  if(!missing(lambda)) lambda <- .vec2mn(lambda)
  # problem dimensions
  PQ1 <- .getPQ(X = y, Lambda = lambda, Psi = V)
  PQ2 <- .getPQ(X = y, Psi = A)
  p <- PQ1[1]
  q <- PQ1[2]
  if(p != 1) stop("y and lambda must be vectors or matrices.")
  # format arguments
  y <- .setDims(y, p = p, q = q)
  lambda <- .setDims(lambda, p = p, q = q)
  if(anyNA(lambda)) stop("lambda and y have incompatible dimensions.")
  V <- .setDims(V, q = q)
  if(anyNA(V)) stop("V and y have incompatible dimensions.")
  A <- .setDims(A, q = q)
  if(anyNA(A)) stop("A and y have incompatible dimensions.")
  # check lengths
  N1 <- .getN(p = p, q = q, X = y, Lambda = lambda, Psi = V)
  N2 <- .getN(p = p, q = q, X = y, Psi = A)
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
  Mu <- GenerateRandomEffectsNormal(N, lambda, y, V, A)
  Mu
}

#' Bayesian Inference for Multivariate Normal Random-Effects Model
#'
#' Gibbs sampler for posterior distribution of parameters and hyperparameters of the Multivariate Normal Random Effects (mNormRE) model.
#'
#' @param nsamples Number of posterior samples to draw.
#' @param Y \code{N x q} matrix of responses.
#' @param V Either a \code{q x q} variance matrix or an \code{q x q x N} array of such matrices.
#' @param X \code{N x p} matrix of covariates.
#' @param prior Parameters of the prior MNIW distribution on the hyperparameters.  See Details.
#' @param burn Integer number of burn-in samples, or fraction of \code{nsamples} to prepend as burn-in.
#' @param init Optional list with elements \code{Beta}, \code{Sigma}, and \code{Mu} providing the initial values for these.  Default values are \code{Beta = matrix(0, p, q)}, \code{Sigma = diag(q)}, and \code{Mu = Y}.
#' @param updateHyp,storeHyp Whether or not to update/store the hyperparameter draws.
#' @param updateRE,storeRE Whether or not to update/store the random-effects draws.
#' @details The mNormRE model is given by
#' \deqn{
#' y_i \mid \mu_i \sim_iid N(\mu_i, V_i)
#' }
#' \deqn{
#' \mu_i \mid \beta, \Sigma ~sim_ind N(x_i' \beta, \Sigma)
#' }
#' \deqn{
#' \beta, \Sigma \sim MNIW(\Lambda, \Omega^{-1}, \Psi, \nu),
#' }
#' where \eqn{y_i} and \eqn{\mu_i} are response and random-effects vectors of length \eqn{q}, \eqn{x_i} are covariate vectors of length \eqn{p}, and \eqn{(\beta, \Sigma)} are hyperparameter matrices of size \eqn{p \times q} and \eqn{\q \times q}.
#'
#' The MNIW prior distribution is given by a list with elements \code{Lambda}, \code{Omega}, \code{Psi}, and \code{nu}.  If any of these is \code{NULL} or missing, the default value is 0.  Note that \code{Omega == 0} gives a Lebesgue prior to \eqn{\beta}.
#' @return A list with elements:
#' \describe{
#'   \item{\code{Beta}}{An \code{}
#' }
#' @export
mNormRE.post <- function(nsamples, Y, V, X, prior = NULL, init, burn,
			 updateHyp = TRUE, storeHyp = TRUE,
                         updateRE = TRUE, storeRE = FALSE) {
  # argument check
  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  if(nrow(X) != n) stop("Y and X have incompatible dimensions.")
  if(is.matrix(V)) V <- array(c(V), dim = c(dim(V), n))
  if(!identical(dim(V), c(q,q,n))) stop("Y and V have incompatible dimensions.")
  # prior
  if(is.null(prior)) prior <- list()
  Lambda <- .setDims(prior$Lambda, p, q)
  if(anyNA(Lambda)) stop("Lambda and (Y, X) have incompatible dimensions.")
  Omega <- .setDims(prior$Omega, p, p)
  if(anyNA(Omega)) stop("Omega and X have incompatible dimensions.")
  Psi <- .setDims(prior$Psi, q, q)
  if(anyNA(Psi)) stop("Psi and Y have incompatible dimensions.")
  nu <- if(is.null(prior$nu)) 0 else prior$nu
  if(length(nu) != 1) stop("nu must be a scalar.")
  # initial values
  if(is.null(init)) init <- list()
  Beta0 <- .setDims(init$Beta, p, q)
  if(anyNA(Beta0)) {
    stop("init$Beta and (Y, X) have incompatible dimensions.")
  }
  Sigma0 <- .setDims(init$Sigma, q)
  if(anyNA(Sigma0)) {
    stop("init$Sigma and Y have incompatible dimensions.")
  }
  Mu0 <- if(is.null(init$Mu)) Y else init$Mu
  if(!identical(dim(Mu0), c(n, q))) {
    stop("init$Mu and Y have incompatible dimensions.")
  }
  # burn-in
  if (missing(burn)) burn <- max(0.1, 1000)
  if (burn < 1) burn <- nsamples * burn
  burn <- floor(burn)
  # don't store if don't update
  if(!updateHyp) storeHyp <- FALSE
  if(!updateRE) storeRE <- FALSE
  # MCMC
  post <- HierUneqVModelGibbs(nSamples = as.integer(nsamples),
                              nBurn = as.integer(burn),
                              Y = Y, X = X, V = V,
                              Lambda = Lamda, Omega = Omega,
                              Psi = Psi, nu = nu,
                              Beta0 = Beta0, iSigma0 = .solveV(Sigma0),
                              Mu0 = Mu0,
                              updateBetaSigma = as.logical(updateHyp),
                              updateMu = as.logical(updateRE),
                              storeBetaSigma = as.logical(storeHyp),
                              storeMu = as.logical(storeRE))
  # output format
  out <- NULL
  if(storeHyp) {
    out <- list(Beta = array(post$Beta, dim = c(p, q, nsamples)),
                Sigma = array(post$Sigma, dim = c(q, q, nsamples)))
  }
  if(storeRE) {
    out <- c(out,
             list(Mu = aperm(array(post$Mu, dim = c(q, n, nsamples)),
                             perm = c(2,1,3))))
  }
  out
}

# solve for variance matrices
.solveV <- function(V, x) {
  C <- chol(V)
  if(missing(x)) x <- diag(nrow(V))
  backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
}

