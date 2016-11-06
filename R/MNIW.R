#--- MNIW distribution ----------------------------------------------------------

#' @title Generate Samples from the Matrix-Normal Inverse-Wishart Distribution.
#' @param n Number of samples.
#' @param Lambda \code{p x q} mean matrix of \code{X} (see details).
#' @param Sigma \code{p x p} conditional row variance matrix of \code{X}.
#' @param Psi \code{q x q} Precision matrix of \code{V}.
#' @param nu Scalar degrees of freedom of \code{V}.
#' @param prec Logical, whether or not to return \code{V} or \code{C = V^{-1}}.
#' @details the Matrix-Normal Inverse-Wishart (MNIW) distribution on \code{p x q} matrix \code{X} and \code{q x q} matrix \code{V} is:
#' \deqn{X | V \sim MN(Lambda, Sigma, V), \qquad V \sim iWish(Psi, nu).}
#' @export
rMNIW <- function(n, Lambda, Sigma, Psi, nu, prec = FALSE, debug = FALSE) {
  # get dimensions
  PQ <- .getPQ(Lambda = Lambda, Sigma = Sigma, Psi = Psi)
  p <- PQ[1]
  q <- PQ[2]
  if(debug) browser()
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  # format arguments
  Lambda <- .setDims(Lambda, p = p, q = q)
  if(anyNA(Lambda)) stop("Something went wrong.  Please report bug.")
  Sigma <- .setDims(Sigma, p = p)
  if(anyNA(Sigma)) stop("Sigma and Lambda have incompatible dimensions.")
  Psi <- .setDims(Psi, q = q)
  if(anyNA(Psi)) stop("Psi and Lambda have incompatible dimensions.")
  nu <- c(nu)
  # check lengths
  N <- .getN(p = p, q = q, Lambda = Lambda, Sigma = Sigma,
             Psi = Psi, nu = nu)
  if(length(N) > 2 || (length(N) == 2 && N[2] != n)) {
    stop("Arguments don't all have length n.")
  }
  XV <- .Call('mniw_GenerateMatrixNIW', PACKAGE = 'mniw',
              n, Lambda, Sigma, Psi, nu, prec)
  if(n > 1) {
    XV$X <- array(XV$X, dim = c(p,q,n))
    XV$V <- array(XV$V, dim = c(q,q,n))
  }
  XV
}
