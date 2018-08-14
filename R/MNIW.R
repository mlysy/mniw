#--- MNIW distribution -----------------------------------------------------

#' Generate samples from the Matrix-Normal Inverse-Wishart distribution.
#'
#' @param n number of samples.
#' @param Lambda \code{p x q} mean matrix of \code{X} (see details).
#' @param Sigma \code{p x p} conditional row variance matrix of \code{X}.
#' @param Psi \code{q x q} scale matrix of \code{V}.
#' @param nu scalar degrees of freedom of \code{V}.
#' @param prec logical. Whether or not to return \code{V} or \code{C = V^{-1}}.
#' @details The Matrix-Normal Inverse-Wishart (MNIW) distribution on \code{p x q} matrix \code{X} and \code{q x q} matrix \code{V} is:
#' \deqn{X | V ~ MN(Lambda, Sigma, V),}
#' \deqn{V ~ iWish(Psi, nu).}
#'
#' @return A list with the folllwing elements:
#' \describe{
#' \item{X}{Random samples for X ~ MN(Lambda, Sigma, V).}
#' \item{V}{Random samples for V ~ iWish(Psi, nu).}
#' }
#'
#' @examples
#' ## Sample from MNIW distribution
#' n = 100
#' p = 2
#' q = 2
#' Lambda = matrix(c(1,-0.5,-1,1),p,q)
#' Sigma = matrix(c(1,0.5,0.5,1),p,p)
#' Psi = matrix(c(1,0.1,0.1,1),p,q)
#' nu = q + 1
#' rMNIW(n, Lambda, Sigma, Psi, nu)
#'
#' @export
rMNIW <- function(n, Lambda, Sigma, Psi, nu, prec = FALSE) {
  # get dimensions
  PQ <- .getPQ(Lambda = Lambda, Sigma = Sigma, Psi = Psi)
  p <- PQ[1]
  q <- PQ[2]
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
  XV <- GenerateMatrixNIW(n, Lambda, Sigma, Psi, nu, prec)
  if(n > 1) {
    XV$X <- array(XV$X, dim = c(p,q,n))
    XV$V <- array(XV$V, dim = c(q,q,n))
  }
  XV
}
