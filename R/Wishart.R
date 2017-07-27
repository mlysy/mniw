#---------------------------------------------------------------------------
#
# Wishart and Inverse-Wishart distributions
#
#---------------------------------------------------------------------------

#' @name Wishart
#' @title Wishart and Inverse-Wishart Distributions.
#' @description Densities and Random Sampling for the Wishart and Inverse-Wishart distributions.
#' @aliases dwish rwish diwish riwish dwishart rwishart
#' @param X Argument to the density function.  Either a \code{q x q} matrix or a \code{q x q x n} array.
#' @param n Number of random matrices to generate.
#' @param Psi Scale parameter.  Either a \code{q x q} matrix or a \code{q x q x n} array.
#' @param nu Degrees-of-freedom parameter.  A scalar or vector.
#' @param inverse Whether or not to use the Inverse-Wishart distribution.
#' @param log Whether or not to compute the log-density.
#' @details \code{dwish} and \code{diwish} are convenience wrappers for \code{dwishart}, and similarly \code{rwish} and \code{riwish} are wrappers for \code{rwishart}.
#'
#' The distribution of a \code{q x q} Wishart random matrix is
#' \deqn{\code{f(X) = const |X|^{(\nu-q-1)/2} * \exp{-0.5 * Tr(\Psi^{-1} X)}},}
#' where \code{Psi} is a symmetric positive-definite matrix and \code{nu > q-1}.
#'
#' The distribution of an Inverse-Wishart random matrix is
#' \deqn{\code{f(X) \propto |X|^{-(\nu+q+1)/2} * \exp{-0.5 * Tr(\Psi X^{-1})}}.}
#' If \code{X ~ Wish(Psi, nu)}, then \code{X^{-1} ~ iWish(Psi^{-1}, nu)}.
#'
#' Quadratic forms involving a nonrandom vector \code{a} and a Wishart random matrix \code{X} have a chi-square distribution.  That is, for fixed \code{a} and \code{X ~ Wish(Psi, nu)},
#' \deqn{\code{(a' X a) / (a' Psi a) ~ chisq(df = nu)}.}
#' @return A vector for densities, or a \code{q x q x n} array for random sampling.

#--- convenience wrappers --------------------------------------------------

# wishart density
#' @rdname Wishart
#' @export
dwish <- function(X, Psi, nu, log = FALSE) {
  dwishart(X, Psi, nu, inverse = FALSE, log)
}

# wishart simulation
#' @rdname Wishart
#' @export
rwish <- function(n, Psi, nu) {
  rwishart(n, Psi, nu, inverse = FALSE)
}

# inverse wishart density
#' @rdname Wishart
#' @export
diwish <- function(X, Psi, nu, log = FALSE) {
  dwishart(X, Psi, nu, inverse = TRUE, log)
}

# inverse wishart simulation
#' @rdname Wishart
#' @export
riwish <- function(n, Psi, nu) {
  rwishart(n, Psi, nu, inverse = TRUE)
}

#--- lower level functions -------------------------------------------------

# density of wishart and inverse wishart

#' @rdname Wishart
#' @export
dwishart <- function(X, Psi, nu, inverse = FALSE, log = FALSE, debug = FALSE) {
  if(debug) browser()
  # get dimensions
  PQ <- .getPQ(X = X, Psi = Psi)
  if(is.na(PQ[2])) stop("Undetermined problem dimensions.")
  q <- PQ[2]
  # format arguments
  X <- .setDims(X, q = q)
  Psi <- .setDims(Psi, q = q)
  nu <- c(nu)
  if(anyNA(X)) stop("Something went wrong.  Please report bug.")
  if(anyNA(Psi)) stop("Psi and X have incompatible dimensions.")
  if(length(X) == 1) X <- as.matrix(X)
  if(length(Psi) == 1) Psi <- as.matrix(Psi)
  # check lengths
  N <- .getN(q = q, X = X, Psi = Psi, nu = nu)
  ## check X
  #q <- dim(X)[1]
  #if(dim(X)[2] != q) stop("X has non-square dimensions.")
  #X <- matrix(X,q)
  ## check Psi
  #if(missing(Psi)) Psi <- diag(q)
  #if(!all(dim(Psi)[1:2] == q)) stop("X and Psi have incompatible dimensions.")
  #Psi <- matrix(Psi,q)
  ## check nu
  #nu <- c(nu)
  ## check lengths
  #N <- unique(sort(c(ncol(X)/q, ncol(Psi)/q, length(nu))))
  #N <- c(1, N[N>1])
  if(length(N) > 2) stop("Arguments have different lengths.")
  ans <- LogDensityWishart(X, Psi, nu, inverse)
  if(!log) ans <- exp(ans)
  ans
}

# random sampling for wishart and inverse wishart
#' @rdname Wishart
#' @export
rwishart <- function(n, Psi, nu, inverse = FALSE) {
  # get problem dimensions
  PQ <- .getPQ(Psi = Psi)
  if(is.na(PQ[2])) stop("Undetermined problem dimensions.")
  q <- PQ[2]
  # format arguments
  Psi <- .setDims(Psi, q = q)
  nu <- c(nu)
  # check lengths
  N <- .getN(q = q, Psi = Psi, nu = nu)
  ## determine the dimensions
  #q <- dim(Psi)[1]
  #if(dim(Psi)[2] != q) stop("Psi has non-square dimensions.")
  #Psi <- matrix(Psi,q)
  #nu <- c(nu)
  ## determine lengths
  #N <- sort(unique(c(ncol(Psi)/q, length(nu))))
  #N <- c(1,N[N>1])
  if(length(N) > 2 || (length(N) == 2 && N[2] != n))
    stop("Arguments don't all have length n.")
  X <- GenerateWishart(n, Psi, nu, inverse)
  if(n > 1) X <- array(X, dim = c(q,q,n))
  X
}

