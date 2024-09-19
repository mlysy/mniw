#' Format the MNIW inputs correctly.
#'
#' Get dimensions of problem from the inputs.
#'
#' @details For density evaluation, `X` and `V` are given and thus the dimensions of the problem are known immediately.
#'
#' For simulation, need to gather these values from the other inputs.
#'
#' What to do when arguments are vectors?
#'
#' Always default to dimensions `1 x 1`.  This is the obvious choice for variance arguments.  For mean arguments, the other choice is eg., `length(arg) x 1`, but this should be handled with the simpler functions `mNorm` and `mNIW`.
#' @noRd
.getPQ <- function(X, V, Lambda, Sigma, Psi) {
  p <- NA
  q <- NA
  if(!.minu(X)) {
    if(is.vector(X)) {
      p <- 1
      q <- 1
    } else {
      p <- dim(X)[1]
      q <- dim(X)[2]
    }
  }
  if(!.minu(V)) {
    if(is.na(q)) {
      q <- ifelse(is.vector(V), 1, dim(V)[1])
    }
  }
  if(!.minu(Lambda)) {
    if(is.na(p)) p <- ifelse(is.vector(Lambda), 1, dim(Lambda)[1])
    if(is.na(q)) q <- ifelse(is.vector(Lambda), 1, dim(Lambda)[2])
  }
  if(!.minu(Sigma)) {
    if(is.na(p)) p <- ifelse(is.vector(Sigma), 1, dim(Sigma)[1])
  }
  if(!.minu(Psi)) {
    if(is.na(q)) q <- ifelse(is.vector(Psi), 1, dim(Psi)[1])
  }
  as.numeric(c(p, q))
}

#' Is `x` missing or NULL?
#'
#' @noRd
.minu <- function(x) missing(x) || is.null(x)


#' Format a vector, matrix, or array to a matrix with `p` rows and `q` columns.
#'
#' @details If dimensions are incompatible return `NA`, otherwise the formated matrix.
#' @noRd
.setDims <- function(X, p, q) {
  var.X <- missing(p) || missing(q)
  if(var.X) {
    if(missing(p)) p <- q
    if(missing(q)) q <- p
  }
  if(.minu(X)) {
    if(!var.X) {
      X <- matrix(0,p,q)
    } else {
      X <- diag(p)
    }
  }
  if(is.vector(X)) X <- array(X, dim = c(1,1,length(X)))
  if(!all(dim(X)[1:2] == c(p,q))) {
    X <- NA
  } else {
    X <- matrix(X,p)
  }
  X
}

#' Get the sample size from arguments for random sampling.
#'
#' @return A vector of integers concatenating `1` with all sample sizes > 1 detected.
#' @noRd
.getN <- function(p, q, X, V, Lambda, Sigma, Psi, nu) {
  N <- NULL
  if(!.minu(X)) N <- c(N, ncol(X)/q)
  if(!.minu(V)) N <- c(N, ncol(V)/q)
  if(!.minu(Lambda)) N <- c(N, ncol(Lambda)/q)
  if(!.minu(Sigma)) N <- c(N, ncol(Sigma)/p)
  if(!.minu(Psi)) N <- c(N, ncol(Psi)/q)
  if(!.minu(nu)) N <- c(N, length(nu))
  N <- unique(sort(N))
  c(1, N[N>1])
}

#' Convert a vector or matrix to MN format.
#'
#' @details Promote to 2-d or 3-d array with (second dimension) `q = 1`.
#' @noRd
.vec2mn <- function(x) {
  if(is.vector(x)) {
    x <- matrix(x, ncol = 1)
  } else {
    x <- t(x)
    x <- array(x, dim = c(dim(x)[1], 1, dim(x)[2]))
  }
  x
}
