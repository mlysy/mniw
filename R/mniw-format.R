# format the MNIW inputs correctly
# that is, get dimensions of problem from the inputs
# for density evaluation, X and V are given and thus the dimensions of the problem are known immediately.
# for simulation, need to gather these values from the other inputs.
# what to do when arguments are vectors?
# always default to dimensions 1 x 1.  This is the obvious choice for variance arguments.  For mean arguments, the other choice is eg., length(arg) x 1, but this should be handled with the simpler functions mNorm and mNIW.
.getPQ <- function(X, V, Lambda, Sigma, Psi) {
  p <- NA
  q <- NA
  if(!missing(X)) {
    if(is.vector(X)) {
      p <- 1
      q <- 1
    } else {
      p <- dim(X)[1]
      q <- dim(X)[2]
    }
  }
  if(!missing(V)) {
    if(is.na(q)) {
      q <- ifelse(is.vector(V), 1, dim(V)[1])
    }
  }
  if(!missing(Lambda)) {
    if(is.na(p)) p <- ifelse(is.vector(Lambda), 1, dim(Lambda)[1])
    if(is.na(q)) q <- ifelse(is.vector(Lambda), 1, dim(Lambda)[2])
  }
  if(!missing(Sigma)) {
    if(is.na(p)) p <- ifelse(is.vector(Sigma), 1, dim(Sigma)[1])
  }
  if(!missing(Psi)) {
    if(is.na(q)) q <- ifelse(is.vector(Psi), 1, dim(Psi)[1])
  }
  c(p = p, q = q)
}

# format a vector, matrix, or array to a matrix with p rows.
# if dimensions are incompatible return NA, otherwise the formated matrix.
# if only
.setDims <- function(X, p, q) {
  var.X <- missing(p) || missing(q)
  if(var.X) {
    if(missing(p)) p <- q
    if(missing(q)) q <- p
  }
  if(missing(X)) {
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

# get the sample size from arguments for random sampling
# returns 1 and possibly all sample sizes > 1 detected.
.getN <- function(p, q, X, V, Lambda, Sigma, Psi, nu) {
  N <- NULL
  if(!missing(X)) N <- c(N, ncol(X)/q)
  if(!missing(V)) N <- c(N, ncol(V)/q)
  if(!missing(Lambda)) N <- c(N, ncol(Lambda)/q)
  if(!missing(Sigma)) N <- c(N, ncol(Sigma)/p)
  if(!missing(Psi)) N <- c(N, ncol(Psi)/q)
  if(!missing(nu)) N <- c(N, length(nu))
  N <- unique(sort(N))
  c(1, N[N>1])
}
