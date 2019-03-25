# lower and upper triangular elements of a matrix
triL <- function(M) {
  M[upper.tri(M)] <- 0
  M
}
triU <- function(M) {
  M[lower.tri(M)] <- 0
  M
}

# reverse cholesky decomposition x = L'L
revchol <- function(x) solve(t(chol(solve(x)))) # L

# convert a list of matrices into an array.
# optionally keep only first element (sArg = TRUE)
# optionally drop the last dimension if it is == 1.
unlistM <- function(X, sArg, drop) {
  n <- length(X)
  PQ <- dim(X[[1]])
  X <- array(unlist(X), dim = c(PQ, length(X)))
  if(sArg) X <- X[,,1,drop=FALSE]
  if(drop) {
    if(dim(X)[3] == 1) X <- array(X, dim = PQ)
  }
  X
}

# convert a list of vectors into a matrix.
# optionally keep only first element (sArg = TRUE)
unlistV <- function(X, sArg) {
  X <- do.call(rbind, X)
  if(sArg) X <- X[1,]
  X
}

# drop the last dimension of an array if it is == 1.
dropN <- function(X) {
  dd <- dim(X)
  nd <- length(dd)
  if(dd[nd] == 1) X <- array(X, dim = dd[1:(nd-1)])
  X
}

# create n random p x q matrices.
# If only one of p or q provided the random matrix is a variance.
# if noArg = TRUE, defaults to matrix of zeros or identity matrix (for variance).
rMM <- function(n, p, q, noArg) {
  var.X <- missing(p) || missing(q)
  if(missing(p)) p <- q
  if(missing(q)) q <- p
  if(noArg) {
    if(!var.X) {
      X <- replicate(n, matrix(0,p,q), simplify = FALSE)
    } else {
      X <- replicate(n, diag(p), simplify = FALSE)
    }
  } else {
    if(!var.X) {
      X <- replicate(n, rMnorm(p,q), simplify = FALSE)
    } else {
      X <- replicate(n, crossprod(rMnorm(p)), simplify = FALSE)
    }
  }
  X
}
