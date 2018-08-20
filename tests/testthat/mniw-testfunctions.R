########################

# R functions for the mniw package

############################

# matrix of iid N(0,1)
rMnorm <- function(p, q) {
  if(missing(q)) q <- p
  matrix(rnorm(p*q), p, q)
}

# lower and upper triangular elements of a matrix
triL <- function(M) {
  M[upper.tri(M)] <- 0
  M
}
triU <- function(M) {
  M[lower.tri(M)] <- 0
  M
}

# log-determinant
ldet <- function(V) {
  determinant(V, log = TRUE)$mod[1]
}

# inner product t(X) %*% {V, V^{-1}} %*% X
# written in C++ now.
#innerprod <- function(X, V, inv = TRUE) {
#  if(inv) return(crossprod(X, solve(V, X)))
#  else return(crossprod(X, V %*% X))
#}

# log multi-gamma function
lmgamma <- function(x, p) p*(p-1)/4 * log(pi) + sum(lgamma(x + (1-1:p)/2))

## density of inverse wishart
#diwish <- function(V, Psi, nu, log = FALSE) {
#  d <- nrow(V)
#  ans <- (nu+d+1) * determinant(V, logarithm = TRUE)$modulus[1]
#  if(!all(Psi == 0)) {
#    ans <- ans + sum(diag(solve(V, Psi)))
#    ans <- ans - nu * determinant(Psi, logarithm = TRUE)$modulus[1]
#    ans <- ans + nu*d * log(2) + 2*lmgamma(nu/2, d)
#  }
#  ans <- -0.5 * ans
#  if(!log) ans <- exp(ans)
#  ans
#}

# density of wishart and inverse wishart
dwishR <- function(X, Psi, nu, inverse = FALSE, log = FALSE,
                   debug = FALSE) {
  if(length(X) == 1) {
    X <- as.matrix(X)
    Psi <- as.matrix(Psi)
  }
  d <- nrow(X)
  if(debug) browser()
  if(!inverse) {
    ans <- (nu-d-1) * ldet(X)
    ans <- ans - nu * ldet(Psi)
    ans <- ans - sum(diag(solve(Psi, X)))
  } else {
    ans <- -(nu+d+1) * ldet(X)
    ans <- ans + nu * ldet(Psi)
    ans <- ans - sum(diag(solve(X, Psi)))
  }
  ans <- ans - nu*d * log(2)
  ans <- 0.5 * ans - lmgamma(nu/2, d)
  if(!log) ans <- exp(ans)
  ans
}

# simulation of wishart/inverse wishart
rwishR <- function(Psi, nu, inverse = FALSE) {
  q <- nrow(Psi)
  XL <- matrix(0, q,q)
  # populate normals
  for(ii in 1:q) {
    XL[ii,ii] <- sqrt(rchisq(1, df = nu-ii+1))
    if(ii > 1) {
      for(jj in 1:(ii-1)) {
        XL[ii,jj] <- rnorm(1)
      }
    }
  }
  if(!inverse) {
    PsiL <- t(chol(Psi))
    XL <- PsiL %*% XL
    X <- XL %*% t(XL)
  } else {
    PsiL <- solve(t(chol(solve(Psi))))
    XL <- solve(PsiL, XL)
    X <- solve(XL %*% t(XL))
  }
  X
}

# density of matrix normal
dMNormR <- function(X, Lambda, SigmaU, SigmaV, log = FALSE,
                    debug = FALSE) {
  X <- as.matrix(X)
  Lambda <- as.matrix(Lambda)
  SigmaU <- as.matrix(SigmaU)
  SigmaV <- as.matrix(SigmaV)
  if(debug) browser()
  n <- nrow(X)
  p <- ncol(X)
  ans <- n*p*log(2*pi) + sum(diag(solve(SigmaV, t(X-Lambda)) %*% solve(SigmaU, X-Lambda)))
  ans <- ans + n*determinant(SigmaV, log = TRUE)$mod[1] + p*determinant(SigmaU, log = TRUE)$mod[1]
  ans <- -ans/2
  if(!log) ans <- exp(ans)
  ans
}

# simulation of matrix normal
rMNormR <- function(Lambda, SigmaU, SigmaV) {
  p <- nrow(Lambda)
  q <- ncol(Lambda)
  Z <- t(rMnorm(q, p))
  X <- Z %*% chol(SigmaV)
  X <- t(chol(SigmaU)) %*% X
  X + Lambda
}

# density of mniw distribution
dmniwR <- function(X, V, Lambda, Sigma, Psi, nu, log = FALSE) {
  ans <- diwish(V = V, Psi = Psi, nu = nu, log = TRUE)
  if(!all(Sigma == 0)) ans <- ans + dMnorm(X, Lambda, Sigma, V, log = TRUE)
  if(!log) ans <- exp(ans)
  ans
}

# simulation of mniw. row variance can be on precision scale by seting Omega.
rmniwR <- function(Lambda, Sigma, Psi, nu, prec = FALSE, debug = FALSE) {
  p <- nrow(Lambda)
  q <- ncol(Lambda)
  V <- rwishR(Psi = Psi, nu = nu, inverse = TRUE)
  CL <- t(chol(solve(V)))
  Z <- t(rMnorm(q, p))
  Z <- Z %*% solve(CL)
  if(!prec) {
    X <- t(chol(Sigma)) %*% Z + Lambda
  } else {
    X <- solve(chol(Sigma)) %*% Z + Lambda
  }
  list(X = X, V = V)
  #PsiL <- solve(t(chol(solve(Psi))))
  #q <- nrow(PsiL)
  #for(ii in 1:q) {
  #  XL[ii,ii] <- sqrt(rchisq(1, df = nu-ii+1))
  #  if(ii > 1) {
  #    for(jj in 1:(ii-1)) {
  #      XL[ii,jj] <- rnorm(1)
  #    }
  #  }
  #}
  #XL <- solve(XiL, XL)
}

# random effects normal simulation in R
rmNormRER <- function(y, V, lambda, A) {
  C <- solveV(V)
  Q <- solveV(A)
  G <- C + Q
  GU <- chol(G)
  z <- rnorm(length(y))
  mu <- backsolve(r = GU, x = Q %*% (lambda - y), transpose = TRUE)
  drop(backsolve(r = GU, x = z + mu) + y)
}

# reverse cholesky decomposition x = L'L
revchol <- function(x) solve(t(chol(solve(x)))) # L

# convert a list of matrices into an array.
# optionally keep only first element (sArg = TRUE)
# optionally drop the last dimension if it is == 1.
unlistM <- function(X, sArg, dropLast) {
  n <- length(X)
  PQ <- dim(X[[1]])
  X <- array(unlist(X), dim = c(PQ, length(X)))
  if(sArg) X <- X[,,1,drop=FALSE]
  if(dropLast) {
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
# if noArg defaults to a matrix of zeros or identity.
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

# lower of maximum abs/rel error
arDiff <- function(x0, xhat) {
  mx <- abs(xhat - x0)
  min(max(mx), max(mx/abs(x0)))
}

# solve for variance matrices
solveV <- function(V, x) {
  C <- chol(V)
  if(missing(x)) x <- diag(nrow(V))
  backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
}

# check R vs cpp code
expect_Rcpp_equal <- function(fun, icase, mx, ...) {
  Rcpp_comp <- function(fun, icase) mx
  eval(bquote({
    expect_equal(object = Rcpp_comp(.(fun), .(icase)),
                 expected = .(rep(0, length(mx))), ...)
  }))
}
