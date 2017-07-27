#require(testthat)
library(mniw)
source("mniw-testfunctions.R")
context("Cross-Products")

test_that("Cross-product X'VX is same in C++ as R", {
  case.par <- expand.grid(p = 1:4, q = 1:4,
                          single.X = c(TRUE, FALSE), single.V = c(TRUE, FALSE),
                          inverse = c(TRUE, FALSE))
  single.ind <- function(is.single) {
    ind <- if(is.single) 1 else 1:n
    ind
  }
  n <- 10
  for(ii in 1:nrow(case.par)) {
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    X <- replicate(n, rMnorm(p,q), simplify = FALSE)
    V <- replicate(n, crossprod(rMnorm(p)), simplify = FALSE)
    # R test
    ipR <- array(NA, dim = c(q,q,n))
    for(jj in 1:n) {
      if(cp$inverse) {
        VX <- solve(V[[ifelse(cp$single.V,1,jj)]], X[[ifelse(cp$single.X,1,jj)]])
      } else {
        VX <- V[[ifelse(cp$single.V,1,jj)]] %*% X[[ifelse(cp$single.X,1,jj)]]
      }
      ipR[,,jj] <- crossprod(X[[ifelse(cp$single.X,1,jj)]], VX)
    }
    # C++ test
    X <- array(unlist(X), dim = c(p,q,n))
    V <- array(unlist(V), dim = c(p,p,n))
    ipcpp <- crossprodV(X = X[,,single.ind(cp$single.X),drop=FALSE],
                       V = V[,,single.ind(cp$single.V),drop=FALSE],
                       inverse = cp$inverse)
    if(all(unlist(cp[c("single.X", "single.V")]))) {
      ipcpp <- array(ipcpp, c(q,q,n))
    }
    expect_equal(ipR, ipcpp)
  }
})

test_that("Cross-product X'VY is same in C++ as R", {
  case.par <- expand.grid(p = 1:4, q = 1:4, r = 1:4,
                          single.X = c(TRUE, FALSE),
                          single.V = c(TRUE, FALSE),
                          single.Y = c(TRUE, FALSE),
                          inverse = c(TRUE, FALSE))
  single.ind <- function(is.single) {
    ind <- if(is.single) 1 else 1:n
    ind
  }
  n <- 10
  for(ii in 1:nrow(case.par)) {
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    r <- cp$r
    sX <- cp$single.X
    sY <- cp$single.Y
    sV <- cp$single.V
    X <- replicate(n, rMnorm(p,q), simplify = FALSE)
    V <- replicate(n, crossprod(rMnorm(p)), simplify = FALSE)
    Y <- replicate(n, rMnorm(p,r), simplify = FALSE)
    # R test
    ipR <- array(NA, dim = c(q,r,n))
    for(jj in 1:n) {
      if(cp$inverse) {
        VY <- solve(V[[ifelse(sV,1,jj)]], Y[[ifelse(sY,1,jj)]])
      } else {
        VY <- V[[ifelse(sV,1,jj)]] %*% Y[[ifelse(sY,1,jj)]]
      }
      ipR[,,jj] <- crossprod(X[[ifelse(sX,1,jj)]], VY)
    }
    # C++ test
    X <- array(unlist(X), dim = c(p,q,n))
    V <- array(unlist(V), dim = c(p,p,n))
    Y <- array(unlist(Y), dim = c(p,r,n))
    ipcpp <- crossprodV(X = X[,,single.ind(sX),drop=FALSE],
                        V = V[,,single.ind(sV),drop=FALSE],
                        Y = Y[,,single.ind(sY),drop=FALSE],
                       inverse = cp$inverse)
    if(all(unlist(cp[c("single.X", "single.V", "single.Y")]))) {
      ipcpp <- array(ipcpp, c(q,r,n))
    }
    #message("ii = ", ii, ", max(abs) = ", signif(max(abs(ipR - ipcpp)),2))
    expect_equal(ipR, ipcpp)
  }
})
