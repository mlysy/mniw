#require(testthat)
library(mniw)
source("mniw-testfunctions.R")
context("Inner Products")

test_that("Inner product is same in C++ as R", {
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
