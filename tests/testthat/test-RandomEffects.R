library(mniw)
source("mniw-testfunctions.R")
context("Random-Effects Normal Distribution")

tol <- 1e-6

test_that("Random-Effects Normal sampling is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(q = c(1,2,4),
                          y = c("single", "multi"),
                          V = c("single", "multi"),
                          lambda = c("single", "multi"),
                          A = c("single", "multi"),
                          drop = c(TRUE, FALSE))
  ncases <- nrow(case.par)
  n <- 12 # number of random draws
  test.seed <- sample(1e6, ncases)
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    q <- cp$q
    sy <- cp$y == "single"
    sV <- cp$V == "single"
    slambda <- cp$lambda == "single"
    sA <- cp$A == "single"
    arr.drop <- cp$drop
    y <- rMM(n = n, p = 1, q = q, noArg = FALSE)
    y <- lapply(y, drop)
    V <- rMM(n = n, q = q, noArg = FALSE)
    lambda <- rMM(n = n, p = 1, q = q, noArg = FALSE)
    lambda <- lapply(lambda, drop)
    A <- rMM(n = n, q = q, noArg = FALSE)
    # R test
    muR <- matrix(NA, n, q)
    set.seed(test.seed[ii])
    for(jj in 1:n) {
      muR[jj,] <- rmNormRER(y = y[[ifelse(sy, 1, jj)]],
                            V = V[[ifelse(sV, 1, jj)]],
                            lambda = lambda[[ifelse(slambda, 1, jj)]],
                            A = A[[ifelse(sA, 1, jj)]])
    }
    # C++ test
    y <- unlistV(y, sy)
    ## y <- drop(unlistM(y, sy, arr.drop))
    ## if(!sy) y <- t(y)
    V <- unlistM(V, sV, arr.drop)
    lambda <- unlistV(lambda, slambda)
    ## lambda <- drop(unlistM(lambda, sy, arr.drop))
    ## if(!slambda) lambda <- t(lambda)
    A <- unlistM(A, sA, arr.drop)
    set.seed(test.seed[ii])
    mucpp <- rmNormRE(n = n, y = y, V = V, lambda = lambda, A = A)
    mx <- arDiff(muR, mucpp)
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      ## expect_equal(mx, 0, tolerance = tol)
      expect_Rcpp_equal("rmNormRE", ii, mx, tolerance = tol)
    }
  }
})
