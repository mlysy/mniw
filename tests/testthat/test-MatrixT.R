## require(testthat)
library(mniw)
source("mniw-testfunctions.R")
context("Matrix-T Distribution")

tol <- 1e-6

test_that("Matrix-T density is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(1,2,4), q = c(1,2,3),
                          X = c("single", "multi"),
                          Mu = c("none", "single", "multi"),
                          RowV = c("none", "single", "multi"),
                          ColV = c("none", "single", "multi"),
                          nu = c("single", "multi"),
                          drop = c(TRUE, FALSE))
  ncases <- nrow(case.par)
  n <- 10
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    sX <- cp$X != "multi"
    noX <- cp$X == "none"
    sMu <- cp$Mu != "multi"
    noMu <- cp$Mu == "none"
    sRowV <- cp$RowV != "multi"
    noRowV <- cp$RowV == "none"
    sColV <- cp$ColV != "multi"
    noColV <- cp$ColV == "none"
    snu <- cp$nu != "multi"
    arr.drop <- cp$drop
    X <- rMM(n, p = p, q = q, noArg = noX)
    Mu <- rMM(n, p = p, q = q, noArg = noMu)
    RowV <- rMM(n, p = p, noArg = noRowV)
    ColV <- rMM(n, q = q, noArg = noColV)
    nu <- sapply(runif(n, 2*q, 3*q),c,simplify = FALSE)
    # R test
    llR <- rep(NA, n)
    for(jj in 1:n) {
      X_R <- X[[ifelse(sX, 1, jj)]]
      Mu_R <- Mu[[ifelse(sMu, 1, jj)]]
      RowV_R <- RowV[[ifelse(sRowV, 1, jj)]]
      ColV_R <- ColV[[ifelse(sColV, 1, jj)]]
      nu_R <- nu[[ifelse(snu, 1, jj)]]
      llR[jj] <- dMTR(X = X_R, Lambda = Mu_R, SigmaU = RowV_R,
                      SigmaV = ColV_R, nu = nu_R, log = TRUE)
    }
    # C++ test
    X_cpp <- unlistM(X, sX, arr.drop)
    Mu_cpp <- unlistM(Mu, sMu, arr.drop)
    RowV_cpp <- unlistM(RowV, sRowV, arr.drop)
    ColV_cpp <- unlistM(ColV, sColV, arr.drop)
    nu_cpp <- unlist(nu)
    if(snu) nu_cpp <- nu_cpp[1]
    arg.list <- list(X = X_cpp, Mu = Mu_cpp,
                     RowV = RowV_cpp, ColV = ColV_cpp)
    arg.list <- arg.list[!c(noX, noMu, noRowV, noColV)]
    arg.list <- c(arg.list, list(nu = nu_cpp, log = TRUE))
    llcpp <- do.call(dMT, args = arg.list)
    # if all inputs to c++ are single it returns only one value
    if(all(c(sX, sMu, sRowV, sColV))) {
      llcpp <- rep(llcpp, n)
    }
    mx <- abs(llR-llcpp)
    mx <- min(max(mx), max(mx/abs(llR)))
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      ## expect_equal(mx, 0, tolerance = tol)
      expect_Rcpp_equal("dMT", ii, mx, tolerance = tol)
    }
  }
})

test_that("Matrix-T sampling is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(1,2,4), q = c(1,2,3),
                          Mu = c("none", "single", "multi"),
                          RowV = c("none", "single", "multi"),
                          ColV = c("none", "single", "multi"),
                          nu = c("single", "multi"),
                          drop = c(TRUE, FALSE))
  case.par <- case.par[with(case.par, {
    !(Mu == "none" & ((RowV == "none") | (ColV == "none")))
  }),]
  ncases <- nrow(case.par)
  rownames(case.par) <- 1:ncases
  n <- 10
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  TestSeed <- sample(1e6, ncases)
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    sMu <- cp$Mu != "multi"
    noMu <- cp$Mu == "none"
    sRowV <- cp$RowV != "multi"
    noRowV <- cp$RowV == "none"
    sColV <- cp$ColV != "multi"
    noColV <- cp$ColV == "none"
    snu <- cp$nu != "multi"
    arr.drop <- cp$drop
    Mu <- rMM(n, p = p, q = q, noArg = noMu)
    RowV <- rMM(n, p = p, noArg = noRowV)
    ColV <- rMM(n, q = q, noArg = noColV)
    nu <- sapply(runif(n, 2*q, 3*q),c,simplify = FALSE)
    # R test
    XR <- array(NA, dim = c(p,q,n))
    set.seed(TestSeed[ii]) # seed
    for(jj in 1:n) {
      Mu_R <- Mu[[ifelse(sMu, 1, jj)]]
      RowV_R <- RowV[[ifelse(sRowV, 1, jj)]]
      ColV_R <- ColV[[ifelse(sColV, 1, jj)]]
      nu_R <- nu[[ifelse(snu, 1, jj)]]
      XR[,,jj] <- rMTR(Lambda = Mu_R, SigmaU = RowV_R,
                       SigmaV = ColV_R, nu = nu_R, prec = FALSE)
    }
    # C++ test
    Mu_cpp <- unlistM(Mu, sMu, arr.drop)
    RowV_cpp <- unlistM(RowV, sRowV, arr.drop)
    ColV_cpp <- unlistM(ColV, sColV, arr.drop)
    nu_cpp <- unlist(nu)
    if(snu) nu_cpp <- nu_cpp[1]
    arg.list <- list(Mu = Mu_cpp,
                     RowV = RowV_cpp, ColV = ColV_cpp)
    arg.list <- arg.list[!c(noMu, noRowV, noColV)]
    arg.list <- c(arg.list, list(n = n, nu = nu_cpp))
    set.seed(TestSeed[ii]) # seed
    Xcpp <- do.call(rMT, args = arg.list)
    mx <- arDiff(XR, Xcpp)
    mx <- min(max(mx), max(mx/abs(XR)))
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      ## expect_equal(mx, 0, tolerance = tol)
      expect_Rcpp_equal("rMT", ii, mx, tolerance = tol)
    }
  }
})

#--- scratch --------------------------------------------------------------

## p <- 3
## q <- 4

## Z <- rMM(n = 1, p = p, q = q, noArg = FALSE)[[1]]
## RowV <- rMM(n = 1, p = p, noArg = FALSE)[[1]]
## ColV <- rMM(n = 1, q = q, noArg = FALSE)[[1]]

## ldet(diag(p) + solve(RowV, Z) %*% solve(ColV, t(Z)))
## ldet(ColV + crossprod(Z, solve(RowV, Z))) - ldet(ColV)
## ldet(RowV + Z %*% solve(ColV, t(Z))) - ldet(RowV)
