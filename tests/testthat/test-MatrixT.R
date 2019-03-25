## require(testthat)
## library(mniw)
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
                          drop = c(TRUE, FALSE), stringsAsFactors = FALSE)
  ncases <- nrow(case.par)
  n <- 5
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    args <- list(X = list(p = p, q = q, rtype = cp$X, vtype = "matrix"),
                 Mu = list(p = p, q = q, rtype = cp$Mu, vtype = "matrix"),
                 RowV = list(p = p, rtype = cp$RowV, vtype = "matrix"),
                 ColV = list(q = q, rtype = cp$ColV, vtype = "matrix"),
                 nu = list(q = q, rtype = cp$nu, vtype = "scalar"))
    args <- get_args(n = n, args = args, drop = cp$drop)
    # R test
    llR <- rep(NA, n)
    for(jj in 1:n) {
      llR[jj] <- dMTR(X = args$R$X[[jj]],
                      Lambda = args$R$Mu[[jj]],
                      SigmaU = args$R$RowV[[jj]],
                      SigmaV = args$R$ColV[[jj]],
                      nu = args$R$nu[[jj]], log = TRUE)
    }
    # C++ test
    llcpp <- do.call(dMT, args = c(args$cpp, list(log = TRUE)))
    # if all inputs to c++ are single it returns only one value
    if(all_single(cp)) {
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
                          drop = c(TRUE, FALSE), stringsAsFactors = FALSE)
  case.par <- case.par[with(case.par, {
    !(Mu == "none" & ((RowV == "none") | (ColV == "none")))
  }),]
  ncases <- nrow(case.par)
  rownames(case.par) <- 1:ncases
  n <- 5
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  TestSeed <- sample(1e6, ncases)
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    args <- list(Mu = list(p = p, q = q, rtype = cp$Mu, vtype = "matrix"),
                 RowV = list(p = p, rtype = cp$RowV, vtype = "matrix"),
                 ColV = list(q = q, rtype = cp$ColV, vtype = "matrix"),
                 nu = list(q = q, rtype = cp$nu, vtype = "scalar"))
    args <- get_args(n = n, args = args, drop = cp$drop)
    # R test
    XR <- array(NA, dim = c(p,q,n))
    set.seed(TestSeed[ii]) # seed
    for(jj in 1:n) {
      XR[,,jj] <- rMTR(Lambda = args$R$Mu[[jj]],
                       SigmaU = args$R$RowV[[jj]],
                       SigmaV = args$R$ColV[[jj]],
                       nu = args$R$nu[[jj]], prec = FALSE)
    }
    # C++ test
    set.seed(TestSeed[ii]) # seed
    ## Xcpp <- do.call(rMT, args = arg.list)
    Xcpp <- do.call(rMT, args = c(args$cpp, list(n = n)))
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
