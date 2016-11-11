#require(testthat)
library(mniw)
source("mniw-testfunctions.R")
context("Matrix Normal Distribution")

test_that("Matrix Normal density is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(1,2,4), q = c(1,2,3),
                          X = c("single", "multi"),
                          Mu = c("none", "single", "multi"),
                          RowV = c("none", "single", "multi"),
                          ColV = c("none", "single", "multi"),
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
    arr.drop <- cp$drop
    X <- rMM(n, p = p, q = q, noArg = noX)
    Mu <- rMM(n, p = p, q = q, noArg = noMu)
    RowV <- rMM(n, p = p, noArg = noRowV)
    ColV <- rMM(n, q = q, noArg = noColV)
    # R test
    llR <- rep(NA, n)
    for(jj in 1:n) {
      llR[jj] <- dMNormR(X = X[[ifelse(sX, 1, jj)]],
                         Lambda = Mu[[ifelse(sMu, 1, jj)]],
                         SigmaU = RowV[[ifelse(sRowV, 1, jj)]],
                         SigmaV = ColV[[ifelse(sColV, 1, jj)]],
                         log = TRUE,
                         debug = FALSE)
    }
    # C++ test
    X <- unlistM(X, sX, arr.drop)
    Mu <- unlistM(Mu, sMu, arr.drop)
    RowV <- unlistM(RowV, sRowV, arr.drop)
    ColV <- unlistM(ColV, sColV, arr.drop)
    arg.list <- list(X = X, Mu = Mu, RowV = RowV, ColV = ColV)
    arg.list <- arg.list[!c(noX, noMu, noRowV, noColV)]
    llcpp <- do.call(dMNorm, args = c(arg.list, list(log = TRUE)))
    # if all inputs to c++ are single it returns only one value
    if(all(c(sX, sMu, sRowV, sColV))) {
      llcpp <- rep(llcpp, n)
    }
    mx <- abs(llR-llcpp)
    mx <- min(max(mx), max(mx/abs(llR)))
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      expect_equal(mx, 0)
    }
  }
})

test_that("Matrix Normal simulation is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(1,2,4), q = c(1,2,3),
                          Mu = c("none", "single", "multi"),
                          RowV = c("none", "single", "multi"),
                          ColV = c("none", "single", "multi"),
                          drop = c(TRUE, FALSE))
  # remove cases where dimensions can't be identified
  case.par <- case.par[!with(case.par, {
    Mu == "none" & ((RowV == "none") | (ColV == "none"))}),]
  ncases <- nrow(case.par)
  rownames(case.par) <- 1:ncases
  n <- 10
  TestSeed <- sample(1e6, ncases)
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:nrow(case.par)) {
    set.seed(TestSeed[ii]) # seed
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    sMu <- cp$Mu != "multi"
    noMu <- cp$Mu == "none"
    sRowV <- cp$RowV != "multi"
    noRowV <- cp$RowV == "none"
    sColV <- cp$ColV != "multi"
    noColV <- cp$ColV == "none"
    arr.drop <- cp$drop
    Mu <- rMM(n, p = p, q = q, noArg = noMu)
    RowV <- rMM(n, p = p, noArg = noRowV)
    ColV <- rMM(n, q = q, noArg = noColV)
    # R test
    set.seed(TestSeed[ii]) # seed
    llR <- array(NA, dim = c(p,q,n))
    for(jj in 1:n) {
      llR[,,jj] <- rMNormR(Lambda = Mu[[ifelse(sMu, 1, jj)]],
                           SigmaU = RowV[[ifelse(sRowV, 1, jj)]],
                           SigmaV = ColV[[ifelse(sColV, 1, jj)]])
    }
    # C++ test
    Mu <- unlistM(Mu, sMu, arr.drop)
    RowV <- unlistM(RowV, sRowV, arr.drop)
    ColV <- unlistM(ColV, sColV, arr.drop)
    arg.list <- list(Mu = Mu, RowV = RowV, ColV = ColV)
    arg.list <- arg.list[!c(noMu, noRowV, noColV)]
    set.seed(TestSeed[ii]) # seed
    llcpp <- do.call(rMNorm, args = c(arg.list, list(n = n)))
    #llcpp <- rMNorm(n = n, Mu = Mu, RowV = RowV, ColV = ColV)
    mx <- abs(range(llR - llcpp))
    mx <- min(max(mx), max(mx/abs(llR)))
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      expect_equal(mx, 0)
    }
  }
})
