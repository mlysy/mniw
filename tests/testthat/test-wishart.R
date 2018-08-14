#require(testthat)
library(mniw)
source("mniw-testfunctions.R")
context("Wishart and Inverse Wishart Distributions")

tol <- 1e-6

test_that("Wishart density is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(q = c(1,2,4),
                          X = c("single", "multi"),
                          Psi = c("none", "single", "multi"),
                          nu = c("single", "multi"),
                          inverse = c(TRUE, FALSE),
                          drop = c(TRUE, FALSE))
  ncases <- nrow(case.par)
  n <- 10
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    q <- cp$q
    sX <- cp$X != "multi"
    noX <- cp$X == "none"
    sPsi <- cp$Psi != "multi"
    noPsi <- cp$Psi == "none"
    snu <- cp$nu != "multi"
    arr.drop <- cp$drop
    X <- rMM(n = n, q = q, noArg = noX)
    Psi <- rMM(n = n, q = q, noArg = noPsi)
    nu <- sapply(runif(n, q, 2*q),c,simplify = FALSE)
    # R test
    llR <- rep(NA, n)
    for(jj in 1:n) {
      llR[jj] <- dwishR(X = X[[ifelse(sX, 1, jj)]],
                        Psi = Psi[[ifelse(sPsi, 1, jj)]],
                        nu = nu[[ifelse(snu, 1, jj)]],
                        inverse = cp$inverse, log = TRUE)
    }
    # C++ test
    X <- unlistM(X, sX, arr.drop)
    Psi <- unlistM(Psi, sPsi, arr.drop)
    nu <- unlist(nu)
    if(snu) nu <- nu[1]
    arg.list <- list(X = X, Psi = Psi)
    arg.list <- arg.list[!c(noX, noPsi)]
    llcpp <- do.call(dwishart,
                     args = c(arg.list,
                       list(inverse = cp$inverse, nu = nu, log = TRUE)))
    #llcpp <- dwishart(X = X, Psi = Psi, nu = nu,
    #                  inverse = cp$inverse, log = TRUE)
    # C++ produces single output if all inputs are single
    if(all(c(sX, sPsi, snu))) {
      llcpp <- rep(llcpp, n)
    }
    mx <- abs(llR - llcpp)
    mx <- min(max(mx), max(mx/abs(llR)))
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      expect_equal(mx, 0, tolerance = tol)
    }
  }
})

test_that("Wishart sampling is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(q = c(1,2,4),
                          Psi = c("single", "multi"),
                          nu = c("single", "multi"),
                          inverse = c(TRUE, FALSE),
                          drop = c(TRUE, FALSE))
  ncases <- nrow(case.par)
  n <- 10
  if(calc.diff) {
    MaxDiff <- rep(NA, ncases)
  }
  for(ii in 1:ncases) {
    cp <- case.par[ii,]
    q <- cp$q
    sPsi <- cp$Psi != "multi"
    noPsi <- cp$Psi == "none"
    snu <- cp$nu != "multi"
    arr.drop <- cp$drop
    Psi <- rMM(n = n, q = q, noArg = noPsi)
    nu <- sapply(runif(n, q, 2*q),c,simplify = FALSE)
    #Psi <- replicate(n, crossprod(rMnorm(q)), simplify = FALSE)
    #nu <- sapply(runif(n, q, 2*q), c, simplify = FALSE)
    #sPsi <- cp$single.Psi
    #snu <- cp$single.nu
    #arr.drop <- cp$drop
    # R test
    test.seed <- sample(1e6, 1)
    set.seed(test.seed)
    llR <- array(NA, dim = c(q,q,n))
    for(jj in 1:n) {
      llR[,,jj] <- rwishR(Psi = Psi[[ifelse(sPsi, 1, jj)]],
                          nu = nu[[ifelse(snu, 1, jj)]],
                          inverse = cp$inverse)
    }
    # C++ test
    Psi <- unlistM(Psi, sPsi, arr.drop)
    nu <- unlist(nu)
    if(snu) nu <- nu[1]
    set.seed(test.seed)
    #Psi <- array(unlist(Psi), dim = c(q,q,n))
    #Psi <- Psi[,,single.ind(sPsi),drop=FALSE]
    #nu <- unlist(nu)[single.ind(snu)]
    #if(arr.drop) {
    #  Psi <- dropN(Psi)
    #}
    llcpp <- rwishart(n = n, Psi = Psi, nu = nu,
                      inverse = cp$inverse)
    mx <- abs(llR - llcpp)
    mx <- min(max(mx), max(mx/abs(llR)))
    if(calc.diff) {
      MaxDiff[ii] <- mx
    } else {
      expect_equal(mx, 0, tolerance = tol)
    }
  }
})
