library(mniw)
source("mniw-testfunctions.R")
context("Matrix-Normal Inverse-Wishart Distribution")

tol <- 1e-6

test_that("MNIW sampling is same in C++ as R", {
  calc.diff <- FALSE
  case.par <- expand.grid(p = c(1,2,4), q = c(1,2,3),
                          Lambda = c("none", "single", "multi"),
                          Sigma = c("none", "single", "multi"),
                          Psi = c("none", "single", "multi"),
                          nu = c("single", "multi"),
                          prec = c(TRUE, FALSE),
                          drop = c(TRUE, FALSE))
  case.par <- case.par[with(case.par, {
    !(Lambda == "none" & ((Sigma == "none") | (Psi == "none")))
  }),]
  ncases <- nrow(case.par)
  rownames(case.par) <- 1:ncases
  n <- 10
  TestSeed <- sample(1e6, ncases)
  if(calc.diff) {
    MaxDiff <- matrix(NA, ncases, 2)
  }
  for(ii in 1:ncases) {
    set.seed(TestSeed[ii])
    cp <- case.par[ii,]
    p <- cp$p
    q <- cp$q
    sLambda <- cp$Lambda != "multi"
    noLambda <- cp$Lambda == "none"
    sSigma <- cp$Sigma != "multi"
    noSigma <- cp$Sigma == "none"
    sPsi <- cp$Psi != "multi"
    noPsi <- cp$Psi == "none"
    snu <- cp$nu != "multi"
    prec <- cp$prec
    arrDrop <- cp$drop
    Lambda <- rMM(n = n, p = p, q = q, noArg = noLambda)
    Sigma <- rMM(n = n, p = p, noArg = noSigma)
    Psi <- rMM(n = n, q = q, noArg = noPsi)
    nu <- sapply(runif(n, 2*q, 3*q),c,simplify = FALSE)
    # R test
    XR <- array(NA, dim = c(p,q,n))
    VR <- array(NA, dim = c(q,q,n))
    set.seed(TestSeed[ii]) # seed
    for(jj in 1:n) {
      XVR <- rmniwR(Lambda = Lambda[[ifelse(sLambda, 1, jj)]],
                    Sigma = Sigma[[ifelse(sSigma, 1, jj)]],
                    Psi = Psi[[ifelse(sPsi, 1, jj)]],
                    nu = nu[[ifelse(snu, 1, jj)]],
                    prec = prec)
      #XVR <- rmniwR(Lambda=Lambda,Sigma=Sigma,Psi=Psi,nu=nu,prec=prec)
      XR[,,jj] <- XVR$X
      VR[,,jj] <- XVR$V
    }
    # C++ test
    Lambda <- unlistM(Lambda, sLambda, arrDrop)
    Sigma <- unlistM(Sigma, sSigma, arrDrop)
    Psi <- unlistM(Psi, sPsi, arrDrop)
    nu <- unlist(nu)
    if(snu) nu <- nu[1]
    arg.list <- list(Lambda = Lambda, Sigma = Sigma, Psi = Psi)
    arg.list <- arg.list[!c(noLambda, noSigma, noPsi)]
    set.seed(TestSeed[ii]) # seed
    XVcpp <- do.call(rMNIW,
                     args = c(arg.list, list(n = n, nu = nu, prec = prec)))
    mx <- c(arDiff(XR, XVcpp$X), arDiff(VR, XVcpp$V))
    if(calc.diff) {
      MaxDiff[ii,] <- mx
    } else {
      expect_equal(mx, c(0,0), tolerance = tol)
    }
  }
})


if(FALSE) {
  range(MaxDiff)

  cp <- case.par[ii,]
  p <- cp$p
  q <- cp$q
  Lambda <- if(cp$Lambda == "none") matrix(0,p,q) else rMnorm(p,q)
  Sigma <- if(cp$Sigma == "none") diag(p) else crossprod(rMnorm(p))
  Psi <- if(cp$Psi == "none") crossprod(rMnorm(q)) else diag(q)
  nu <- runif(1, q, 2*q)
  prec <- cp$prec

  n <- 1
  set.seed(test.seed)
  XVR2 <- rmniwR(Lambda = Lambda,
                 Sigma = Sigma,
                 Psi = Psi,
                 nu = nu,
                 prec = prec)
  set.seed(test.seed)
  XVcpp2 <- rMNIW(n = 1,
                  Lambda = Lambda,
                  Sigma = Sigma,
                  Psi = Psi,
                  nu = nu,
                  prec = prec)
  tmp <- expect_equal(XVR2, XVcpp2, tolerance = tol)
  tmp
}
