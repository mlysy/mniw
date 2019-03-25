# graphical tests for the hierarchical model

require(mniw)

#--- check that new and old version give same outputs --------------------------

case.par <- expand.grid(p = 1:3, q = 1:3, flat.prior = c(TRUE, FALSE),
                        updateHyp = c(TRUE, FALSE), storeHyp = c(TRUE, FALSE),
                        updateRE = c(TRUE, FALSE), storeRE = c(TRUE, FALSE))
ncases <- nrow(case.par)

seed <- 123
test_eq <- rep(NA, ncases)
for(ii in 1:ncases) {
  case <- case.par[ii,]
  # problem dimensions
  nsamples <- sample(8:12, 1)
  N <- sample(15:25, 1)
  p <- case$p
  q <- case$q
  # prior specification
  Lambda0 <- rMNorm(1, Mu = matrix(0, p, q))
  Omega0 <- crossprod(rMNorm(1, Mu = matrix(0, p, p)))
  Psi0 <- crossprod(rMNorm(1, Mu = matrix(0, q, q)))
  nu0 <- rexp(1) + (q+1)
  prior_list <- list(Psi = Psi0, Lambda = Lambda0,
                     Omega = solve(Omega0), nu = nu0)
  if(case$flat.prior) prior_list$Omega <- matrix(0, p, p)
  # simulate data
  X <- rMNorm(1, matrix(0, N, p))
  V <- array(rwish(N, diag(q), q+1),dim=c(q,q,N))
  Sigma00 = rwish(1, Psi0, nu0)
  Beta00 = rMNorm(1, Lambda0, Omega0, Sigma00)
  Mu00 = rmNorm(N, X %*% Beta00, Sigma00) # random effects
  Y00 = rmNorm(N, Mu00, V) # Data
  # mcmc initialization
  init_list <- list(Beta = Beta00, Sigma = Sigma00, Mu = Mu00)
  # Gibbs sampling to get posterior
  set.seed(seed)
  r_fit <- mNormRE.post(nsamples, Y=Y00, V=V, X=X,
                        prior = prior_list,
                        init = init_list,
                        burn = ceiling(nsamples/2),
                        updateHyp = case$updateHyp, updateRE = case$updateRE,
                        storeHyp = case$storeHyp, storeRE = case$storeRE)
  set.seed(seed)
  r_fit2 <- mNormRE.post2(nsamples, Y=Y00, V=V, X=X,
                          prior = prior_list,
                          init = init_list,
                          burn = ceiling(nsamples/2),
                          updateHyp = case$updateHyp, updateRE = case$updateRE,
                          storeHyp = case$storeHyp, storeRE = case$storeRE)
  test_eq[ii] <- identical(r_fit, r_fit2)
}

all(test_eq)
