# Input Specifications
nsamples <- 10
N <- 20
p <- 1
q <- 2
# prior specification
Lambda0 <- rMNorm(1, Lambda = matrix(0, p, q))
Omega0 <- crossprod(rMNorm(1, Lambda = matrix(0, p, p)))
Psi0 <- crossprod(rMNorm(1, Lambda = matrix(0, q, q)))
nu0 <- rexp(1) + (q+1)
# covariate specification
X <- rMNorm(1, matrix(0, N, p))
V <- array(rwish(N, diag(q), q+1),dim=c(q,q,N))

# generate response data
Sigma00 <- rwish(1, Psi0, nu0)
Beta00 <- rMNorm(1, Lambda0, Omega0, Sigma00)
Mu00 <- rmNorm(N, X %*% Beta00, Sigma00) # random effects
Y00 <- rmNorm(N, Mu00, V) # Data

# Gibbs sampling to get posterior
prior_list <- list(Psi=Psi0, Lambda=Lambda0, Omega=solve(Omega0), nu=nu0) # Prior
init_list <- list(Beta=Beta00,Sigma=Sigma00,Mu=Mu00) # Initialize MCMC
r_fit <- RxNormLM(nsamples, Y=Y00, V=V, X=X,
                  prior = prior_list,
                  init = init_list,
                  burn = ceiling(nsamples/2))
