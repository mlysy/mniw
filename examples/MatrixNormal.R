# problem dimensions
p <- 4
q <- 2
n <- 10 # number of observations
# parameter values
Lambda <- matrix(rnorm(p*q),p,q) # mean matrix
# row-wise variance matrix (positive definite)
SigmaR <- crossprod(matrix(rnorm(p*p), p, p))
SigmaC <- rwish(n, Psi = diag(q), nu = q + 1) # column-wise variance (vectorized)

# random sample
X <- rMNorm(n, Lambda = Lambda, SigmaR = SigmaR, SigmaC = SigmaC)

# log-density at each sampled value
dMNorm(X, Lambda = Lambda, SigmaR = SigmaR, SigmaC = SigmaC, log = TRUE)
