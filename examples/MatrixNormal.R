# problem dimensions
p <- 4
q <- 2
n <- 10 # number of observations
# parameter values
Mu <- matrix(rnorm(p*q),p,q) # mean matrix
# row-wise variance matrix (positive definite)
RowV <- crossprod(matrix(rnorm(p*p), p, p))
ColV <- rwish(n, Psi = diag(q), nu = q + 1) # column-wise variance (vectorized)

# random sample
X <- rMNorm(n, Mu = Mu, RowV = RowV, ColV = ColV)

# log-density at each sampled value
dMNorm(X, Mu = Mu, RowV = RowV, ColV = ColV, log = TRUE)
