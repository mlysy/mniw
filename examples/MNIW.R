# problem dimensions
p <- 2
q <- 3
n <- 10 # number of samples
# parameter specification
Lambda <- matrix(rnorm(p*q),p,q) # single argument
Sigma <- rwish(n, Psi = diag(p), nu = p + rexp(1)) # vectorized argument
Psi <- rwish(n = 1, Psi = diag(q), nu = q + rexp(1)) # single argument
nu <- q + rexp(1)
# simulate n draws
rMNIW(n, Lambda = Lambda, Sigma = Sigma, Psi = Psi, nu = nu)
