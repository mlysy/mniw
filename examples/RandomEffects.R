# data specification
q <- 5
y <- rnorm(q)
V <- rwish(1, diag(q), q+1)
# prior specification
lambda <- rep(0,q)
A <- diag(q)
n <- 10
# random sampling
rmNormRE(n, y, V, lambda, A)
