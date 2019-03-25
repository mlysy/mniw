# test the GenerateOL method

require(mniw)

ans <- replicate(n = 100, {
  q <- sample(1:20,1) # number of dimensions
  x <- rnorm(q) # observation vector
  V <- rwish(1, diag(q), q+1) # observation variance
  A <- rwish(1, diag(q), q+1) # prior variance
  C <- mniw:::.solveV(V)
  Omega <- mniw:::.solveV(A)
  OmegaL <- t(chol(Omega))
  mniw:::GenerateRandomEffectsOL(n = 1, x = x, C = C, OmegaL = OmegaL)
})

# check the matrix product

p <- 4
V <- crossprod(matrix(rnorm(p^2), p))
A <- crossprod(matrix(rnorm(p^2), p))
C <- solve(V)
Q <- solve(A)

B1 <- solve(C+Q, Q)
B2 <- V %*% solve(V+A)
