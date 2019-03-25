#--- compare speed to mvnfast ----------------------------------------------

require(mvnfast)
require(mniw)

n <- 1e4
d <- 500

mu <- rnorm(d)
V <- rWishart(n = 1, df = d+5, Sigma = diag(d))[,,1]

system.time({
  X1 <- rmNorm(n = n, mu = mu, V = V)
})
system.time({
  X2 <- rmvn(n = n, mu = mu, sigma = V)
})
