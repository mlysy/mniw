#--- random testing -------------------------------------------------------------

require(mniw)

q <- 3
N <- 5

Sigma <- replicate(N,
                   crossprod(rMNorm(n = 1, Mu = matrix(0,q,q))))
Lambda <- rMNorm(n = 1, Mu = matrix(0,q,N))
#Lambda[,-1] <- Lambda[,1]

set.seed(1)
Z <- matrix(rnorm(q*N),q,N)
X <- matrix(NA, q,N)
for(ii in 1:N) {
  X[,ii] <- t(chol(Sigma[,,ii])) %*% Z[,ii] + Lambda[,ii]
}
X <- t(X)
set.seed(1)
X2 <- rmNorm(n = N, mu = t(Lambda), V = Sigma)
range(X-X2)

