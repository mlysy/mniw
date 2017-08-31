#--- basic testing ---------------------------------------------------------

require(Rcpp)
require(RcppEigen)
sourceCpp("CrossProdVXY.cpp")

# cross-products
rMN <- function(n, p = n) {
  matrix(rnorm(n*p), n, p)
}

p <- 3
q <- 5
r <- 2

X <- rMN(p, q)
V <- crossprod(rMN(p))
Y <- rMN(p,r)

ipR <- t(X) %*% solve(V) %*% Y
ipC <- crossprodV(X = X, Y = Y, V = V, inverse = TRUE)[,,1]
range(ipR - ipC)

#--- random testing --------------------------------------------------------

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

#--- rmNormRE ------------------------------------------------------------------


require(mniw)

nSamples = 10
q = 3
y = rnorm(q)
V = diag(q) #rwish(1, diag(q), q+1)
rmNormRE(n=nSamples, y=y, V=V)
