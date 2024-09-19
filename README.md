<!-- comment -->

# mniw: The Matrix-Normal Inverse-Wishart Distribution

*Martin Lysy, Bryan Yates*

<!-- badges: start -->

[![R-CMD-check](https://github.com/mlysy/mniw/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mlysy/mniw/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->


---

### Description

Density evaluation and random number generation for the Matrix-Normal Inverse-Wishart (MNIW) distribution, as well as the the Matrix-Normal, Matrix-T, Wishart, and Inverse-Wishart distributions.  Core calculations are implemented in a portable (header-only) C++ library, with matrix manipulations using the [**Eigen**](http://eigen.tuxfamily.org/index.php?title=Main_Page) library for linear algebra.  Also provided is a Gibbs sampler for Bayesian inference on a random-effects model with Matrix-Normal observations.

### Installation

To install the CRAN version (1.0.1):
``` r
install.packages("mniw", INSTALL_opts = "--install-tests")
```

To install the latest development version: first install the [**devtools**](https://CRAN.R-project.org/package=devtools), then:

```r
devtools::install_github("mlysy/mniw", INSTALL_opts = "--install-tests")
```

### Usage

The primary advantage of the **mniw** package is that it "vectorizes" over its input arguments.  Take for example the simulation of a Wishart distribution, which can be done with the built-in R function `stats::rWishart()`: 

```r
n <- 10
p <- 3
nu <- 6
# produces an array of size p x p x n
Psi <- stats::rWishart(n = n, df = nu, Sigma = diag(p))
```

Now suppose we want to generate Wishart random variables each with a different `Sigma`:

```r
# Vectorizing over the 'Sigma' argument
X <- apply(Psi, 3, stats::rWishart, n = 1, df = nu)
X <- array(X, dim = c(p, p, n))
```

However, the code above is both slow for large `n`, and inconvenient due to the reshaping of the `apply()` output.  The equivalent code using **mniw** is:

```r
X <- rwish(n, df = nu, Psi = Psi) # produces an array of size p x p x n
```

It is both simpler, and much faster for large `n` and `p`.

The other functions in **mniw** behave much the same way.  A complete description of the distributions provided by the package is available [here](https://mlysy.github.io/mniw/articles/mniw-distributions.html).
