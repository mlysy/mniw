# mniw: Tools for the Matrix-Normal Inverse-Wishart Distribution

*Martin Lysy, Bryan Yates*

---

### Description

Density evaluation and random number generation for the Matrix-Normal Inverse-Wishart (MNIW) distribution, as well as the the Matrix-Normal, Matrix-T, Wishart, and Inverse-Wishart distributions.  Core calculations are implemented in a portable (header-only) C++ library, with matrix manipulations using the [**Eigen**](http://eigen.tuxfamily.org/index.php?title=Main_Page) library for linear algebra.  Also provided is a Gibbs sampler for Bayesian inference with a random-effects model with Matrix-Normal observations.

### Installation

Install the R package [**devtools**](https://CRAN.R-project.org/package=devtools) and run
```r
devtools::install_github("mlysy/mniw")
```

### Usage

Please see function examples, e.g, `?rMNIW`.  A complete description of the distributions provided by the package is available in `vignette("mniw-distributions")`.
