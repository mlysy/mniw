// Posterior inference for mNormRE model

data {
  int<lower=1> q; // number of responses
  int<lower=1> p; // number of covariates
  int<lower=1> n; // number of observations

}
