// temporary exports of C++ functions for debugging purposes

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;
#include "RandomEffects.h"
#include <iostream>

// [[Rcpp::export]]
Eigen::MatrixXd GenerateRandomEffectsOL(int n,
					Eigen::VectorXd x,
					Eigen::MatrixXd C,
					Eigen::MatrixXd OmegaL) {
  int q = x.size();
  MatrixXd Mu(q, n);
  RanEffNormal renorm(q);
  // MatrixXd Uq = MatrixXd::Zero(q,q);
  // MatrixXd Omega(q,q);
  // CrossProdLLt(Omega, OmegaL, Uq);
  // return Omega;
  // std::cout << Omega << std::endl;
  for(int ii=0; ii<n; ii++) {
    renorm.GenerateOL(Mu.col(ii), x, C, OmegaL);
  }
  return Mu;
}
