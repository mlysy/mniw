///////////////////////////////////////////////////////////////////

// Exported Matrix-T functions

//////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "MatrixT.h"

//[[Rcpp::export]]
Eigen::VectorXd LogDensityMatrixT(Eigen::MatrixXd X,
				  Eigen::MatrixXd Mu,
				  Eigen::MatrixXd RowV,
				  Eigen::MatrixXd ColV,
				  Eigen::VectorXd nu) {
  // dimensions of the problem
  int p = RowV.rows();
  int q = ColV.rows();
  int N = X.cols()/q;
  N = std::max<int>(N, Mu.cols()/q);
  N = std::max<int>(N, RowV.cols()/p);
  N = std::max<int>(N, ColV.cols()/q);
  N = std::max<int>(N, nu.size());
  // output variables
  VectorXd logDens(N);
  // internal variables
  LLT<MatrixXd> cholRowV(p);
  LLT<MatrixXd> cholColV(q);
  double ldRowV, ldColV;
  bool singleX = (X.cols() == q);
  bool singleMu = (Mu.cols() == q);
  bool singleRowV = (RowV.cols() == p);
  bool singleColV = (ColV.cols() == q);
  bool singleNu = (nu.size() == 1);
  MatrixT matT(p,q);
  // precompute Cholesky factors and log-determinants
  if(singleRowV) {
    cholRowV.compute(RowV);
    ldRowV = logDetCholV(cholRowV);
  }
  if(singleColV) {
    cholColV.compute(ColV);
    ldColV = logDetCholV(cholColV);
  }
  // main loop
  for(int ii=0; ii<N; ii++) {
    // compute Chol/log-det if not precomputed
    if(!singleRowV) {
      cholRowV.compute(RowV.block(0,p*ii,p,p));
      ldRowV = logDetCholV(cholRowV);
    }
    if(!singleColV) {
      cholColV.compute(ColV.block(0,q*ii,q,q));
      ldColV = logDetCholV(cholColV);
    }
    // density evaluation
    logDens(ii) = matT.LogDens(X.block(0,q*ii*(!singleX),p,q),
			       Mu.block(0,q*ii*(!singleMu),p,q),
			       RowV.block(0,p*ii*(!singleRowV),p,p),
			       cholRowV, ldRowV,
			       ColV.block(0,q*ii*(!singleColV),q,q),
			       cholColV, ldColV,
			       nu(ii*(!singleNu)));
  }
  return logDens;
}
