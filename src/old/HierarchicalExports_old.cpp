// #include <Rcpp.h>
// using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
//#include <iostream>
#include "mniwWishart.h"
#include "mniwMatNorm.h"
#include "mniwRandomEffects.h"
#include "TriUtils.h"
// #include "RandomEffects.h"

// Gibbs sampler for the posterior distribution of the parameters of the "uneqV"
// hierarchical model described by Everson and Morris (2000):
//
// y_i | mu_i ~ind N(mu_i, V_i)
// mu_i | beta, Sigma ~ind N(x_i' beta, Sigma)
// beta, Sigma ~ MNIW(Lambda, Omega^{-1}, Psi, nu)
//
// the model parameters are: beta, Sigma, Mu = (mu_1, ..., mu_n)
// the data is Y = (y_1, ..., y_n), X = (x_1, ..., x_n), V = (V_1, ..., V_n)
// the hyperparameters are Lambda, Omega, Psi, nu.
// setting all(Omega == 0) = TRUE is equivalent to a Lebesgue prior on beta.
//[[Rcpp::export]]
List HierUneqVModelGibbs_old(int nSamples, int nBurn,
			     Eigen::MatrixXd Y, Eigen::MatrixXd X,
			     Eigen::MatrixXd V,
			     Eigen::MatrixXd Lambda, Eigen::MatrixXd Omega,
			     Eigen::MatrixXd Psi, double nu,
			     Eigen::MatrixXd Beta0, Eigen::MatrixXd iSigma0,
			     Eigen::MatrixXd Mu0,
			     bool updateBetaSigma, bool updateMu,
			     bool storeBetaSigma, bool storeMu) {
  // dimensions of the problem
  int N = Y.rows();
  int p = X.cols();
  int q = Y.cols();
  int ii, jj;
  // // don't store if don't update
  // if(!updateBetaSigma) {
  //   storeBetaSigma = false;
  // }
  // if(!updateMu) {
  //   storeMu = false;
  // }
  // whether we have the flat prior on Beta
  bool flatBeta = Omega.squaredNorm() == 0.0;
  // output variables
  MatrixXd BetaOut(p, (storeBetaSigma ? nSamples*q : q));
  MatrixXd SigmaOut(q, (storeBetaSigma ? nSamples*q: q));
  MatrixXd MuOut(q, (storeMu ? nSamples*N : N));
  if(!storeBetaSigma) {
    BetaOut.setZero();
    SigmaOut.setZero();
  }
  if(!storeMu) {
    MuOut.setZero();
  }
  // internal variables
  TempPQ *tmp = new TempPQ(p,q);
  double nuHat = nu + N - (flatBeta ? p : 0);
  MatrixXd Yt = Y.adjoint();
  // initialize Gibbs sampler
  MatrixXd Mu = Mu0.adjoint();
  MatrixXd Beta = Beta0;
  MatrixXd OmegaL = iSigma0.llt().matrixL();
  // MatrixXd Mu = Yt;
  // MatrixXd Beta = MatrixXd::Zero(p,q);
  // MatrixXd OmegaL = MatrixXd::Zero(q,q);
  // internal values
  MatrixXd iXtXO = (X.adjoint() * X + Omega).eval().llt().solve(MatrixXd::Identity(p,p));
  MatrixXd RowiOmegaHatL = iXtXO.llt().matrixL();
  MatrixXd IH = MatrixXd::Identity(N,N) -
    (X * (X.adjoint() * X).llt().solve(X.adjoint())).eval();
  MatrixXd LtO = Lambda.adjoint() * Omega;
  MatrixXd LtOL = LtO * Lambda;
  MatrixXd PsiHat(q,q);
  MatrixXd XiHatL(q,q);
  MatrixXd MuXLtO(q,p);
  MatrixXd IHMut(N,q);
  MatrixXd LambdaHat(p,q);
  MatrixXd x(q,N);
  // inverses of the variances
  MatrixXd iV(q,N*q);
  for(jj=0; jj<N; jj++) {
    iV.block(0,jj*q,q,q) = V.block(0,jj*q,q,q).llt().solve(tmp->Iq);
  }
  for(ii=-nBurn; ii<nSamples; ii++) {
    if(updateBetaSigma) {
      // update Beta, Sigma
      // new MNIW parameters
      MuXLtO.noalias() = Mu * X + LtO;
      LambdaHat.noalias() = iXtXO * MuXLtO.adjoint();
      if(!flatBeta) {
	PsiHat.noalias() = Mu * Mu.adjoint();
	PsiHat.noalias() -= MuXLtO * LambdaHat;
	PsiHat += LtOL + Psi;
      }
      else {
	IHMut.noalias() = IH * Mu.adjoint();
	PsiHat.noalias() = Mu * IHMut;
	PsiHat += Psi;
      }
      // std::cout << "PsiHat = \n" << PsiHat << "\nnuHat = \n" << nuHat << std::endl;
      // reverse cholesky decomposition PsiHat = XiHatU * XiHatL
      ReverseCholesky(XiHatL, PsiHat, tmp->lltq);
      // generate Beta and OmegaL, where Sigma^{-1} = OmegaL * OmegaL'
      GenerateWishartLowerTriXi(OmegaL, XiHatL, nuHat);
      GenerateMatrixNormalRowSColO(Beta, LambdaHat, RowiOmegaHatL, OmegaL, tmp->Xpq);
    }
    if(updateMu) {
      // update Mu
      x = (X * Beta).adjoint();
      x -= Yt;
      for(jj=0; jj<N; jj++) {
	GenerateRandomEffectsOL(Mu.col(jj), x.col(jj), iV.block(0,jj*q,q,q),
				OmegaL, tmp);
      }
      Mu += Yt;
    }
    // storage
    if(ii >= 0) {
      if(storeBetaSigma) {
	BetaOut.block(0,ii*q,p,q) = Beta;
	InverseLLt(SigmaOut.block(0,ii*q,q,q), OmegaL, tmp->Lq, tmp->Uq, tmp->Iq);
      }
      if(storeMu) {
	MuOut.block(0,ii*N,q,N) = Mu;
      }
    }
  }
  delete tmp;
  return List::create(_["Beta"] = BetaOut,
		      _["Sigma"] = SigmaOut,
		      _["Mu"] = MuOut);
}
