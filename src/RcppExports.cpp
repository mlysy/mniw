// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// CrossProdVXX
Eigen::MatrixXd CrossProdVXX(Eigen::MatrixXd X, Eigen::MatrixXd V, int p, int q, bool inverse);
RcppExport SEXP _mniw_CrossProdVXX(SEXP XSEXP, SEXP VSEXP, SEXP pSEXP, SEXP qSEXP, SEXP inverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< bool >::type inverse(inverseSEXP);
    rcpp_result_gen = Rcpp::wrap(CrossProdVXX(X, V, p, q, inverse));
    return rcpp_result_gen;
END_RCPP
}
// CrossProdVXY
Eigen::MatrixXd CrossProdVXY(Eigen::MatrixXd X, Eigen::MatrixXd Y, Eigen::MatrixXd V, int p, int q, int r, bool inverse);
RcppExport SEXP _mniw_CrossProdVXY(SEXP XSEXP, SEXP YSEXP, SEXP VSEXP, SEXP pSEXP, SEXP qSEXP, SEXP rSEXP, SEXP inverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< bool >::type inverse(inverseSEXP);
    rcpp_result_gen = Rcpp::wrap(CrossProdVXY(X, Y, V, p, q, r, inverse));
    return rcpp_result_gen;
END_RCPP
}
// LogDensityMultivariateNormal
Eigen::VectorXd LogDensityMultivariateNormal(Eigen::MatrixXd X, Eigen::MatrixXd Mu, Eigen::MatrixXd V);
RcppExport SEXP _mniw_LogDensityMultivariateNormal(SEXP XSEXP, SEXP MuSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Mu(MuSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(LogDensityMultivariateNormal(X, Mu, V));
    return rcpp_result_gen;
END_RCPP
}
// LogDensityMatrixNormal
Eigen::VectorXd LogDensityMatrixNormal(Eigen::MatrixXd X, Eigen::MatrixXd Mu, Eigen::MatrixXd RowV, Eigen::MatrixXd ColV);
RcppExport SEXP _mniw_LogDensityMatrixNormal(SEXP XSEXP, SEXP MuSEXP, SEXP RowVSEXP, SEXP ColVSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Mu(MuSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type RowV(RowVSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ColV(ColVSEXP);
    rcpp_result_gen = Rcpp::wrap(LogDensityMatrixNormal(X, Mu, RowV, ColV));
    return rcpp_result_gen;
END_RCPP
}
// GenerateMultivariateNormal
Eigen::MatrixXd GenerateMultivariateNormal(int N, Eigen::MatrixXd Lambda, Eigen::MatrixXd Sigma);
RcppExport SEXP _mniw_GenerateMultivariateNormal(SEXP NSEXP, SEXP LambdaSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(GenerateMultivariateNormal(N, Lambda, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// GenerateMatrixNormal
Eigen::MatrixXd GenerateMatrixNormal(int N, Eigen::MatrixXd Lambda, Eigen::MatrixXd RowSigma, Eigen::MatrixXd ColSigma);
RcppExport SEXP _mniw_GenerateMatrixNormal(SEXP NSEXP, SEXP LambdaSEXP, SEXP RowSigmaSEXP, SEXP ColSigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type RowSigma(RowSigmaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ColSigma(ColSigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(GenerateMatrixNormal(N, Lambda, RowSigma, ColSigma));
    return rcpp_result_gen;
END_RCPP
}
// GenerateMatrixNIW
List GenerateMatrixNIW(int N, Eigen::MatrixXd Lambda, Eigen::MatrixXd Sigma, Eigen::MatrixXd Psi, Eigen::VectorXd nu, bool inverse);
RcppExport SEXP _mniw_GenerateMatrixNIW(SEXP NSEXP, SEXP LambdaSEXP, SEXP SigmaSEXP, SEXP PsiSEXP, SEXP nuSEXP, SEXP inverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type inverse(inverseSEXP);
    rcpp_result_gen = Rcpp::wrap(GenerateMatrixNIW(N, Lambda, Sigma, Psi, nu, inverse));
    return rcpp_result_gen;
END_RCPP
}
// GenerateRandomEffectsNormal
Eigen::MatrixXd GenerateRandomEffectsNormal(int N, Eigen::MatrixXd lambda, Eigen::MatrixXd y, Eigen::MatrixXd V, Eigen::MatrixXd A);
RcppExport SEXP _mniw_GenerateRandomEffectsNormal(SEXP NSEXP, SEXP lambdaSEXP, SEXP ySEXP, SEXP VSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type V(VSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(GenerateRandomEffectsNormal(N, lambda, y, V, A));
    return rcpp_result_gen;
END_RCPP
}
// HierUneqVModelGibbs
List HierUneqVModelGibbs(int nSamples, int nBurn, Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd V, Eigen::MatrixXd Lambda, Eigen::MatrixXd Omega, Eigen::MatrixXd Psi, double nu, Eigen::MatrixXd Beta0, Eigen::MatrixXd iSigma0, Eigen::MatrixXd Mu0, bool updateBetaSigma, bool updateMu, bool storeBetaSigma, bool storeMu);
RcppExport SEXP _mniw_HierUneqVModelGibbs(SEXP nSamplesSEXP, SEXP nBurnSEXP, SEXP YSEXP, SEXP XSEXP, SEXP VSEXP, SEXP LambdaSEXP, SEXP OmegaSEXP, SEXP PsiSEXP, SEXP nuSEXP, SEXP Beta0SEXP, SEXP iSigma0SEXP, SEXP Mu0SEXP, SEXP updateBetaSigmaSEXP, SEXP updateMuSEXP, SEXP storeBetaSigmaSEXP, SEXP storeMuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nSamples(nSamplesSEXP);
    Rcpp::traits::input_parameter< int >::type nBurn(nBurnSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type V(VSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Beta0(Beta0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type iSigma0(iSigma0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Mu0(Mu0SEXP);
    Rcpp::traits::input_parameter< bool >::type updateBetaSigma(updateBetaSigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type updateMu(updateMuSEXP);
    Rcpp::traits::input_parameter< bool >::type storeBetaSigma(storeBetaSigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type storeMu(storeMuSEXP);
    rcpp_result_gen = Rcpp::wrap(HierUneqVModelGibbs(nSamples, nBurn, Y, X, V, Lambda, Omega, Psi, nu, Beta0, iSigma0, Mu0, updateBetaSigma, updateMu, storeBetaSigma, storeMu));
    return rcpp_result_gen;
END_RCPP
}
// LogDensityWishart
Eigen::VectorXd LogDensityWishart(Eigen::MatrixXd X, Eigen::MatrixXd Psi, Eigen::VectorXd nu, bool inverse);
RcppExport SEXP _mniw_LogDensityWishart(SEXP XSEXP, SEXP PsiSEXP, SEXP nuSEXP, SEXP inverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type inverse(inverseSEXP);
    rcpp_result_gen = Rcpp::wrap(LogDensityWishart(X, Psi, nu, inverse));
    return rcpp_result_gen;
END_RCPP
}
// GenerateWishart
Eigen::MatrixXd GenerateWishart(int N, Eigen::MatrixXd Psi, Eigen::VectorXd nu, bool inverse);
RcppExport SEXP _mniw_GenerateWishart(SEXP NSEXP, SEXP PsiSEXP, SEXP nuSEXP, SEXP inverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type inverse(inverseSEXP);
    rcpp_result_gen = Rcpp::wrap(GenerateWishart(N, Psi, nu, inverse));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mniw_CrossProdVXX", (DL_FUNC) &_mniw_CrossProdVXX, 5},
    {"_mniw_CrossProdVXY", (DL_FUNC) &_mniw_CrossProdVXY, 7},
    {"_mniw_LogDensityMultivariateNormal", (DL_FUNC) &_mniw_LogDensityMultivariateNormal, 3},
    {"_mniw_LogDensityMatrixNormal", (DL_FUNC) &_mniw_LogDensityMatrixNormal, 4},
    {"_mniw_GenerateMultivariateNormal", (DL_FUNC) &_mniw_GenerateMultivariateNormal, 3},
    {"_mniw_GenerateMatrixNormal", (DL_FUNC) &_mniw_GenerateMatrixNormal, 4},
    {"_mniw_GenerateMatrixNIW", (DL_FUNC) &_mniw_GenerateMatrixNIW, 6},
    {"_mniw_GenerateRandomEffectsNormal", (DL_FUNC) &_mniw_GenerateRandomEffectsNormal, 5},
    {"_mniw_HierUneqVModelGibbs", (DL_FUNC) &_mniw_HierUneqVModelGibbs, 16},
    {"_mniw_LogDensityWishart", (DL_FUNC) &_mniw_LogDensityWishart, 4},
    {"_mniw_GenerateWishart", (DL_FUNC) &_mniw_GenerateWishart, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_mniw(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
