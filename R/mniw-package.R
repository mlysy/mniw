#' Density and Simulation for the Matrix Normal Inverse Wishart (MNIW) Distribution
#'
#' @details The Matrix-Normal Inverse-Wishart (MNIW) distribution on \code{(X,V)} is given by
#' \deqn{X | V ~ MN(Lambda, Sigma, V), \qquad V ~ iW(Psi, nu).}
#' The default parameters are \code{Lambda = matrix(0,p,q)}, \code{Sigma = diag(p)}, \code{Psi = diag(q)}.
#'
#' @docType package
#' @name mniw
#' @importFrom Rcpp evalCpp
#' @useDynLib mniw
NULL
