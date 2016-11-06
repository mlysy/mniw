#' @title Matrix Inner Product
#' @description Vectorized matrix inner products \code{Y = X' V X} or \code{Y = X' V^{-1} X}.
#'
#' @param X Either a matrix of size \code{p x q}, or an array of size \code{p x q x n}.
#' @param V Either a matrix of size \code{p x q}, or an array of size \code{p x q x n}.
#' @param inverse Whether or not the inner product should be calculated with \code{V} (\code{inverse = FALSE}) or \code{V^{-1}}.
#' @return An array of size \code{q x q x n}.
#' @export
crossprodV <- function(X, V, inverse = FALSE) {
  # dimensions of X and V
  if(is.vector(X)) X <- as.matrix(X)
  p <- dim(X)[1]
  q <- dim(X)[2]
  X <- matrix(X,p)
  if(!all(dim(V)[1:2] == p)) stop("X and V have incompatible dimensions.")
  V <- matrix(V,p)
  # check lengths
  n <- unique(sort(c(ncol(X)/q, ncol(V)/p)))
  n <- c(1, n[n>1])
  if(length(n) > 2) stop("X and V have different lengths.")
  n <- max(n)
  # finally call C code
  Y <- .Call('mniw_CrossProdV', PACKAGE = 'mniw', X, V, p, q, inverse)
  array(Y, dim = c(q,q,n))
}
