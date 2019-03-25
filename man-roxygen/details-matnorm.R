#' @details The Matrix-Normal distribution \eqn{\boldsymbol{X} \sim \textrm{Matrix-Normal}(\boldsymbol{\Lambda}, \boldsymbol{\Sigma}_R, \boldsymbol{\Sigma}_C)}{X ~ Matrix-Normal(\Lambda, \Sigma_R, \Sigma_C)} on the random matrix \eqn{\boldsymbol{X}_{p \times q}}{X_(p x q)} is defined as
#' \deqn{
#' \textrm{vec}(\boldsymbol{X}) \sim \mathcal{N}(\textrm{vec}(\boldsymbol{\Lambda}), \boldsymbol{\Sigma}_C \otimes \boldsymbol{\Sigma}_R),
#' }{
#' vec(X) ~ N(vec(\Lambda), \Sigma_C \%x\% \Sigma_R),
#' }
#' where \eqn{\textrm{vec}(\boldsymbol{X})}{vec(X)} is a vector stacking the columns of \eqn{\boldsymbol{X}}{X}, and \eqn{\boldsymbol{\Sigma}_C \otimes \boldsymbol{\Sigma}_R}{\Sigma_C \%x\% \Sigma_R} denotes the Kronecker product.
