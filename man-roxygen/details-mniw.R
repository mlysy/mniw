#' @details The Matrix-Normal Inverse-Wishart (MNIW) distribution \eqn{(\boldsymbol{X}, \boldsymbol{V}) \sim \textrm{MNIW}(\boldsymbol{\Lambda}, \boldsymbol{\Sigma}, \boldsymbol{\Psi}, \nu)}{(X, V) ~ MNIW(\Lambda, \Sigma, \Psi, \nu)} on random matrices \eqn{\boldsymbol{X}_{p \times q}}{X_(p x q)} and symmetric positive-definite \eqn{\boldsymbol{V}_{q\times q}}{V_(q x q)} is defined as
#' \deqn{
#' \arraycolsep=1.4pt
#' \begin{array}{rcl}
#' \boldsymbol{V} & \sim & \textrm{Inverse-Wishart}(\boldsymbol{\Psi}, \nu) \\
#' \boldsymbol{X} \mid \boldsymbol{V} & \sim & \textrm{Matrix-Normal}(\boldsymbol{\Lambda}, \boldsymbol{\Sigma}, \boldsymbol{V}),
#' \end{array}
#' }{
#' V ~ Inverse-Wishart(\Psi, \nu)
#' }
#' \deqn{
#' \vspace{-1em}
#' }{
#' X | V ~ Matrix-Normal(\Lambda, \Sigma, V),
#' }
#' where the Matrix-Normal distribution is defined as the multivariate normal
#' \deqn{
#' \textrm{vec}(\boldsymbol{X}) \sim \mathcal{N}(\textrm{vec}(\boldsymbol{\Lambda}), \boldsymbol{V} \otimes \boldsymbol{\Sigma}),
#' }{
#' vec(X) ~ N(vec(\Lambda), V \%x\% \Sigma),
#' }
#' where \eqn{\textrm{vec}(\boldsymbol{X})}{vec(X)} is a vector stacking the columns of \eqn{\boldsymbol{X}}{X}, and \eqn{\boldsymbol{V} \otimes \boldsymbol{\Sigma}}{V \%x\% \Sigma} denotes the Kronecker product.
