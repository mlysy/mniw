#' @details The Wishart distribution \eqn{\boldsymbol{X} \sim \textrm{Wishart}(\boldsymbol{\Psi}, \nu)}{X ~ Wishart(\Psi, \nu)} on a symmetric positive-definite random matrix \eqn{\boldsymbol{X}}{X} of size \eqn{q \times q}{q x q} has PDF
#' \deqn{
#' f(\boldsymbol{X} \mid \boldsymbol{\Psi}, \nu) = \frac{|\boldsymbol{X}|^{(\nu-q-1)/2}\exp\big\{-\textrm{tr}(\boldsymbol{\Psi}^{-1}\boldsymbol{X})/2\big\}}{2^{q\nu/2}|\boldsymbol{\Psi}|^{\nu/2} \Gamma_q(\frac \nu 2)},
#' }{
#' f(X | \Psi, \nu) = [ |X|^((\nu-q-1)/2) * e^(-trace(\Psi^{-1} X)/2) ] / [ 2^(q\nu/2) * |\Psi|^(\nu/2) * \Gamma_q(\nu/2) ],
#' }
#' where \eqn{\Gamma_q(\alpha)} is the multivariate gamma function,
#' \deqn{
#' \Gamma_q(\alpha) = \pi^{q(q-1)/4} \prod_{i=1}^q \Gamma(\alpha + (1-i)/2).
#' }{
#' \Gamma_q(\alpha) = \pi^(q(q-1)/4) \prod_{i=1}^q \Gamma(\alpha + (1-i)/2).
#' }
#' The Inverse-Wishart distribution \eqn{\boldsymbol{X} \sim \textrm{Inverse-Wishart}(\boldsymbol{\Psi}, \nu)}{X ~ Inverse-Wishart(\Psi, \nu)} is defined as \eqn{\boldsymbol{X}^{-1} \sim \textrm{Wishart}(\boldsymbol{\Psi}^{-1}, \nu)}{X^{-1} ~ Wishart(\Psi^{-1}, \nu)}.
