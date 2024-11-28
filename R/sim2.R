#' Simulated Group Survival Data: Multiple-source Transfer Learning
#'
#' A dataset containing survival data for 5 groups, where part of the effects
#' and baseline hazards are heterogeneous. Each group consists of 100
#' individuals. The survival time \eqn{T} for group \eqn{i} is generated
#' according to the following model:
#' \deqn{
#'     \lambda^{(i)}(t)=\lambda_{0,i}(t)\exp\left(x^\top\beta^{(i)}\right),
#' }
#' where \eqn{\lambda_{0,i}(t)} represents the baseline hazard for group
#' \eqn{i}, \eqn{\beta^{(i)}} represents the effects for group \eqn{i}, and
#' \eqn{x} represents the covariates. The covariates
#' \eqn{x=(x_{1},\ldots,x_{20})} are generated from a multivariate normal
#' distribution with mean 0 and covariance matrix \eqn{\Sigma=I_{20}}. The
#' baseline hazard function is defined as:
#' \deqn{
#'     \lambda_{0,i}(t)=\left\{\begin{array}{ll}
#'         t^{2}, & \text{if } i=1,3,5 \\
#'         t,     & \text{if } i=2,4.
#'     \end{array}\right.
#' }
#' The effects are defined as:
#' \deqn{
#'     \beta^{(i)}=\left\{\begin{array}{ll}
#'         (0.3,0.3,0.3,0.3,0,\ldots,0),   & \text{if } i=1, \\
#'         (0.9,0.9,0.3,0.3,0,\ldots,0), & \text{if } i=2,4, \\
#'         (-0.3,-0.3,0.3,0.3,0,\ldots,0), & \text{if } i=3,5.
#'     \end{array}\right.
#' }
#' The maximum censoring time is fixed at 2, with an approximate censoring rate
#' of 20%.
#'
#' @name sim2
#' @docType data
#' @format A data frame with 500 rows and 24 variables:
#' \describe{
#' \item{id}{Individual identifier, 1-500}
#' \item{group}{Group indicator, 1-5}
#' \item{time}{Survival time}
#' \item{status}{Status indicator, 0=censored, 1=event}
#' \item{X1}{Covariate 1}
#' \item{X2}{Covariate 2}
#' \item{X3}{Covariate 3}
#' \item{X4}{Covariate 4}
#' \item{X5}{Covariate 5}
#' \item{X6}{Covariate 6}
#' \item{X7}{Covariate 7}
#' \item{X8}{Covariate 8}
#' \item{X9}{Covariate 9}
#' \item{X10}{Covariate 10}
#' \item{X11}{Covariate 11}
#' \item{X12}{Covariate 12}
#' \item{X13}{Covariate 13}
#' \item{X14}{Covariate 14}
#' \item{X15}{Covariate 15}
#' \item{X16}{Covariate 16}
#' \item{X17}{Covariate 17}
#' \item{X18}{Covariate 18}
#' \item{X19}{Covariate 19}
#' \item{X20}{Covariate 20}
#' }
#' @source \url{../articles/simulate-data.html}
NULL
