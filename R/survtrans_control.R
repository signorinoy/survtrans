#' Ancillary arguments for controlling survtrans fitting
#' @param abstol the absolute tolerance for the proposed algorithm.
#' Default is 1e-4.
#' @param reltol the relative tolerance for the proposed algorithm.
#' Default is 1e-3.
#' @param maxit the maximum number of iterations for the proposed algorithm.
#'  Default is 300.
#' @param verbose a logical value indicating whether to print messages
#'  during the fitting process. Default is \code{FALSE}.
survtrans_control <- function(
    abstol = 1e-4, reltol = 1e-3, maxit = 300, verbose = FALSE) {
  if (abstol <= 0) stop("Invalid absolute tolerance")
  if (reltol <= 0) stop("Invalid relative tolerance")
  if (maxit < 0) stop("Invalid value for iterations")
  if (!is.logical(verbose)) stop("Invalid value for verbose")
  list(
    abstol = abstol, reltol = reltol, maxit = as.integer(maxit),
    verbose = verbose
  )
}
