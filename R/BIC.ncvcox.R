#' BIC for a \code{ncvcox} object
#'
#' @param object An object of class \code{ncvcox}.
#' @param type A character string specifying the type of BIC to compute.
#' "traditional" corresponds to Cn=1, and "modified" corresponds to
#' Cn=log(log(d)).
#' @param ... Additional arguments (not unused).
#' @return A numeric value representing the BIC of the fitted \code{ncvcox}
#' object.
#' @export
BIC.ncvcox <- function(object, type = c("traditional", "modified"), ...) {
  type <- match.arg(type)

  # Properties of the ncvcox object
  coefficients <- object$coefficients
  n_samples <- nrow(object$x)
  n_features <- ncol(object$x)

  # The number of parameters should minus the number of active constraints
  n_parameters <- sum(coefficients != 0)

  # Log-likelihood of the model
  loglik <- logLik(object)

  # The parameter of the BIC
  c_n <- ifelse(type == "traditional", 1, log(log(n_features)))

  return(-2 * loglik + c_n * n_parameters * log(n_samples))
}
