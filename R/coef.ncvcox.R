#' Extract the coefficients from a \code{ncvcox} object
#' @param object An object of class \code{ncvcox}.
#' @param ... Additional arguments (not unused).
#' @return A named numeric vector containing the coefficients of the fitted
#' \code{ncvcox} object. Zero coefficients are removed.
#' @export
coef.ncvcox <- function(object, ...) {
  coefficients <- object$coefficients
  names(coefficients) <- colnames(object$x)
  coefficients <- coefficients[coefficients != 0]
  return(coefficients)
}
