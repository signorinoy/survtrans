#' Extract coefficients from a \code{cv.ncvcox} object
#' @param object An object of class \code{cv.ncvcox}.
#' @param s A character string specifying the rule to select the coefficients.
#' The default is "lambda.min". Other options are "lambda.1se".
#' @param ... Unused.
#' @return A named numeric vector containing the coefficients of the fitted
#' \code{cv.ncvcox} object.
#' @export
coef.cv.ncvcox <- function(object, s = "lambda.min", ...) {
  if (!s %in% c("lambda.min", "lambda.1se")) {
    stop("s must be one of 'lambda.min' or 'lambda.1se'")
  }
  lambda_opt <- switch(s,
    "lambda.min" = object$lambda.min,
    "lambda.1se" = object$lambda.1se
  )
  coefficients <- object$coefficients[which(object$lambdas == lambda_opt), ]
  names(coefficients) <- colnames(object$coefficients)
  return(coefficients)
}
