#' Prediction method for \code{coxtrans} objects.
#' @param object An object of class \code{coxtrans}.
#' @param newdata Optional new data for making predictions. If omitted,
#'   predictions are made using the data used for fitting the model.
#' @param newgroup Optional new group for making predictions. If omitted,
#'   predictions are made using the groups from the original data.
#' @param type The type of prediction to perform. Options include:
#'   \describe{
#'     \item{\code{"lp"}}{The linear predictor.}
#'     \item{\code{"terms"}}{The components of the linear predictor.}
#'     \item{\code{"risk"}}{The risk score \eqn{\exp(\text{lp})}.}
#'     \item{\code{"expected"}}{The expected number of events, given the
#'              covariates and follow-up time.}
#'     \item{\code{"survival"}}{The survival probability, given the covariates
#'             and follow-up time.}
#'   }
#' @param ... Additional arguments (not unused).
#' @return A numeric vector of predictions.
#' @export
predict.ncvcox <- function(
    object, newdata = NULL, newgroup = NULL,
    type = c("lp", "terms", "risk", "expected", "survival"), ...) {
  type <- match.arg(type)
  x <- stats::model.matrix(object$formula, newdata)[, -1]

  # Properties of the coxtrans object
  coefficients <- object$coefficients * attr(x, "scale")

  lp <- x %*% coefficients

  if (type == "lp") {
    return(lp)
  } else if (type == "risk") {
    risk <- exp(lp)
    return(risk)
  } else {
    stop("type must be one of 'lp' or 'risk'")
  }
}
