#' Log-likelihood for a \code{coxtrans} object
#'
#' @param object An object of class \code{coxstream}.
#' @param ... Additional arguments (not unused).
#' @return A numeric value representing the log-likelihood of the fitted
#' \code{coxtrans} object.
#' @export
logLik.coxtrans <- function(object, ...) {
  # Properties of the coxtrans object
  time <- object$time
  status <- object$status
  group <- object$group
  x <- object$x
  n_groups <- length(unique(group))
  group_levels <- levels(group)
  group_idxs <- lapply(group_levels, function(x) which(group == x))
  coefficients <- object$coefficients
  coefficients <- sweep(coefficients, 1, attr(x, "scale"), "*")

  beta <- coefficients[, 1:n_groups] + coefficients[, (n_groups + 1)]
  # Calculate the log-likelihood
  offset <- numeric(nrow(x))
  for (k in seq_len(n_groups)) {
    idx <- group_idxs[[k]]
    offset[idx] <- x[idx, ] %*% beta[, k]
  }
  hazard <- exp(offset)
  risk_set <- stats::ave(hazard, group, FUN = cumsum)
  for (k in seq_len(n_groups)) {
    idx <- group_idxs[[k]]
    risk_set[idx] <- stats::ave(risk_set[idx], time[idx], FUN = max)
  } # Update the risk set for each group based on unique time points
  return(sum(status * (offset - log(risk_set))))
}
