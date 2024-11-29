#' @title Log-likelihood for a \code{ncvcox} object
#'
#' @param object An object of class \code{ncvcox}.
#' @param ... Additional arguments (not unused).
#' @return A numeric value representing the log-likelihood of the fitted
#' \code{ncvcox} object.
#' @export
logLik.ncvcox <- function(object, ...) {
  # Properties of the coxens object
  time <- object$time
  status <- object$status
  group <- object$group
  x <- object$x
  n_groups <- length(unique(group))
  group_levels <- levels(group)
  group_idxs <- lapply(group_levels, function(g) which(group == g))
  coefficients <- object$coefficients

  # Calculate the log-likelihood
  offset <- x %*% coefficients
  hazard <- exp(offset)
  risk_set <- stats::ave(hazard, group, FUN = cumsum)
  risk_set <- unlist(lapply(seq_len(n_groups), function(k) {
    idx <- group_idxs[[k]]
    stats::ave(risk_set[idx], time[idx], FUN = max)
  })) # Update the risk set for each group based on unique time points
  return(sum(status * (offset - log(risk_set))))
}
