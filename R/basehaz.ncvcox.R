#' Predict the cumulative baseline hazard function for \code{ncvcox} objects
#'
#' @param object An object of class \code{ncvcox}.
#' @param newdata A numeric vector of time points at which to predict the
#' baseline hazard function. If \code{NULL}, the function will predict the
#' baseline hazard function at the unique event times in the fitted data.
#' @param ... Additional arguments (not Additional arguments (not unused).).
#'
#' @return A \code{data.frame} with one row for each time point, and columns
#' containing the event time, the cumulative baseline hazard function, and the
#' strata.
#' @export
basehaz.ncvcox <- function(object, newdata, ...) {
  # Properties of the ncvcox object
  time <- object$time
  status <- object$status
  group <- object$group
  x <- object$x
  n_groups <- length(unique(group))
  group_levels <- levels(group)
  group_idxs <- lapply(group_levels, function(x) which(group == x))
  coefficients <- object$coefficients * attr(x, "scale")

  offset <- x %*% coefficients
  hazard <- exp(offset)
  risk_set <- stats::ave(hazard, group, FUN = cumsum)
  for (k in seq_len(n_groups)) {
    idx <- group_idxs[[k]]
    risk_set[idx] <- stats::ave(risk_set[idx], time[idx], FUN = max)
  }

  basehaz_df <- data.frame()
  for (k in seq_len(n_groups)) {
    idx <- group_idxs[[k]]
    time_rev <- rev(time[idx])
    status_rev <- rev(status[idx])
    risk_set_rev <- rev(risk_set[idx])
    basehaz <- cumsum(status_rev / risk_set_rev)
    basehaz_df <- rbind(basehaz_df, data.frame(
      time = time_rev[status_rev == 1],
      basehaz = basehaz[status_rev == 1],
      strata = group_levels[k]
    ))
  }
  if (n_groups == 1) {
    basehaz_df$strata <- NULL
  }
  return(basehaz_df)
}
