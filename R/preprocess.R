preprocess <- function(formula, data, group, offset) {
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  time <- y[, 1]
  status <- y[, 2]
  x <- stats::model.matrix(formula, data)[, -1]

  # Properties of the data
  n_samples <- nrow(x)

  # Standardize the covariates
  x <- scale(x)
  x_center <- attr(x, "scaled:center")
  x_scale <- attr(x, "scaled:scale")

  # Check the offset and group arguments
  if (missing(offset)) offset <- rep(0.0, n_samples)
  if (missing(group)) group <- rep(0, n_samples)
  if (!is.factor(group)) group <- factor(group)

  # Sort the data by time
  sorted <- order(time, decreasing = TRUE)
  time <- time[sorted]
  status <- status[sorted]
  x <- x[sorted, , drop = FALSE]
  attr(x, "center") <- x_center
  attr(x, "scale") <- x_scale

  offset <- offset[sorted]
  group <- group[sorted]

  list(
    x = x, time = time, status = status, group = group, offset = offset
  )
}
