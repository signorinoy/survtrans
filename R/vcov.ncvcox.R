#' Variance-covariance matrix for a \code{ncvcox} object.
#' @param object An object of class \code{ncvcox}.
#' @param ... Additional arguments (unused).
#' @return A matrix representing the variance-covariance matrix of the
#' coefficients.
#' @export
vcov.ncvcox <- function(object, ...) {
  # Properties of the coxens object
  time <- object$time
  status <- object$status
  group <- object$group
  x <- object$x
  n_samples <- nrow(x)
  group_levels <- levels(group)
  n_groups <- length(group_levels)
  group_idxs <- lapply(group_levels, function(g) which(group == g))
  coefficients <- object$coefficients * attr(x, "scale")

  # Select the non-zero coefficients and corresponding variables
  coefs <- coefficients[coefficients != 0]
  x <- x[, coefficients != 0]
  n_nonzero <- length(coefs)

  # Calculate the gradient and Hessian for the non-zero coefficients
  lp <- x %*% coefs
  gradients <- matrix(0, nrow = n_samples, ncol = n_nonzero)
  hessians <- matrix(0, nrow = n_samples, ncol = n_nonzero^2)
  for (k in 1:n_groups) {
    idx <- group_idxs[[k]]
    ghs <- calc_grad_hess(lp[idx], x[idx, ], time[idx], status[idx])
    gradients[idx, ] <- ghs$grad
    hessians[idx, ] <- ghs$hess
  }
  hess <- matrix(colSums(hessians), n_nonzero, n_nonzero)
  cov_grad <- stats::cov(gradients) * n_samples
  hess_inv <- solve(hess)

  vcov <- hess_inv %*% cov_grad %*% hess_inv
  return(vcov)
}
