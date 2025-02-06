#' Variance-covariance matrix for a \code{coxtrans} object.
#' @param object An object of class \code{coxstream}.
#' @param ... Additional arguments (unused).
#' @return A matrix representing the variance-covariance matrix of the
#' coefficients.
#' @export
vcov.coxtrans <- function(object, ...) {
  # Properties of the coxtrans object
  time <- object$time
  status <- object$status
  group <- object$group
  x <- object$x
  n_samples_total <- nrow(x)
  n_features <- ncol(x)
  group_levels <- levels(group)
  n_groups <- length(group_levels)
  group_idxs <- lapply(group_levels, function(g) which(group == g))
  n_samples_group <- sapply(group_idxs, length)
  coefficients <- object$coefficients
  coefficients <- sweep(coefficients, 1, attr(x, "scale"), "*")

  # Extract the coefficients from the object
  eta <- coefficients[, 1:n_groups]
  beta <- coefficients[, (n_groups + 1)]

  # Reassign the coefficients' group and track expanded coefficients
  eta_expanded <- numeric()
  coefs_processed <- numeric()
  eta_idx <- matrix(0, n_features, n_groups)
  n_total_groups <- 0
  for (j in seq_len(n_features)) {
    feature_groups <- as.factor(eta[j, ])
    feature_levels <- unique(as.character(feature_groups))
    for (k in seq_along(feature_levels)) {
      if (feature_levels[k] == "0") next
      idx <- feature_groups == feature_levels[k]
      eta_idx[j, idx] <- k + n_total_groups
    }
    eta_expanded <- c(eta_expanded, as.numeric(feature_levels))
    coefs_processed <- c(coefs_processed, as.numeric(feature_levels) + beta[j])
    n_total_groups <- n_total_groups + length(feature_levels)
  }
  coefs_expanded <- c(eta_expanded, beta)

  # Check if there are non-zero coefficients
  n_expanded_nonzero <- sum(coefs_expanded != 0)
  n_processed_nonzero <- sum(coefs_processed != 0)
  if (n_expanded_nonzero == 0 || n_processed_nonzero == 0) {
    stop("No non-zero coefficients to compute the variance-covariance matrix")
  }

  # Group the samples according to the estimated coefficients' group
  x <- cbind(
    do.call(rbind, lapply(seq_len(n_groups), function(k) {
      # Construct the feature map for each group
      feature_map <- matrix(0, nrow = n_features, ncol = n_total_groups)
      for (j in seq_len(n_features)) {
        feature_map[j, eta_idx[j, k]] <- 1
      }
      x[group_idxs[[k]], ] %*% feature_map
    })),
    do.call(rbind, lapply(group_idxs, function(idx) x[idx, ]))
  )
  time <- unlist(lapply(group_idxs, function(idx) time[idx]))
  status <- unlist(lapply(group_idxs, function(idx) status[idx]))

  # Select the non-zero coefficients and corresponding variables
  coefs <- coefs_expanded[coefs_expanded != 0]
  x <- x[, coefs_expanded != 0]
  from <- sort(unique(c(eta_idx)))
  to <- seq_len(length(from)) - 1
  eta_idx <- apply(eta_idx, c(1, 2), function(e) to[match(e, from)])

  # Calculate the gradient and Hessian for the non-zero coefficients
  lp <- x %*% coefs
  gradients <- matrix(0, nrow = n_samples_total, ncol = n_expanded_nonzero)
  hessians <- matrix(0, nrow = n_samples_total, ncol = n_expanded_nonzero^2)
  n_passes <- 0
  for (k in seq_len(n_groups)) {
    idx <- n_passes + seq_len(length(group_idxs[[k]]))
    n_passes <- n_passes + length(idx)
    ghs <- calc_grad_hess(lp[idx], x[idx, ], time[idx], status[idx])
    gradients[idx, ] <- ghs$grad
    hessians[idx, ] <- ghs$hess
  }
  hess <- matrix(colSums(hessians), n_expanded_nonzero, n_expanded_nonzero)
  cov_grad <- stats::cov(gradients) * n_samples_total

  # Construct the Null space of the constraints
  n_constr <- sum(rowSums(eta_idx != 0) > 0)
  if (n_constr > 0) {
    constr_idx <- which(rowSums(eta_idx != 0) > 0)
    constr <- matrix(0, nrow = n_constr, ncol = sum(coefs_expanded != 0))
    for (i in seq_len(n_constr)) {
      idx <- constr_idx[i]
      group_levels <- unique(eta_idx[idx, ])
      group_levels <- group_levels[group_levels != 0]
      constr[i, group_levels] <- sapply(
        group_levels,
        function(level) sum((eta_idx[idx, ] == level) * n_samples_group)
      ) / n_samples_total
    }
    null_constr <- MASS::Null(t(constr))
    hess_inv <- null_constr %*%
      solve(t(null_constr) %*% hess %*% null_constr) %*% t(null_constr)
    vcov_constr <- hess_inv %*% cov_grad %*% hess_inv
  } else {
    hess_inv <- solve(hess)
    vcov_constr <- hess_inv %*% cov_grad %*% hess_inv
  }

  # Construct the variance-covariance matrix for the original coefficients
  prox_origional <- matrix(0, n_expanded_nonzero, n_processed_nonzero)
  beta_idx <- cumsum(beta != 0) + sum(eta_expanded != 0)
  beta_idx[beta == 0] <- 0
  j <- 1
  for (i in seq_len(n_features)) {
    idx1 <- unique(eta_idx[i, ])
    idx2 <- beta_idx[i]
    if (any(idx1 != 0) || idx2 != 0) {
      for (idx in idx1) {
        prox_origional[idx, j] <- 1
        prox_origional[idx2, j] <- 1
        j <- j + 1
      }
    }
  }
  vcov_processed <- t(prox_origional) %*% vcov_constr %*% prox_origional

  return(vcov_processed)
}
