#' Transfer Learning for Cox Model with Global and Local Shrinkage
#'
#' @param formula A formula object, with the response on the left of a \code{~}
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function.
#' @param data A data frame containing the variables in the model.
#' @param group A factor variable indicating the group of each observation.
#' @param lambda1 A non-negative value specifying the sparse penalty parameter.
#' The default is 0.
#' @param lambda2 A non-negative value specifying the global bias penalty
#' parameter. The default is 0.
#' @param lambda3 A non-negative value specifying the local bias penalty
#' parameter. The default is 0.
#' @param penalty A character string specifying the penalty function. The
#' default is "lasso". Other options are "MCP" and "SCAD".
#' @param gamma A non-negative value specifying the penalty parameter. The
#' default is 3.7 for SCAD and 3.0 for MCP.
#' @param rho A value larger than 1 specifying the increase/decrease factor
#' for the augmented Lagrangian's penalty parameter. The default is 2.0.
#' @param tau A value larger than 1 specifying the tolerance for the
#' trade-off between the primal and dual residuals. The default is 10.0.
#' @param init A numeric vector of initial values for the coefficients. The
#' default is a zero vector.
#' @param control An object of class \link{survtrans_control} containing control
#' parameters for the fitting algorithm. Default is
#' \code{survtrans_control(...)}.
#' @param ... Additional arguments to be passed to the fitting algorithm.
#'
#' @return An object of class \code{coxtrans}.
#' @export
#'
#' @examples
#' formula <- survival::Surv(time, status) ~ . - group - id
#' fit <- coxtrans(
#'   formula, sim2, sim2$group,
#'   lambda1 = 0.03, lambda2 = 0.01, lambda3 = 0.01, penalty = "SCAD"
#' )
#' summary(fit)
coxtrans <- function(
    formula, data, group, lambda1 = 0, lambda2 = 0, lambda3 = 0,
    penalty = c("lasso", "MCP", "SCAD"),
    gamma = switch(penalty,
      SCAD = 3.7,
      MCP = 3,
      1
    ), rho = 2.0, tau = 10.0, init,
    control, ...) {
  # Load the data
  data <- preprocess(formula, data, group = group)
  x <- data$x
  x_scale <- attr(x, "scale")
  time <- data$time
  status <- data$status
  group <- data$group

  # Properties of the data
  n_samples_total <- nrow(x)
  n_features <- ncol(x)
  n_groups <- length(unique(group))
  group_levels <- levels(group)
  group_idxs <- lapply(group_levels, function(x) which(group == x))
  n_samples_group <- sapply(group_idxs, length)

  risk_set_size <- stats::ave(rep(1, n_samples_total), group, FUN = cumsum)
  risk_set_size <- unlist(lapply(1:n_groups, function(k) {
    idx <- group_idxs[[k]]
    stats::ave(risk_set_size[idx], time[idx], FUN = max)
  }))
  null_deviance <- -sum(status * log(risk_set_size))

  ## Check the penalty argument
  penalty <- match.arg(penalty, choices = c("lasso", "MCP", "SCAD"))

  ## Check the init argument
  if (!missing(init) && length(init) > 0) {
    if (length(init) != n_features * (n_groups + 1)) {
      stop("Wrong length for inital values")
    }
  } else {
    init <- rep(0, n_features * (n_groups + 1))
  }
  ## Check the control argument
  if (missing(control)) control <- survtrans_control(...)

  # Extract the coefficients from init vector
  init <- sweep(matrix(init, nrow = n_features), 1, x_scale, `*`)
  theta <- matrix(
    data = init[, 1:(n_groups + 1), drop = FALSE],
    nrow = n_features * (n_groups + 1), ncol = 1
  )

  # Initialize the training process
  n_iterations <- 0
  message <- ""
  convergence <- FALSE

  offset <- numeric(n_samples_total)
  w <- numeric(n_samples_total)
  z <- numeric(n_samples_total)

  idx <- which(lower.tri(matrix(1, n_groups, n_groups)), arr.ind = TRUE)
  e <- matrix(0, nrow = nrow(idx), ncol = n_groups)
  e[cbind(seq_len(nrow(idx)), idx[, 1])] <- 1
  e[cbind(seq_len(nrow(idx)), idx[, 2])] <- -1

  contr_sum <- Matrix::Matrix(
    cbind(
      kronecker(
        matrix(n_samples_group / n_samples_total, nrow = 1),
        diag(n_features)
      ),
      0 * diag(n_features)
    ),
    sparse = TRUE
  )
  contr_penalty <- rbind(
    cbind(
      Matrix::Diagonal(n_features * n_groups),
      kronecker(matrix(1, n_groups, 1), diag(n_features))
    ),
    cbind(
      Matrix::Diagonal(n_features * n_groups),
      Matrix::Matrix(0, n_features * n_groups, n_features, sparse = TRUE)
    ),
    cbind(
      kronecker(e, diag(n_features)),
      Matrix::Matrix(
        0, n_groups * (n_groups - 1) * n_features / 2, n_features,
        sparse = TRUE
      )
    )
  )

  n_constraints_sum <- nrow(contr_sum)
  n_constraints_penalty <- nrow(contr_penalty)
  n_constraints <- n_constraints_sum + n_constraints_penalty
  n_parameters <- n_features * (n_groups + 1)

  contr_sum2 <- Matrix::crossprod(contr_sum)
  contr_penalty2 <- Matrix::crossprod(contr_penalty)
  contr2 <- contr_sum2 + contr_penalty2

  sparse_idx <- 1:(n_groups * n_features)
  global_idx <- (n_groups * n_features + 1):(2 * n_groups * n_features)
  local_idx <- (2 * n_groups * n_features + 1):n_constraints_penalty

  x_tilde <- cbind(
    Matrix::bdiag(lapply(group_idxs, function(idx) x[idx, ])),
    do.call(rbind, lapply(group_idxs, function(idx) x[idx, ]))
  )
  time_tilde <- unlist(lapply(group_idxs, function(idx) time[idx]))
  status_tilde <- unlist(lapply(group_idxs, function(idx) status[idx]))

  alpha <- Matrix::Matrix(0, n_constraints_penalty, 1, sparse = TRUE)
  mu <- Matrix::Matrix(0, n_features, 1, sparse = TRUE)
  nu <- Matrix::Matrix(0, n_constraints_penalty, 1, sparse = TRUE)
  vartheta <- 1

  history <- c()

  repeat {
    n_iterations <- n_iterations + 1

    # Calculate the weights and residuals
    offset <- x_tilde %*% theta
    n_passes <- 0
    for (k in 1:n_groups) {
      idx <- n_passes + seq_len(length(group_idxs[[k]]))
      n_passes <- n_passes + length(group_idxs[[k]])
      wls <- approx_likelihood(
        offset = offset[idx], time = time_tilde[idx], status = status_tilde[idx]
      )
      w[idx] <- wls$weights
      z[idx] <- wls$residuals + offset[idx]
    }

    # Update the coefficients
    xwx <- Matrix::crossprod(x_tilde, w * x_tilde) / n_samples_total
    xwz <- Matrix::crossprod(x_tilde, w * z) / n_samples_total
    lhs <- xwx + vartheta * contr2
    rhs <- xwz - Matrix::crossprod(contr_sum, mu) +
      vartheta * Matrix::crossprod(contr_penalty, alpha - nu / vartheta)
    theta <- Matrix::solve(lhs, rhs, sparse = TRUE, tol = 1e-6)

    # Update the auxiliary variables
    alpha_old <- alpha
    alpha <- contr_penalty %*% theta + nu / vartheta

    alpha[sparse_idx] <- threshold_prox(
      alpha[sparse_idx], vartheta, penalty, lambda1, gamma
    )
    alpha[global_idx] <- threshold_prox(
      alpha[global_idx], vartheta, penalty, lambda2, gamma
    )
    alpha[local_idx] <- threshold_prox(
      alpha[local_idx], vartheta, penalty, lambda3, gamma
    )
    mu <- mu + vartheta * contr_sum %*% theta
    nu <- nu + vartheta * (contr_penalty %*% theta - alpha)

    r_norm <- norm(
      as.matrix(
        rbind(contr_sum %*% theta, contr_penalty %*% theta - alpha)
      ), "2"
    )
    s_norm <- norm(
      as.matrix(Matrix::crossprod(contr_penalty, alpha - alpha_old)), "2"
    ) * vartheta

    eps_pri <- sqrt(n_constraints) * control$abstol +
      control$reltol * pmax(
        norm(
          as.matrix(rbind(contr_sum %*% theta, contr_penalty %*% theta)), "2"
        ),
        norm(as.matrix(alpha), "2")
      )
    eps_dual <- sqrt(n_parameters) * control$abstol +
      control$reltol * norm(
        as.matrix(
          Matrix::crossprod(contr_sum, mu) +
            Matrix::crossprod(contr_penalty, nu)
        ),
        "2"
      )

    # Check the convergence
    if (n_iterations >= control$maxit) {
      convergence <- TRUE
      message <- paste0(
        "Maximum number of iterations reached (", control$maxit, ")."
      )
    }
    if (r_norm < eps_pri && s_norm < eps_dual) {
      convergence <- TRUE
      message <- paste0(
        "Convergence reached at iteration ", n_iterations, "."
      )
    }

    offset <- x_tilde %*% theta
    hazard <- exp(offset)
    risk_set <- numeric(n_samples_total)
    n_passes <- 0
    for (k in 1:n_groups) {
      idx <- n_passes + seq_len(length(group_idxs[[k]]))
      n_passes <- n_passes + length(group_idxs[[k]])
      risk_set[idx] <- cumsum(hazard[idx])
      risk_set[idx] <- stats::ave(risk_set[idx], time_tilde[idx], FUN = max)
    }
    loss <- sum(status_tilde * (offset - log(risk_set)))
    if (loss / null_deviance < 0.01) {
      convergence <- TRUE
      message <- paste0(
        "The log-likelihood is too small (", loss / null_deviance,
        "). Stopping the algorithm."
      )
    }

    alpha_ <- contr_penalty %*% theta
    loss_penalty <- penalty(alpha_[sparse_idx], penalty, lambda1, gamma) +
      penalty(alpha_[global_idx], penalty, lambda2, gamma) +
      penalty(alpha_[local_idx], penalty, lambda3, gamma)
    loss_penalty <- loss_penalty * n_samples_total
    loss_total <- loss - loss_penalty

    if (control$verbose) {
      cat(
        "========================================\n",
        sprintf("Iteration Number       : %d", n_iterations), "\n",
        sprintf("Residuals (pri, dual)  : %.4f, %.4f", r_norm, s_norm), "\n",
        sprintf("Epsilon (pri, dual)    : %.4f, %.4f", eps_pri, eps_dual), "\n",
        sprintf("Augmented Parameter    : %.4f", vartheta), "\n",
        sprintf("Total Loss             : %.4f", loss_total), "\n",
        sprintf(" - Log Likelihood      : %.4f", loss), "\n",
        sprintf(" - Penalty             : %.4f", loss_penalty), "\n",
        "========================================\n"
      )
    }
    history <- rbind(
      history,
      c(
        n_iterations, r_norm, s_norm, eps_pri, eps_dual, vartheta,
        loss, loss_penalty, loss_total
      )
    )

    if (convergence) break

    # Update the penalty parameter
    if (r_norm > tau * s_norm) vartheta <- vartheta * rho
    if (s_norm > tau * r_norm) vartheta <- vartheta / rho
    vartheta <- min(max(vartheta, 1e-3), 2.118034)
  }

  theta_ <- theta
  eps_theta <- Matrix::norm(
    Matrix::crossprod(contr_penalty, alpha - alpha_old), "I"
  )
  eps_alpha <- Matrix::norm(alpha - alpha_old, "I")
  alpha_local <- matrix(alpha[local_idx, 1], nrow = n_features)
  theta <- matrix(theta, nrow = n_features, ncol = n_groups + 1)
  eta <- matrix(0, nrow = n_features, ncol = n_groups)
  eta_idx <- matrix(0, nrow = n_features, ncol = n_groups)
  for (i in 1:n_features) {
    is_processed <- rep(FALSE, n_groups)
    for (j in 1:n_groups) {
      if (is_processed[j]) next
      is_processed[j] <- TRUE
      eta_idx[i, j] <- j
      for (k in 1:n_groups) {
        if (is_processed[k]) next
        pos <- pair_index(j, k, n_groups)
        if (abs(alpha_local[i, pos]) < eps_alpha) {
          eta_idx[i, k] <- j
          is_processed[k] <- TRUE
        }
      }
    }
    for (j in unique(eta_idx[i, ])) {
      idx <- which(eta_idx[i, ] == j)
      eta[i, idx] <- mean(theta[i, idx])
    }
  }

  # Handling the extreme small values to zero
  eta[abs(eta) < eps_theta] <- 0
  beta <- theta[, n_groups + 1]
  beta[abs(beta) < eps_theta] <- 0

  # Forcing beta satisfying sum(eta) = 0
  for (j in 1:n_features) {
    idx <- which(eta[j, ] != 0)
    if (length(idx) == 0) next
    eta_mean <- mean(eta[j, idx])
    eta[j, idx] <- eta[j, idx] - eta_mean
    beta[j] <- beta[j] + eta_mean
  }

  # Unscale the coefficients
  coefficients <- sweep(cbind(eta, beta), 1, x_scale, `/`)
  colnames(coefficients) <- c(group_levels, "Center")
  rownames(coefficients) <- colnames(x)

  colnames(history) <- c(
    "Iteration", "Primal Residual", "Dual Residual", "Primal Epsilon",
    "Dual Epsilon", "Augmented Parameter", "Negative Log Likelihood",
    "Penalty Term", "Objective"
  )
  history <- as.data.frame(history)

  # Return the fitted model
  fit <- list(
    coefficients = coefficients,
    logLik = loss,
    iter = n_iterations,
    message = message,
    history = history,
    penalty = penalty,
    lambda1 = lambda1,
    lambda2 = lambda2,
    lambda3 = lambda3,
    gamma = gamma,
    formula = formula,
    call = match.call(),
    time = time,
    status = status,
    group = group,
    x = x,
    theta = theta_,
    alpha = alpha
  )
  class(fit) <- "coxtrans"
  return(fit)
}
