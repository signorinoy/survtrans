#' Transfer Learning for Cox Model with Global and Local Shrinkage
#'
#' @param formula A formula object, with the response on the left of a \code{~}
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function.
#' @param data A data frame containing the variables in the model.
#' @param group A factor variable indicating the group of each observation.
#' @param target A character string specifying the target group.
#' @param lambda1 A non-negative value specifying the sparse penalty parameter.
#' The default is 0.
#' @param lambda2 A non-negative value specifying the biased penalty
#' parameter. The default is 0.
#' @param alpha A value between 0 and 1 specifying the degree of the local and
#' global shrinkage. The default is 1, which means no global shrinkage.
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
#'   formula, sim2, sim2$group, 1,
#'   lambda1 = 0.075, lambda2 = 0.08, alpha = 0.5, penalty = "SCAD"
#' )
#' summary(fit)
coxtrans <- function(
    formula, data, group, target,
    lambda1 = 0.0, lambda2 = 0.0, alpha = 1.0,
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

  # Keep the target group at the first position
  target_level <- as.character(target)
  group_levels <- c(target_level, group_levels[group_levels != target_level])

  group_idxs <- lapply(group_levels, function(x) which(group == x))
  n_samples_group <- sapply(group_idxs, length)

  risk_set_size <- stats::ave(rep(1, n_samples_total), group, FUN = cumsum)
  risk_set_size <- unlist(lapply(1:n_groups, function(k) {
    idx <- group_idxs[[k]]
    stats::ave(risk_set_size[idx], time[idx], FUN = max)
  }))
  null_deviance <- -sum(status * log(risk_set_size))

  ## Check the input arguments
  if (lambda1 < 0 || lambda2 < 0) stop("Lambda parameters must be non-negative")
  if (alpha < 0 || alpha > 1) stop("Alpha must be in [0,1]")
  penalty <- match.arg(penalty, choices = c("lasso", "MCP", "SCAD"))
  if (!missing(init) && length(init) > 0) {
    if (length(init) != n_features * (n_groups + 1)) {
      stop("Wrong length for inital values")
    }
  } else {
    init <- rep(0, n_features * (n_groups + 1))
  }
  if (missing(control)) control <- survtrans_control(...)

  # Extract the coefficients from init vector
  init <- sweep(matrix(init, nrow = n_features), 1, x_scale, `*`)
  theta <- matrix(
    data = init[, 1:(n_groups + 1), drop = FALSE],
    nrow = n_features * (n_groups + 1), ncol = 1
  )

  # Construct pre-computed matrices
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
  contr_penalty <- Matrix::sparseMatrix(
    i = c(
      seq_len(n_features), seq_len(n_features),
      n_features + seq_len(n_features * n_groups),
      n_features + seq_len(n_features * (n_groups - 1))
    ),
    j = c(
      seq_len(n_features), (n_features * n_groups) + seq_len(n_features),
      rep(1:n_features, times = n_groups),
      n_features + seq_len(n_features * (n_groups - 1))
    ),
    x = c(
      rep(1, n_features * (n_groups + 2)), rep(-1, n_features * (n_groups - 1))
    ),
    dims = c(n_features * (n_groups + 1), n_features * (n_groups + 1))
  )
  n_constraints_sum <- nrow(contr_sum)
  n_constraints_penalty <- nrow(contr_penalty)
  n_constraints <- n_constraints_sum + n_constraints_penalty
  n_parameters <- n_features * (n_groups + 1)
  contr_sum2 <- Matrix::crossprod(contr_sum)
  contr_penalty2 <- Matrix::crossprod(contr_penalty)
  contr2 <- contr_sum2 + contr_penalty2

  sparse_idx <- seq_len(n_features)
  local_idx <- n_features + seq_len(n_features * (n_groups - 1))
  global_idx <- n_features * n_groups + seq_len(n_features)

  x_tilde <- cbind(
    Matrix::bdiag(lapply(group_idxs, function(idx) x[idx, ])),
    do.call(rbind, lapply(group_idxs, function(idx) x[idx, ]))
  )
  time_tilde <- unlist(lapply(group_idxs, function(idx) time[idx]))
  status_tilde <- unlist(lapply(group_idxs, function(idx) status[idx]))

  # Initialize the training process
  n_iterations <- 0
  message <- ""
  convergence <- FALSE
  history <- matrix(NA_real_, nrow = control$maxit, ncol = 9)
  colnames(history) <- c(
    "Iteration", "Primal.Residual", "Dual.Residual", "Primal.Epsilon",
    "Dual.Epsilon", "Augmented.Parameter", "Negative.Log.Likelihood",
    "Penalty.Term", "Objective"
  )

  offset <- numeric(n_samples_total)
  w <- numeric(n_samples_total)
  z <- numeric(n_samples_total)

  eta <- Matrix::Matrix(0, n_constraints_penalty, 1, sparse = TRUE)
  mu <- Matrix::Matrix(0, n_features, 1, sparse = TRUE)
  nu <- Matrix::Matrix(0, n_constraints_penalty, 1, sparse = TRUE)
  vartheta <- 1

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
    lhs_chol <- Matrix::Cholesky(lhs, perm = TRUE)
    rhs <- xwz - Matrix::crossprod(contr_sum, mu) +
      vartheta * Matrix::crossprod(contr_penalty, eta - nu / vartheta)
    theta <- Matrix::solve(lhs_chol, rhs, system = "A")

    # Update the auxiliary variables
    eta_old <- eta
    eta <- contr_penalty %*% theta + nu / vartheta

    eta[sparse_idx] <- threshold_prox(
      eta[sparse_idx], vartheta, penalty, lambda1, gamma
    )
    eta[global_idx] <- threshold_prox(
      eta[global_idx], vartheta, penalty, lambda2 * (1 - alpha), gamma
    )
    eta[local_idx] <- threshold_prox(
      eta[local_idx], vartheta, penalty, lambda2 * alpha, gamma
    )
    mu <- mu + vartheta * contr_sum %*% theta
    nu <- nu + vartheta * (contr_penalty %*% theta - eta)

    r_norm <- norm(
      as.matrix(
        rbind(contr_sum %*% theta, contr_penalty %*% theta - eta)
      ), "2"
    )
    s_norm <- norm(
      as.matrix(Matrix::crossprod(contr_penalty, eta - eta_old)), "2"
    ) * vartheta

    eps_pri <- sqrt(n_constraints) * control$abstol +
      control$reltol * pmax(
        norm(
          as.matrix(rbind(contr_sum %*% theta, contr_penalty %*% theta)), "2"
        ),
        norm(as.matrix(eta), "2")
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

    eta_ <- contr_penalty %*% theta
    loss_penalty <- penalty(eta_[sparse_idx], penalty, lambda1, gamma) +
      penalty(eta_[global_idx], penalty, lambda2 * (1 - alpha), gamma) +
      penalty(eta_[local_idx], penalty, lambda2 * alpha, gamma)
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
    history[n_iterations, ] <- c(
      n_iterations, r_norm, s_norm, eps_pri, eps_dual, vartheta,
      loss, loss_penalty, loss_total
    )

    if (convergence) break

    # Update the penalty parameter
    if (r_norm > tau * s_norm) vartheta <- vartheta * rho
    if (s_norm > tau * r_norm) vartheta <- vartheta / rho
    vartheta <- min(max(vartheta, 1e-3), 2.118034)
  }

  eps <- norm(as.matrix(Matrix::crossprod(contr_penalty, eta - eta_old)), "I")
  eps_local <- norm(as.matrix(eta[local_idx] - eta_old[local_idx]), "I")
  flag_local <- matrix(abs(eta[local_idx]) <= eps_local, nrow = n_features)
  flag_local <- cbind(rep(TRUE, n_features), flag_local)

  theta <- matrix(theta, nrow = n_features, ncol = n_groups + 1)
  delta <- theta[, seq_len(n_groups)]
  w <- theta[, (n_groups + 1)]
  beta <- matrix(NA, nrow = n_features, ncol = n_groups)

  for (i in seq_len(n_features)) {
    delta_global <- mean(delta[i, ])
    idx <- which(flag_local[i, ])
    delta_local <- mean(delta[i, idx])
    # Biased Penalty & Constraints
    delta[i, idx] <- ifelse(
      abs(delta_local - delta_global) < eps, 0, delta_local - delta_global
    )
    # Sparse Penalty
    beta[i, ] <- ifelse(
      abs(delta[i, ] + w[i]) < eps, 0, delta[i, ] + w[i]
    )
  }
  w <- rowMeans(beta)

  coefficients <- cbind(beta - w, w)
  coefficients[abs(coefficients) < eps] <- 0
  colnames(coefficients) <- c(group_levels, "Center")
  rownames(coefficients) <- colnames(x)

  colnames(history) <- c(
    "Iteration", "Primal.Residual", "Dual.Residual", "Primal.Epsilon",
    "Dual.Epsilon", "Augmented.Parameter", "Log.Likelihood",
    "Penalty.Term", "Total.Loss"
  )
  history <- history[seq_len(n_iterations), ]

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
    alpha = alpha,
    gamma = gamma,
    formula = formula,
    call = match.call(),
    time = time,
    status = status,
    group = group,
    x = x,
    control = control
  )
  class(fit) <- "coxtrans"
  return(fit)
}

#' Diagnose Cox Transfer Model's Optimization Process
#'
#' @param object An object of class \code{coxstream}.
#' @param ... Additional arguments (currently unused).
#' @details This function produces two plots:
#' - Residuals Convergence: Plots the evolution of primal and dual residuals
#' along with their tolerance levels.
#' - Loss Decomposition: Plots the negative log-likelihood, total loss, and
#' penalty term.
#' @export
diagnose.coxtrans <- function(object, ...) {
  history <- as.list(as.data.frame(object$history))
  iter_range <- range(history$Iteration)
  colors <- c(
    Augmented = "#1B9E77", Primal = "#D95F02", Dual = "#7570B3",
    NLL = "#E7298A", Penalty = "#66A61E", Total = "#E6AB02"
  )
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(
    mfrow = c(1, 2), mar = c(5, 4.5, 4, 4) + 0.1, mgp = c(2.5, 0.8, 0),
    cex.axis = 0.9, cex.lab = 1.1, cex.main = 1.2
  )

  plot_residuals <- function() {
    y_lim <- range(
      0, history$Primal.Residual, history$Dual.Residual,
      history$Primal.Epsilon, history$Dual.Epsilon
    ) * 1.05
    plot(NA,
      xlim = iter_range, ylim = y_lim, xlab = "Iteration", ylab = "Residuals",
      main = "Residuals Convergence"
    )
    grid(col = "gray90", lty = 2)
    lines(
      history$Iteration, history$Primal.Residual,
      col = colors["Primal"], lwd = 2
    )
    lines(
      history$Iteration, history$Dual.Residual,
      col = colors["Dual"], lwd = 2
    )
    lines(
      history$Iteration, history$Primal.Epsilon,
      col = colors["Primal"], lty = 2, lwd = 1.5
    )
    lines(
      history$Iteration, history$Dual.Epsilon,
      col = colors["Dual"], lty = 2, lwd = 1.5
    )
    par(new = TRUE)
    plot(
      history$Iteration, history$Augmented.Parameter,
      type = "l", col = colors["Augmented"], lwd = 2, axes = FALSE,
      xlab = "", ylab = "",
      ylim = range(history$Augmented.Parameter) + c(-0.05, 0.05)
    )
    axis(4, col.axis = colors["Augmented"], col = colors["Augmented"])
    mtext("Augmented Parameter",
      side = 4, line = 2.5, col = colors["Augmented"], cex = 0.9
    )
    legend(
      "top",
      legend = c("Primal", "Dual", "Residual", "Tolerance", expression(rho)),
      col = c(
        colors["Primal"], colors["Dual"], "black", "black",
        colors["Augmented"]
      ),
      lty = c(1, 1, 1, 2, 1),
      lwd = c(2, 2, 1.5, 1.5, 2),
      horiz = FALSE,
      ncol = 3,
      xpd = NA,
      inset = c(0, 0.01),
      cex = 0.75
    )
  }

  plot_objective <- function() {
    y_lim_obj <- range(
      -history$Log.Likelihood, -history$Total.Loss
    ) + c(-0.05, 0.05)
    plot(NA,
      xlim = iter_range, ylim = y_lim_obj, xlab = "Iteration", ylab = "Loss",
      main = "Loss Decomposition"
    )
    grid(col = "gray90", lty = 2)
    lines(
      history$Iteration, -history$Log.Likelihood,
      col = colors["NLL"], lwd = 2
    )
    lines(
      history$Iteration, -history$Total.Loss,
      col = colors["Total"], lwd = 2
    )
    par(new = TRUE)
    plot(
      history$Iteration, history$Penalty.Term,
      type = "l", col = colors["Penalty"], lwd = 2, axes = FALSE,
      xlab = "", ylab = "", ylim = range(history$Penalty.Term) + c(-0.05, 0.05)
    )
    axis(4, col.axis = colors["Penalty"], col = colors["Penalty"])
    mtext("Penalty Term",
      side = 4, line = 2.5, col = colors["Penalty"], cex = 0.9
    )
    legend(
      "top",
      legend = c("Total", "NLL", "Penalty"),
      col = colors[c("Total", "NLL", "Penalty")],
      lty = 1,
      lwd = 2,
      horiz = FALSE,
      ncol = 3,
      xpd = NA,
      inset = c(0, 0.01),
      cex = 0.75
    )
  }

  plot_residuals()
  plot_objective()
}

#' Extract the coefficients from a \code{coxtrans} object
#' @param object An object of class \code{coxstream}.
#' @param ... Additional arguments (not unused).
#' @return A named numeric vector containing the coefficients of the fitted
#' \code{coxtrans} object. The names indicate the group(s) to which the
#' coefficients belong. Zero coefficients are removed.
#' @export
coef.coxtrans <- function(object, ...) {
  # Properties of the coxtrans object
  coefficients <- object$coefficients
  n_features <- nrow(coefficients)
  n_groups <- ncol(coefficients) - 1

  eta <- coefficients[, 1:n_groups]
  beta <- coefficients[, (n_groups + 1)]

  group_names <- colnames(eta)
  variable_names <- colnames(object$x)

  # Generate the coefficients and names for each feature
  coef <- numeric()
  coef_names <- character()
  for (j in seq_len(n_features)) {
    feature_groups <- as.factor(eta[j, ])
    feature_levels <- unique(as.character(feature_groups))
    coef_postfix <- sapply(feature_levels, function(level) {
      idx <- feature_groups == level
      ifelse(
        sum(idx) == n_groups,
        "ALL", paste0(group_names[idx], collapse = ", ")
      )
    })
    coef_names <- c(
      coef_names, paste0(variable_names[j], " (", coef_postfix, ")")
    )
    coef <- c(coef, as.numeric(feature_levels) + beta[j])
  }
  names(coef) <- coef_names
  coef <- coef[abs(coef) > object$control$abstol]
  coef
}

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
  n_processed_nonzero <- sum(abs(coefs_processed) > object$control$abstol)
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
        if (
          (idx != 0) && (abs(coefs[idx] + coefs[idx2]) < object$control$abstol)
        ) {
          next
        }
        prox_origional[idx, j] <- 1
        prox_origional[idx2, j] <- 1
        j <- j + 1
      }
    }
  }
  t(prox_origional) %*% vcov_constr %*% prox_origional
}


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

  beta <- coefficients[, 1:n_groups] + coefficients[, (n_groups + 1)]
  offset <- numeric(nrow(x))
  for (k in seq_len(n_groups)) {
    idx <- group_idxs[[k]]
    offset[idx] <- x[idx, ] %*% beta[, k]
  }
  hazard <- exp(offset)
  risk_set <- stats::ave(hazard, group, FUN = cumsum)

  # Update the risk set for each group based on unique time points
  for (k in seq_len(n_groups)) {
    idx <- group_idxs[[k]]
    risk_set[idx] <- stats::ave(risk_set[idx], time[idx], FUN = max)
  }
  sum(status * (offset - log(risk_set)))
}

#' BIC for a \code{coxtrans} object
#'
#' @param object An object of class \code{coxstream}.
#' @param type A character string specifying the type of BIC to compute.
#' "traditional" corresponds to Cn=1, and "modified" corresponds to
#' Cn=log(log(d)).
#' @param ... Additional arguments (not unused).
#' @return A numeric value representing the BIC of the fitted \code{coxtrans}
#' object.
#' @export
BIC.coxtrans <- function(object, type = c("traditional", "modified"), ...) {
  type <- match.arg(type)

  # Properties of the coxtrans object
  coefficients <- object$coefficients
  n_samples <- nrow(object$x)
  n_features <- nrow(coefficients)
  n_groups <- ncol(coefficients) - 1

  # Number of active constraints
  eta_group <- matrix(0, n_features, n_groups)
  for (j in seq_len(n_features)) {
    feature_values <- coefficients[j, 1:n_groups]
    feature_levels <- unique(feature_values)
    for (k in seq_along(feature_levels)) {
      eta_group[j, feature_values == feature_levels[k]] <- k
    }
  }
  n_active_constraints <- sum(apply(eta_group, 1, function(row) {
    length(unique(row)) > 1
  }))

  # The number of parameters should minus the number of active constraints
  n_parameters <- sum(apply(coefficients, 1, function(coef_row) {
    length(unique(coef_row[coef_row != 0]))
  })) - n_active_constraints

  # Log-likelihood of the model
  loglik <- logLik(object)

  # The parameter of the BIC
  c_n <- ifelse(type == "traditional", 1, log(log(n_features * n_groups)))

  return(-2 * loglik + c_n * n_parameters * log(n_samples))
}

#' Summary method for a \code{coxtrans} object
#'
#' @param object An object of class \code{coxtrans}.
#' @param conf.int A numeric value between 0 and 1 indicating the confidence
#' level of the confidence interval. Default is 0.95.
#' @param compressed Logical; if \code{TRUE}, the summary is compressed and
#' only includes the group-level coefficients. Default is \code{TRUE}.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class \code{summary.coxtrans}, with the following
#' components:
#' \item{\code{n}, \code{nevent}}{Number of observations and number of events,
#' respectively, in the fit.}
#' \item{\code{logLik}}{The log partial likelihood at the final value.}
#' \item{\code{BIC}}{The Bayesian Information Criterion at the final value.}
#' \item{\code{coefficients}}{A matrix with one row for each coefficient, and
#' columns containing the coefficient, the hazard ratio exp(coef), standard
#' error, Wald statistic, and P value.}
#' \item{\code{conf.int}}{A matrix with one row for each coefficient, containing
#' the confidence limits for exp(coef).}
#'
#' @export
summary.coxtrans <- function(object, conf.int = 0.95, compressed = TRUE, ...) {
  # Extract necessary components from the object
  n_samples <- nrow(object$x)
  n_events <- sum(object$status)
  loglik <- logLik(object)
  bic_value <- BIC(object)
  group_levels <- levels(object$group)
  variable_names <- colnames(object$x)

  # Standard errors
  vcov_matrix <- vcov(object)
  if (is.null(vcov_matrix)) {
    stop("Variance-covariance matrix is not available.")
  }
  se <- sqrt(diag(vcov_matrix))

  coefficients <- coef(object)
  z_scores <- coefficients / se
  p_values <- stats::pchisq(z_scores^2, 1, lower.tail = FALSE)
  coef_matrix <- cbind(
    coefficients, exp(coefficients), se, z_scores, p_values
  )
  dimnames(coef_matrix) <- list(
    names(coefficients), c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
  )

  z <- stats::qnorm((1 + conf.int) / 2)
  conf_int_matrix <- cbind(
    exp(coefficients), exp(-coefficients),
    exp(coefficients - z * se), exp(coefficients + z * se)
  )
  dimnames(conf_int_matrix) <- list(
    names(coefficients), c(
      "exp(coef)", "exp(-coef)",
      paste("lower .", round(100 * conf.int, 2), sep = ""),
      paste("upper .", round(100 * conf.int, 2), sep = "")
    )
  )

  if (!compressed) {
    coef_matrix_extract <- matrix(nrow = 0, ncol = ncol(coef_matrix))
    conf_int_matrix_extract <- matrix(nrow = 0, ncol = ncol(conf_int_matrix))
    coef_groups <- c()
    coef_variables <- c()
    coef_names <- rownames(coef_matrix)

    extract_var_group <- function(coef_name) {
      variable_name <- trimws(stringr::str_extract(coef_name, "^[^\\(]+"))
      group_names <- stringr::str_extract(coef_name, "\\(([^)]+)\\)")
      group_names <- gsub("[()]", "", trimws(group_names))
      group_names_split <- stringr::str_split(group_names, ", ")[[1]]
      if ("ALL" %in% group_names_split) {
        group_names_split <- group_levels
      }
      list(variable = variable_name, groups = group_names_split)
    }

    for (i in seq_along(coef_names)) {
      var_group_info <- extract_var_group(coef_names[i])
      variable_name <- var_group_info$variable
      group_names <- var_group_info$groups
      for (group_name in group_names) {
        coef_groups <- c(coef_groups, trimws(group_name))
        coef_variables <- c(coef_variables, variable_name)
        coef_matrix_extract <- rbind(
          coef_matrix_extract, coef_matrix[i, , drop = FALSE]
        )
        conf_int_matrix_extract <- rbind(
          conf_int_matrix_extract, conf_int_matrix[i, , drop = FALSE]
        )
      }
    }

    coef_names <- paste0(coef_groups, " ", coef_variables)
    all_combinations <- expand.grid(group_levels, variable_names)
    all_combinations <- apply(all_combinations, 1, paste, collapse = " ")
    missing_names <- setdiff(all_combinations, coef_names)

    if (length(missing_names) > 0) {
      for (name in missing_names) {
        coef_matrix_extract <- rbind(coef_matrix_extract, c(0, 1, NA, NA, NA))
        conf_int_matrix_extract <- rbind(
          conf_int_matrix_extract, c(1, 1, NA, NA)
        )
        coef_names <- c(coef_names, name)
      }
    }

    coef_names_split <- do.call(rbind, stringr::str_split(coef_names, " "))
    coef_names_df <- as.data.frame(coef_names_split)
    coef_names_df[, 1] <- factor(coef_names_df[, 1], levels = group_levels)
    coef_order <- order(coef_names_df[, 1])
    coef_names <- apply(coef_names_df[coef_order, ], 1, paste, collapse = " ")

    coef_matrix <- coef_matrix_extract[coef_order, ]
    conf_int_matrix <- conf_int_matrix_extract[coef_order, ]
    rownames(coef_matrix) <- coef_names
    rownames(conf_int_matrix) <- coef_names
  }

  # Create a summary list
  summary_list <- list(
    n = n_samples,
    nevent = n_events,
    logLik = loglik,
    call = object$call,
    BIC = bic_value,
    coefficients = coef_matrix,
    conf.int = conf_int_matrix
  )

  class(summary_list) <- "summary.coxtrans"
  return(summary_list)
}

#' Print method for a \code{summary.coxtrans} object
#' @description This function prints a summary of the results of a
#' \code{coxtrans} model in a formatted and user-friendly manner, including the
#' model call, number of samples, number of events, coefficients, and confidence
#' intervals. It also includes
#' significance stars for p-values.
#' @param x A summary object produced from a fitted \code{coxtrans} model. This
#' object contains information such as model coefficients and confidence
#' intervals.
#' @param digits An integer controlling the number of significant digits to
#' print for numeric values.
#' @param signif.stars Logical; if \code{TRUE}, significance stars are printed
#' along with the p-values.
#' @param ... Additional arguments (not unused).
#' @return The function prints the summary of the \code{coxtrans} model and
#' returns the object \code{x} invisibly.
#' @details The function provides a formatted output that includes:
#' \itemize{
#'   \item \strong{Call:} The original function call that produced the model.
#'   \item \strong{n and events:} The total number of samples and the number of
#'         events (e.g., deaths).
#'   \item \strong{Coefficients:} The regression coefficients, their standard
#'         errors, z-values, and p-values, formatted in a table. Significance
#'         stars are shown next to p-values if \code{signif.stars} is
#'         \code{TRUE}.
#'   \item \strong{Confidence intervals:} The exponentiated coefficients along
#'         with their confidence intervals.
#' }
#' @export
print.summary.coxtrans <- function(
    x, digits = max(getOption("digits") - 3, 3),
    signif.stars = getOption("show.signif.stars"), ...) {
  # Print call
  cat("Call:\n")
  print(x$call)
  cat("\n")

  # Print number of samples and events
  cat("  n=", x$n, ", number of events=", x$nevent, "\n\n", sep = "")

  # Print coefficients with formatted output
  stats::printCoefmat(x$coefficients,
    digits = digits, signif.stars = signif.stars,
    cs.ind = 1:3, tst.ind = 4, P.values = TRUE, has.Pvalue = TRUE
  )

  # Print confidence intervals
  print(format(x$conf.int, digits = digits), quote = FALSE)

  invisible(x)
}

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
predict.coxtrans <- function(
    object, newdata = NULL, newgroup = NULL,
    type = c("lp", "terms", "risk", "expected", "survival"), ...) {
  type <- match.arg(type)
  x <- stats::model.matrix(object$formula, newdata)[, -1]
  group <- factor(newgroup, levels = levels(object$group))

  # Properties of the coxtrans object
  group_levels <- levels(object$group)
  n_groups <- length(unique(object$group))
  coefficients <- object$coefficients
  coefficients <- sweep(coefficients, 1, attr(object$x, "scale"), "*")

  beta <- coefficients[, 1:n_groups] + coefficients[, (n_groups + 1)]

  lp <- numeric(nrow(x))
  for (k in seq_len(n_groups)) {
    idx <- which(group == group_levels[k])
    if (length(idx) > 0) {
      lp[idx] <- x[idx, ] %*% beta[, k]
    }
  }

  if (type == "lp") {
    return(lp)
  } else if (type == "risk") {
    risk <- exp(lp)
    return(risk)
  } else {
    stop("type must be one of 'lp' or 'risk'")
  }
}

#' Predict the cumulative baseline hazard function for \code{coxtrans} objects
#'
#' @param object An object of class \code{coxtrans}.
#' @param newdata A numeric vector of time points at which to predict the
#' baseline hazard function. If \code{NULL}, the function will predict the
#' baseline hazard function at the unique event times in the fitted data.
#' @param ... Additional arguments (not Additional arguments (not unused).).
#'
#' @return A \code{data.frame} with one row for each time point, and columns
#' containing the event time, the cumulative baseline hazard function, and the
#' strata.
#' @export
basehaz.coxtrans <- function(object, newdata, ...) {
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
  basehaz_df
}
