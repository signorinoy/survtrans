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
#' @param lambda2 A non-negative value specifying the global biased penalty
#' parameter. The default is 0.
#' @param lambda3 A non-negative value specifying the local biased penalty
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
#'   formula, sim2, sim2$group, 1,
#'   lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
#' )
#' summary(fit)
coxtrans <- function(
    formula, data, group, target,
    lambda1 = 0.0, lambda2 = 0.0, lambda3 = 0.0,
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
  if (lambda1 < 0 || lambda2 < 0 || lambda3 < 0) {
    stop("Lambda parameters must be non-negative")
  }
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
      eta[global_idx], vartheta, penalty, lambda2, gamma
    )
    eta[local_idx] <- threshold_prox(
      eta[local_idx], vartheta, penalty, lambda3, gamma
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
      message <- stringr::str_glue(
        "Maximum number of iterations reached ({control$maxit})."
      )
    }
    if (r_norm < eps_pri && s_norm < eps_dual) {
      convergence <- TRUE
      message <- stringr::str_glue(
        "Convergence reached at iteration {n_iterations}."
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
      message <- stringr::str_glue(
        "The log-likelihood is too small ({loss / null_deviance}). ",
        "Stopping the algorithm."
      )
    }

    eta_ <- contr_penalty %*% theta
    loss_penalty <- penalty(eta_[sparse_idx], penalty, lambda1, gamma) +
      penalty(eta_[global_idx], penalty, lambda2, gamma) +
      penalty(eta_[local_idx], penalty, lambda3, gamma)
    loss_penalty <- loss_penalty * n_samples_total
    loss_total <- loss - loss_penalty

    if (control$verbose) {
      cli::cli_h2("Iteration Info")
      cli::cli_text("Iteration       : {.val {n_iterations}}")

      cli::cli_h2("Residuals")
      cli::cli_text(stringr::str_glue(
        "Primal Residual : {r_norm}  (Tol: {eps_pri})"
      ))
      cli::cli_text(stringr::str_glue(
        "Dual Residual   : {s_norm}  (Tol: {eps_dual})"
      ))

      cli::cli_h2("Optimization Parameters")
      cli::cli_text(stringr::str_glue("Augmented Param : {vartheta}"))

      cli::cli_h2("Loss Summary")
      cli::cli_text(stringr::str_glue("Total Loss      : {loss_total}"))
      cli::cli_text(stringr::str_glue("     - LogLik       : {loss}"))
      cli::cli_text(stringr::str_glue("     - Penalty      : {loss_penalty}"))

      cli::cli_rule()
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

  coefficients <- sweep(coefficients, 1, x_scale, "/")
  x <- sweep(x, 2, x_scale, "*")

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
  coefficients <- object$coefficients
  n_groups <- ncol(coefficients) - 1L
  n_features <- nrow(coefficients)
  beta <- coefficients[, 1:n_groups, drop = FALSE] +
    coefficients[, n_groups + 1L]

  phi_list <- lapply(
    seq_len(n_features),
    function(j) unique(beta[j, ])
  )

  is_global <- coefficients[, 1] == 0
  psi_list <- lapply(seq_len(n_features), function(j) {
    vals <- phi_list[[j]]
    if (is_global[j] && length(vals) > 1) vals[-1] else vals
  })
  psi <- unlist(psi_list)
  names(psi) <- unlist(lapply(seq_along(psi_list), function(j) {
    vals <- psi_list[[j]]
    if (is_global[j] && length(vals) > 1) {
      idx <- seq_along(vals) + 1
    } else {
      idx <- seq_along(vals)
    }
    stringr::str_c(rownames(coefficients)[j], idx, sep = ".")
  }))
  psi
}

#' Variance-covariance matrix for a \code{coxtrans} object.
#' @param object An object of class \code{coxstream}.
#' @param ... Additional arguments (unused).
#' @return A matrix representing the variance-covariance matrix of the
#' coefficients.
#' @export
vcov.coxtrans <- function(object, ...) {
  time <- object$time
  status <- object$status
  group <- object$group
  x <- object$x

  n_samples <- nrow(x)
  group_levels <- levels(group)
  n_groups <- length(group_levels)
  group_idxs <- lapply(group_levels, function(g) which(group == g))

  coefficients <- object$coefficients

  psi <- coef(object)
  link_matrix <- build_link_matrix(coefficients)

  z <- Matrix::bdiag(lapply(group_idxs, function(idx) x[idx, ])) %*% link_matrix
  time <- unlist(lapply(group_idxs, function(idx) time[idx]))
  status <- unlist(lapply(group_idxs, function(idx) status[idx]))

  is_nonzero <- as.vector(psi) != 0
  n_nonzero <- sum(is_nonzero)
  z1 <- as.matrix(z[, is_nonzero, drop = FALSE])
  psi1 <- psi[is_nonzero]
  lp <- z1 %*% psi1

  gradients <- matrix(0, nrow = n_samples, ncol = n_nonzero)
  hessians <- matrix(0, nrow = n_samples, ncol = n_nonzero^2)
  n_passes <- 0
  for (k in seq_len(n_groups)) {
    idx <- n_passes + seq_len(length(group_idxs[[k]]))
    n_passes <- n_passes + length(idx)
    ghs <- calc_grad_hess(lp[idx], z1[idx, ], time[idx], status[idx])
    gradients[idx, ] <- ghs$grad
    hessians[idx, ] <- ghs$hess
  }
  hess <- matrix(colSums(hessians), n_nonzero, n_nonzero)
  hess_inv <- solve(hess)
  grad_cov <- crossprod(gradients)
  vcov <- hess_inv %*% grad_cov %*% hess_inv
  dimnames(vcov) <- list(
    names(psi1), names(psi1)
  )
  vcov
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
  n_parameters <- length(coef(object))

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
#' @param target_only Logical; if \code{TRUE}, only the coefficients for the
#' target group are shown in the summary. Default is \code{TRUE}.
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
summary.coxtrans <- function(object, conf.int = 0.95, target_only = TRUE, ...) {
  # Extract necessary components from the object
  n_samples <- nrow(object$x)
  n_events <- sum(object$status)
  loglik <- logLik(object)
  bic_value <- BIC(object)

  # Standard errors
  vcov_matrix <- vcov(object)
  if (is.null(vcov_matrix)) {
    stop("Variance-covariance matrix is not available.")
  }

  coefficients <- coef(object)
  is_nonzero <- coefficients != 0
  coefficients <- coefficients[is_nonzero]
  if (target_only) {
    n_features <- nrow(object$coefficients)
    link_matrix <- as.matrix(build_link_matrix(object$coefficients)[
      seq_len(n_features), is_nonzero,
      drop = FALSE
    ])
    coefficients <- as.vector(link_matrix %*% coefficients)
    vcov_matrix <- link_matrix %*% vcov_matrix %*% t(link_matrix)
    names(coefficients) <- rownames(object$coefficients)
  }

  se <- sqrt(diag(vcov_matrix))
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
      stringr::str_c("lower .", round(100 * conf.int, 2)),
      stringr::str_c("upper .", round(100 * conf.int, 2), sep = "")
    )
  )

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

  is_nonzero <- x$coefficients[, "coef"] != 0

  # Print coefficients with formatted output
  stats::printCoefmat(x$coefficients[is_nonzero, , drop = FALSE],
    digits = digits, signif.stars = signif.stars,
    cs.ind = 1:3, tst.ind = 4, P.values = TRUE, has.Pvalue = TRUE
  )

  # Print confidence intervals
  print(
    format(x$conf.int[is_nonzero, , drop = FALSE], digits = digits),
    quote = FALSE
  )

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
