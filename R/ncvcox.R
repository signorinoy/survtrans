#' Non-convex penalized Cox proportional hazards model
#'
#' @param formula A formula object, with the response on the left of a \code{~}
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function.
#' @param data A data frame containing the variables in the model.
#' @param group A factor variable indicating the group of each observation.
#' @param lambda A non-negative value specifying the penalty parameter. The
#' default is 0.
#' @param penalty A character string specifying the penalty function. The
#' default is "lasso". Other options are "MCP" and "SCAD".
#' @param gamma A non-negative value specifying the penalty parameter. The
#' default is 3.7 for SCAD and 3.0 for MCP.
#' @param init A numeric vector of initial values for the coefficients. The
#' default is a zero vector.
#' @param control An object of class \link{survtrans_control} containing control
#' parameters for the fitting algorithm. Default is
#' \code{survtrans_control(...)}.
#' @param ... Additional arguments to be passed to the fitting algorithm.
#' @return An object of class \code{ncvcox}.
#' @export
#' @examples
#' formula <- survival::Surv(time, status) ~ . - group - id
#' df <- sim2[sim2$group == 2 | sim2$group == 4, ]
#' fit <- ncvcox(formula, df, df$group, lambda = 0.1, penalty = "SCAD")
ncvcox <- function(
    formula, data, group, lambda = 0,
    penalty = c("lasso", "MCP", "SCAD"),
    gamma = switch(penalty,
      SCAD = 3.7,
      MCP = 3,
      1
    ), init, control, ...) {
  # Load the data
  data <- preprocess(formula, data, group)
  x <- data$x
  x_scale <- attr(x, "scale")
  time <- data$time
  status <- data$status
  group <- data$group

  # Properties of the data
  n_samples <- nrow(x)
  n_features <- ncol(x)
  n_groups <- length(unique(group))
  group_levels <- levels(group)
  group_idxs <- lapply(group_levels, function(x) which(group == x))

  risk_set_size <- stats::ave(rep(1, n_samples), group, FUN = cumsum)
  risk_set_size <- unlist(lapply(1:n_groups, function(k) {
    idx <- group_idxs[[k]]
    stats::ave(risk_set_size[idx], time[idx], FUN = max)
  }))
  null_deviance <- -sum(status * log(risk_set_size))

  # Check the penalty argument
  penalty <- match.arg(penalty, choices = c("lasso", "MCP", "SCAD"))

  # Check the init argument
  if (!missing(init) && length(init) > 0) {
    if (length(init) != n_features) {
      stop("Wrong length for initial values")
    }
  } else {
    init <- numeric(n_features)
  }

  # Check the control argument
  if (missing(control)) control <- survtrans_control(...)

  # Initialize the training process
  n_iterations <- 0
  message <- ""
  convergence <- FALSE
  beta <- init

  offset <- x %*% beta
  w <- numeric(n_samples)
  z <- numeric(n_samples)

  x2 <- x**2

  repeat {
    n_iterations <- n_iterations + 1
    beta_prev <- beta

    # Calculate the weights and residuals
    for (k in 1:n_groups) {
      idx <- group_idxs[[k]]
      wls <- approx_likelihood(
        offset = offset[idx], time = time[idx], status = status[idx]
      )
      w[idx] <- wls$weights
      z[idx] <- wls$residuals + offset[idx]
    }
    xw <- x * w
    xwx <- colMeans(x2 * w)

    # Update beta by cyclic coordinate descent
    repeat {
      max_diff <- 0
      for (j in seq_len(n_features)) {
        beta_j <- beta[j]
        r <- z - x[, -j, drop = FALSE] %*% beta[-j]
        beta[j] <- close_update(
          mean(xw[, j] * r), xwx[j], penalty, lambda, gamma
        )
        delta_beta <- beta[j] - beta_j
        if (delta_beta != 0) {
          max_diff <- max(max_diff, abs(delta_beta))
        }
      }
      if (max_diff <= control$abstol) break
    }

    # Calculate the log-likelihood
    offset <- x %*% beta
    hazard <- exp(offset)
    risk_set <- stats::ave(hazard, group, FUN = cumsum)
    for (k in 1:n_groups) {
      idx <- group_idxs[[k]]
      risk_set[idx] <- stats::ave(risk_set[idx], time[idx], FUN = max)
    }
    loss <- sum(status * (offset - log(risk_set)))
    loss_penalty <- penalty(beta, penalty, lambda, gamma) * n_samples
    loss_total <- loss - loss_penalty
    if (control$verbose) {
      cat(
        "========================================\n",
        sprintf("Iteration Number       : %d", n_iterations), "\n",
        sprintf("Total Loss             : %.4f", loss_total), "\n",
        sprintf(" - Log Likelihood      : %.4f", loss), "\n",
        sprintf(" - Penalty             : %.4f", loss_penalty), "\n",
        "========================================\n"
      )
    }

    # Check the convergence
    if (is.infinite(loss) || is.nan(loss)) {
      stop("The log-likelihood is not finite. Stopping the algorithm.")
    }
    if (n_iterations >= control$maxit) {
      convergence <- TRUE
      message <- paste0(
        "Maximum number of iterations reached (", control$maxit, ")."
      )
    }
    if (max(abs(beta - beta_prev)) <= control$abstol) {
      convergence <- TRUE
      message <- paste0(
        "Convergence reached at iteration ", n_iterations, "."
      )
    }
    if (loss / null_deviance < 0.01) {
      convergence <- TRUE
      message <- paste0(
        "The log-likelihood is too small (", loss / null_deviance,
        "). Stopping the algorithm."
      )
    }
    if (convergence) break
  }

  # Unstandardize the beta
  coefficients <- beta / x_scale
  names(coefficients) <- colnames(x)

  # Return the fit
  fit <- list(
    coefficients = beta,
    logLik = loss,
    iter = n_iterations,
    message = message,
    penalty = penalty,
    lambda = lambda,
    gamma = gamma,
    formula = formula,
    call = match.call(),
    time = time,
    status = status,
    group = group,
    x = x
  )
  class(fit) <- "ncvcox"
  return(fit)
}
