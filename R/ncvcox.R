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
#' summary(fit)
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
  r <- numeric(n_samples)

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
      r[idx] <- wls$residuals
    }
    xw <- x * w
    xwx <- colMeans(x2 * w)

    # Update beta by cyclic coordinate descent
    repeat {
      max_diff <- 0
      for (j in seq_len(n_features)) {
        beta_j <- beta[j]
        u <- mean(xw[, j] * r) + xwx[j] * beta_j
        v <- xwx[j]
        beta[j] <- close_update(u, v, penalty, lambda, gamma)
        shift <- beta[j] - beta_j
        if (shift != 0) {
          r <- r - x[, j] * shift
          max_diff <- max(max_diff, abs(shift) * sqrt(v))
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

#' Extract the coefficients from a \code{ncvcox} object
#' @param object An object of class \code{ncvcox}.
#' @param ... Additional arguments (not unused).
#' @return A named numeric vector containing the coefficients of the fitted
#' \code{ncvcox} object. Zero coefficients are removed.
#' @export
coef.ncvcox <- function(object, ...) {
  coefficients <- object$coefficients
  names(coefficients) <- colnames(object$x)
  coefficients[coefficients != 0]
}

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
  x <- as.matrix(x[, coefficients != 0, drop = FALSE])
  n_nonzero <- length(coefs)

  # Calculate the gradient and Hessian for the non-zero coefficients
  lp <- x %*% coefs
  gradients <- matrix(0, nrow = n_samples, ncol = n_nonzero)
  hessians <- matrix(0, nrow = n_samples, ncol = n_nonzero^2)
  for (k in 1:n_groups) {
    idx <- group_idxs[[k]]
    ghs <- calc_grad_hess(
      lp[idx], x[idx, , drop = FALSE], time[idx], status[idx]
    )
    gradients[idx, ] <- ghs$grad
    hessians[idx, ] <- ghs$hess
  }
  hess <- matrix(colSums(hessians), n_nonzero, n_nonzero)
  cov_grad <- stats::cov(gradients) * n_samples
  hess_inv <- solve(hess)

  hess_inv %*% cov_grad %*% hess_inv
}

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
  coefficients <- object$coefficients * attr(x, "scale")

  # Calculate the log-likelihood
  offset <- x %*% coefficients
  hazard <- exp(offset)
  risk_set <- stats::ave(hazard, group, FUN = cumsum)
  risk_set <- unlist(lapply(seq_len(n_groups), function(k) {
    idx <- group_idxs[[k]]
    stats::ave(risk_set[idx], time[idx], FUN = max)
  })) # Update the risk set for each group based on unique time points
  sum(status * (offset - log(risk_set)))
}

#' BIC for a \code{ncvcox} object
#'
#' @param object An object of class \code{ncvcox}.
#' @param type A character string specifying the type of BIC to compute.
#' "traditional" corresponds to Cn=1, and "modified" corresponds to
#' Cn=log(log(d)).
#' @param ... Additional arguments (not unused).
#' @return A numeric value representing the BIC of the fitted \code{ncvcox}
#' object.
#' @export
BIC.ncvcox <- function(object, type = c("traditional", "modified"), ...) {
  type <- match.arg(type)

  # Properties of the ncvcox object
  coefficients <- object$coefficients
  n_samples <- nrow(object$x)
  n_features <- ncol(object$x)

  # The number of parameters should minus the number of active constraints
  n_parameters <- sum(coefficients != 0)

  # Log-likelihood of the model
  loglik <- logLik(object)

  # The parameter of the BIC
  c_n <- ifelse(type == "traditional", 1, log(log(n_features)))

  -2 * loglik + c_n * n_parameters * log(n_samples)
}

#' Summary method for a \code{ncvcox} object
#'
#' @param object An object of class \code{ncvcox}.
#' @param conf.int A numeric value between 0 and 1 indicating the confidence
#' level of the confidence interval. Default is 0.95.
#' @param compressed Logical; if \code{TRUE}, the summary is compressed and
#' only includes the group-level coefficients. Default is \code{TRUE}.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class \code{summary.ncvcox}, with the following
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
summary.ncvcox <- function(object, conf.int = 0.95, compressed = TRUE, ...) {
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
    coef_matrix_extract <- coef_matrix
    conf_int_matrix_extract <- conf_int_matrix

    coef_names <- names(coefficients)
    variable_names <- colnames(object$x)
    missing_names <- setdiff(variable_names, coef_names)

    if (length(missing_names) > 0) {
      for (name in missing_names) {
        coef_matrix_extract <- rbind(coef_matrix_extract, c(NA, NA, NA, NA, NA))
        conf_int_matrix_extract <- rbind(
          conf_int_matrix_extract, c(NA, NA, NA, NA)
        )
        coef_names <- c(coef_names, name)
      }
    }
    coef_order <- match(variable_names, coef_names)
    coef_matrix <- coef_matrix_extract[coef_order, ]
    conf_int_matrix <- conf_int_matrix_extract[coef_order, ]
    rownames(coef_matrix) <- variable_names
    rownames(conf_int_matrix) <- variable_names
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

  class(summary_list) <- "summary.ncvcox"
  summary_list
}

#' Print method for a \code{summary.ncvcox} object
#' @description This function prints a summary of the results of a
#' \code{ncvcox} model in a formatted and user-friendly manner, including the
#' model call, number of samples, number of events, coefficients, and confidence
#' intervals. It also includes
#' significance stars for p-values.
#' @param x A summary object produced from a fitted \code{ncvcox} model. This
#' object contains information such as model coefficients and confidence
#' intervals.
#' @param digits An integer controlling the number of significant digits to
#' print for numeric values.
#' @param signif.stars Logical; if \code{TRUE}, significance stars are printed
#' along with the p-values.
#' @param ... Additional arguments (not unused).
#' @return The function prints the summary of the \code{ncvcox} model and
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
print.summary.ncvcox <- function(
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
predict.ncvcox <- function(
    object, newdata = NULL, newgroup = NULL,
    type = c("lp", "terms", "risk", "expected", "survival"), ...) {
  type <- match.arg(type)
  x <- stats::model.matrix(object$formula, newdata)[, -1]

  # Properties of the coxtrans object
  coefficients <- object$coefficients * attr(x, "scale")

  lp <- x %*% coefficients

  if (type == "lp") {
    return(lp)
  } else if (type == "risk") {
    risk <- exp(lp)
    return(risk)
  } else {
    stop("type must be one of 'lp' or 'risk'")
  }
}
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
  basehaz_df
}
