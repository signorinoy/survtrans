#' Generic function for basehaz
#'
#' @param object Any object.
#' @param ... Additional arguments.
#'
#' @return A numeric vector of baseline hazard.
#' @export
basehaz <- function(object, ...) {
  UseMethod("basehaz")
}

#' Generic function for diagnose
#'
#' @param object Any object.
#' @param ... Additional arguments.
#' @return A ggplot object.
#' @export
diagnose <- function(object, ...) {
  UseMethod("diagnose")
}

#' Preprocess Survival Data
#'
#' This function preprocesses survival data for a Cox transfer analysis. It
#' performs several steps including extracting the response and covariates from
#' a model frame, standardizing the covariates, ensuring the validity of offsets
#' and grouping, and sorting the data by survival time in descending order.
#'
#' @param formula A formula specifying the survival model
#' (e.g., Surv(time, status) ~ covariates).
#' @param data A data frame containing the variables referenced in the formula.
#' @param group Optional grouping vector.
#' @param offset Optional offset vector.
#'
#' @details
#' The function first creates a model frame from the provided formula and data,
#' then extracts the response variable, which should contain the survival time
#' and censoring status. The covariates are extracted using a model matrix and
#' standardized using the scale function. The function also checks the
#' \code{group} and \code{offset} arguments, assigning default values if they
#' are not provided. Finally, all components (time, status, covariates, group,
#' and offset) are sorted in descending order based on the survival time.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{x}{A matrix of standardized covariates.}
#'   \item{time}{A vector of survival times, sorted in descending order.}
#'   \item{status}{A vector of censoring indicators, corresponding to the sorted
#'                 survival times.}
#'   \item{group}{A factor vector representing the group classification for each
#'                sample, sorted by time.}
#'   \item{offset}{A numeric vector representing offsets for each sample, sorted
#'                 by time.}
#' }
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
  if (missing(offset) || is.null(offset)) offset <- rep(0.0, n_samples)
  if (missing(group) || is.null(group)) group <- rep(0, n_samples)
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


#' Compute penalty for a numeric vector using specified regularization type.
#'
#' This function computes the penalty for a given parameter vector `x` based on
#' the specified penalty approach.
#'
#' @param x A numeric vector for which the penalty is computed.
#' @param penalty A character string specifying the penalty type. Valid options
#' are "lasso", "MCP", or "SCAD".
#' @param lambda A numeric value representing the regularization parameter.
#' @param gamma A numeric value used in the penalty for MCP and SCAD.
#'
#' @return A numeric value representing the computed penalty.
#'
#' @details The computation differs according to the penalty type:
#'   \describe{
#'     \item{lasso}{Returns the L1 norm of `x` scaled by `lambda`.}
#'     \item{MCP}{Applies a minimax concave penalty where the penalty function
#'                behaves linearly when the absolute values are below
#'                `lambda * gamma` and quadratically otherwise.}
#'     \item{SCAD}{Applies the smoothly clipped absolute deviation penalty;
#'                 different expressions are used when `abs(x)` is below
#'                 `lambda`, between `lambda` and `gamma * lambda`, or above
#'                 these thresholds.}
#'   }
penalty <- function(x, penalty, lambda, gamma) {
  x_abs <- abs(x)
  switch(penalty,
    lasso = lambda * sum(x_abs),
    MCP = {
      condition <- x_abs <= lambda * gamma
      sum(lambda * x_abs[condition] - 0.5 * x_abs[condition]^2 / gamma) +
        sum(0.5 * lambda * gamma^2 * (!condition))
    },
    SCAD = {
      condition1 <- x_abs <= lambda
      condition2 <- x_abs > lambda & x_abs <= gamma * lambda
      sum(lambda * x_abs[condition1]) +
        sum(
          (2 * gamma * lambda * x_abs[condition2]
            - x_abs[condition2]^2 - lambda^2) / (2 * (gamma - 1))
        ) +
        sum(0.5 * (gamma + 1) * lambda^2 * (!condition1 & !condition2))
    },
    stop("Invalid penalty type. Please choose 'lasso', 'MCP', or 'SCAD'.")
  )
}

#' Calculate the maximum value of the penalty parameter lambda
#'
#' @param formula A formula expression for regression models, in the form
#' \code{response ~ predictors}. The response must be a survival object as
#' returned by the \code{\link{Surv}} function.
#' @param data A data frame containing the variables in the model.
#' @param group A factor specifying the group of each sample.
#' @param offset A numeric vector specifying the offset.
#'
#' @return The maximum value of the penalty parameter lambda, which shrinks all
#' the coefficients to zero.
#' @export
calc_lambda_max <- function(formula, data, group, offset) {
  # Load the data
  data <- preprocess(formula, data, group, offset)
  x <- data$x
  time <- data$time
  status <- data$status
  group <- data$group
  offset <- data$offset

  # Properties of the data
  group_levels <- levels(group)

  # Calculate the lambda_max
  lambdas_max <- numeric(length(group_levels))
  for (i in seq_along(group_levels)) {
    idx <- which(group == group_levels[i])
    wls <- approx_likelihood(offset[idx], time[idx], status[idx])
    if (length(idx) > 1) {
      xwr <- colMeans(sweep(x[idx, ], 1, wls$residuals * wls$weights, `*`))
    } else {
      xwr <- 0
    }
    lambdas_max[i] <- max(abs(xwr), na.rm = TRUE)
  }
  lambda_max <- max(lambdas_max, na.rm = TRUE)
  lambda_max
}
#' Simulate survival data for a multi-source Cox model
#'
#' @param beta A vector of length p representing the common coefficients,
#' where p is the number of features.
#' @param eta A matrix of size p x K representing the group-specific
#' coefficients, where K is the number of groups.
#' @param lambda A vector of length K representing the baseline hazard's scale
#' parameters.
#' @param gamma A vector of length K representing the baseline hazard's shape
#' parameters.
#' @param dist A string specifying the distribution of the baseline hazard,
#' either "exponential", "weibull", or "gompertz".
#' @param maxt A positive number specifying the maximum time to simulate.
#' @param n_samples A vector of length K specifying the number of samples per
#' group, or a single number specifying the total number of samples.
#' @param seed An integer specifying the random seed, with a default value of 0.
#'
#' @return A data frame with columns "id", "group", "X1", "X2", ..., "Xp",
#' "time", and "status".
#'
#' @export
#'
#' @examples
#' beta <- c(1, 1)
#' eta <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
#' lambda <- c(1, 2)
#' gamma <- c(2, 1)
#' dist <- c("gompertz", "weibull")
#' maxt <- 3
#' n_samples <- 100
#' df <- simsurv_tl(beta, eta, lambda, gamma, dist, maxt, n_samples)
#' head(df)
simsurv_tl <- function(
    beta, eta, lambda, gamma, dist, maxt, n_samples, seed = 0) {
  set.seed(seed)
  n_groups <- ncol(eta)
  n_features <- nrow(eta)
  names(beta) <- paste0("X", 1:n_features)
  if (length(n_samples) == 1) n_samples <- rep(n_samples, n_groups)
  mu <- rep(0, n_features)
  sigma <- diag(n_features)
  covs <- c()
  times <- c()
  for (k in 1:n_groups) {
    cov <- data.frame(
      id = seq_len(n_samples[k]),
      MASS::mvrnorm(n_samples[k], mu = mu, Sigma = sigma)
    )
    time <- simsurv::simsurv(
      lambdas = lambda[k], gammas = gamma[k], dist = dist[k],
      x = cov, betas = beta + eta[, k], maxt = maxt
    )
    cov$group <- k
    time$group <- k
    covs <- rbind(covs, cov)
    times <- rbind(times, time)
  }
  df <- merge(covs, times, by = c("id", "group"), all.x = TRUE)
  names(df)[names(df) == "eventtime"] <- "time"
  df$id <- seq_len(nrow(df))
  df
}
