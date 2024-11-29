#' Cross-validation for \code{ncvcox}
#'
#' @param formula A formula object, with the response on the left of a \code{~}
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function.
#' @param data A data frame containing the variables in the model.
#' @param group A factor variable indicating the group of each observation.
#' @param penalty A character string specifying the penalty function. The
#' default is "lasso". Other options are "MCP" and "SCAD".
#' @param gamma A non-negative value specifying the penalty parameter. The
#' default is 3.7 for SCAD and 3.0 for MCP.
#' @param nfolds An integer specifying the number of folds.
#' @param nlambdas An integer specifying the number of lambda values.
#' @param lambda_min_ratio a numeric value specifying the minimum lambda value
#' as a fraction of lambda_max.
#' @param seed An integer specifying the random seed. Default is 0.
#' @param control An object of class \link{survtrans_control} containing control
#' parameters for the fitting algorithm. Default is
#' \code{survtrans_control(...)}.
#' @param ... Additional arguments to be passed to the fitting algorithm.
#' @return An object of class \code{cv.ncvcox}.
#' @export
#' @examples
#' formula <- survival::Surv(time, status) ~ . - group - id
#' df <- sim2[sim2$group == 2 | sim2$group == 4, ]
#' cvfit <- cv.ncvcox(formula, df, df$group, penalty = "SCAD")
#' coef(cvfit)
cv.ncvcox <- function(
    formula, data, group, penalty = c("lasso", "MCP", "SCAD"),
    gamma = switch(penalty,
      SCAD = 3.7,
      MCP = 3,
      1
    ), nfolds = 10, nlambdas = 100, lambda_min_ratio = NULL,
    seed = 0, control, ...) {
  # Load the data
  data_ <- data
  group_ <- group
  data <- preprocess(formula, data, group)
  x <- data$x
  x_scale <- attr(x, "scale")
  x <- sweep(x, 2, x_scale, "/")
  time <- data$time
  status <- data$status
  group <- data$group

  # Properties of the data
  n_samples <- nrow(x)
  n_features <- ncol(x)
  n_groups <- length(unique(group))
  group_levels <- levels(group)
  group_idxs <- lapply(group_levels, function(x) which(group == x))

  # Check the arguments
  penalty <- match.arg(penalty, choices = c("lasso", "MCP", "SCAD"))
  if (is.null(lambda_min_ratio)) {
    lambda_min_ratio <- ifelse(n_samples < n_features, 0.01, 1e-04)
  }
  if (missing(control)) control <- survtrans_control(...)

  # Determmine the lambda sequence
  lambda_max <- calc_lambda_max(formula, data_, group_)
  lambda_min <- lambda_max * lambda_min_ratio
  lambdas <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambdas))

  coefficients <- matrix(0, nrow = nlambdas, ncol = n_features)
  set.seed(seed)
  idx <- sample(1:nfolds, n_samples, replace = TRUE)

  criterions <- matrix(0, nrow = nlambdas, ncol = nfolds)
  for (i in seq_along(lambdas)) {
    fit <- ncvcox(
      formula, data_, group_,
      lambda = lambdas[i], penalty = penalty,
      gamma = gamma, control = control, ...
    )
    coefficients[i, ] <- fit$coefficients
    for (k in 1:nfolds) {
      fit <- ncvcox(
        formula = formula, data = data_[idx != k, ], group = group_[idx != k],
        lambda = lambdas[i], penalty = penalty, gamma = gamma,
        control = control, ...
      )
      offset <- x %*% fit$coefficients
      hazard <- exp(offset)
      risk_set <- stats::ave(hazard, group, FUN = cumsum)
      for (j in 1:n_groups) {
        risk_set[group_idxs[[j]]] <- stats::ave(
          risk_set[group_idxs[[j]]], time[group_idxs[[j]]],
          FUN = max
        )
      }
      loss <- sum(status * (offset - log(risk_set)))
      criterions[i, k] <- loss - fit$logLik
    }
  }
  colnames(coefficients) <- colnames(x)
  colnames(criterions) <- paste0("Fold", 1:nfolds)
  cvm <- rowMeans(criterions)
  cvsd <- apply(criterions, 1, stats::sd)
  lambda.min <- lambdas[which.max(cvm)]
  lambda.1se <- lambdas[which.max(cvm + cvsd)]

  fit <- list(
    lambdas = lambdas, coefficients = coefficients, criterions = criterions,
    name = "Cross-validated log-likelihood", cvm = cvm, cvsd = cvsd,
    lambda.min = lambda.min, lambda.1se = lambda.1se, call = match.call()
  )
  class(fit) <- "cv.ncvcox"
  return(fit)
}
