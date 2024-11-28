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
