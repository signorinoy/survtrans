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
  coef <- coef[coef != 0]
  return(coef)
}
