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
