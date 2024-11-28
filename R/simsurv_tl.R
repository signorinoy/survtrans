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
  return(df)
}
