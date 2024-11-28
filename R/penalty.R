penalty <- function(x, penalty, lambda, gamma) {
  x_abs <- abs(x)
  switch(penalty,
    lasso = lambda * sum(x_abs),
    MCP = sum(
      ifelse(
        x_abs <= lambda * gamma,
        lambda * x_abs - 0.5 * x_abs^2 / gamma,
        0.5 * lambda * gamma^2
      )
    ),
    SCAD = sum(
      ifelse(
        x_abs <= lambda,
        lambda * x_abs,
        ifelse(
          x_abs <= gamma * lambda,
          (2 * gamma * lambda * x_abs - x_abs^2 - lambda^2) / (2 * (gamma - 1)),
          0.5 * (gamma + 1) * lambda^2
        )
      )
    ),
    cli::cli_alert_danger("Invalid penalty")
  )
}
