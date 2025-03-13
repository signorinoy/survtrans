pair_index <- function(i, j, s) {
  if (is.na(i) || is.na(j) || i > s || j > s || i == j) {
    cli::cli_alert_danger("Invalid input")
    return(NULL)
  }

  if (i > j) {
    tmp <- i
    i <- j
    j <- tmp
  }

  if (i == 1) {
    pos <- j - 1
  } else {
    pos <- (i - 1) * (2 * s - i) / 2 + (j - i)
  }
  pos
}
