# helper-prepare_functions.R

# Simple valid saemvsData object
create_simple_data <- function() {
  y <- list(rnorm(10))
  t <- list(seq_len(10))
  x_candidates <- matrix(1:10, nrow = 1)
  saemvsData(y = y, t = t, x_candidates = x_candidates)
}

# Simple valid saemvsInit object
create_simple_init <- function() {
  saemvsInit(
    intercept = c(0),
    beta_forced = NULL,
    beta_candidates = NULL,
    cov_re = diag(1),
    sigma2 = 1
  )
}
