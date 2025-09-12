# Helpers for testing saemvsInit ----

# Valid intercept
.intercept_valid <- c(0.5, -1)

# Valid beta matrices (match intercept length = 2)
.beta_forced_valid <- matrix(1:4, nrow = 2, ncol = 2)
.beta_candidates_valid <- matrix(5:8, nrow = 2, ncol = 2)

# Valid covariance matrix (2x2, symmetric, positive definite)
.cov_re_valid <- diag(2)

# Invalid covariance matrices
.cov_re_not_square <- matrix(1:6, nrow = 2, ncol = 3)
.cov_re_not_symmetric <- matrix(c(1, 2, 0, 1), nrow = 2)
.cov_re_not_posdef <- matrix(c(1, 0, 0, 0), nrow = 2)

# Valid sigma2
.sigma2_valid <- 0.8

# Invalid sigma2
.sigma2_negative <- -1
.sigma2_vector <- c(1, 2)
