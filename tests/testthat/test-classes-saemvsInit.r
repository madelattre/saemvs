# =========================================================
# test-classes-saemvsInit.r
#
# Unit tests for saemvsInit class and saemvsInit constructor
# (see classes.r) using testthat
# =========================================================

# Valid intercept
intercept_2 <- c(1, 2)
intercept_empty <- numeric(0)

# Valid beta matrices
beta_forced_valid <- matrix(1:4, 2, 2)
beta_candidates_valid <- matrix(5:8, 2, 2)

# Invalid beta matrices
beta_forced_bad <- matrix(1:3, 1, 3)
beta_candidates_bad <- matrix(1:3, 1, 3)

# Covariance matrices
cov_re_valid <- diag(2)
cov_re_non_square <- matrix(1:6, 2, 3)
cov_re_neg_diag <- diag(c(-1, 2))
cov_re_zero_diag_nonzero <- matrix(c(0, 1, 0, 2), 2, 2)

# sigma2
sigma2_valid <- 1
sigma2_zero <- 0
sigma2_vector <- c(1, 2)

test_that("Valid creation with default = FALSE", {
  obj <- saemvsInit(
    intercept = intercept_2,
    beta_forced = beta_forced_valid,
    beta_candidates = beta_candidates_valid,
    cov_re = cov_re_valid,
    sigma2 = sigma2_valid
  )
  expect_s4_class(obj, "saemvsInit")
})

test_that("Valid creation with default = TRUE", {
  obj <- saemvsInit(intercept = intercept_2, default = TRUE)
  expect_s4_class(obj, "saemvsInit")
  expect_null(obj@beta_forced)
  expect_null(obj@beta_candidates)
})

test_that("Empty intercept raises error", {
  expect_error(saemvsInit(intercept = intercept_empty))
})

test_that("Beta_candidates column mismatch raises error", {
  expect_error(saemvsInit(
    intercept = intercept_2,
    beta_candidates = beta_candidates_bad
  ))
})

test_that("Beta_forced column mismatch raises error", {
  expect_error(saemvsInit(
    intercept = intercept_2,
    beta_forced = beta_forced_bad
  ))
})

test_that("cov_re not square raises error", {
  expect_error(saemvsInit(intercept = intercept_2, cov_re = cov_re_non_square))
})

test_that("cov_re negative diagonal raises error", {
  expect_error(saemvsInit(intercept = intercept_2, cov_re = cov_re_neg_diag))
})

test_that("cov_re zero diagonal with non-zero row/col raises error", {
  expect_error(saemvsInit(
    intercept = intercept_2,
    cov_re = cov_re_zero_diag_nonzero
  ))
})

test_that("sigma2 non-positive or length>1 raises error", {
  expect_error(saemvsInit(intercept = intercept_2, sigma2 = sigma2_zero))
  expect_error(saemvsInit(intercept = intercept_2, sigma2 = sigma2_vector))
})
