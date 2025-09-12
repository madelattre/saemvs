test_that("saemvsInit accepts valid inputs", {
  init <- saemvsInit(
    intercept = .intercept_valid,
    beta_forced = .beta_forced_valid,
    beta_candidates = .beta_candidates_valid,
    cov_re = .cov_re_valid,
    sigma2 = .sigma2_valid
  )

  expect_s4_class(init, "saemvsInit")
  expect_equal(init@intercept, .intercept_valid)
  expect_equal(init@sigma2, .sigma2_valid)
  expect_equal(ncol(init@beta_forced), length(.intercept_valid))
})

test_that("beta_forced columns must match intercept length", {
  bad_beta <- matrix(1:4, nrow = 2, ncol = 2) # 2 cols, intercept length will be 3
  expect_error(
    saemvsInit(
      intercept = c(0.1, 0.2, 0.3),
      beta_forced = bad_beta,
      cov_re = diag(3)
    ),
    "Number of columns in 'beta_forced'"
  )
})

test_that("beta_candidates columns must match intercept length", {
  bad_beta <- matrix(1:9, nrow = 3, ncol = 3)
  expect_error(
    saemvsInit(
      intercept = c(0.1, 0.2),
      beta_candidates = bad_beta,
      cov_re = diag(2)
    ),
    "Number of columns in 'beta_candidates'"
  )
})

test_that("cov_re must be a valid covariance matrix", {
  intercept <- c(0.1, 0.2)

  expect_error(
    saemvsInit(intercept, cov_re = .cov_re_not_square),
    "cov_re"
  )
  expect_error(
    saemvsInit(intercept, cov_re = .cov_re_not_symmetric),
    "cov_re"
  )
  expect_error(
    saemvsInit(intercept, cov_re = .cov_re_not_posdef),
    "cov_re"
  )
})

test_that("sigma2 must be strictly positive scalar", {
  intercept <- c(0.1, 0.2)

  expect_error(
    saemvsInit(intercept, cov_re = diag(2), sigma2 = .sigma2_negative),
    "sigma2"
  )
  expect_error(
    saemvsInit(intercept, cov_re = diag(2), sigma2 = .sigma2_vector),
    "sigma2"
  )

  ok <- saemvsInit(intercept, cov_re = diag(2), sigma2 = .sigma2_valid)
  expect_equal(ok@sigma2, .sigma2_valid)
})
