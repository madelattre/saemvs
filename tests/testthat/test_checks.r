# =========================================================
# test-checks.r
#
# Unit tests for checks.r using testthat
# =========================================================

cov_valid <- matrix(c(2, 1, 1, 2), nrow = 2) # positive definite, symmetric
cov_non_numeric <- matrix(c("a", "b", "c", "d"), 2, 2)
cov_non_square <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
cov_non_symmetric <- matrix(c(1, 2, 3, 4), nrow = 2)
cov_not_pos_def <- matrix(c(1, 2, 2, 1), nrow = 2) # eigenvalues less or equal to zero

# ---------------------------------------------
# --- Tests for check_covariance() ---
# ---------------------------------------------

test_that("check_covariance identifies numeric matrices", {
  expect_null(check_covariance(cov_valid, "PD matrix"))
  expect_match(
    check_covariance(cov_non_numeric, "Non-numeric"),
    "must be a numeric matrix"
  )
})

test_that("check_covariance identifies square matrices", {
  expect_match(
    check_covariance(cov_non_square, "Non-square"),
    "must be a square matrix"
  )
})

test_that("check_covariance identifies symmetric matrices", {
  expect_match(
    check_covariance(cov_non_symmetric, "Non-symmetric"),
    "must be symmetric"
  )
})

test_that("check_covariance identifies positive definite matrices", {
  expect_null(check_covariance(cov_valid, "PD matrix"))
  expect_match(
    check_covariance(cov_not_pos_def, "Non-PD matrix"),
    "must be positive definite"
  )
})

# ---------------------------------------------
# --- Tests for check_beta_gamma() ---
# ---------------------------------------------

test_that("check_beta_gamma validates non-empty matrices", {
  expect_error(check_beta_gamma(NULL, gamma_ns, q = 2, case = "ns"))
  expect_error(check_beta_gamma(beta_ns, NULL, q = 2, case = "ns"))
})

# ---------------------------------------------
# --- Tests for check_beta_support() ---
# ---------------------------------------------

test_that("check_beta_support checks support consistency", {
  expect_error(check_beta_support(beta_ns, support_valid))
  expect_error(check_beta_support(matrix(c(1, 2, 3, 4), nrow = 2), NULL))
})

# ---------------------------------------------
# --- Tests for numeric slot checks ---
# ---------------------------------------------

pos_val <- 5
neg_val <- -3
pos_int_val <- as.integer(7)
pos_non_int_val <- 7.0 # numeric but not integer

test_that("check_positive_slot works", {
  expect_null(check_positive_slot(pos_val, "slot"))
  expect_match(
    check_positive_slot(neg_val, "slot"),
    "must be strictly positive"
  )
})

test_that("check_positive_or_null_slot works", {
  expect_null(check_positive_or_null_slot(NULL, "slot"))
  expect_null(check_positive_or_null_slot(pos_val, "slot"))
  expect_match(
    check_positive_or_null_slot(neg_val, "slot"),
    "must be NULL or strictly positive"
  )
})

test_that("check_positive_integer_slot works", {
  expect_null(check_positive_integer_slot(pos_int_val, "slot"))
  expect_match(
    check_positive_integer_slot(pos_non_int_val, "slot"),
    "must be a strictly positive integer"
  )
  expect_match(
    check_positive_integer_slot(0, "slot"),
    "must be a strictly positive integer"
  )
})
