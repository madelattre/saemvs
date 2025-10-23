# =========================================================
# test-classes-saemvsHyperSlab.r
#
# Unit tests for saemvsHyperSlab class and saemvsHyperSlab
# constructor (see classes.r) using testthat
# =========================================================

# Valid matrices
cov_scale_2 <- diag(2)
cov_scale_3 <- diag(3)

# Invalid matrices
cov_non_square <- matrix(1:6, 2, 3)
cov_non_numeric <- matrix(letters[1:4], 2, 2)

# Valid slab parameters
slab_default <- 12000
slab_custom <- 5000

# Invalid slab parameters
slab_negative <- -100
slab_non_numeric <- "abc"

# Valid degrees of freedom
df_valid <- 3

# Invalid degrees of freedom
df_zero <- 0
df_non_numeric <- "abc"


test_that("Valid creation with defaults", {
  obj <- saemvsHyperSlab(
    slab_parameter = slab_default,
    cov_re_prior_scale = cov_scale_2,
    cov_re_prior_df = df_valid
  )
  expect_s4_class(obj, "saemvsHyperSlab")
  expect_equal(obj@slab_parameter, slab_default)
})

test_that("Custom positive slab works", {
  obj <- saemvsHyperSlab(
    slab_parameter = slab_custom,
    cov_re_prior_scale = cov_scale_2,
    cov_re_prior_df = df_valid
  )
  expect_equal(obj@slab_parameter, slab_custom)
})

test_that("Non-positive slab raises error", {
  expect_error(
    saemvsHyperSlab(
      slab_parameter = slab_negative,
      cov_re_prior_scale = cov_scale_2,
      cov_re_prior_df = df_valid
    ),
    "slab_parameter must be NULL or strictly positive"
  )
})

test_that("Non-numeric slab raises error", {
  expect_error(
    saemvsHyperSlab(
      slab_parameter = slab_non_numeric,
      cov_re_prior_scale = cov_scale_2,
      cov_re_prior_df = df_valid
    )
  )
})

test_that("Invalid cov_re_prior_scale raises error", {
  expect_error(
    saemvsHyperSlab(
      slab_parameter = slab_default,
      cov_re_prior_scale = cov_non_square,
      cov_re_prior_df = df_valid
    )
  )
  expect_error(
    saemvsHyperSlab(
      slab_parameter = slab_default,
      cov_re_prior_scale = cov_non_numeric,
      cov_re_prior_df = df_valid
    )
  )
})

test_that("Non-positive cov_re_prior_df raises error", {
  expect_error(
    saemvsHyperSlab(
      slab_parameter = slab_default,
      cov_re_prior_scale = cov_scale_2,
      cov_re_prior_df = df_zero
    )
  )
  expect_error(
    saemvsHyperSlab(
      slab_parameter = slab_default,
      cov_re_prior_scale = cov_scale_2,
      cov_re_prior_df = df_non_numeric
    )
  )
})

test_that("NULL cov_re_prior_scale is allowed", {
  obj <- saemvsHyperSlab(
    slab_parameter = slab_default,
    cov_re_prior_scale = NULL,
    cov_re_prior_df = df_valid
  )
  expect_s4_class(obj, "saemvsHyperSlab")
})
