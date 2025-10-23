# =========================================================
# test-methods-check.r
#
# Unit tests for methods-check.r using testthat
# =========================================================

# --- Prepare useful objects

y_series_list <- list(1:3, 4:6)
t_series_list <- list(1:3, 1:3)
x_candidates_ok <- matrix(1:4, nrow = 2)
x_candidates_bad <- matrix(1:2, nrow = 1) # mauvaise nrow
x_forced_ok <- matrix(c(1, 0), ncol = 1)
x_forced_bad <- matrix(letters[1:4], 2, 2)

data_base <- saemvsData(
  y = y_series_list,
  t = t_series_list,
  x_candidates = x_candidates_ok,
  x_forced = x_forced_ok
)

model_base <- saemvsModel(
  g = function(phi, t) phi[1] + phi[2] * t,
  phi_dim = 2,
  phi_to_select_idx = 1,
  phi_fixed_idx = 2,
  x_forced_support = matrix(c(1, 1), ncol = 2)
)

model_test <- saemvsModel(
  g = function(phi, t) phi,
  phi_dim = 3,
  phi_to_select_idx = c(1, 2),
  phi_fixed_idx = 3
)

data_test <- saemvsData(
  y = list(1:5, 1:5, 1:5),
  t = list(1:5, 1:5, 1:5),
  x_candidates = NULL,
  x_forced = NULL
)

init_correct <- saemvsInit(
  intercept = c(0.1, 0.2, 0.3),
  beta_forced = matrix(0, nrow = 2, ncol = 3),
  beta_candidates = matrix(0, nrow = 2, ncol = 3),
  cov_re = diag(c(1, 1, 0)),
  sigma2 = 1,
  default = FALSE
)

x_forced_support_test <- matrix(
  c(1, 0, 1, 0, 1, 0),
  nrow = 2, ncol = 3, byrow = TRUE
)


hyper_slab_base <- saemvsHyperSlab(
  slab_parameter = 1000,
  cov_re_prior_scale = diag(1),
  cov_re_prior_df = 3
)

hyper_base <- saemvsHyperSpikeAndSlab(
  spike_parameter = 0.1,
  hyper_slab = hyper_slab_base
)

tuning_base <- saemvsTuning(
  niter = 5,
  nburnin = 2,
  niter_mh = 2,
  kernel_mh = "random_walk",
  covariance_decay = 0.98,
  mh_proposal_scale = 1.0,
  spike_values_grid = c(0.05, 0.08),
  n_is_samples = 1000,
  seed = 123,
  nb_workers = 1
)

# ---------------------------------------------
# --- Tests for check_data() ---
# ---------------------------------------------

test_that("check_data error handling", {
  model <- model_base
  data <- data_base

  # phi_to_select_idx non-empty but x_candidates NULL
  model@phi_to_select_idx <- 1
  data@x_candidates <- NULL
  expect_error(
    check_data(data, model),
    "x_candidates.*missing"
  )

  # x_candidates present but wrong number of rows
  data <- data_base
  data@x_candidates <- x_candidates_bad
  expect_error(
    check_data(data, model),
    "must have .* rows"
  )

  # x_forced present, incompatible number of columns
  data <- data_base
  model <- model_base
  model@x_forced_support <- matrix(0, nrow = 2, ncol = 2)
  data@x_forced <- matrix(1:3, nrow = 1, ncol = 3)
  expect_error(
    check_data(data, model),
    "Number of columns in 'x_forced'"
  )

  # x_forced present, non-numeric
  data <- data_base
  data@x_forced <- matrix(letters[1:4], nrow = 2, ncol = 2)
  expect_error(
    check_data(data, model),
    "'x_forced' must be numeric"
  )


  # x_forced present, incompatible number of columns
  data <- data_base
  model <- model_base
  x_forced_support_bad <- matrix(0, nrow = 3, ncol = 2)
  model@x_forced_support <- x_forced_support_bad
  data@x_forced <- matrix(1:4, nrow = 2, ncol = 2)
  expect_error(
    check_data(data, model),
    "Number of columns in 'x_forced'"
  )

  # Everything correct
  model <- model_base
  data <- data_base
  expect_silent(check_data(data, model))
})

# ---------------------------------------------
# --- Tests for check_init() ---
# ---------------------------------------------

test_that("check_init errors", {
  # intercept length mismatch
  init1 <- init_correct
  init1@intercept <- c(0.1, 0.2)
  expect_error(
    check_init(init1, data_test, model_test),
    "Length of 'intercept'"
  )

  # beta_candidates wrong ncol
  init2 <- init_correct
  init2@beta_candidates <- matrix(0, nrow = 2, ncol = 2)
  expect_error(
    check_init(init2, data_test, model_test),
    "beta_candidates.*must be a matrix"
  )

  # beta_forced wrong ncol
  init3 <- init_correct
  init3@beta_forced <- matrix(0, nrow = 2, ncol = 2)
  expect_error(
    check_init(init3, data_test, model_test),
    "beta_forced.*must be a matrix"
  )

  # cov_re not square or wrong size
  init4 <- init_correct
  init4@cov_re <- matrix(0, nrow = 2, ncol = 3)
  expect_error(
    check_init(init4, data_test, model_test),
    "cov_re.*numeric matrix"
  )

  # cov_re diag=0 mismatch phi_fixed_idx
  init5 <- init_correct
  init5@cov_re <- diag(c(1, 0, 1))
  model_test@phi_fixed_idx <- 1
  expect_error(
    check_init(init5, data_test, model_test),
    "Mismatch between zero diagonals"
  )

  # sigma2 less than or equal to zero
  init6 <- init_correct
  init6@sigma2 <- -1
  model_test@phi_fixed_idx <- 3
  expect_error(
    check_init(init6, data_test, model_test),
    "sigma2.*strictly positive"
  )

  # beta_forced inconsistent with x_forced_support
  init7 <- init_correct
  init7@beta_forced <- matrix(1, nrow = 2, ncol = 3)
  model_test@x_forced_support <- x_forced_support_test
  expect_error(
    check_init(init7, data_test, model_test),
    "Inconsistent forced values"
  )

  # all correct
  init8 <- init_correct
  model_test@phi_fixed_idx <- 3
  model_test@x_forced_support <- NULL
  expect_silent(
    check_init(init8, data_test, model_test)
  )
})

# ---------------------------------------------
# --- Tests for check_hyper() ---
# ---------------------------------------------


test_that("check_hyper errors", {
  # inclusion_prob_prior_a wrong length
  hyper_a <- hyper_base
  hyper_a@inclusion_prob_prior_a <- c(0.5, 0.5)
  expect_error(
    check_hyper(hyper_a, model_base, tuning_base),
    "inclusion_prob_prior_a"
  )

  # inclusion_prob_prior_b wrong length
  hyper_b <- hyper_base
  hyper_b@inclusion_prob_prior_b <- c(0.5, 0.5)
  expect_error(
    check_hyper(hyper_b, model_base, tuning_base),
    "inclusion_prob_prior_b"
  )

  # cov_re_prior_scale wrong dimension
  hyper_cov <- hyper_base
  hyper_cov@cov_re_prior_scale <- diag(3)
  expect_error(
    check_hyper(hyper_cov, model_base, tuning_base),
    "cov_re_prior_scale"
  )

  # spike_values_grid greater than or equal to slab_parameter
  tuning_spike <- tuning_base
  tuning_spike@spike_values_grid <- c(0.05, 1500)
  expect_error(
    check_hyper(hyper_base, model_base, tuning_spike),
    "spike_values_grid"
  )
})

test_that("check_hyper succeeds for correct inputs", {
  expect_null(
    check_hyper(hyper_base, model_base, tuning_base)
  )
})
