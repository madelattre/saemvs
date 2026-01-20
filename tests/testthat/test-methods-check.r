# =========================================================
# test-methods-check.r
#
# Unit tests for methods-check.r using testthat
# =========================================================

# --- Prepare useful objects

# y_series_list <- list(1:3, 4:6)
# t_series_list <- list(1:3, 1:3)
# x_candidates_ok <- matrix(1:4, nrow = 2)
# x_candidates_bad <- matrix(1:2, nrow = 1) # mauvaise nrow
# x_forced_ok <- matrix(c(1, 0), ncol = 1)
# x_forced_bad <- matrix(letters[1:4], 2, 2)

# data_base <- saemvsData(
#   y = y_series_list,
#   t = t_series_list,
#   x_candidates = x_candidates_ok,
#   x_forced = x_forced_ok
# )

data_base <- data.frame(
  id    = rep(1:2, each = 3),
  y     = c(1, 2, 3, 4,  5, 6),
  t     = rep(1:3, 2),
  cov1  = c(10, 10, 10, 20, 20, 20),
  cov2  = c(5, 5, 5, 7, 7, 7),
  forced_cov = c(0, 0, 0, 1, 1, 1)
)

formula_example <- y ~ . + group(id) + repeated(t) + forced_cov

# Construction de data_base
data_base <- saemvsData_from_df(
  formula = formula_example,
  data    = data_base
)

model_base <- saemvsModel(
  g = function(t, a, b) a + b * t,
  phi_to_select = c("a"),
  phi_fixed = c("b"),
  x_forced_support = list(
    a = c("forced_cov"),
    b = c("forced_cov")
  )
)

model_proc <- prepare_model(data_base, model_base)

init_correct <- saemvsInit(
  intercept = c(0.1, 0.2, 0.3),
  beta_forced = matrix(0, nrow = 2, ncol = 3),
  beta_candidates = matrix(0, nrow = 2, ncol = 3),
  cov_re = diag(c(1, 1, 0)),
  sigma2 = 1,
  default = FALSE
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

test_that("check_data_and_model error handling", {
  model <- model_base
  data <- data_base

  # phi_to_select non-empty but x_candidates NULL
  model@phi_to_select <- c("a")
  data@x_candidates <- NULL
  expect_error(
    check_data_and_model(data, model),
    "x_candidates.*missing"
  )

  # x_candidates present but wrong number of rows
  data <- data_base
  data@x_candidates <- matrix(1:2, nrow = 1) # mauvaise nrow
  expect_error(
    check_data_and_model(data, model),
    "must have .* rows"
  )

  # x_forced present but non-numeric
  data <- data_base
  data@x_forced <- matrix(letters[1:4], nrow = 2, ncol = 2)
  expect_error(
    check_data_and_model(data, model),
    "'x_forced' must be numeric"
  )

  # x_forced present but missing column names
  data <- data_base
  data@x_forced <- matrix(1:4, nrow = 2, ncol = 2)
  colnames(data@x_forced) <- NULL
  expect_error(
    check_data_and_model(data, model),
    "'x_forced' must have column names"
  )

  # x_forced present but missing covariates declared in model
  data <- data_base
  data@x_forced <- matrix(0, nrow = 2, ncol = 1)
  colnames(data@x_forced) <- "wrong_cov"
  model@x_forced_support <- list(a = "forced_cov")
  expect_error(
    check_data_and_model(data, model),
    "Some covariates were declared for forced inclusion.*not declared"
  )

  # x_forced present but extra covariates not used in model
  data <- data_base
  data@x_forced <- matrix(0, nrow = 2, ncol = 2)
  colnames(data@x_forced) <- c("forced_cov", "extra_cov")
  model@x_forced_support <- list(a = "forced_cov")
  expect_error(
    check_data_and_model(data, model),
    "Some covariates were declared as forced in the data.*not used"
  )

  # Everything correct
  data <- data_base
  model <- model_base
  expect_silent(check_data_and_model(data, model))
})

# ---------------------------------------------
# --- Tests for check_hyper() ---
# ---------------------------------------------


test_that("check_hyper errors", {
  model_proc <- prepare_model(data_base, model_base)
  # inclusion_prob_prior_a wrong length
  hyper_a <- hyper_base
  hyper_a@inclusion_prob_prior_a <- c(0.5, 0.5)
  expect_error(
    check_hyper(hyper_a, model_proc, tuning_base),
    "inclusion_prob_prior_a"
  )

  # inclusion_prob_prior_b wrong length
  hyper_b <- hyper_base
  hyper_b@inclusion_prob_prior_b <- c(0.5, 0.5)
  expect_error(
    check_hyper(hyper_b, model_proc, tuning_base),
    "inclusion_prob_prior_b"
  )

  # cov_re_prior_scale wrong dimension
  hyper_cov <- hyper_base
  hyper_cov@cov_re_prior_scale <- diag(3)
  expect_error(
    check_hyper(hyper_cov, model_proc, tuning_base),
    "cov_re_prior_scale"
  )

  # spike_values_grid greater than or equal to slab_parameter
  tuning_spike <- tuning_base
  tuning_spike@spike_values_grid <- c(0.05, 1500)
  expect_error(
    check_hyper(hyper_base, model_proc, tuning_spike),
    "spike_values_grid"
  )
})

test_that("check_hyper succeeds for correct inputs", {
  expect_null(
    check_hyper(hyper_base, model_proc, tuning_base)
  )
})
