# =========================================================
# test-prepare_init.r
#
# Unit tests for prepare_init method (see
# methods-prepare-objects.r) using testthat
# =========================================================


# ------------------------------ #
# --- Prepare useful objects --- #
# ------------------------------ #

initprep_model_simple <- saemvsModel(
  g = function(phi, t) phi[2] + phi[1] / (1 + exp(-(t - phi[3]))),
  phi_dim = 3,
  phi_to_select_idx = 2,
  phi_fixed_idx = c(1, 3),
  x_forced_support = matrix(
    c(1, 1, 1),
    nrow = 1
  )
)

initprep_data_simple <- saemvsData(
  y = list(1:5, 6:10, 11:15),
  t = list(1:5, 1:5, 1:5),
  x_candidates = matrix(seq(1, 6), nrow = 3, ncol = 2),
  x_forced = matrix(seq(7, 9), ncol = 1)
)

initprep_data_processed <- prepare_data(
  initprep_data_simple, initprep_model_simple
)

initprep_init_simple <- saemvsInit(
  intercept = rep(0.1, 3),
  beta_candidates = matrix(0.1, nrow = 2, ncol = 3),
  beta_forced = matrix(0.5, nrow = 1, ncol = 3),
  cov_re = diag(3),
  sigma2 = 1,
  default = FALSE
)

# --- No parameter to select ---
initprep_model_no_select <- saemvsModel(
  g = function(phi, t) phi[2] + phi[1] / (1 + exp(-(t - phi[3]))),
  phi_dim = 3,
  phi_to_select_idx = integer(0),
  phi_fixed_idx = 1:3,
  x_forced_support = matrix(1, nrow = 1, ncol = 3)
)

initprep_init_no_select <- saemvsInit(
  intercept = rep(0.1, 3),
  beta_candidates = matrix(0.2, nrow = 2, ncol = 3),
  beta_forced = matrix(0.3, nrow = 1, ncol = 3),
  cov_re = diag(0.5, 3),
  sigma2 = 1,
  default = FALSE
)

initprep_init_base <- saemvsInit(
  intercept = rep(0.1, initprep_model_simple@phi_dim),
  beta_candidates = matrix(0.2, nrow = 2, ncol = initprep_model_simple@phi_dim),
  beta_forced = matrix(0.3,
    nrow = nrow(initprep_model_simple@x_forced_support),
    ncol = initprep_model_simple@phi_dim
  ),
  cov_re = diag(0.5, initprep_model_simple@phi_dim),
  sigma2 = 1,
  default = FALSE
)


initprep_model_no_forced <- saemvsModel(
  g = function(phi, t) phi[2] + phi[1] / (1 + exp(-(t - phi[3]))),
  phi_dim = 3,
  phi_to_select_idx = 1:3,
  phi_fixed_idx = integer(0),
  x_forced_support = matrix(0, nrow = 1, ncol = 3)
)

initprep_init_no_forced <- saemvsInit(
  intercept = rep(0.1, 3),
  beta_candidates = matrix(0.2, nrow = 2, ncol = 3),
  beta_forced = matrix(0, nrow = 1, ncol = 3),
  cov_re = diag(0.5, 3),
  sigma2 = 1,
  default = FALSE
)

# ------------------ #
# --- Unit tests --- #
# ------------------ #

test_that("prepare_init nominal case works (minimal)", {
  result <- prepare_init(
    initprep_init_simple,
    initprep_model_simple,
    initprep_data_processed
  )

  expect_s4_class(result, "saemvsProcessedInit")
  expect_true(is.matrix(result@beta_to_select))
  expect_true(is.matrix(result@gamma_to_select))
  expect_true(is.numeric(result@inclusion_prob))
  expect_equal(
    length(result@inclusion_prob),
    length(initprep_model_simple@phi_to_select_idx)
  )
})

test_that("prepare_init output structure is consistent (minimal)", {
  result <- prepare_init(
    initprep_init_simple,
    initprep_model_simple,
    initprep_data_processed
  )

  expect_s4_class(result, "saemvsProcessedInit")
  expected_slots <- c(
    "beta_to_select", "beta_not_to_select",
    "gamma_to_select", "gamma_not_to_select",
    "sigma2", "inclusion_prob",
    "intercept", "beta_forced", "beta_candidates", "cov_re", "default"
  )
  expect_true(all(expected_slots %in% slotNames(result)))
})

test_that("prepare_init applies fixed-effect variance correction (minimal)", {
  model_fixed <- initprep_model_simple
  model_fixed@phi_fixed_idx <- 3
  init_small <- initprep_init_simple
  init_small@intercept[3] <- 0.01

  result <- prepare_init(init_small, model_fixed, initprep_data_processed)

  diag_vals <- diag(result@gamma_not_to_select)
  expect_true(all(diag_vals > 0))
})


test_that("prepare_init nominal case works", {
  result <- prepare_init(
    initprep_init_simple,
    initprep_model_simple,
    initprep_data_processed
  )

  expect_s4_class(result, "saemvsProcessedInit")
  expect_true(is.matrix(result@beta_to_select))
  expect_true(is.matrix(result@gamma_to_select))
  expect_true(is.numeric(result@inclusion_prob))
  expect_equal(
    length(result@inclusion_prob),
    length(initprep_model_simple@phi_to_select_idx)
  )
})

test_that("prepare_init handles no parameter to select", {
  data_proc <- prepare_data(initprep_data_simple, initprep_model_no_select)
  result <- prepare_init(
    initprep_init_no_select,
    initprep_model_no_select,
    data_proc
  )

  expect_null(result@beta_to_select)
  expect_null(result@gamma_to_select)
  expect_null(result@inclusion_prob)
  expect_true(is.matrix(result@beta_not_to_select))
})

test_that("prepare_init handles no forced covariates", {
  data_proc <- prepare_data(initprep_data_simple, initprep_model_no_forced)
  result <- prepare_init(
    initprep_init_no_forced,
    initprep_model_no_forced,
    data_proc
  )

  expect_true(is.matrix(result@beta_to_select))
  expect_true(is.matrix(result@gamma_to_select))
})

test_that("prepare_init handles all parameters forced", {
  data_proc <- prepare_data(initprep_data_simple, initprep_model_no_select)
  result <- prepare_init(
    initprep_init_no_select,
    initprep_model_no_select,
    data_proc
  )

  expect_null(result@beta_to_select)
  expect_true(is.matrix(result@beta_not_to_select))
})

test_that("prepare_init applies fixed-effect variance correction", {
  model_fixed <- initprep_model_simple
  model_fixed@phi_fixed_idx <- 2:3
  init_small <- initprep_init_base
  init_small@intercept[2:3] <- 0.01

  data_proc <- prepare_data(initprep_data_simple, model_fixed)
  result <- prepare_init(init_small, model_fixed, data_proc)

  diag_vals <- diag(result@gamma_not_to_select)
  expect_true(all(diag_vals > 0))
})

test_that("prepare_init output structure is consistent", {
  result <- prepare_init(
    initprep_init_base,
    initprep_model_simple,
    initprep_data_processed
  )

  expect_s4_class(result, "saemvsProcessedInit")
  expected_slots <- c(
    "beta_to_select", "beta_not_to_select",
    "gamma_to_select", "gamma_not_to_select",
    "sigma2", "inclusion_prob"
  )
  expect_true(all(expected_slots %in% slotNames(result)))
})
