# =========================================================
# test-calc-utils.r
#
# Unit tests for calc-utils.r using testthat
# =========================================================

# ---------------------------------------------
# --- Tests for threshold() ---
# ---------------------------------------------

nu1_valid <- 2
nu0_valid <- 0.5
alpha_valid <- 0.2
alpha_invalid <- -0.1 # out of 0-1 range

test_that("threshold returns correct numeric value", {
  expect_type(threshold(nu1_valid, nu0_valid, alpha_valid), "double")
  expect_gt(threshold(nu1_valid, nu0_valid, alpha_valid), 0)
})

# ---------------------------------------------
# --- Tests for p_star() ---
# ---------------------------------------------

beta_matrix <- matrix(c(0.5, 1.0, -0.5, 0), nrow = 2, ncol = 2)
alpha_vec <- c(0.2, 0.8)
nu0 <- 0.1
nu1 <- 1.0

test_that("p_star returns numeric matrix with entries in [0,1]", {
  ps <- p_star(beta_matrix, alpha_vec, nu0, nu1)
  expect_equal(dim(ps), dim(beta_matrix))
  expect_true(all(ps >= 0 & ps <= 1))
})

# ---------------------------------------------
# --- Tests for get_case() ---
# ---------------------------------------------

config_map_full <- list(
  method_type = "map",
  num_parameters_not_to_select = 0,
  fixed_parameters_indices = integer(0)
)
config_map_part_nofixed <- list(
  method_type = "map",
  num_parameters_not_to_select = 2,
  fixed_parameters_indices = integer(0)
)
config_map_part_fixed <- list(
  method_type = "map",
  num_parameters_not_to_select = 2,
  fixed_parameters_indices = c(1)
)
config_mle_nofixed <- list(
  method_type = "mle",
  num_parameters_not_to_select = 2,
  fixed_parameters_indices = integer(0)
)
config_mle_fixed <- list(
  method_type = "mle",
  num_parameters_not_to_select = 2,
  fixed_parameters_indices = c(1)
)

test_that("get_case returns correct string", {
  expect_equal(get_case(config_map_full), "map_full_select")
  expect_equal(get_case(config_map_part_nofixed), "map_part_select_nofixed")
  expect_equal(get_case(config_map_part_fixed), "map_part_select_fixed")
  expect_equal(get_case(config_mle_nofixed), "mle_nofixed")
  expect_equal(get_case(config_mle_fixed), "mle_fixed")
})

# ---------------------------------------------
# --- Tests for estim_phi_individuals() ---
# ---------------------------------------------

data_example <- saemvsData(
  y = list(c(1, 2, 3), c(2, 3, 4)),
  t = list(c(1, 2, 3), c(1, 2, 3))
)

model_example <- saemvsModel(
  phi_dim = 2,
  g = function(phi, t) phi[1] + phi[2] * t
)

init_example <- saemvsInit(
  intercept = c(1, 1)
)

test_that("estimate_phi_individuals returns numeric matrix of correct size", {
  res <- estimate_phi_individuals(data_example, model_example,
    init_example,
    maxit = 10
  )
  expect_equal(dim(res), c(
    length(data_example@y_series),
    model_example@phi_dim
  ))
  expect_true(is.numeric(res))
})

# ---------------------------------------------
# --- Tests for build_init_from_phi_lasso() ---
# ---------------------------------------------

est_indiv <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)

data_processed <- methods::new("saemvsProcessedData",
  y_series = list(c(1, 2, 3), c(2, 3, 4)),
  t_series = list(c(1, 2, 3), c(1, 2, 3)),
  x_phi_to_select = matrix(c(1, 2, 2, 1, 3, 2), nrow = 2),
  x_phi_not_to_select = matrix(c(1, 2, 3, 4), nrow = 2)
)

model_build <- saemvsModel(
  phi_dim = 2,
  g = function(phi, t) phi[1] + phi[2] * t,
  phi_to_select_idx = c(1, 2),
  x_forced_support = matrix(0, nrow = 1, ncol = 2)
)

init_build <- saemvsInit(
  intercept = c(1, 1),
  default = TRUE,
  cov_re = diag(2),
  sigma2 = 1
)

test_that("build_init_from_phi_lasso returns correct saemvsInit object", {
  init_obj <- build_init_from_phi_lasso(
    est_indiv, data_processed,
    model_build, init_build
  )
  expect_s4_class(init_obj, "saemvsInit")
  expect_equal(
    ncol(init_obj@beta_candidates),
    model_build@phi_to_select_idx %>% length()
  )
})
