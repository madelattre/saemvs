# =========================================================
# test-prepare_data.r
#
# Unit tests for prepare_data method (see
# methods-prepare-objects.r) using testthat
# =========================================================

y_list_prep <- list(1:5, 1:5, 1:5)
t_list_prep <- list(1:5, 1:5, 1:5)

x_candidates_mat_prep <- matrix(1:15, nrow = 3)
x_forced_mat_prep <- matrix(16:21, nrow = 3, ncol = 2)


model_base_prep <- saemvsModel(
  g = function(phi, t) t %*% phi,
  phi_dim = 5,
  phi_to_select_idx = 1:3,
  phi_fixed_idx = 4:5,
  x_forced_support = matrix(
    c(1, 0, 1, 1, 0, 0, 1, 0, 1, 0),
    ncol = 5, byrow = TRUE
  )
)
data_base_prep <- saemvsData(
  y = y_list_prep,
  t = t_list_prep,
  x_candidates = x_candidates_mat_prep,
  x_forced = x_forced_mat_prep
)

model_all_forced_prep <- saemvsModel(
  g = function(phi, t) t %*% phi,
  phi_dim = 3,
  phi_to_select_idx = integer(0),
  phi_fixed_idx = 1:3,
  x_forced_support = matrix(1, nrow = 1, ncol = 3)
)

data_no_forced_prep <- saemvsData(
  y = list(1:3),
  t = list(1:3),
  x_candidates = matrix(1:3, nrow = 1),
  x_forced = NULL
)

model_minimal_prep <- saemvsModel(
  g = function(phi, t) t %*% phi,
  phi_dim = 1,
  phi_to_select_idx = 1,
  phi_fixed_idx = integer(0),
  x_forced_support = matrix(numeric(0), nrow = 0, ncol = 0)
)

test_that("prepare_data nominal case", {
  result <- prepare_data(data_base_prep, model_base_prep)

  expect_s4_class(result, "saemvsProcessedData")
  expect_length(result@y_series, length(data_base_prep@y_series))
  expect_true(is.matrix(result@x_forced))
  expect_true(is.matrix(result@x_phi_to_select))
  expect_equal(nrow(result@x_forced), nrow(data_base_prep@x_forced))
})

test_that("prepare_data handles active forced covariates correctly", {
  result <- prepare_data(data_base_prep, model_base_prep)
  expect_true(ncol(result@x_forced) > 0)
  expect_true(all(dim(result@x_forced) <= dim(data_base_prep@x_forced)))
})

test_that("prepare_data handles all parameters forced", {
  result <- prepare_data(data_base_prep, model_all_forced_prep)
  expect_true(
    is.null(result@x_phi_to_select) || ncol(result@x_phi_to_select) == 0
  )
  expect_true(is.matrix(result@x_forced))
})

test_that("prepare_data handles absence of forced covariates", {
  result <- prepare_data(data_no_forced_prep, model_minimal_prep)
  expect_true(
    is.null(result@x_forced) || ncol(result@x_forced) == 0
  )
  expect_true(is.matrix(result@x_phi_to_select))
})


test_that("prepare_data output structure is consistent", {
  result <- prepare_data(data_base_prep, model_base_prep)
  expect_equal(length(result@y_series), length(data_base_prep@y_series))
  expect_equal(nrow(result@x_candidates), nrow(data_base_prep@x_candidates))
  expect_true(all(sapply(result@y_series, is.numeric)))
})
