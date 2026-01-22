# =========================================================
# test-prepare_model.r
#
# Unit tests for prepare_model method (see
# methods-prepare-objects.r) using testthat
# =========================================================


model_func <- function(t, phi) phi

model <- new(
  "saemvsModel",
  model_func = model_func,
  phi_names = c("phi1", "phi2"),
  phi_to_select = "phi1",
  phi_fixed = "phi2",
  x_forced_support = list(phi1 = "z1")
)


data_no_forced <- new(
  "saemvsData",
  y_series = list(c(1, 2)),
  t_series = list(c(0, 1)),
  x_candidates = matrix(1, nrow = 1, ncol = 1),
  x_forced = NULL
)

test_that("prepare_model dispatches and returns saemvsProcessedModel", {
  res <- prepare_model(data_no_forced, model)

  expect_s4_class(res, "saemvsProcessedModel")
})

test_that("prepare_model handles NULL x_forced correctly", {
  res <- prepare_model(data_no_forced, model)

  expect_equal(res@phi_dim, 2)
  expect_equal(res@phi_to_select_idx, 1)
  expect_equal(res@phi_fixed_idx, 2)

  # matrice vide mais bien dimensionnée
  expect_true(is.matrix(res@x_forced_support))
  expect_equal(dim(res@x_forced_support), c(0, 2))
  expect_equal(colnames(res@x_forced_support), c("phi1", "phi2"))
})


data_forced <- new(
  "saemvsData",
  y_series = list(c(1, 2)),
  t_series = list(c(0, 1)),
  x_candidates = matrix(3, nrow = 1, ncol = 1),
  x_forced = matrix(
    c(1, 2),
    nrow = 1,
    dimnames = list(NULL, c("z1", "z2"))
  )
)

model_no_support <- model
model_no_support@x_forced_support <- NULL

test_that("prepare_model creates zero forced support when none specified", {
  res <- prepare_model(data_forced, model_no_support)

  expect_equal(dim(res@x_forced_support), c(2, 2))
  expect_equal(rownames(res@x_forced_support), c("z1", "z2"))
  expect_equal(colnames(res@x_forced_support), c("phi1", "phi2"))

  expect_true(all(res@x_forced_support == 0L))
})

test_that("prepare_model fails if x_forced has no column names", {
  data_bad <- data_forced
  colnames(data_bad@x_forced) <- NULL

  expect_error(
    prepare_model(data_bad, model),
    "Forced covariates in data must have column names"
  )
})

test_that("prepare_model correctly computes phi indices", {
  res <- prepare_model(data_no_forced, model)

  expect_equal(res@phi_to_select_idx, 1)
  expect_equal(res@phi_fixed_idx, 2)
})

test_that("prepare_model handles empty phi_to_select", {
  model2 <- model
  model2@phi_to_select <- character(0)

  res <- prepare_model(data_no_forced, model2)

  expect_null(res@phi_to_select_idx)
})


test_that("prepare_model fills forced support matrix correctly", {
  res <- prepare_model(data_forced, model)

  expect_equal(res@x_forced_support["z1", "phi1"], 1L)
  expect_equal(res@x_forced_support["z2", "phi1"], 0L)
  expect_equal(
    res@x_forced_support[, "phi2"],
    setNames(c(0L, 0L), c("z1", "z2"))
  )
})
