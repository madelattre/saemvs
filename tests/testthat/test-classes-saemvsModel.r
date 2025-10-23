# =========================================================
# test-classes-saemvsModel.r
#
# Unit tests for saemvsModel class and saemvsModel constructor
# (see classes.r) using testthat
# =========================================================


g_valid <- function(phi, t) phi * t
g_invalid <- 3 # Not a function

x_forced_valid <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2)
x_forced_empty <- matrix(numeric(0), nrow = 0, ncol = 0)
x_forced_bad <- matrix(1:4, nrow = 2, ncol = 2) # phi_dim mismatch example

phi_to_select_ok <- c(1)
phi_fixed_ok <- c(2)
phi_to_select_bad <- c(3) # out of range
phi_fixed_overlap <- c(1) # overlaps with phi_to_select_ok


invalid_support <- matrix(c(0, 1, 2, 0), nrow = 2, ncol = 2)

test_that("saemvsModel: valid object creation works", {
  model <- saemvsModel(
    g = g_valid,
    phi_dim = 2,
    phi_to_select_idx = phi_to_select_ok,
    phi_fixed_idx = phi_fixed_ok,
    x_forced_support = x_forced_valid
  )
  expect_s4_class(model, "saemvsModel")
  expect_true(is.function(model@model_func))
  expect_equal(model@phi_dim, 2L)
  expect_equal(model@phi_to_select_idx, phi_to_select_ok)
  expect_equal(model@phi_fixed_idx, phi_fixed_ok)
  expect_equal(model@x_forced_support, x_forced_valid)
})

test_that("saemvsModel: invalid g (non-function) raises error", {
  expect_error(saemvsModel(g_invalid, phi_dim = 2), "function")
})

test_that("saemvsModel: negative or non-integer phi_dim raises error", {
  expect_error(saemvsModel(g_valid, phi_dim = -1))
  expect_error(saemvsModel(g_valid, phi_dim = 2.5))
})

test_that("saemvsModel: invalid index vectors raise error", {
  expect_error(saemvsModel(g_valid,
    phi_dim = 2,
    phi_to_select_idx = phi_to_select_bad
  ))
  expect_error(saemvsModel(
    g_valid,
    phi_dim = 2,
    phi_to_select_idx = phi_to_select_ok,
    phi_fixed_idx = phi_fixed_overlap
  ))
})

test_that("saemvsModel: invalid x_forced_support dimension raises error", {
  expect_error(saemvsModel(g_valid,
    phi_dim = 3,
    x_forced_support = x_forced_bad
  ))
})

test_that("saemvsModel: empty x_forced_support is accepted", {
  model <- saemvsModel(g_valid, phi_dim = 2, x_forced_support = x_forced_empty)
  expect_s4_class(model, "saemvsModel")
})

test_that("x_forced_support contains only 0 or 1", {
  expect_error(
    saemvsModel(
      g = g_valid,
      phi_dim = 2,
      x_forced_support = invalid_support
    ),
    "'x_forced_support' must only contain 0 or 1"
  )
})
