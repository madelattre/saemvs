# =========================================================
# test-classes-saemvsModel.r
#
# Unit tests for saemvsModel class and saemvsModel constructor
# (see classes.r) using testthat
# =========================================================


# ----- valid and invalid model functions -----
g_valid <- function(t, phi1, phi2) phi1 * t + phi2
g_invalid <- 3 # Not a function
g_bad_return <- function(t, phi1, phi2) c(1, 2) # returns vector instead of scalar # nolint : line_length_linter.

# ----- x_forced_support examples -----
x_forced_valid <- list(phi1 = c("cov1"), phi2 = c("cov2"))
x_forced_empty <- list()
x_forced_not_named <- list("cov1", "cov2") # unnamed
x_forced_invalid_type <- list(phi1 = 1:2)   # not character vector
x_forced_unknown_name <- list(phi3 = c("cov1")) # phi3 not in phi_names

# ----- phi_to_select / phi_fixed examples -----
phi_to_select_ok <- "phi1"
phi_fixed_ok <- "phi2"
phi_to_select_bad <- "phi3" # unknown name
phi_fixed_overlap <- "phi1" # overlaps with phi_to_select_ok

# ------------------------------
test_that("saemvsModel: valid object creation works", {
  model <- saemvsModel(
    g = g_valid,
    phi_to_select = phi_to_select_ok,
    phi_fixed = phi_fixed_ok,
    x_forced_support = x_forced_valid
  )
  expect_s4_class(model, "saemvsModel")
  expect_true(is.function(model@model_func))
  expect_equal(model@phi_to_select, phi_to_select_ok)
  expect_equal(model@phi_fixed, phi_fixed_ok)
  expect_equal(model@x_forced_support, x_forced_valid)
  expect_equal(model@phi_names, c("phi1", "phi2"))
})

test_that("saemvsModel: invalid g (non-function) raises error", {
  expect_error(saemvsModel(g_invalid), "'g' must be a function")
})

test_that("saemvsModel: model_func must return numeric scalar", {
  expect_error(saemvsModel(g_bad_return), "must return a numeric scalar")
})

test_that("saemvsModel: unknown phi_to_select or phi_fixed names raise error", {
  expect_error(
    saemvsModel(
      g = g_valid,
      phi_to_select = phi_to_select_bad
    ),
    regexp = "phi_to_select"
  )
  # expect_error(
  #   saemvsModel(
  #     g = g_valid,
  #     phi_to_select = phi_to_select_ok,
  #     phi_fixed = phi_fixed_overlap
  #   ),
  #   regexp = "Unknown parameter(s) in 'phi_fixed'"
  # )
})

test_that("saemvsModel: x_forced_support validation", {
  # empty is accepted
  model <- saemvsModel(g = g_valid, x_forced_support = x_forced_empty)
  expect_s4_class(model, "saemvsModel")

  # not named
  expect_error(
    saemvsModel(g = g_valid, x_forced_support = x_forced_not_named),
    "'x_forced_support' must be a named list"
  )

  # element not character
  expect_error(
    saemvsModel(g = g_valid, x_forced_support = x_forced_invalid_type),
    "must be a character vector of covariate names"
  )

  # unknown names
  expect_error(
    saemvsModel(g = g_valid, x_forced_support = x_forced_unknown_name),
    "Unknown parameter\\(s\\) in 'x_forced_support'"
  )
})

test_that("saemvsModel: works without phi_to_select or phi_fixed", {
  model <- saemvsModel(g = g_valid)
  expect_s4_class(model, "saemvsModel")
  expect_null(model@phi_to_select)
  expect_null(model@phi_fixed)
})
