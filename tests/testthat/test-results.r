# =========================================================
# test-results.r
#
# Unit tests for summary_saemvs method (see
# methods-summary.r) and convergence_plot method
# (see methods-plots.r) using testthat
# =========================================================

beta_est <- matrix(c(0.1, 0.2), nrow = 2, ncol = 1)
gamma_est <- matrix(0.5, nrow = 1, ncol = 1)

saemvs_res <- new(
  "saemvsResults",
  criterion = "BIC",
  criterion_values = 0.5,
  thresholds = list(c(0.05)),
  beta_map = list(matrix(0.1, nrow = 2, ncol = 1)),
  mle_estimates = list(list(beta = beta_est, gamma = gamma_est)),
  support = list(matrix(TRUE, nrow = 2, ncol = 1)),
  unique_support = list(matrix(TRUE, nrow = 2, ncol = 1)),
  support_mapping = 1,
  spike_values_grid = 0.1,
  phi_fixed_idx = numeric(0),
  phi_to_select_idx = 1,
  forced_variables_idx = list(integer(0)),
  selected_variables_idx = list(1)
)


beta_to_select <- list(matrix(c(0.1, 0.2), nrow = 2, ncol = 1))
beta_not_to_select <- list(matrix(0.3, nrow = 1, ncol = 1))
gamma_to_select <- list(matrix(0.5, nrow = 1, ncol = 1))
gamma_not_to_select <- list(matrix(0.7, nrow = 1, ncol = 1))
sigma2 <- c(0.25, 0.2, 0.15)

saem_res <- methods::new(
  "saemResults",
  beta_to_select = beta_to_select,
  beta_not_to_select = beta_not_to_select,
  gamma_to_select = gamma_to_select,
  gamma_not_to_select = gamma_not_to_select,
  sigma2 = sigma2,
  phi_to_select_idx = 1,
  phi_not_to_select_idx = 2
)

# ---------------------------------------------
# --- Tests for summary ---
# ---------------------------------------------
test_that("summary_saemvs runs without error and produces output", {
  expect_output(summary(saemvs_res))
})


# ---------------------------------------------
# --- Tests for print ---
# ---------------------------------------------
# Fake phi draws: 3 iterations, 2 individus, 2 paramètres
phi_iter <- list(
  matrix(c(1, 2, 3, 4), nrow = 2,
         dimnames = list(c("id1", "id2"), c("phi1", "phi2"))),
  matrix(c(2, 3, 4, 5), nrow = 2,
         dimnames = list(c("id1", "id2"), c("phi1", "phi2"))),
  matrix(c(3, 4, 5, 6), nrow = 2, 
         dimnames = list(c("id1", "id2"), c("phi1", "phi2")))
)

saemvs_res@phi <- list(phi_iter)
saemvs_res@phi_names <- c("phi1", "phi2")

test_that("print(saemvsResults) runs without error", {
  expect_output(print(saemvs_res))
})

test_that("print(saemvsResults) displays key information", {
  output <- capture.output(print(saemvs_res))

  expect_true(any(grepl("Object of class 'saemvsResults'", output)))
  expect_true(any(grepl("Selection criterion:", output)))
  expect_true(any(grepl("Best model:", output)))
  expect_true(any(grepl("Number of selected covariates:", output)))
})


# ---------------------------------------------
# --- Tests for predict ---
# ---------------------------------------------

test_that("predict(saemvsResults) returns a matrix", {
  res <- predict(saemvs_res, k = 2)

  expect_true(is.matrix(res))
  expect_equal(dim(res), c(2, 2))
  expect_equal(colnames(res), c("phi1", "phi2"))
  expect_equal(rownames(res), c("id1", "id2"))
})

test_that("predict(saemvsResults) computes correct averages", {
  res <- predict(saemvs_res, k = 2)

  # Expected means
  # phi1: (2+3)/2 = 2.5, (3+4)/2 = 3.5
  # phi2: (4+5)/2 = 4.5, (5+6)/2 = 5.5
  expected <- matrix(
    c(2.5, 3.5, 4.5, 5.5),
    nrow = 2,
    dimnames = list(c("id1", "id2"), c("phi1", "phi2"))
  )

  expect_equal(res, expected)
})

test_that("predict_support() works for valid support index", {
  res <- predict_support(saemvs_res, support_idx = 1, k = 2)

  expect_true(is.matrix(res))
  expect_equal(dim(res), c(2, 2))
})

test_that("predict_support() fails with invalid support index", {
  expect_error(
    predict_support(saemvs_res, support_idx = 2),
    "support_idx must be between 1 and"
  )
})

test_that("predict() fails when k is larger than available iterations", {
  expect_error(
    predict(saemvs_res, k = 10),
    "Not enough iterations"
  )
})

# ---------------------------------------------
# --- Tests for summary(saemvsResults) ---
# ---------------------------------------------

mock_saemvsResults <- function() {
  new(
    "saemvsResults",
    criterion = "BIC",
    criterion_values = 0.5,
    thresholds = list(c(0.05)),
    beta_map = list(matrix(0.1, nrow = 2, ncol = 1)),
    mle_estimates = list(
      list(
        beta = matrix(c(0.1, 0.3), nrow = 2),
        gamma = diag(1)
      )
    ),
    support = list(matrix(TRUE, nrow = 2, ncol = 1)),
    unique_support = list(matrix(TRUE, nrow = 2, ncol = 1)),
    support_mapping = 1,
    spike_values_grid = 0.1,
    phi_fixed_idx = numeric(0),
    phi_to_select_idx = 1,
    forced_variables_idx = list(integer(0)),
    selected_variables_idx = list(1),
    phi_names = "phi1"
  )
}


test_that("summary(saemvsResults) runs without error", {
  saemvs_res <- mock_saemvsResults()
  expect_output(summary(saemvs_res))
})

test_that("summary(saemvsResults) displays expected sections", {
  saemvs_res <- mock_saemvsResults()
  output <- capture.output(summary(saemvs_res))

  expect_true(any(grepl("Best model", output)))
  expect_true(any(grepl("Selected Variables", output)))
  expect_true(any(grepl("Estimated Parameters", output)))
  expect_true(any(grepl("Coefficients", output)))
  expect_true(any(grepl("Covariance Matrix", output)))
})

test_that("summary() respects digits argument", {
  saemvs_res <- mock_saemvsResults()
  output <- capture.output(summary(saemvs_res, digits = 1))
  expect_true(any(grepl("0.1", output)))
})

test_that("summary_support() works for a valid support", {
  saemvs_res <- mock_saemvsResults()
  expect_output(summary_support(saemvs_res, support_idx = 1))
})

test_that("summary_support() fails with invalid support index", {
  saemvs_res <- mock_saemvsResults()
  expect_error(
    summary_support(saemvs_res, support_idx = 2),
    "support_idx must be between 1 and"
  )
})

test_that("summary_all() runs without error", {
  saemvs_res <- mock_saemvsResults()
  expect_output(summary_all(saemvs_res))
})

# ---------------------------------------------
# --- Tests for coef(saemvsResults) ---
# ---------------------------------------------

beta_est <- matrix(
  c(0.1, 0.2,   # mu
    0.3, 0.4),  # x1
  nrow = 2,
  ncol = 2,
  byrow = TRUE
)

saemvs_res@mle_estimates <- list(
  list(
    beta = beta_est,
    gamma = diag(2)
  )
)

saemvs_res@phi_names <- c("phi1", "phi2")

test_that("coef(saemvsResults) returns a labeled matrix", {
  beta <- coef(saemvs_res)

  expect_true(is.matrix(beta))
  expect_equal(colnames(beta), c("phi1", "phi2"))
  expect_true("μ" %in% rownames(beta))
})

test_that("coef() fails if beta and phi_names are inconsistent", {
  saemvs_bad <- saemvs_res
  saemvs_bad@phi_names <- "phi1"

  expect_error(
    coef(saemvs_bad),
    "dimnames"
  )
})

test_that("coef() returns correct beta estimates", {
  beta <- coef(saemvs_res)

  expect_equal(beta["μ", "phi1"], 0.1)
  expect_equal(beta[2, "phi1"], 0.3)
})

test_that("coef_support() works for a valid support", {
  beta <- coef_support(saemvs_res, support_idx = 1)

  expect_true(is.matrix(beta))
  expect_equal(dim(beta), c(2, 2))
})

test_that("coef_support() fails with invalid support index", {
  expect_error(
    coef_support(saemvs_res, support_idx = 2),
    "support_idx must be between 1 and"
  )
})
