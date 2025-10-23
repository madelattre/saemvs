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

saem_res <- new(
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
# --- Tests for summary_saemvs() ---
# ---------------------------------------------
test_that("summary_saemvs runs without error and produces output", {
  expect_output(summary_saemvs(saemvs_res))
})

# ---------------------------------------------
# --- Tests for convergence_plot() ---
# ---------------------------------------------
test_that("convergence_plot runs without error for sigma2", {
  expect_silent(convergence_plot(saem_res, component = "sigma2"))
})

test_that("convergence_plot runs without error for beta_to_select", {
  expect_silent(suppressMessages(convergence_plot(saem_res,
    component = "coef_phi_sel"
  )))
})
