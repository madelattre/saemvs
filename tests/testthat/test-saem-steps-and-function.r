# =========================================================
# test-saem-steps-and-function.r
#
# Unit tests for s-step.r, sa-steps.r, m-steps.r and
# saem-functions.r using testthat
# =========================================================


# Need to use one model with all parameters subject to selection
# (_full_select), another without parameters subject to selection
# (_no_select), and a third with both types of parameters to test
# the three variants of each step function


# ------------------------------ #
# --- Prepare useful objects --- #
# ------------------------------ #


# ---------------------------------------------
# --- Tests for s-step.r
# ---------------------------------------------

# --- Model with both selected and unselected parameters

initprep_model_simple <- saemvsModel(
  g = function(t, a, b, c) b + a / (1 + exp(-(t - c))),
  phi_to_select = c("b"),
  phi_fixed = c("a", "c"),
  x_forced_support = list(
    a = c("F1"),
    b = c("F1"),
    c = c("F1")
  )
)

initprep_data_simple <- saemvsData(
  y = list(1:5, 6:10, 11:15),
  t = list(1:5, 1:5, 1:5),
  x_candidates = matrix(seq(1, 6), nrow = 3, ncol = 2),
  x_forced = matrix(seq(7, 9), ncol = 1)
)

initprep_model_processed <- prepare_model(
  initprep_data_simple,
  initprep_model_simple
)

initprep_data_processed <- prepare_data(
  initprep_data_simple, initprep_model_processed
)

initprep_init_simple_bis <- saemvsInit(
  intercept = rep(0.1, 3),
  beta_candidates = matrix(0.1, nrow = 2, ncol = 3),
  beta_forced = matrix(0.5, nrow = 1, ncol = 3),
  cov_re = diag(c(0, 3, 0)),
  sigma2 = 1,
  default = FALSE
)

initprep_init_processed <- prepare_init(
  initprep_init_simple_bis,
  initprep_model_processed,
  initprep_data_processed
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

hyper_slab_saem <- saemvsHyperSlab(
  slab_parameter = 1000,
  cov_re_prior_scale = diag(1),
  cov_re_prior_df = 3
)

hyper_saem <- saemvsHyperSpikeAndSlab(
  spike_parameter = 0.1,
  hyper_slab = hyper_slab_saem
)

test_config <- make_config(
  initprep_data_processed,
  initprep_model_processed,
  tuning_base,
  initprep_init_processed,
  hyper_saem
)

test_state <- init_state(test_config)

# --- Model with selected parameters only

initprep_model_full_select <- saemvsModel(
  g = function(t, a, b, c) b + a / (1 + exp(-(t - c))),
  phi_to_select = c("a", "b", "c"),
  x_forced_support = list(
    a = c("F1"),
    b = c("F1"),
    c = c("F1")
  )
)

initprep_model_full_select_processed <- prepare_model( # nolint : object_length_linter
  initprep_data_simple,
  initprep_model_full_select
)

initprep_data_processed_full_select <- prepare_data( # nolint : object_length_linter
  initprep_data_simple, initprep_model_full_select_processed
)

initprep_init_full_select <- saemvsInit(
  intercept = rep(0.1, 3),
  beta_candidates = matrix(0.2, nrow = 2, ncol = 3),
  beta_forced = matrix(0.3, nrow = 1, ncol = 3),
  cov_re = diag(0.5, 3),
  sigma2 = 1,
  default = FALSE
)

hyper_slab_saem_full_select <- saemvsHyperSlab(
  slab_parameter = 1000,
  cov_re_prior_scale = diag(3),
  cov_re_prior_df = 3
)

hyper_slab_saem_full_select <- prepare_hyper(
  hyper_slab_saem_full_select,
  initprep_data_processed_full_select,
  initprep_model_full_select_processed
)

hyper_saem_full_select <- saemvsHyperSpikeAndSlab(
  spike_parameter = 0.1,
  hyper_slab = hyper_slab_saem_full_select
)

initprep_init_full_select_processed <- prepare_init( # nolint : object_length_linter
  initprep_init_full_select,
  initprep_model_full_select_processed,
  initprep_data_processed_full_select
)


test_config_full_select <- make_config(
  initprep_data_processed_full_select,
  initprep_model_full_select_processed,
  tuning_base,
  initprep_init_full_select_processed,
  hyper_saem_full_select
)

test_state_full_select <- init_state(test_config_full_select)


# --- Model with unselected parameters only

initprep_model_no_select <- saemvsModel(
  g = function(t, a, b, c) b + a / (1 + exp(-(t - c))),
  phi_to_select = c(),
  phi_fixed = c("a", "b", "c"),
  x_forced_support = list(
    a = c("F1"),
    b = c("F1"),
    c = c("F1")
  )
)

initprep_model_no_select_processed <- prepare_model( # nolint : object_length_linter
  initprep_data_simple,
  initprep_model_no_select
)

initprep_init_no_select <- saemvsInit(
  intercept = rep(0.1, 3),
  beta_candidates = matrix(0.2, nrow = 2, ncol = 3),
  beta_forced = matrix(0.3, nrow = 1, ncol = 3),
  cov_re = diag(0.5, 3),
  sigma2 = 1,
  default = FALSE
)


initprep_data_processed_no_select <- prepare_data( # nolint : object_length_linter
  initprep_data_simple, initprep_model_no_select_processed
)

initprep_init_no_select_processed <- prepare_init( # nolint : object_length_linter  
  initprep_init_no_select,
  initprep_model_no_select_processed,
  initprep_data_processed_no_select
)

test_config_no_select <- make_config(
  initprep_data_processed_no_select,
  initprep_model_no_select_processed,
  tuning_base,
  initprep_init_no_select_processed,
  hyper_saem
)

test_state_no_select <- init_state(test_config_no_select)


# ------------------ #
# --- Unit tests --- #
# ------------------ #

fake_backend <- list(
  metropolis_vector = function(...) {
    list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9))
  },
  g_vector = function(t, phi) {
    phi[2] + phi[1] / (1 + exp(-(t - phi[3])))
  }
)

# no distinction between selectable and non-selectable parameters
# at the metropolis-hastings stage
test_that("metropolis_s_step updates phi correctly", {
  state_new <- metropolis_s_step(test_config, 1, test_state, fake_backend)
  expect_true(is.matrix(state_new$phi[[2]]))
  expect_equal(dim(state_new$phi[[2]]), dim(test_state$phi[[1]]))
  expect_true(all(state_new$phi[[2]] == matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    3, 3
    ,
    byrow = TRUE
  )))
})

test_that("update_proposal_mh_all updates mean and variance proposals", {
  state_new <- update_proposal_mh_all(test_config, 1, test_state)
  expect_true(is.matrix(state_new$mprop_mh[[2]]))
  expect_true(is.matrix(state_new$vprop_mh[[2]]))
})

test_that("update_proposal_mh_to_select updates mean and variance proposals", {
  state_new <- update_proposal_mh_to_select(
    test_config_full_select, 1, test_state_full_select
  )
  expect_true(is.matrix(state_new$mprop_mh[[2]]))
  expect_true(is.matrix(state_new$vprop_mh[[2]]))
})

test_that(
  "update_proposal_mh_not_to_select updates mean and variance proposals",
  {
    state_new <- update_proposal_mh_not_to_select(
      test_config_no_select, 1, test_state_no_select
    )
    expect_true(is.matrix(state_new$mprop_mh[[2]]))
    expect_true(is.matrix(state_new$vprop_mh[[2]]))
  }
)


# ---------------------------------------------
# --- Tests for sa-steps.r
# --- verifies that sufficient statistics are of the correct type and have
# --- the correct dimensions (i.e., the same as in the previous iteration)
# ---------------------------------------------

test_that("sa_step_to_select updates sufficient statistics correctly", {
  state_new <- sa_step_to_select(
    test_config_full_select, 1, test_state_full_select
  )
  expect_true(is.matrix(state_new$s2_to_select[[2]]))
  expect_true(
    all(
      dim(state_new$s2_to_select[[2]]) ==
        dim(test_state_full_select$s2_to_select[[1]])
    )
  )
  expect_true(is.matrix(state_new$s3_to_select[[2]]))
  expect_true(
    all(
      dim(state_new$s3_to_select[[2]]) ==
        dim(test_state_full_select$s3_to_select[[1]])
    )
  )
  expect_type(state_new$s1, "double")
})

test_that("sa_step_not_to_select updates sufficient statistics correctly", {
  state_new <- sa_step_not_to_select(
    test_config_no_select, 1, test_state_no_select
  )
  expect_true(is.matrix(state_new$s2_not_to_select[[2]]))
  expect_true(
    all(
      dim(state_new$s2_not_to_select[[2]]) ==
        dim(test_state_no_select$s2_not_to_select[[1]])
    )
  )
  expect_true(is.matrix(state_new$s3_not_to_select[[2]]))
  expect_true(
    all(
      dim(state_new$s3_not_to_select[[2]]) ==
        dim(test_state_no_select$s3_not_to_select[[1]])
    )
  )
  expect_type(state_new$s1, "double")
})


test_that("sa_step_all updates sufficient statistics correctly", {

  state_new <- sa_step_all(test_config, 1, test_state, fake_backend)

  expect_type(state_new$s1, "double")
  expect_true(is.matrix(state_new$s2_to_select[[2]]))
  expect_true(is.matrix(state_new$s2_not_to_select[[2]]))
  expect_true(is.matrix(state_new$s3_to_select[[2]]))
  expect_true(is.matrix(state_new$s3_not_to_select[[2]]))
  expect_true(
    all(
      dim(state_new$s2_to_select[[2]]) ==
        dim(test_state$s2_to_select[[1]])
    )
  )
  expect_true(
    all(
      dim(state_new$s2_not_to_select[[2]]) ==
        dim(test_state$s2_not_to_select[[1]])
    )
  )
  expect_true(
    all(
      dim(state_new$s3_to_select[[2]]) ==
        dim(test_state$s3_to_select[[1]])
    )
  )
  expect_true(
    all(
      dim(state_new$s3_not_to_select[[2]]) ==
        dim(test_state$s3_not_to_select[[1]])
    )
  )
})


# ---------------------------------------------
# --- Tests for m-steps.r
# --- verifies that updated parameters are of the correct type and have
# --- the correct dimensions (i.e., the same as in the previous iteration)
# ---------------------------------------------

test_that("m_step_not_to_select updates beta and gamma correctly", {
  mockery::stub(m_step_not_to_select, "mat_inv", solve)
  mockery::stub(m_step_not_to_select, "solve_linear_syst", solve)

  state_new <- m_step_not_to_select(
    test_config_no_select, 1, test_state_no_select
  )
  expect_true(is.matrix(state_new$beta_not_to_select[[2]]))
  expect_true(is.matrix(state_new$gamma_not_to_select[[2]]))
  expect_true(
    all(
      dim(state_new$beta_not_to_select[[2]]) ==
        dim(test_state_no_select$beta_not_to_select[[1]])
    )
  )
  expect_true(
    all(
      dim(state_new$gamma_not_to_select[[2]]) ==
        dim(test_state_no_select$gamma_not_to_select[[1]])
    )
  )
})

test_that("m_step_to_select updates beta, gamma, and alpha correctly", {
  mockery::stub(m_step_to_select, "mat_inv", solve)
  mockery::stub(m_step_to_select, "solve_linear_syst", solve)

  fake_kronecker_gamma_diag_mult <- function(gamma, d_diag, p) {
    q <- nrow(gamma)
    dim <- p * q
    res <- matrix(0, nrow = dim, ncol = dim)

    for (i in seq_len(q)) {
      for (j in seq_len(q)) {
        gij <- gamma[i, j]
        if (gij == 0) next
        row_offset <- (i - 1) * p
        col_offset <- (j - 1) * p
        for (k in seq_len(p)) {
          row <- row_offset + k
          col <- col_offset + k
          res[row, col] <- gij * d_diag[row]
        }
      }
    }
    res
  }

  mockery::stub(
    where = m_step_to_select,
    what  = "kronecker_gamma_diag_mult",
    how   = fake_kronecker_gamma_diag_mult
  )

  state_new <- m_step_to_select(
    test_config_full_select, 1, test_state_full_select
  )
  expect_true(is.matrix(state_new$beta_to_select[[2]]))
  expect_true(is.matrix(state_new$gamma_to_select[[2]]))
  expect_type(state_new$alpha[[2]], "double")
})


# ---------------------------------------------
# --- Tests for perform_saem_iteration()
# --- saem-functions.r
# --- in the no parameters to select case only
# --- just check that the output is a list with expected entries
# ---------------------------------------------

test_that("perform_saem_iteration runs one iteration end-to-end", {
  result <- perform_saem_iteration(
    1, test_config_no_select, test_state_no_select,
    sa_func = sa_step_not_to_select,
    m_func = m_step_not_to_select,
    mh_update_func = update_proposal_mh_not_to_select,
    backend = fake_backend
  )

  expect_true(is.list(result))
  expect_true("phi" %in% names(result))
  expect_true("beta_not_to_select" %in% names(result))
  expect_true("gamma_not_to_select" %in% names(result))
})



# ---------------------------------------------
# --- Tests for run_saem_full()
# --- saem-functions.r
# --- in the no parameters to select case only
# --- just check that the output is a list with
# --- expected entries after multiple iterations
# ---------------------------------------------

test_that("run_saem_full executes multiple iterations without error", {

  result <- run_saem_full(
    test_config_no_select, test_state_no_select,
    sa_func = sa_step_not_to_select,
    m_func = m_step_not_to_select,
    mh_update_func = update_proposal_mh_not_to_select,
    backend = fake_backend
  )

  expect_true(is.list(result))
  expect_true(length(result$phi) >= test_config$num_iterations)
  expect_true(length(result$beta_not_to_select) >= test_config$num_iterations)
  expect_true(length(result$gamma_not_to_select) >= test_config$num_iterations)
})
