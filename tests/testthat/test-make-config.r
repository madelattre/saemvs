# =========================================================
# test-make-config.r
#
# Unit tests for the make_config method (see
# methods-prepare-objects.r) using testthat
# =========================================================


# --- saemvsData ---
y_series_list <- list(1:3, 4:6)
t_series_list <- list(1:3, 1:3)
x_candidates_ok <- matrix(1:4, nrow = 2)
x_forced_ok <- matrix(c(1, 0), ncol = 1)

data_base <- saemvsData(
  y = y_series_list,
  t = t_series_list,
  x_candidates = x_candidates_ok,
  x_forced = x_forced_ok
)

# --- saemvsModel ---
model_base <- saemvsModel(
  g = function(phi, t) phi[1] + phi[2] * t,
  phi_dim = 2,
  phi_to_select_idx = 1,
  phi_fixed_idx = 2,
  x_forced_support = matrix(c(1, 1), ncol = 2)
)

# --- saemvsInit ---
init_base <- saemvsInit(
  intercept = c(1, 2),
  beta_forced = matrix(c(1, 1), ncol = 2),
  beta_candidates = x_candidates_ok,
  cov_re = diag(2),
  sigma2 = 1
)

# Base slab object
hyper_slab_base <- saemvsHyperSlab(
  slab_parameter = 1000,
  cov_re_prior_scale = diag(2),
  cov_re_prior_df = 3
)

test_that(
  "make_config returns a correctly structured list with minimal safe example",
  {
    data_proc <- prepare_data(data_base, model_base)
    init_proc <- prepare_init(init_base, model_base, data_proc)

    cov_scale <- diag(model_base@phi_dim)

    tuning_algo_base <- saemvsTuning(
      spike_values_grid = c(0.01, 0.1, 1)
    )

    hyper_base <- saemvsHyperSpikeAndSlab(
      spike_parameter = 0.1,
      hyper_slab = hyper_slab_base
    )

    cfg <- make_config(
      data = data_proc,
      model = model_base,
      tuning_algo = tuning_algo_base,
      init = init_proc,
      hyperparam = hyper_base
    )

    expect_type(cfg, "list")
    expect_equal(cfg$total_parameters, model_base@phi_dim)
    expect_equal(
      length(cfg$parameters_to_select_indices),
      length(model_base@phi_to_select_idx)
    )
    expect_equal(cfg$method_type, "map")
  }
)
