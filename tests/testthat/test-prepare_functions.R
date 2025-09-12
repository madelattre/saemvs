# test_that("prepare_data centers and scales candidate covariates", {
#   data <- saemvsData(
#     y = list(rnorm(10)),
#     t = list(seq_len(10)),
#     x_candidates = matrix(rnorm(10), nrow = 1, ncol = 10)
#   )
  
#   model <- saemvsModel(
#     g = function(phi, t) phi[1] * t,
#     phi_dim = 1,
#     phi_to_select_idx = 1L
#   )

#   prepared <- prepare_data(data, model)

#   expect_equal(mean(prepared@x_candidates), 0, tolerance = 1e-8)
#   expect_equal(sd(prepared@x_candidates), 1, tolerance = 1e-8)
# })


test_that("prepare_init fills missing initial values", {
  model <- saemvsModel(
    g = function(phi, t) phi[1] * t,
    phi_dim = 2,
    phi_to_select_idx = 1L
  )

  # Construire un init minimal compatible avec le nouveau saemvsInit
  intercept_vector <- rep(0, model@phi_dim)
  init <- saemvsInit(
    intercept = intercept_vector,
    beta_forced = NULL,
    beta_candidates = NULL,
    cov_re = diag(model@phi_dim),
    sigma2 = 1
  )

  prepared <- prepare_init(init, model)

  expect_true(is.numeric(prepared@beta_not_to_select))  
  expect_true(is.numeric(prepared@gamma_not_to_select))  
  expect_true(is.numeric(prepared@sigma2))
})


test_that("prepare_hyper completes missing fields for MAP models", {
  model <- saemvsModel(
    g = function(phi, t) phi[1] * t,
    phi_dim = 2,
    phi_to_select_idx = 1L
  )

  data <- saemvsData(
    y = list(rnorm(10)),
    t = list(seq_len(10)),
    x_candidates = matrix(rnorm(1), ncol = 1)
  )
  
  data_prepared <- prepare_data(data, model)

  # Utiliser le helper corrigé pour créer un hyper objet valide
  hyper <- create_test_hyper_slab(n_selected = length(model@phi_to_select_idx),
                             n_candidates = 1)

  prepared <- prepare_hyper(hyper, data_prepared, model)

  expect_equal(
    length(prepared@inclusion_prob_prior_a),
    length(model@phi_to_select_idx)
  )
  expect_equal(
    nrow(prepared@cov_re_prior_scale),
    length(model@phi_to_select_idx)
  )
  expect_equal(
    ncol(prepared@cov_re_prior_scale),
    length(model@phi_to_select_idx)
  )
})


# test_that("prepare_hyper is skipped for MLE models", {
#   model <- saemvsModel(
#     g = function(phi, t) phi[1] * t,
#     phi_dim = 1,
#     phi_to_select_idx = integer(0)
#   )

#   data <- saemvsData(
#     y = list(rnorm(10)),
#     t = list(seq_len(10)),
#     x_candidates = matrix(rnorm(10), ncol = 10, nrow = 1)
#   )

#   hyper <- create_test_hyper(n_selected = 0, n_candidates = 1)

#   expect_error(
#     prepare_hyper(hyper, data, model),
#     regexp = "not.*intended|phi_to_select_idx"
#   )
# })
