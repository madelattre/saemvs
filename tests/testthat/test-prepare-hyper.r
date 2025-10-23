# =========================================================
# test-prepare_init.r
#
# Unit tests for prepare_init method (see
# methods-prepare-objects.r) using testthat
# =========================================================



y_list <- list(1:3, 1:3)
t_list <- list(1:3, 1:3)
x_candidates <- matrix(1, nrow = 2, ncol = 2)
x_forced <- matrix(1.1, nrow = 2, ncol = 1)

phi_dim <- 3
cov_scale <- diag(phi_dim)

test_that("prepare_hyper: MLE case returns NULL priors", {
  model_mle <- saemvsModel(
    g = function(phi, t) t %*% phi,
    phi_dim = phi_dim,
    phi_to_select_idx = integer(0)
  )

  data_base <- saemvsData(y_list, t_list, x_candidates, x_forced)
  data_proc <- prepare_data(data_base, model_mle)

  hyper <- saemvsHyperSlab(
    slab_parameter = 1000,
    cov_re_prior_scale = cov_scale
  )
  res <- prepare_hyper(hyper, data_proc, model_mle)

  expect_s4_class(res, "saemvsHyperSlab")
  expect_null(res@inclusion_prob_prior_a)
  expect_null(res@inclusion_prob_prior_b)
})

test_that("prepare_hyper: MAP case fills missing priors", {
  model_map <- saemvsModel(
    g = function(phi, t) t %*% phi,
    phi_dim = phi_dim,
    phi_to_select_idx = 1:2
  )

  data_base <- saemvsData(y_list, t_list, x_candidates, x_forced)
  data_proc <- prepare_data(data_base, model_map)

  hyper <- saemvsHyperSlab(
    slab_parameter = 1000,
    cov_re_prior_scale = cov_scale
  )
  res <- prepare_hyper(hyper, data_proc, model_map)

  expect_equal(
    res@inclusion_prob_prior_a,
    rep(1, length(model_map@phi_to_select_idx))
  )
  expect_equal(
    res@inclusion_prob_prior_b,
    rep(ncol(data_proc@x_candidates), length(model_map@phi_to_select_idx))
  )
})

test_that("prepare_hyper: MAP case keeps existing priors", {
  model_map <- saemvsModel(
    g = function(phi, t) t %*% phi,
    phi_dim = phi_dim,
    phi_to_select_idx = 1:2
  )

  data_base <- saemvsData(y_list, t_list, x_candidates, x_forced)
  data_proc <- prepare_data(data_base, model_map)

  hyper <- saemvsHyperSlab(
    slab_parameter = 1000,
    cov_re_prior_scale = cov_scale
  )
  hyper@inclusion_prob_prior_a <- c(2, 3)
  hyper@inclusion_prob_prior_b <- c(5, 5)

  res <- prepare_hyper(hyper, data_proc, model_map)
  expect_equal(res@inclusion_prob_prior_a, c(2, 3))
  expect_equal(res@inclusion_prob_prior_b, c(5, 5))
})

test_that("prepare_hyper: MAP case errors if priors wrong length", {
  model_map <- saemvsModel(
    g = function(phi, t) t %*% phi,
    phi_dim = phi_dim,
    phi_to_select_idx = 1:2
  )

  data_base <- saemvsData(y_list, t_list, x_candidates, x_forced)
  data_proc <- prepare_data(data_base, model_map)

  hyper <- saemvsHyperSlab(
    slab_parameter = 1000,
    cov_re_prior_scale = cov_scale
  )
  hyper@inclusion_prob_prior_a <- 1
  hyper@inclusion_prob_prior_b <- c(1, 2)

  expect_error(
    prepare_hyper(hyper, data_proc, model_map),
    "Length of 'inclusion_prob_prior_a'.*does not match"
  )
})
