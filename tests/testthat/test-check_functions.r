test_that("check_data catches missing x_candidates", {
  data <- create_test_data()
  model <- create_test_model(phi_to_select_idx = 1:2)
  data@x_candidates <- NULL
  expect_error(check_data(data, model), "x_candidates")
})

test_that("check_data catches wrong number of rows in x_candidates", {
  data <- create_test_data()
  model <- create_test_model(phi_to_select_idx = 1:2)
  data@x_candidates <- matrix(rnorm(10), nrow = 5, ncol = 2)
  expect_error(check_data(data, model), "'x_candidates' must have")
})

test_that("check_init catches incorrect intercept length", {
  data <- create_test_data()
  model <- create_test_model(phi_dim = 3)
  init <- create_test_init(model)
  init@intercept <- c(0, 0)  # Incorrect length
  expect_error(check_init(init, data, model), "intercept")
})

test_that("check_init catches sigma2 <= 0", {
  data <- create_test_data()
  model <- create_test_model()
  init <- create_test_init(model)
  init@sigma2 <- -1
  expect_error(check_init(init, data, model), "sigma2")
})

# test_that("check_hyper catches wrong inclusion_prob_prior length", {
#   data <- create_test_data()
#   model <- create_test_model(phi_to_select_idx = 1:2)
#   tuning <- methods::new("saemvsTuning", niter = 10, nburnin = 5, niter_mh = 1,
#                          kernel_mh = "random_walk", mh_proposal_scale = 0.1,
#                          step = rep(1,10), covariance_decay = 0.95, spike_values_grid = c(0.1, 0.2))
#   hyper <- create_test_hyper(n_selected = 3, n_candidates = 2)

#   expect_error(check_hyper(hyper, model, tuning), "inclusion_prob_prior_a")
# })
