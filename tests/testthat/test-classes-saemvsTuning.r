# =========================================================
# test-classes-saemvsTuning.r
#
# Unit tests for saemvsTuning class and saemvsTuning
# constructor (see classes.r) using testthat
# =========================================================

# Spike values
spike_grid_valid <- c(0.1, 0.5, 1)
spike_grid_empty <- numeric(0)
spike_grid_negative <- c(-0.1, 0.5)

# niter and nburnin
niter_valid <- 500
nburnin_valid <- 350
niter_invalid <- 10.5
nburnin_invalid <- 600

# niter_mh
niter_mh_valid <- 5
niter_mh_invalid <- 2.5

# kernel
kernel_valid <- "random_walk"
kernel_invalid <- "invalid"

# covariance_decay
cov_decay_valid <- 0.98
cov_decay_invalid0 <- 0
cov_decay_invalid1 <- 1

# mh_proposal_scale
mh_scale_valid <- 1
mh_scale_invalid <- 0

# n_is_samples
n_is_samples_valid <- 10000
n_is_samples_invalid <- 3.5

# seed
seed_valid <- 123
seed_invalid <- 2.5

# nb_workers
nb_workers_valid <- 2
nb_workers_invalid <- 0


test_that("Valid creation", {
  obj <- saemvsTuning(
    niter = niter_valid,
    nburnin = nburnin_valid,
    niter_mh = niter_mh_valid,
    kernel_mh = kernel_valid,
    covariance_decay = cov_decay_valid,
    mh_proposal_scale = mh_scale_valid,
    spike_values_grid = spike_grid_valid,
    n_is_samples = n_is_samples_valid,
    seed = seed_valid,
    nb_workers = nb_workers_valid
  )
  expect_s4_class(obj, "saemvsTuning")
})

test_that("niter, nburnin and niter_mh validation", {
  expect_error(saemvsTuning(
    niter = niter_invalid,
    spike_values_grid = spike_grid_valid
  ))
  expect_error(saemvsTuning(nburnin = -1, spike_values_grid = spike_grid_valid))
  expect_error(saemvsTuning(
    nburnin = nburnin_invalid, niter = niter_valid,
    spike_values_grid = spike_grid_valid
  ))
  expect_error(saemvsTuning(
    niter_mh = niter_mh_invalid,
    spike_values_grid = spike_grid_valid
  ))
})

test_that("kernel_mh validation", {
  expect_error(saemvsTuning(
    kernel_mh = kernel_invalid,
    spike_values_grid = spike_grid_valid
  ))
})

test_that("covariance_decay validation", {
  expect_error(saemvsTuning(
    covariance_decay = cov_decay_invalid0,
    spike_values_grid = spike_grid_valid
  ))
  expect_error(saemvsTuning(
    covariance_decay = cov_decay_invalid1,
    spike_values_grid = spike_grid_valid
  ))
})

test_that("mh_proposal_scale validation", {
  expect_error(saemvsTuning(
    mh_proposal_scale = mh_scale_invalid,
    spike_values_grid = spike_grid_valid
  ))
})

test_that("spike_values_grid validation", {
  expect_error(saemvsTuning(spike_values_grid = spike_grid_empty))
  expect_error(saemvsTuning(spike_values_grid = spike_grid_negative))
})

test_that("n_is_samples validation", {
  expect_error(saemvsTuning(
    n_is_samples = n_is_samples_invalid,
    spike_values_grid = spike_grid_valid
  ))
})

test_that("seed validation", {
  expect_error(saemvsTuning(
    seed = seed_invalid,
    spike_values_grid = spike_grid_valid
  ))
})

test_that("nb_workers validation", {
  expect_error(saemvsTuning(
    nb_workers = nb_workers_invalid,
    spike_values_grid = spike_grid_valid
  ))
})
