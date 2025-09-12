# Helpers to create valid or invalid saemvsTuning objects for unit tests

# Valid tuning object
.valid_saemvsTuning <- function() {
  saemvsTuning(
    niter = 500,
    nburnin = 350,
    niter_mh = 5,
    kernel_mh = "random_walk",
    covariance_decay = 0.98,
    mh_proposal_scale = 1.5,
    spike_values_grid = c(0.01, 0.05, 0.1),
    n_is_samples = 10000,
    seed = 123,
    nb_workers = 2
  )
}

# Invalid: nburnin > niter
.invalid_nbigger <- function() {
  saemvsTuning(
    niter = 100,
    nburnin = 150,
    niter_mh = 5,
    kernel_mh = "random_walk",
    covariance_decay = 0.98,
    mh_proposal_scale = 1.5,
    spike_values_grid = c(0.01, 0.05),
    n_is_samples = 10000,
    seed = 123,
    nb_workers = 2
  )
}

# Invalid: spike_values_grid contains zero
.invalid_spike_zero <- function() {
  saemvsTuning(
    niter = 100,
    nburnin = 50,
    niter_mh = 5,
    kernel_mh = "random_walk",
    covariance_decay = 0.98,
    mh_proposal_scale = 1.5,
    spike_values_grid = c(0, 0.05),
    n_is_samples = 10000,
    seed = 123,
    nb_workers = 2
  )
}

# Invalid: kernel_mh wrong value
.invalid_kernel <- function() {
  saemvsTuning(
    niter = 100,
    nburnin = 50,
    niter_mh = 5,
    kernel_mh = "invalid_kernel",
    covariance_decay = 0.98,
    mh_proposal_scale = 1.5,
    spike_values_grid = c(0.01, 0.05),
    n_is_samples = 10000,
    seed = 123,
    nb_workers = 2
  )
}
