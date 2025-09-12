test_that("saemvsTuning valid object works", {
  tuning <- .valid_saemvsTuning()
  expect_s4_class(tuning, "saemvsTuning")
  expect_length(tuning@step, tuning@niter)
  expect_true(all(tuning@spike_values_grid > 0))
})

test_that("saemvsTuning catches invalid nburnin > niter", {
  expect_error(
    .invalid_nbigger(),
    "'nburnin' must be smaller than 'niter'"
  )
})

test_that("saemvsTuning catches invalid spike_values_grid", {
  expect_error(
    .invalid_spike_zero(),
    "All values in 'spike_values_grid' must be strictly positive"
  )
})

test_that("saemvsTuning catches invalid kernel_mh", {
  expect_error(
    .invalid_kernel(),
    "'kernel_mh' must be either 'random_walk' or 'pop'"
  )
})
