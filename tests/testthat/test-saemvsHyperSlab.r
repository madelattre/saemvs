# --------------------------
# Valid tests for saemvsHyperSlab
# --------------------------
test_that("Valid saemvsHyperSlab object is created", {
  obj <- do.call(saemvsHyperSlab, valid_hyper_slab)
  expect_s4_class(obj, "saemvsHyperSlab")
  expect_true(obj@slab_parameter > 0)
  expect_true(obj@cov_re_prior_df > 0)
  expect_true(is.matrix(obj@cov_re_prior_scale))
})

# --------------------------
# Invalid tests for saemvsHyperSlab
# --------------------------
test_that("Invalid saemvsHyperSlab cases fail as expected", {
  for (case_name in names(invalid_cases_slab)) {
    args <- invalid_cases_slab[[case_name]]$args
    expected_msg <- invalid_cases_slab[[case_name]]$msg

    expect_error(
      methods::validObject(do.call(saemvsHyperSlab, args)),
      expected_msg,
      info = paste("Case:", case_name)
    )
  }
})

# --------------------------
# Valid tests for saemvsHyperSpikeAndSlab
# --------------------------
test_that("Valid saemvsHyperSpikeAndSlab object is created", {
  obj <- do.call(saemvsHyperSpikeAndSlab, valid_hyper_spike)
  expect_s4_class(obj, "saemvsHyperSpikeAndSlab")
  expect_true(obj@spike_parameter > 0)
})

# --------------------------
# Invalid tests for saemvsHyperSpikeAndSlab
# --------------------------
test_that("Invalid saemvsHyperSpikeAndSlab cases fail as expected", {
  for (case_name in names(invalid_cases_spike)) {
    args <- invalid_cases_spike[[case_name]]$args
    expected_msg <- invalid_cases_spike[[case_name]]$msg

    expect_error(
      methods::validObject(do.call(saemvsHyperSpikeAndSlab, args)),
      expected_msg,
      info = paste("Case:", case_name)
    )
  }
})
