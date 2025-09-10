testthat::test_that("Valid saemvsModel object is created", {
  obj <- do.call(saemvsModel, valid_args_model)
  testthat::expect_s4_class(obj, "saemvsModel")
  testthat::expect_equal(obj@phi_dim, 3)
  testthat::expect_equal(obj@phi_to_select_idx, c(1, 2))
  testthat::expect_equal(obj@phi_fixed_idx, c(3))
})

testthat::test_that("Invalid saemvsModel cases fail as expected", {
  for (case_name in names(invalid_cases_model)) {
    args <- invalid_cases_model[[case_name]]$args
    expected_msg <- invalid_cases_model[[case_name]]$msg

    testthat::expect_error(
      do.call(saemvsModel, args),
      expected_msg,
      info = paste("Case:", case_name)
    )
  }
})
