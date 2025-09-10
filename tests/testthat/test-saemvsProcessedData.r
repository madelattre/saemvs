# ---- Test valid object ----
testthat::test_that("Valid saemvsProcessedData object is created", {
  obj <- do.call(saemvsProcessedData, valid_args_processed)
  testthat::expect_s4_class(obj, "saemvsProcessedData")

  # Vérifications supplémentaires
  testthat::expect_equal(length(obj@y_series), 2)
  testthat::expect_equal(nrow(obj@x_phi_to_select), 2)
  testthat::expect_equal(nrow(obj@x_phi_not_to_select), 2)
  testthat::expect_length(obj@x_phi_not_to_select_list, 2)
})

# ---- Test invalid objects ----
testthat::test_that("Invalid saemvsProcessedData cases fail as expected", {
  for (case_name in names(invalid_cases_processed)) {
    args <- invalid_cases_processed[[case_name]]$args
    expected_msg <- invalid_cases_processed[[case_name]]$msg

    testthat::expect_error(
      do.call(saemvsProcessedData, args),
      expected_msg,
      info = paste("Case:", case_name)
    )
  }
})
