testthat::test_that("Valid saemvsData object is created", {
  obj <- do.call(saemvsData, valid_args)
  testthat::expect_s4_class(obj, "saemvsData")
})

testthat::test_that("Invalid saemvsData cases fail as expected", {
  for (case_name in names(invalid_cases)) {
    args <- invalid_cases[[case_name]]$args
    expected_msg <- invalid_cases[[case_name]]$msg

    testthat::expect_error(
      do.call(saemvsData, args),
      expected_msg,
      info = paste("Case:", case_name)
    )
  }
})
