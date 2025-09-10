testthat::test_that("saemvsData constructor creates valid object", {
  y <- list(rnorm(5), rnorm(3))
  t <- list(1:5, 1:3)
  x_candidates <- matrix(rnorm(4), nrow = 2)
  x_forced <- matrix(rnorm(2), nrow = 2)

  obj <- saemvsData(y, t, x_candidates, x_forced)

  testthat::expect_s4_class(obj, "saemvsData")
  testthat::expect_equal(length(obj@y_series), 2)
  testthat::expect_equal(length(obj@t_series), 2)
  testthat::expect_equal(nrow(obj@x_candidates), 2)
  testthat::expect_equal(nrow(obj@x_forced), 2)
})

testthat::test_that("saemvsData constructor fails on invalid input", {
  # Mismatched y and t length
  testthat::expect_error(
    saemvsData(list(rnorm(5)), list(1:5, 1:3)),
    "must have the same number of elements"
  )

  # Non-numeric y
  testthat::expect_error(
    saemvsData(list(letters[1:5]), list(1:5)),
    "must be numeric"
  )

  # Wrong number of rows in x_candidates
  testthat::expect_error(
    saemvsData(list(rnorm(5), rnorm(3)), list(1:5, 1:3),
      x_candidates = matrix(rnorm(6), nrow = 3)
    ),
    "must have 2 rows"
  )
})

# ---------------------------------------------------------------------

testthat::test_that("saemvsDataFromDFs builds valid object", {
  # Example longitudinal data
  long_df <- data.frame(
    id = rep(1:2, each = 5),
    t = rep(1:5, 2),
    y = rnorm(10)
  )

  # Example covariate data
  covar_df <- data.frame(
    id = 1:2,
    x1 = rnorm(2),
    x2 = rnorm(2)
  )

  obj <- saemvsDataFromDFs(
    long_df, covar_df,
    id_col = "id", y_col = "y", t_col = "t",
    x_candidates_cols = "x1",
    x_forced_cols = "x2"
  )

  testthat::expect_s4_class(obj, "saemvsData")
  testthat::expect_length(obj@y_series, 2)
  testthat::expect_equal(ncol(obj@x_candidates), 1)
  testthat::expect_equal(ncol(obj@x_forced), 1)
})

testthat::test_that("saemvsDataFromDFs fails on missing columns", {
  long_df <- data.frame(
    id = rep(1:2, each = 5),
    t = rep(1:5, 2),
    y = rnorm(10)
  )
  covar_df <- data.frame(id = 1:2, x1 = rnorm(2))

  # Missing y_col in long_df
  testthat::expect_error(
    saemvsDataFromDFs(long_df, covar_df,
      id_col = "id", y_col = "not_y", t_col = "t"
    ),
    "Column 'not_y' is missing in long_df"
  )

  # Missing id_col in covar_df
  testthat::expect_error(
    saemvsDataFromDFs(long_df, covar_df[, "x1", drop = FALSE],
      id_col = "id", y_col = "y", t_col = "t"
    ),
    "Column 'id' is missing in covar_df"
  )
})

testthat::test_that("saemvsDataFromDFs fails on non-numeric covariates", {
  long_df <- data.frame(
    id = rep(1:2, each = 5),
    t = rep(1:5, 2),
    y = rnorm(10)
  )

  covar_df <- data.frame(
    id = 1:2,
    x1 = c("a", "b"), # non-numeric candidate
    x2 = rnorm(2)
  )

  testthat::expect_error(
    saemvsDataFromDFs(long_df, covar_df,
      id_col = "id", y_col = "y", t_col = "t",
      x_candidates_cols = "x1"
    ),
    "x_candidates must be numeric"
  )

  testthat::expect_error(
    saemvsDataFromDFs(long_df, covar_df,
      id_col = "id", y_col = "y", t_col = "t",
      x_forced_cols = "x1"
    ),
    "x_forced must be numeric"
  )
})

testthat::test_that("saemvsDataFromDFs detects individuals mismatch", {
  long_df <- data.frame(
    id = rep(1:3, each = 5),
    t = rep(1:5, 3),
    y = rnorm(15)
  )

  covar_df <- data.frame(
    id = 1:2, # missing individual 3
    x1 = rnorm(2)
  )

  testthat::expect_error(
    saemvsDataFromDFs(long_df, covar_df,
      id_col = "id", y_col = "y", t_col = "t"
    ),
    "Some individuals in long_df are not present in covar_df"
  )
})
