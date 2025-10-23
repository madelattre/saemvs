# =========================================================
# test-classes-saemvsData.r
#
# Unit tests for saemvsData class and saemvsData constructor
# (see classes.r) using testthat
# =========================================================

test_that("saemvsData: valid object creation", {
  valid_data <- saemvsData(
    y = list(c(1, 2, 3), c(4, 5, 6)),
    t = list(c(0, 1, 2), c(0, 1, 2)),
    x_candidates = matrix(1:6, nrow = 2),
    x_forced = matrix(8:11, nrow = 2)
  )

  expect_s4_class(valid_data, "saemvsData")
  expect_true(validObject(valid_data))
})


test_that("saemvsData: invalid y type", {
  expect_error(
    saemvsData(
      y = list(c("a", "b", "c")),
      t = list(1:3)
    ),
    "numeric"
  )
})

test_that("saemvsData: mismatched lengths between y and t", {
  expect_error(
    saemvsData(
      y = list(1:3, 4:6),
      t = list(1:3)
    ),
    "y_series and t_series must have the same number of elements"
  )
})

test_that("saemvsData: non-numeric t", {
  expect_error(
    saemvsData(
      y = list(1:3),
      t = list(c("a", "b", "c"))
    ),
    "numeric"
  )
})

test_that("saemvsData: x_candidates with wrong number of rows", {
  expect_error(
    saemvsData(
      y = list(1:3, 4:6),
      t = list(1:3, 1:3),
      x_candidates = matrix(1:9, nrow = 3) # should have 2 rows
    ),
    "rows"
  )
})

test_that("saemvsData: x_forced with wrong number of rows", {
  expect_error(
    saemvsData(
      y = list(1:3, 4:6),
      t = list(1:3, 1:3),
      x_forced = matrix(1:9, nrow = 3) # should have 2 rows
    ),
    "rows"
  )
})

test_that("saemvsData: x_candidates non-numeric", {
  expect_error(
    saemvsData(
      y = list(1:3),
      t = list(1:3),
      x_candidates = matrix(c("a", "b"), nrow = 1)
    ),
    "numeric"
  )
})

test_that("saemvsData: x_candidates not a matrix", {
  expect_error(
    saemvsData(
      y = list(1:3),
      t = list(1:3),
      x_forced = c(1, 2, 3)
    ),
    "matrix"
  )
})

test_that(
  "saemvsData: error when x_forced has a column identical to x_candidates",
  {
    expect_error(
      saemvsData(
        y = list(1:2, 1:2),
        t = list(1:2, 1:2),
        x_candidates = matrix(c(1, 2, 3, 4), nrow = 2),
        x_forced = matrix(c(3, 4, 5, 6), nrow = 2)
      ),
      "Column 1 of 'x_forced' is identical to column 2 of 'x_candidates'"
    )
  }
)


# Unit tests for saemvsData class creation and validation
# with saemvsDataFromDFs constructor


long_df_valid <- data.frame(
  id = rep(1:3, each = 3),
  time = rep(1:3, times = 3),
  y = rnorm(9)
)

covar_df_valid <- data.frame(
  id = 1:3,
  age = c(30, 40, 50),
  bmi = c(22, 25, 28),
  group = c(1, 0, 1)
)


test_that("saemvsDataFromDFs: valid construction", {
  obj <- saemvsDataFromDFs(
    long_df = long_df_valid,
    covar_df = covar_df_valid,
    id_col = "id",
    y_col = "y",
    t_col = "time",
    x_candidates_cols = c("age", "bmi"),
    x_forced_cols = "group"
  )

  expect_s4_class(obj, "saemvsData")
  expect_true(validObject(obj))
  expect_equal(length(obj@y_series), 3)
  expect_equal(length(obj@t_series), 3)
})

test_that("saemvsDataFromDFs: missing column in long_df", {
  long_df_missing <- long_df_valid
  long_df_missing$time <- NULL
  expect_error(
    saemvsDataFromDFs(
      long_df = long_df_missing,
      covar_df = covar_df_valid,
      id_col = "id",
      y_col = "y",
      t_col = "time"
    ),
    "Column 'time' is missing in long_df."
  )
})

test_that("saemvsDataFromDFs: missing id in covar_df", {
  covar_df_missing <- covar_df_valid
  covar_df_missing$id <- NULL
  expect_error(
    saemvsDataFromDFs(
      long_df = long_df_valid,
      covar_df = covar_df_missing,
      id_col = "id",
      y_col = "y",
      t_col = "time"
    ),
    "Column 'id' is missing in covar_df."
  )
})

test_that("saemvsDataFromDFs: individual missing in covar_df", {
  long_df_extra <- rbind(long_df_valid, data.frame(id = 4, time = 1, y = 0.5))
  expect_error(
    saemvsDataFromDFs(
      long_df = long_df_extra,
      covar_df = covar_df_valid,
      id_col = "id",
      y_col = "y",
      t_col = "time"
    ),
    "Some individuals in long_df are not present in covar_df."
  )
})

test_that("saemvsDataFromDFs: x_forced not numeric", {
  covar_nonnum <- covar_df_valid
  covar_nonnum$group <- factor(covar_nonnum$group)
  expect_error(
    saemvsDataFromDFs(
      long_df = long_df_valid,
      covar_df = covar_nonnum,
      id_col = "id",
      y_col = "y",
      t_col = "time",
      x_forced_cols = "group"
    ),
    "x_forced must be numeric."
  )
})

test_that("saemvsDataFromDFs: x_candidates not numeric", {
  covar_nonnum <- covar_df_valid
  covar_nonnum$age <- as.character(covar_nonnum$age)
  expect_error(
    saemvsDataFromDFs(
      long_df = long_df_valid,
      covar_df = covar_nonnum,
      id_col = "id",
      y_col = "y",
      t_col = "time",
      x_candidates_cols = "age"
    ),
    "x_candidates must be numeric."
  )
})

test_that("saemvsDataFromDFs: no forced covariates", {
  obj <- saemvsDataFromDFs(
    long_df = long_df_valid,
    covar_df = covar_df_valid,
    id_col = "id",
    y_col = "y",
    t_col = "time",
    x_candidates_cols = c("age", "bmi"),
    x_forced_cols = NULL
  )
  expect_s4_class(obj, "saemvsData")
  expect_null(obj@x_forced)
})

test_that("saemvsDataFromDFs: default candidate columns", {
  obj <- saemvsDataFromDFs(
    long_df = long_df_valid,
    covar_df = covar_df_valid,
    id_col = "id",
    y_col = "y",
    t_col = "time"
  )
  expect_s4_class(obj, "saemvsData")
  expect_true(validObject(obj))
  expect_false(is.null(obj@x_candidates))
})
