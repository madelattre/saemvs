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
      x_candidates = c(1, 2, 3),
      x_forced = NULL
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

full_data_valid <- data.frame(
  id = rep(1:3, each = 2),
  time = rep(1:2, times = 3),
  y = rnorm(6),
  x1 = c(1, 1, 0, 0, 1, 1),
  age = c(30, 30, 40, 40, 50, 50),
  bmi = c(22, 22, 25, 25, 28, 28)
)


test_that("saemvsData_from_df: missing repeated column", {
  data_missing <- full_data_valid
  data_missing$time <- NULL

  expect_error(
    saemvsData_from_df(
      formula = y ~ . + repeated(time) + group(id),
      data = data_missing
    ),
    "Column 'time' is missing in data."
  )
})

test_that("saemvsData_from_df: missing group column", {
  data_missing <- full_data_valid
  data_missing$id <- NULL

  expect_error(
    saemvsData_from_df(
      formula = y ~ . + repeated(time) + group(id),
      data = data_missing
    ),
    "Column 'id' is missing in data."
  )
})

test_that("saemvsData_from_df: missing forced covariate", {
  data_missing <- full_data_valid
  data_missing$x1 <- NULL

  expect_error(
    saemvsData_from_df(
      formula = y ~ . + repeated(time) + group(id) + x1,
      data = data_missing
    ),
    "Column 'x1' is missing in data."
  )
})

test_that("saemvsData_from_df: multiple forced covariates, one missing", {
  data_missing <- full_data_valid
  data_missing$x2 <- NULL

  expect_error(
    saemvsData_from_df(
      formula = y ~ . + repeated(time) + group(id) + x1 + x2,
      data = data_missing
    ),
    "Column 'x2' is missing in data."
  )
})

test_that("saemvsData_from_df: data contains NA", {
  data_na <- full_data_valid
  data_na$bmi[1] <- NA

  expect_error(
    saemvsData_from_df(
      formula = y ~ . + repeated(time) + group(id),
      data = data_na
    ),
    "Data contains missing values"
  )
})

test_that("saemvsData_from_df: forced covariate not numeric", {
  data_nonnum <- full_data_valid
  data_nonnum$x1 <- factor(data_nonnum$x1)

  expect_error(
    saemvsData_from_df(
      formula = y ~ . + repeated(time) + group(id) + x1,
      data = data_nonnum
    ),
    "x_forced must be numeric."
  )
})

test_that("saemvsData_from_df: candidate covariate not numeric", {
  data_nonnum <- full_data_valid
  data_nonnum$age <- as.character(data_nonnum$age)

  expect_error(
    saemvsData_from_df(
      formula = y ~ . + repeated(time) + group(id),
      data = data_nonnum
    ),
    "x_candidates must be numeric."
  )
})

test_that("saemvsData_from_df: no forced covariates", {
  obj <- saemvsData_from_df(
    formula = y ~ . + repeated(time) + group(id),
    data = full_data_valid
  )

  expect_s4_class(obj, "saemvsData")
  expect_null(obj@x_forced)
})

test_that("saemvsData_from_df: default candidate covariates", {
  obj <- saemvsData_from_df(
    formula = y ~ . + repeated(time) + group(id),
    data = full_data_valid
  )

  expect_true(validObject(obj))
  expect_equal(
    colnames(obj@x_candidates),
    c("x1", "age", "bmi")
  )
})

test_that("saemvsData_from_df: forced covariates excluded from candidates", {
  obj <- saemvsData_from_df(
    formula = y ~ . + repeated(time) + group(id) + x1,
    data = full_data_valid
  )

  expect_equal(colnames(obj@x_forced), "x1")
  expect_false("x1" %in% colnames(obj@x_candidates))
})
