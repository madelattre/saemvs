# =========================================================
# test-matrix-utils.r
#
# Unit tests for matrix-utils.r using testthat
# =========================================================

# --- Mock for x_per_indiv_rcpp ---
x_per_indiv_rcpp <- function(x_list, index_list) {
  n <- length(x_list)
  total_cols <- sum(sapply(index_list, length))
  mat <- matrix(0, nrow = n, ncol = total_cols)

  col_offset <- 0
  for (i in seq_len(n)) {
    x <- x_list[[i]]
    idx <- index_list[[i]]
    len <- length(idx)
    mat[i, (col_offset + 1):(col_offset + len)] <- x[idx]
    col_offset <- col_offset + len
  }
  mat
}

# ---------------------------------------------
# --- Tests for expand_to_list() ---
# ---------------------------------------------

test_that("expand_to_list handles basic input correctly", {
  x_mat <- matrix(1:4, nrow = 2, ncol = 2)
  support <- matrix(c(1, 0, 0, 1), nrow = 2)
  q <- 2
  res <- expand_to_list(x_mat, support, q)
  expect_true(is.list(res))

  n_indiv <- nrow(x_mat)
  total_cols <- sum(sapply(
    split(support, col(support)),
    function(x) sum(x == 1)
  ))
  expect_equal(length(res), n_indiv)
})


# ---------------------------------------------
# --- Tests for shrink_covariance_matrix() ---
# ---------------------------------------------

test_that("shrink_covariance_matrix applies shrinkage correctly", {
  sigma_old <- matrix(c(1, 2, 3, 4), 2, 2)
  sigma_temp <- matrix(c(10, 20, 30, 40), 2, 2)
  shrink_indices <- c(1)
  alphas <- c(0.5)

  res <- shrink_covariance_matrix(sigma_old, sigma_temp, shrink_indices, alphas)

  expect_equal(res[1, 1], sigma_old[1, 1] * alphas[1])
  expect_equal(res[1, 2], 0)
  expect_equal(res[2, 1], 0)
  expect_equal(res[2, 2], sigma_temp[2, 2])
})

test_that("shrink_covariance_matrix works with multiple shrink indices", {
  sigma_old <- diag(3)
  sigma_temp <- matrix(1, 3, 3)
  shrink_indices <- c(1, 3)
  alphas <- c(0.5, 0.2)

  res <- shrink_covariance_matrix(sigma_old, sigma_temp, shrink_indices, alphas)

  expect_equal(res[1, 1], 0.5)
  expect_equal(res[3, 3], 0.2)
  expect_equal(res[1, 2], 0)
  expect_equal(res[2, 1], 0)
  expect_equal(res[2, 2], 1)
})

test_that("shrink_covariance_matrix handles edge cases", {
  # No shrink
  sigma_old <- diag(2)
  sigma_temp <- matrix(9, 2, 2)
  res <- shrink_covariance_matrix(sigma_old, sigma_temp, integer(0), numeric(0))
  expect_equal(res, sigma_temp)

  # All indices
  shrink_indices <- 1:2
  alphas <- c(0.1, 0.2)
  res <- shrink_covariance_matrix(sigma_old, sigma_temp, shrink_indices, alphas)
  expect_equal(res[1, 1], 0.1)
  expect_equal(res[2, 2], 0.2)
  expect_equal(res[1, 2], 0)
  expect_equal(res[2, 1], 0)
})

test_that("shrink_covariance_matrix throws errors on invalid input", {
  sigma_old <- diag(2)
  sigma_temp <- diag(3)
  expect_error(shrink_covariance_matrix(sigma_old, sigma_temp, 1, 0.5))

  sigma_temp <- diag(2)
  expect_error(shrink_covariance_matrix(
    sigma_old, sigma_temp, c(0, 3),
    c(0.1, 0.2)
  ))

  expect_error(shrink_covariance_matrix(sigma_old, sigma_temp, 1, c(0.1, 0.2)))
})

# ---------------------------------------------
# --- Tests for zero_out_shrinked() ---
# ---------------------------------------------

test_that("zero_out_shrinked zeros out specified rows and columns", {
  mat <- matrix(1:9, nrow = 3, byrow = TRUE)
  indices <- c(1, 3)

  res <- zero_out_shrinked(mat, indices)

  expect_equal(res[1, ], c(0, 0, 0))
  expect_equal(res[3, ], c(0, 0, 0))
  expect_equal(res[, 1], c(0, 0, 0))
  expect_equal(res[, 3], c(0, 0, 0))
  expect_equal(res[2, 2], mat[2, 2])
})

test_that("zero_out_shrinked returns mat unchanged for empty indices", {
  mat <- matrix(1:4, 2, 2)
  expect_equal(zero_out_shrinked(mat, NULL), mat)
  expect_equal(zero_out_shrinked(mat, integer(0)), mat)
})

test_that("zero_out_shrinked throws errors on invalid input", {
  mat <- matrix(1:4, 2, 2)

  expect_error(zero_out_shrinked(NULL, 1))
  expect_error(zero_out_shrinked(1:4, 1))

  expect_error(zero_out_shrinked(mat, 0))
  expect_error(zero_out_shrinked(mat, 3))
  expect_error(zero_out_shrinked(mat, c(1, 4)))
})


# ---------------------------------------------
# --- Tests for merge_support() ---
# ---------------------------------------------

test_that("merge_support combines fixed and selected supports correctly", {
  # Mimicks a situation with 3 parameters: 2 to select, 1 not to select
  # 2 covariates forced, 2 selectable
  fixed_support <- matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, byrow = TRUE)
  selected_support <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  nb_phi_s <- 2
  nb_phi_ns <- 1
  perm <- c(2, 1, 3)

  res <- merge_support(
    fixed_support, selected_support, nb_phi_s, nb_phi_ns,
    perm
  )

  # Check dimensions are correct
  expect_equal(nrow(res), nrow(fixed_support) + nrow(selected_support))
  expect_equal(ncol(res), nb_phi_s + nb_phi_ns)
})

test_that("merge_support handles empty supports", {
  nb_phi_s <- 2
  nb_phi_ns <- 2
  perm <- c(1, 2, 3, 4)

  # empty selected_support
  res1 <- merge_support(
    matrix(c(1, 0, 0, 1), 2, 2), matrix(0, 2, 2), nb_phi_s,
    nb_phi_ns, perm
  )
  expect_equal(nrow(res1), 2)

  # empty fixed_support
  res2 <- merge_support(
    NULL, matrix(c(0, 1, 1, 0), 2, 2), nb_phi_s, nb_phi_ns,
    perm
  )
  expect_equal(nrow(res2), 2)

  # empty selected_support and empty fixed_support
  res3 <- merge_support(NULL, NULL, nb_phi_s, nb_phi_ns, perm)
  expect_null(res3)
})

test_that("merge_support throws errors for invalid input types", {
  expect_error(merge_support(1:4, NULL, 2, 2, c(1, 2)))
  expect_error(merge_support(NULL, 1:4, 2, 2, c(1, 2)))
  expect_error(merge_support(NULL, NULL, "2", 2, c(1, 2)))
  expect_error(merge_support(NULL, NULL, 2, "2", c(1, 2)))
  expect_error(merge_support(NULL, NULL, 2, 2, "perm"))
})

# ---------------------------------------------
# --- Tests for is_empty_support() ---
# ---------------------------------------------

test_that("is_empty_support returns TRUE for empty inputs", {
  expect_true(is_empty_support(NULL))
  expect_true(is_empty_support(numeric(0)))
  expect_true(is_empty_support(matrix(0, nrow = 2, ncol = 2)))
})

test_that("is_empty_support returns FALSE for non-empty inputs", {
  expect_false(is_empty_support(matrix(c(0, 1, 0, 0), nrow = 2)))
  expect_false(is_empty_support(c(0, 1)))
  expect_false(is_empty_support(matrix(1, nrow = 1, ncol = 1)))
})

test_that("is_empty_support throws error for non-numeric input", {
  expect_error(is_empty_support("string"))
  expect_error(is_empty_support(list(1, 2, 3)))
  expect_error(is_empty_support(TRUE))
})

# ---------------------------------------------
# --- Tests for is_empty_matrix() ---
# ---------------------------------------------

test_that("is_empty_matrix returns TRUE for empty inputs", {
  expect_true(is_empty_matrix(NULL))
  expect_true(is_empty_matrix(matrix(nrow = 0, ncol = 2)))
  expect_true(is_empty_matrix(matrix(nrow = 2, ncol = 0)))
})

test_that("is_empty_matrix returns FALSE for non-empty matrices", {
  expect_false(is_empty_matrix(matrix(0, nrow = 1, ncol = 1)))
  expect_false(is_empty_matrix(matrix(1:4, nrow = 2, ncol = 2)))
})

test_that("is_empty_matrix throws error for non-matrix input", {
  expect_error(is_empty_matrix(1:4))
  expect_error(is_empty_matrix(list(1, 2, 3)))
  expect_error(is_empty_matrix("string"))
  expect_error(is_empty_matrix(TRUE))
})


# ---------------------------------------------
# --- Tests for extract_rows_with_ones() ---
# ---------------------------------------------

test_that("extract_rows_with_ones returns integer(0) for empty matrix", {
  expect_identical(extract_rows_with_ones(NULL), integer(0))
  expect_identical(
    extract_rows_with_ones(matrix(nrow = 0, ncol = 2)),
    integer(0)
  )
  expect_identical(
    extract_rows_with_ones(matrix(nrow = 2, ncol = 0)),
    integer(0)
  )
})

test_that("extract_rows_with_ones returns integer(0) when no 1s present", {
  mat <- matrix(0, nrow = 3, ncol = 3)
  expect_identical(extract_rows_with_ones(mat), integer(0))
})

test_that(
  "extract_rows_with_ones returns correct row indices when 1s are present",
  {
    mat <- matrix(c(
      0, 1, 0,
      0, 0, 0,
      1, 0, 1
    ), nrow = 3, byrow = TRUE)
    expect_identical(extract_rows_with_ones(mat), as.integer(c(1, 3)))
  }
)

test_that(
  "extract_rows_with_ones works for single-row and single-column matrices",
  {
    expect_identical(
      extract_rows_with_ones(matrix(1, nrow = 1, ncol = 1)),
      as.integer(1)
    )
    expect_identical(
      extract_rows_with_ones(matrix(0, nrow = 1, ncol = 1)),
      integer(0)
    )
    expect_identical(extract_rows_with_ones(matrix(c(0, 1, 0),
      nrow = 3,
      ncol = 1
    )), as.integer(2))
  }
)


# ---------------------------------------------
# --- Tests for extract_sub_support() ---
# ---------------------------------------------

test_that("extract_sub_support returns submatrix for valid indices", {
  supp <- matrix(c(
    1, 0, 0,
    0, 1, 1
  ), nrow = 2, byrow = TRUE)

  expect_identical(
    extract_sub_support(supp, 1:2),
    supp[, 1:2, drop = FALSE]
  )

  expect_identical(
    extract_sub_support(supp, 3),
    supp[, 3, drop = FALSE]
  )
})

test_that("extract_sub_support returns NULL for empty support", {
  supp_empty <- matrix(0, nrow = 2, ncol = 3)
  expect_null(extract_sub_support(supp_empty, 1))
  expect_null(extract_sub_support(NULL, 1))
})

test_that("extract_sub_support returns NULL for NULL or empty idx", {
  supp <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 2, byrow = TRUE)
  expect_null(extract_sub_support(supp, NULL))
  expect_null(extract_sub_support(supp, integer(0)))
})

test_that("extract_sub_support errors for invalid idx type", {
  supp <- matrix(1:4, nrow = 2)
  expect_error(
    extract_sub_support(supp, "a"),
    "extract_sub_support: 'idx' must be a numeric or integer vector."
  )
})

test_that("extract_sub_support errors for out-of-bounds indices", {
  supp <- matrix(1:4, nrow = 2)
  expect_error(
    extract_sub_support(supp, 0),
    "extract_sub_support: 'idx' contains invalid column indices."
  )
  expect_error(
    extract_sub_support(supp, 5),
    "extract_sub_support: 'idx' contains invalid column indices."
  )
})

test_that(
  "extract_sub_support returns NULL if extracted submatrix is all zeros",
  {
    supp <- matrix(c(
      1, 0, 0,
      0, 0, 0
    ), nrow = 2, byrow = TRUE)
    expect_null(extract_sub_support(supp, 2))
  }
)
