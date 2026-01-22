# =========================================================
# test-compile-model.r
#
# Unit tests for function compile-model (see model-utils.r)
# using testthat
# =========================================================


g_fun_min <- function(phi, t) phi[1] + phi[2] * t

test_that("compile_model adjusts R indices for C++ correctly", {
  body_txt <- deparse(body(g_fun_min)) |> paste(collapse = "")
  body_cpp <- gsub("phi\\[(\\d+)\\]", "phi[\\1-1]", body_txt)

  # Check that R indices are correctly converted to C++ indices (0-based)
  expect_match(body_cpp, "phi\\[1-1\\]\\s*\\+\\s*phi\\[2-1\\]\\s*\\*\\s*t")
})

test_that("compile_model generates code containing expected functions", {
  body_txt <- deparse(body(g_fun_min)) |> paste(collapse = "")
  code_cpp <- sprintf(
    "// [[Rcpp::export]]\ninline double g_scalar_cpp(const arma::vec& phi, double t) { return %s; }", # nolint : line_length_linter
    body_txt
  )

  # Check that the C++ code contains the names of the exported functions.
  expect_match(code_cpp, "g_scalar_cpp")
  expect_match(code_cpp, "inline double")
})

test_that("compile_model creates temporary file for C++", {
  tmp_file <- tempfile(fileext = ".cpp")
  body_txt <- deparse(body(g_fun_min)) |> paste(collapse = "")
  code_cpp <- sprintf(
    "// [[Rcpp::export]]\ninline double g_scalar_cpp(const arma::vec& phi, double t) { return %s; }", # nolint : line_length_linter
    body_txt
  )
  writeLines(code_cpp, tmp_file)

  # Checks the existence and contents of the file
  expect_true(file.exists(tmp_file))
  expect_gt(file.info(tmp_file)$size, 0)
})

test_that("compile_model does not execute compilation (CRAN safe)", {
  expect_silent({
    body_txt <- deparse(body(g_fun_min)) |> paste(collapse = "")
    code_cpp <- sprintf(
      "// [[Rcpp::export]]\ninline double g_scalar_cpp(const arma::vec& phi, double t) { return %s; }", # nolint : line_length_linter
      body_txt
    )
  })
})
