# -----------------------------
# Helper and test functions
# -----------------------------
g_helper <- function(x) x^2  # external helper, not a model

test_functions <- list(
  f1 = function(t, a) a * t + 1,
  f2 = function(time, p1, p2) p1 + p2 / time,
  f3 = function(t, k) exp(-k * t),
  f4 = function(t, a, b) log(a + b * t) / sqrt(t),
  f5 = function(t, a, b) if (t > 0) a / t else b,
  f6 = function(t, x) { y <- x + t; z <- y^2; z }, # nolint : brace_linter
  f7 = function(t, n) { s <- 0; for (i in 1:n) s <- s + i * t; s }, # nolint : brace_linter
  f8 = function(t, max_val) { i <- 0; sum <- 0; while (i < max_val) { sum <- sum + t*i; i <- i + 1 }; sum }, # nolint : line_length_linter
  f9 = function(t, a) g_helper(a + t),
  f10 = function(t, a, vec) vec[1] + vec[2] / a + t,
  f11 = function(t, a) { if (is.na(a)) return(NA_real_); if (is.infinite(a)) return(Inf); a / t }, # nolint : line_length_linter
  f12 = function(t, a) a + NaN,
  f13 = function(t, a, b, c) (a + b) / (1 + exp(-(t - c))) - sqrt(a * b),
  f14 = function(t, x) x^2 + abs(x) - floor(x) + cos(t),
  f15 = function(t, p1, p2, p3) { phi_vec <- c(p1, p2, p3); sum(phi_vec * t) } # nolint : brace_linter
)

# -----------------------------
# 2️⃣ Helper function for testing
# -----------------------------
run_make_phi_tests <- function(f_name, f) {

  # Skip helpers
  if (grepl("helper", f_name)) return(NULL)

  # apply make_phi_fn
  f_phi <- make_phi_fn(f)
  environment(f_phi)$g_helper <- g_helper

  # apply transpile_to_cpp
  cpp_code <- transpile_to_cpp(f_phi, function_name = paste0(f_name, "_cpp"))

  # R evaluation with arbitrary values
  param_names <- names(formals(f_phi))[-1] # skip 't'
  phi_test <- seq_along(param_names)
  t_test <- 1
  args_list <- c(list(t = t_test), setNames(as.list(phi_test), param_names))

  eval_result <- do.call(f_phi, args_list)

  list(f_phi = f_phi, cpp_code = cpp_code, eval_result = eval_result)
}

# -----------------------------
# Create testthat tests
# -----------------------------
test_that("make_phi_fn and transpile_to_cpp produce valid output", {
  for (name in names(test_functions)) {

    # Skip helper functions
    if (grepl("helper", name)) next

    f <- test_functions[[name]]

    res <- run_make_phi_tests(name, f)

    # Ensure f_phi is a function
    expect_type(res$f_phi, "closure")

    # Ensure cpp_code is a string
    expect_type(res$cpp_code, "character")

    # Ensure eval_result does not throw an error
    expect_false(inherits(res$eval_result, "error"))

    # Optionally: 
    # check result is numeric (most functions produce numeric output)
    expect_type(res$eval_result, "double")
  }
})
