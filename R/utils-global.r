# Declare global variable names to suppress R CMD check notes
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("iteration", "value", "crit") # variables used in ggplot or data.table
  )
}

# Create stubs for C++ functions dynamically generated at runtime, only for NOT_CRAN
if (getRversion() >= "2.15.1") {
  if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    g_vector_cpp <- function(...) stop("stub")
    metropolis_vector_cpp <- function(...) stop("stub")
    g_scalar_cpp <- function(...) stop("stub")
    rmvnorm_cpp <- function(...) stop("stub")
    logdmvnorm_cpp <- function(...) stop("stub")
  }
}

# Global variables
utils::globalVariables(c("j", "facet_label", "var"))

utils::globalVariables(c(
  "g_scalar_cpp",
  "g_vector_cpp",
  "rmvnorm_cpp",
  "logdmvnorm_cpp",
  "metropolis_vector_cpp",
  "loglik_cpp",
  "rmvnorm_mat_cpp"
))