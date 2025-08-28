#' Check that a matrix is square, symmetric, and positive definite
#'
#' @param mat A numeric matrix to check.
#' @param name_mat Name of the matrix, used in messages.
#' @return NULL if valid; a character string with an error description otherwise.
#' @keywords internal
check_covariance <- function(mat, name_mat) {
  if (ncol(mat) != nrow(mat)) {
    return(sprintf("%s must be a square matrix.", name_mat))
  }
  if (!all.equal(mat, t(mat), tolerance = 1e-8)) {
    return(sprintf("%s must be symmetric.", name_mat))
  }
  eigenvalues <- eigen(mat, symmetric = TRUE)$values
  if (min(eigenvalues) < 1e-8) {
    return(sprintf("%s must be positive definite.", name_mat))
  }
  return(NULL)
}

#' Check consistency of beta and gamma matrices
#'
#' @param beta Beta coefficient matrix.
#' @param gamma Covariance matrix of random effects.
#' @param q Number of parameters (subject or not subject to selection).
#' @param case Character; either "ns" (not to select) or "s" (to select).
#' @keywords internal
check_beta_gamma <- function(beta, gamma, q, case = "ns") {
  start_case <- if (case == "ns") {
    sprintf(
      "As %d parameters are not subject to selection (%d missing parameter indices in 'index_select'),",
      q, q
    )
  } else {
    sprintf(
      "As %d parameters are subject to selection (%d indices in 'index_select'),",
      q, q
    )
  }

  if (q != 0) {
    if (is.null(beta) || is.null(gamma)) {
      stop(sprintf(
        "%s 'beta_%s' and 'gamma_%s' should both be non-empty matrices.",
        start_case, case, case
      ))
    }
  } else {
    if (ncol(beta) != ncol(gamma)) {
      stop(sprintf(
        "%s 'beta_%s' and 'gamma_%s' must both have %d columns.",
        start_case, case, case, q
      ))
    }
  }
}

#' Check beta matrix dimensions against support
#'
#' @param beta Beta coefficient matrix.
#' @param support Covariate support matrix (can be NULL).
#' @keywords internal
check_beta_support <- function(beta, support) {
  if (is.null(support) && nrow(beta) > 1) {
    stop("beta_ns must have only one row (intercepts) as no 'covariate_support' is provided in the model.")
  }
  if (!is.null(support) && nrow(support) != (nrow(beta) - 1)) {
    stop(sprintf(
      "beta_ns must have %d rows (intercept + %d covariates) as 'support' is designed for %d covariates (%d rows).",
      nrow(support) + 1, nrow(support), nrow(support), nrow(support)
    ))
  }
}

#' Check beta matrix for high-dimensional parameters
#'
#' @param beta Beta coefficient matrix.
#' @param pv Number of parameters to select (including intercept).
#' @keywords internal
check_beta_hdim <- function(beta, pv) {
  if (pv != nrow(beta)) {
    stop(sprintf(
      "beta_s must have %d rows (intercept + %d covariates) as matrix 'v' contains %d covariates.",
      pv, pv - 1, pv - 1
    ))
  }
}

#' Check that a numeric slot is strictly positive
#'
#' @param slot Numeric value.
#' @param name_slot Name of the slot.
#' @keywords internal
check_positive_slot <- function(slot, name_slot) {
  if (slot <= 0) {
    return(sprintf("%s must be strictly positive.", name_slot))
  }
  return(NULL)
}


#' Check that a numeric slot is either NULL or strictly positive
#'
#' @param slot Numeric value or NULL.
#' @param name_slot Name of the slot.
#' @keywords internal
check_positive_or_null_slot <- function(slot, name_slot) {
  if (!is.null(slot) && slot <= 0) {
    return(sprintf("%s must be NULL or strictly positive.", name_slot))
  }
  return(NULL)
}

#' Check that a numeric slot is a strictly positive integer
#'
#' @param slot Numeric value.
#' @param name_slot Name of the slot.
#' @keywords internal
check_positive_integer_slot <- function(slot, name_slot) {
  if ((slot <= 0) || !is.integer(slot)) {
    return(sprintf("%s must be a strictly positive integer.", name_slot))
  }
  return(NULL)
}
