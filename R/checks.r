#' Check that a matrix is square, symmetric, and positive definite
#'
#' @param mat A numeric matrix to check.
#' @param name_mat Name of the matrix, used in messages.
#' @return NULL if valid; a character string with an error description
#' otherwise.
#' @keywords internal
#' @noRd
check_covariance <- function(mat, name_mat) {
  if (!is.numeric(mat)) {
    return(sprintf("%s must be a numeric matrix.", name_mat))
  }
  if (ncol(mat) != nrow(mat)) {
    return(sprintf("%s must be a square matrix.", name_mat))
  }
  if (!isTRUE(all.equal(mat, t(mat), tolerance = 1e-8))) {
    return(sprintf("%s must be symmetric.", name_mat))
  }
  eigenvalues <- eigen(mat, symmetric = TRUE)$values
  if (min(eigenvalues) < 1e-8) {
    return(sprintf("%s must be positive definite.", name_mat))
  }
  return(NULL)
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
#' @noRd
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
#' @noRd
check_positive_integer_slot <- function(slot, name_slot) {
  if ((slot <= 0) || !is.integer(slot)) {
    return(sprintf("%s must be a strictly positive integer.", name_slot))
  }
  return(NULL)
}
