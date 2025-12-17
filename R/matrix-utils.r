#' Expand a design matrix into a list of individual-specific matrices
#'
#' Internal utility function that transforms a design matrix \code{x_mat} into
#' a list of matrices, one per individual, based on a binary support matrix.
#' Each individual matrix is replicated \code{q} times and processed via
#' \code{x_per_indiv_rcpp}.
#'
#' @param x_mat A numeric matrix of covariates (rows = individuals,
#' columns = features).
#' @param support A binary matrix indicating which features are active (1) or
#' inactive (0).
#' @param q Number of replications for each individual matrix.
#' @return A list of length equal to the number of rows in \code{x_mat}, each
#' element being a list of \code{q} matrices.
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item Each column of \code{support} is split and the indices of ones are
#'   extracted.
#'   \item Each row of \code{x_mat} is converted to a 1-row matrix.
#'   \item Each individual matrix is replicated \code{q} times.
#'   \item The replicated matrices are passed to \code{x_per_indiv_rcpp}
#'   together with the support indices.
#' }
#'
#' @examples
#' \dontrun{
#' x <- matrix(1:6, nrow = 2, ncol = 3)
#' supp <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2, byrow = TRUE)
#' expand_to_list(x, supp, q = 2)
#' }
#' @keywords internal
expand_to_list <- function(x_mat, support, q) {
  fun_name <- "expand_to_list"

  if (is.null(x_mat) || length(x_mat) == 0) {
    stop(sprintf("%s: x_mat must be a non-empty matrix", fun_name))
  }
  if (!is.matrix(x_mat)) {
    stop(sprintf("%s: x_mat must be a matrix", fun_name))
  }
  if (is.null(support) || length(support) == 0) {
    stop(sprintf("%s: support must be a non-empty matrix", fun_name))
  }
  if (!is.matrix(support)) {
    stop(sprintf("%s: support must be a matrix", fun_name))
  }
  if (ncol(x_mat) != nrow(support)) {
    stop(sprintf(
      "%s: x_mat must have as many columns as support has rows.",
      fun_name
    ))
  }
  if (!is.numeric(q) || length(q) != 1 || q < 1) {
    stop(sprintf("%s: q must be a single positive integer", fun_name))
  }

  split_cols <- split(support, col(support))
  supp_index <- lapply(split_cols, function(x) which(x == 1))
  x_list <- lapply(split(x_mat, row(x_mat)), matrix, nrow = 1)
  x_q_list <- lapply(x_list, function(x) replicate(q, x, simplify = FALSE))
  x <- lapply(x_q_list, x_per_indiv_rcpp, supp_index)
  return(x) # nolint: return-linter
}

#' Apply shrinkage to a covariance matrix
#'
#' This internal function applies shrinkage to a covariance matrix `sigma_old`
#' for specified indices. The diagonal elements corresponding to the shrink
#' indices are multiplied by `alphas`, while all other entries involving these
#' indices are set to zero. Elements not affected by shrinkage are copied from
#' `sigma_temp`.
#'
#' @param sigma_old Square matrix (\eqn{d x d}) representing the original
#' covariance. The diagonal elements of the shrink indices will be reduced.
#' @param sigma_temp square matrix (\eqn{d x d}) serving as reference for
#' non-shrinked elements. Must have the same dimensions as `sigma_old`.
#' @param shrink_indices Integer vector indicating the rows/columns to shrink.
#' The diagonals of these indices are multiplied by `alphas`, and the other
#' entries involving these indices are set to zero.
#' @param alphas Numeric vector of the same length as `shrink_indices`
#' specifying the shrinkage factors to apply to the diagonal elements of the
#' shrink indices.
#'
#' @return A square \eqn{d x d} matrix where:
#'   \itemize{
#'     \item Diagonal elements of the shrink indices are reduced according to
#'     `alphas`.
#'     \item Off-diagonal elements involving any shrink index are set to zero.
#'     \item Other elements are copied from `sigma_temp`.
#'   }
#'
#' @details
#' - This function is internal and should not be exported.
#' - Input matrices must be numeric, non-NULL, and of the same dimension.
#' - `shrink_indices` must contain valid integers between 1 and the matrix
#' dimension.
#' - The function validates inputs and throws informative errors if they are
#' incorrect (NULL matrices, dimension mismatch, invalid indices, or mismatched
#' `alphas` length).
#'
#' @examples
#' \dontrun{
#' sigma_old <- diag(3)
#' sigma_temp <- matrix(1, 3, 3)
#' shrink_indices <- c(1, 3)
#' alphas <- c(0.5, 0.2)
#' shrink_covariance_matrix(sigma_old, sigma_temp, shrink_indices, alphas)
#' }
#' @keywords internal
shrink_covariance_matrix <- function(
    sigma_old, sigma_temp, shrink_indices, alphas) {
  fun_name <- "shrink_covariance_matrix"

  if (is.null(sigma_old) || !is.matrix(sigma_old)) {
    stop(sprintf("%s: sigma_old must be a non-null matrix", fun_name))
  }
  if (is.null(sigma_temp) || !is.matrix(sigma_temp)) {
    stop(sprintf("%s: sigma_temp must be a non-null matrix", fun_name))
  }
  if (!all(dim(sigma_old) == dim(sigma_temp))) {
    stop(sprintf(
      "%s: sigma_old and sigma_temp must have the same dimensions",
      fun_name
    ))
  }
  d <- nrow(sigma_old)
  if (!is.numeric(shrink_indices) || any(shrink_indices < 1) ||
        any(shrink_indices > d)) {
    stop(sprintf(
      "%s: shrink_indices must be integers between 1 and %d",
      fun_name, d
    ))
  }
  if (!is.numeric(alphas) || length(alphas) != length(shrink_indices)) {
    stop(sprintf(
      "%s: alphas must be a numeric vector of same length as shrink_indices",
      fun_name
    ))
  }


  sigma_new <- matrix(0, nrow = d, ncol = d)
  shrink_mask <- rep(FALSE, d)
  shrink_mask[shrink_indices] <- TRUE

  for (i in 1:d) {
    for (j in 1:d) {
      if (i == j && shrink_mask[i]) {
        sigma_new[i, j] <- sigma_old[i, j] * alphas[which(shrink_indices == i)]
      } else if (shrink_mask[i] || shrink_mask[j]) {
        sigma_new[i, j] <- 0
      } else {
        sigma_new[i, j] <- sigma_temp[i, j]
      }
    }
  }

  return(sigma_new) # nolint: return-linter
}

#' Zero out specified rows and columns of a matrix
#'
#' Internal utility function that sets to zero all rows and columns of a numeric
#' matrix corresponding to the provided indices. All other elements remain
#' unchanged.
#'
#' @param mat A numeric matrix. Must be non-NULL.
#' @param indices Integer vector specifying the rows and columns to zero out.
#'   Must contain valid indices within the dimensions of `mat`.
#'
#' @return A numeric matrix of the same dimensions as `mat`, with the specified
#' rows and columns zeroed out.
#'
#' @details
#' - This is an internal function and should **not** be exported.
#' - If `indices` is NULL or empty, the function returns `mat` unchanged.
#' - If `mat` is NULL or not a matrix, or if any `indices` are invalid, the
#' function throws an informative error.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(1:9, nrow = 3, byrow = TRUE)
#' indices <- c(1, 3)
#' zero_out_shrinked(mat, indices)
#' }
#' @keywords internal
zero_out_shrinked <- function(mat, indices) {
  if (is.null(mat) || !is.matrix(mat)) {
    stop("zero_out_shrinked: 'mat' must be a non-NULL matrix.")
  }

  if (is.null(indices) || length(indices) == 0) {
    return(mat)
  }

  if (!all(indices %in% seq_len(nrow(mat))) ||
        !all(indices %in% seq_len(ncol(mat)))) {
    stop(
      paste0(
        "zero_out_shrinked: 'indices' contain values outside ",
        "the matrix dimensions."
      )
    )
  }

  mat_zeroed <- mat
  mat_zeroed[indices, ] <- 0
  mat_zeroed[, indices] <- 0
  return(mat_zeroed) # nolint: return-linter
}


#' Merge fixed and selected covariate supports
#'
#' Internal utility function that merges two covariate support matrices
#' (fixed and selected) into a single support matrix. The resulting matrix
#' does **not** include an intercept column, which simplifies subsequent
#' computations in SAEM algorithms.
#'
#' @param fixed_support Numeric matrix or NULL. Support for covariates that are
#' forced in the model. Must not include an intercept.
#' @param selected_support Numeric matrix or NULL. Support for covariates
#' selected by the model. Must not include an intercept.
#' @param nb_phi_s Integer. Number of selected parameters (columns) in the
#' model.
#' @param nb_phi_ns Integer. Number of non-selected parameters (columns) in the
#' model.
#' @param perm Integer vector. Permutation to reorder columns of the combined
#' support matrix.
#'
#' @return A numeric matrix combining `fixed_support` and `selected_support`,
#' reordered by `perm`. Returns `NULL` if both inputs are empty or NULL.
#'
#' @details
#' - Internal function; should **not** be exported.
#' - If `fixed_support` or `selected_support` is empty (NULL or all zeros),
#'   it is ignored in the combination.
#' - The output matrix does not contain an intercept column.
#'
#' @examples
#' \dontrun{
#' fixed_support <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
#' selected_support <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
#' perm <- c(2, 1)
#' merge_support(fixed_support, selected_support,
#'   nb_phi_s = 2, nb_phi_ns = 2,
#'   perm
#' )
#' }
#'
#' @keywords internal
merge_support <- function(
    fixed_support, selected_support, nb_phi_s, nb_phi_ns, perm) {
  if (!is.null(fixed_support) && !is.matrix(fixed_support)) {
    stop("merge_support: 'fixed_support' must be a matrix or NULL.")
  }
  if (!is.null(selected_support) && !is.matrix(selected_support)) {
    stop("merge_support: 'selected_support' must be a matrix or NULL.")
  }
  if (!is.integer(nb_phi_s) && !is.numeric(nb_phi_s)) {
    stop("merge_support: 'nb_phi_s' must be an integer.")
  }
  if (!is.integer(nb_phi_ns) && !is.numeric(nb_phi_ns)) {
    stop("merge_support: 'nb_phi_ns' must be an integer.")
  }
  if (!is.integer(perm) && !is.numeric(perm)) {
    stop("merge_support: 'inv_perm' must be an integer vector.")
  }

  if (!is_empty_support(selected_support)) {
    n_selected <- nrow(selected_support)
    selected_extended <- cbind(
      selected_support,
      matrix(0, nrow = n_selected, ncol = nb_phi_ns)
    )
    selected_ordered <- selected_extended[, perm, drop = FALSE]
  } else {
    selected_ordered <- NULL
  }

  if (!is_empty_support(fixed_support)) {
    combined_support <- rbind(fixed_support, selected_ordered)
  } else {
    combined_support <- selected_ordered
  }


  return(combined_support) # nolint: return-linter
}



#' Check if a support matrix is empty
#'
#' Internal utility function that tests whether a support matrix is considered
#' empty. A support is empty if it is `NULL`, has zero length, or contains only
#' zeros.
#'
#' @param support Numeric matrix, vector, or NULL. Support to be checked.
#'
#' @return Logical value:
#'   - `TRUE` if the support is empty,
#'   - `FALSE` otherwise.
#'
#' @details
#' - Internal function; should **not** be exported.
#' - Accepts both matrices and vectors.
#' - Returns `TRUE` if the input is `NULL`, has zero elements, or is entirely
#' zero.
#'
#' @examples
#' \dontrun{
#' is_empty_support(NULL) # TRUE
#' is_empty_support(matrix(0, nrow = 2, ncol = 2)) # TRUE
#' is_empty_support(matrix(c(0, 1), nrow = 1)) # FALSE
#' }
#' @keywords internal
is_empty_support <- function(support) {
  if (is.null(support)) {
    return(TRUE)
  }
  if (!is.numeric(support)) {
    stop("is_empty_support: 'support' must be numeric or NULL.")
  }
  if (length(support) == 0) {
    return(TRUE)
  }
  if (all(support == 0)) {
    return(TRUE)
  }
  return(FALSE) # nolint: return-linter
}

#' Check if a matrix is empty
#'
#' This internal utility function determines whether a given matrix is
#' considered "empty". A matrix is treated as empty if:
#' \itemize{
#'   \item it is \code{NULL},
#'   \item it has zero rows,
#'   \item or it has zero columns.
#' }
#' If the input is not \code{NULL} and not a matrix, the function raises an
#' error.
#'
#' @details
#' The definition of "empty" here refers only to structural properties
#' (i.e. dimensions). A matrix filled with zeros but with nonzero dimensions
#' is \emph{not} considered empty. This function is primarily intended to
#' simplify checks in algorithms where empty matrices must be handled
#' explicitly.
#'
#' @param mat A numeric matrix, or \code{NULL}.
#'
#' @return A logical scalar:
#' \code{TRUE} if \code{mat} is \code{NULL}, has zero rows, or zero columns;
#' \code{FALSE} otherwise.
#'
#' @examples
#' \dontrun{
#' is_empty_matrix(NULL) # TRUE
#' is_empty_matrix(matrix(nrow = 0)) # TRUE
#' is_empty_matrix(matrix(, 2, 0)) # TRUE
#' is_empty_matrix(matrix(0, 2, 2)) # FALSE
#' }
#' @keywords internal
is_empty_matrix <- function(mat) {
  if (is.null(mat)) {
    return(TRUE)
  }
  if (!is.matrix(mat)) {
    stop("is_empty_matrix: 'mat' must be a matrix or NULL.")
  }
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(TRUE)
  }
  return(FALSE) # nolint: return-linter
}


#' Extract indices of rows containing at least one "1"
#'
#' This internal utility function identifies the rows of a matrix that contain
#' at least one element equal to \code{1}. It returns the indices of such rows.
#'
#' @details
#' If the input matrix is considered empty according to
#' \code{\link{is_empty_matrix}}, the function returns an empty integer vector.
#' Rows are scanned using \code{\link[base]{apply}} over the first dimension,
#' checking for the presence of \code{1} in each row.
#'
#' @param mat A numeric matrix, or \code{NULL}. Must not be empty in the sense
#'   of dimensions (see \code{\link{is_empty_matrix}}), otherwise an empty
#'   result is returned.
#'
#' @return An integer vector of row indices where at least one entry equals
#'   \code{1}. Returns \code{integer(0)} if no such rows exist or if the matrix
#'   is empty.
#'
#' @examples
#' \dontrun{
#' m <- matrix(c(0, 1, 0, 0, 1, 0), nrow = 3, byrow = TRUE)
#' extract_rows_with_ones(m)
#' # Returns: 1 2
#'
#' extract_rows_with_ones(matrix(0, 3, 3))
#' # Returns: integer(0)
#'
#' extract_rows_with_ones(NULL)
#' # Returns: integer(0)
#' }
#' @seealso \code{\link{is_empty_matrix}}
#' @keywords internal
extract_rows_with_ones <- function(mat) {
  if (is_empty_matrix(mat)) {
    return(integer(0))
  }
  res <- which(apply(mat, 1, function(x) any(x == 1)))
  return(res) # nolint: return-linter
}

#' Extract a sub-support matrix by column indices
#'
#' This internal utility function extracts a subset of columns from a
#' covariate-support matrix. If the resulting submatrix is empty, invalid,
#' or contains only zeros, \code{NULL} is returned.
#'
#' @details
#' The function is primarily designed to manipulate support matrices
#' representing covariate selection patterns. It ensures consistent
#' handling of empty inputs:
#' \itemize{
#'   \item If \code{support} is empty according to
#'   \code{\link{is_empty_support}}, the result is \code{NULL}.
#'   \item If \code{idx} is \code{NULL} or has length zero, the result is
#'   \code{NULL}.
#'   \item If the extracted submatrix contains only zeros, the result is
#'   \code{NULL}.
#' }
#'
#' @param support A numeric matrix or \code{NULL}, representing a covariate
#'   support structure. Must have nonzero dimensions to be non-empty.
#' @param idx An integer or numeric vector of column indices to extract.
#'
#' @return A numeric matrix containing the selected columns of \code{support},
#'   or \code{NULL} if the input is empty, indices are invalid, or the extracted
#'   submatrix contains only zeros.
#'
#' @examples
#' \dontrun{
#' supp <- matrix(c(1, 0, 0, 1, 0, 0), nrow = 2, byrow = TRUE)
#'
#' extract_sub_support(supp, 1) # Returns matrix with first column
#' extract_sub_support(supp, 2:3) # Returns NULL if only zeros
#' extract_sub_support(NULL, 1) # Returns NULL
#' }
#' @seealso \code{\link{is_empty_support}}
#' @keywords internal
extract_sub_support <- function(support, idx) {
  if (is_empty_support(support)) {
    return(NULL)
  }
  if (is.null(idx) || length(idx) == 0) {
    return(NULL)
  }
  if (!is.numeric(idx)) {
    stop("extract_sub_support: 'idx' must be a numeric or integer vector.")
  }
  if (any(idx < 1 | idx > ncol(support))) {
    stop("extract_sub_support: 'idx' contains invalid column indices.")
  }

  submat <- support[, idx, drop = FALSE]
  if (all(submat == 0)) {
    return(NULL)
  }
  return(submat) # nolint: return-linter
}
