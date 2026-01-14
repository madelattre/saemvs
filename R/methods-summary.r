# Internal helper for summary() to print one support
# @keywords internal
.summary_one_support <- function(
  object,
  support_idx,
  digits = 3
) {
  fixed_param_idx <- object@phi_fixed_idx
  phi_to_select_idx <- object@phi_to_select_idx

  x_candidates_names <- object@x_candidates_names %||% character(0)
  x_forced_names <- object@x_forced_names %||% character(0)

  forced_variables_idx <- object@forced_variables_idx[[support_idx]]
  selected_variables_idx <- object@selected_variables_idx[[support_idx]]

  support_matrix <- object@unique_support[[support_idx]]
  num_phi_sel <- ncol(support_matrix)

  phi_names <- object@phi_names

  cat("\n---- Selected Variables ----\n")
  cat("----------------------------\n\n")

  for (j in seq_len(num_phi_sel)) {
    candidate_rows <- setdiff(
      seq_len(nrow(support_matrix))[-1],
      forced_variables_idx + 1
    )

    selected_rel <- which(support_matrix[candidate_rows, j])
    selected_abs <- intersect(selected_rel, selected_variables_idx)

    selected_names <- if (
      length(selected_abs) > 0 && length(x_candidates_names) > 0
    ) {
      x_candidates_names[selected_abs]
    } else {
      character(0)
    }

    cat(sprintf(
      "  - %s : %s\n",
      phi_names[phi_to_select_idx[j]],
      if (length(selected_names) > 0)
        paste(selected_names, collapse = ", ")
      else
        "(none)"
    ))
  }

  beta_est <- object@mle_estimates[[support_idx]]$beta
  gamma_est <- object@mle_estimates[[support_idx]]$gamma

  colnames(beta_est) <- phi_names

  row_labels <- character(nrow(beta_est))
  row_labels[1] <- "\u03BC"

  if (length(forced_variables_idx) > 0) {
    forced_rows <- 1 + seq_along(forced_variables_idx)
    for (i in seq_along(forced_rows)) {
      idx <- forced_variables_idx[i]
      row_labels[forced_rows[i]] <- x_forced_names[idx] %||%
        paste0("forced_", idx)
    }
  } else {
    forced_rows <- integer(0)
  }

  if (length(selected_variables_idx) > 0) {
    candidate_rows <- setdiff(2:nrow(beta_est), forced_rows)
    for (i in seq_along(selected_variables_idx)) {
      idx <- selected_variables_idx[i]
      row_labels[candidate_rows[i]] <- x_candidates_names[idx] %||%
        paste0("selected_", idx)
    }
  }

  rownames(beta_est) <- row_labels
  colnames(gamma_est) <- phi_names
  rownames(gamma_est) <- phi_names

  cat("\n---- Estimated Parameters ----\n")
  cat("------------------------------\n\n")

  cat("---- Coefficients ----\n")
  print(round(beta_est, digits))

  cat("\n---- Covariance Matrix ----\n\n")
  if (length(fixed_param_idx) == 0) {
    print(round(gamma_est, digits))
  } else {
    print(round(zero_out_shrinked(gamma_est, fixed_param_idx), digits))
  }

  invisible(NULL)
}


#' Summary of SAEM Variable Selection Results
#'
#' Provides a human-readable summary of the results from a SAEM variable
#' selection run.
#' For the best model (chosen by BIC/e-BIC), it displays:
#' \itemize{
#'   \item Selected covariates for each \eqn{\phi} parameter.
#'   \item Estimated regression coefficients (\eqn{\beta}) with labeled rows
#'         for intercept, forced covariates, and selected covariates.
#'   \item Estimated covariance matrix (\eqn{\Gamma}), where fixed-effect
#'   parameters are zeroed out for clarity.
#' }
#'
#' @param object An object of class \code{saemvsResults} containing
#'   the output of a SAEMVS run.
#' @param ... Additional arguments passed to methods or generic functions.
#' Currently supported arguments:
#'   \describe{
#'     \item{digits}{Numeric scalar (default = 3). Number of decimal places
#'       used when printing coefficients and the covariance matrix.}
#'   }
#'
#' @details
#' The best model is identified as the one minimizing the selection criterion.
#' The function separates forced covariates (always included), selected
#' covariates (chosen by the algorithm), and the intercept. Coefficients are
#' displayed in a labeled table, followed by the estimated covariance matrix.
#'
#' @examples
#' \dontrun{
#' # Assuming 'res' is a saemvsResults object from a SAEMVS run
#' summary(res)
#' }
#'

#' @rdname summary
#' @export
setMethod(
  "summary",
  signature(object = "saemvsResults"),
  function(object, ...) {
    digits <- list(...)$digits %||% 3
    best_idx <- which.min(object@criterion_values)

    cat("\n==== Best model (",
        object@criterion,
        ") ====\n", sep = "")
    cat("Criterion value:",
        object@criterion_values[best_idx], "\n")

    .summary_one_support(object, best_idx, digits)
  }
)

#' Summary of a specific support
#' 
#' Displays:
#' \itemize{
#'   \item Selected covariates for each \eqn{\phi} parameter
#'   \item Estimated regression coefficients (\eqn{\beta})
#'   \item Estimated covariance matrix (\eqn{\Gamma})
#' }
#'
#' @param object A \code{saemvsResults} object
#' @param support_idx Integer. Index of the support to display in
#'  \code{unique_support}.
#' @param digits Number of digits for printing numeric values (default = 3)
#'
#' @export
#' @rdname summary
summary_support <- function(object, support_idx, digits = 3) {
  stopifnot(inherits(object, "saemvsResults"))

  n_sup <- length(object@unique_support)
  if (support_idx < 1 || support_idx > n_sup) {
    stop("support_idx must be between 1 and ", n_sup)
  }

  cat("\n==== Support", support_idx, "====\n")
  cat("Criterion value:",
      object@criterion_values[support_idx], "\n")

  .summary_one_support(object, support_idx, digits)
}

#' Summary of all supports
#'
#' Iterates through all unique supports in a SAEMVS run and prints a
#' summary for each.
#'
#' @param object A \code{saemvsResults} object
#' @param digits Number of digits for printing numeric values (default = 3)
#'
#' @export
#' @rdname summary
summary_all <- function(object, digits = 3) {
  stopifnot(inherits(object, "saemvsResults"))

  for (k in seq_along(object@unique_support)) {
    summary_support(object, k, digits)
    if (k < length(object@unique_support)) {
      cat("\n\n")
    }
  }
}
