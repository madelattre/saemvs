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
#' @param saem_results An object of class \code{saemvsResults} containing
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
#' summary_saemvs(res)
#' }
#'
#' @export
setGeneric(
  "summary_saemvs",
  function(saem_results, ...) {
    standardGeneric("summary_saemvs")
  }
)

#' @rdname summary_saemvs
#' @exportMethod summary_saemvs
setMethod(
  "summary_saemvs",
  signature(
    saem_results = "saemvsResults"
  ),
  function(saem_results, ...) {
    digits <- list(...)$digits %||% 3
    fixed_param_idx <- saem_results@phi_fixed_idx
    phi_to_select_idx <- saem_results@phi_to_select_idx

    x_candidates_names <- saem_results@x_candidates_names %||% character(0)
    x_forced_names <- saem_results@x_forced_names %||% character(0)


    best_model_idx <- which.min(saem_results@criterion_values)

    forced_variables_idx <- saem_results@forced_variables_idx[[best_model_idx]]
    selected_variables_idx <- saem_results@selected_variables_idx[[best_model_idx]] # nolint: line_length_linter

    best_support_matrix <- saem_results@unique_support[[best_model_idx]]
    support_dim <- dim(best_support_matrix)
    num_phi_sel <- support_dim[2]

    # Récupération des noms de phi
    phi_names <- saem_results@phi_names


    cat("\n---- Selected Variables ----\n")
    cat("----------------------------\n\n")
    for (j in seq_len(num_phi_sel)) {
      # Candidate rows excluent l'intercept (row 1) et les covariables forcées
      candidate_rows <- setdiff(
        seq_len(support_dim[1])[-1], forced_variables_idx + 1
      )
      selected_rel <- which(best_support_matrix[candidate_rows, j] == TRUE)
      if (length(selected_rel) > 0) {
        selected_abs <- intersect(selected_rel, selected_variables_idx)
        selected_abs_names <- character(0)
        if (length(selected_abs) > 0 && length(x_candidates_names) > 0) {
          selected_abs_names <- x_candidates_names[selected_abs]
        }
        cat(sprintf(
          "  - %s : %s\n", phi_names[phi_to_select_idx[j]],
          if (length(selected_abs_names) > 0) {
            paste(selected_abs_names, collapse = ", ")
          } else {
            "(none)"
          }
        ))
      } else {
        cat(sprintf("  - %s : (none)\n", phi_names[phi_to_select_idx[j]]))
      }
    }

    beta_est <- saem_results@mle_estimates[[best_model_idx]]$beta
    gamma_est <- saem_results@mle_estimates[[best_model_idx]]$gamma


    colnames(beta_est) <- phi_names

    row_labels <- character(nrow(beta_est))
    row_labels[1] <- "\u03BC"

    # Forced variables
    if (length(forced_variables_idx) > 0) {
      forced_rows <- 1 + seq_along(forced_variables_idx) # 1 = intercept # nolint: commented_code_linter
      for (i in seq_along(forced_rows)) {
        idx <- forced_variables_idx[i]
        if (length(x_forced_names) >= idx) {
          name <- x_forced_names[idx]
        } else {
          name <- paste0("forced_", idx)
        }
        row_labels[forced_rows[i]] <- name
      }
    } else {
      forced_rows <- integer(0)
    }

    # Selected variables
    if (length(selected_variables_idx) > 0) {
      candidate_rows <- setdiff(2:nrow(beta_est), forced_rows)
      if (length(candidate_rows) != length(selected_variables_idx)) {
        warning("Mismatch between selected_variables_idx and available beta rows.") # nolint: line_length_linter
      }
      for (i in seq_along(selected_variables_idx)) {
        idx <- selected_variables_idx[i]
        if (length(x_candidates_names) >= idx) {
          name <- x_candidates_names[idx]
        } else {
          name <- paste0("selected_", idx)
        }
        row_labels[candidate_rows[i]] <- name
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
    if (is.null(fixed_param_idx) || length(fixed_param_idx) == 0) {
      print(round(gamma_est, digits))
    } else {
      print(round(zero_out_shrinked(gamma_est, fixed_param_idx), digits))
    }

    cat("\n")
  }
)
