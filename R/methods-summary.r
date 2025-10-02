#' Summary of SAEM Variable Selection Results
#'
#' Provides a human-readable summary of the results from a SAEM variable selection run.
#' For the best model (chosen by BIC/e-BIC), it displays:
#' \itemize{
#'   \item Selected covariates for each \eqn{\phi} parameter.
#'   \item Estimated regression coefficients (\eqn{\beta}) with labeled rows
#'         for intercept, forced covariates, and selected covariates.
#'   \item Estimated covariance matrix (\eqn{\Gamma}), where fixed-effect parameters
#'         are zeroed out for clarity.
#' }
#'
#' @param saem_results An object of class \code{saemvsResults} containing
#'   the output of a SAEMVS run.
#' @param digits Numeric scalar (default = 3). Number of decimal places to use
#'   when printing coefficients and covariance matrix.
#'
#' @details
#' The best model is identified as the one minimizing the selection criterion.
#' The function separates forced covariates (always included), selected covariates
#' (chosen by the algorithm), and the intercept. Coefficients are displayed in a
#' labeled table, followed by the estimated covariance matrix.
#'
#' @examples
#' \dontrun{
#' # Assuming 'res' is a saemvsResults object from a SAEMVS run
#' summary_saemvs(res, digits = 4)
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
  function(saem_results, digits = 3) {
    # --------------------------
    # Fixed-effect indices (used to zero out covariance display if needed)
    fixed_param_idx <- saem_results@phi_fixed_idx

    # --------------------------
    # Identify best model index according to selection criterion (BIC/e-BIC)
    best_model_idx <- which.min(saem_results@criterion_values)

    # Retrieve indices of forced and selected covariates
    forced_variables_idx <- saem_results@forced_variables_idx[[best_model_idx]]
    selected_variables_idx <- saem_results@selected_variables_idx[[best_model_idx]]

    # Support matrix for the best model
    # Rows = intercept + covariates, columns = phi parameters
    best_support_matrix <- saem_results@unique_support[[best_model_idx]]
    support_dim <- dim(best_support_matrix)
    num_phi_sel <- support_dim[2]

    # --------------------------
    # Display selected covariates by φ parameter
    cat("\n---- Selected Variables ----\n")
    cat("----------------------------\n\n")
    for (j in seq_len(num_phi_sel)) {
      # Candidate rows exclude intercept (row 1) and forced covariates
      candidate_rows <- setdiff(seq_len(support_dim[1])[-1], forced_variables_idx + 1)
      selected_rel <- which(best_support_matrix[candidate_rows, j] == TRUE)
      if (length(selected_rel) > 0) {
        selected_abs <- intersect(selected_rel, selected_variables_idx) # selected_variables_idx[selected_rel]
        cat(sprintf("  - \u03C6%d : %s\n", j, paste(selected_abs, collapse = ", ")))
      } else {
        cat(sprintf("  - \u03C6%d : (none)\n", j))
      }
    }

    # --------------------------
    # Extract MLE estimates: coefficients (beta) and covariance matrix (gamma)
    beta_est <- saem_results@mle_estimates[[best_model_idx]]$beta
    gamma_est <- saem_results@mle_estimates[[best_model_idx]]$gamma
    num_phi <- dim(beta_est)[2]

    # Column labels for φ parameters
    colnames(beta_est) <- paste0("\u03C6", seq_len(num_phi))

    # --------------------------
    # Row labels for beta coefficients
    row_labels <- character(nrow(beta_est))
    row_labels[1] <- "\u03BC" # Intercept always in first row

    # Forced covariates (if any)
    if (length(forced_variables_idx) > 0) {
      forced_rows <- 1 + seq_along(forced_variables_idx) # offset by 1 (intercept)
      for (i in seq_along(forced_rows)) {
        row_labels[forced_rows[i]] <- paste0("forced_", forced_variables_idx[i])
      }
    } else {
      forced_rows <- integer(0)
    }

    # Selected covariates (if any)
    if (length(selected_variables_idx) > 0) {
      # Candidate rows = all rows except intercept + forced covariates
      candidate_rows <- setdiff(2:nrow(beta_est), forced_rows)
      if (length(candidate_rows) != length(selected_variables_idx)) {
        warning("Mismatch between selected_variables_idx and available beta rows.")
      }
      for (i in seq_along(selected_variables_idx)) {
        row_labels[candidate_rows[i]] <- paste0("selected_", selected_variables_idx[i])
      }
    }

    # Assign row labels
    rownames(beta_est) <- row_labels

    # Covariance matrix labels
    colnames(gamma_est) <- paste0("\u03C6", seq_len(num_phi))
    rownames(gamma_est) <- paste0("\u03C6", seq_len(num_phi))

    # --------------------------
    # Print results
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
