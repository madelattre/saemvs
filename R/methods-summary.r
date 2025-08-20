#' @export
setGeneric(
  "summary",
  function(res, dec) {
    standardGeneric("summary")
  }
)

#' @exportMethod summary
setMethod(
  "summary",
  signature(
    res = "resSAEMVS",
    dec = "numeric"
  ),
  function(res, dec = 3) {
    # Vérifier que res n'est pas vide

    selected_model <- which.min(res@crit_values)

    selected_support <- res@unique_support[[selected_model]]

    dim_s <- dim(selected_support)


    cat("\n")
    cat("---- Selected variables ----\n")
    cat("----------------------------\n")
    cat("\n")

    all_select_var <- c()

    for (j in seq(dim_s[2])) {
      select_var <- which(selected_support[-1, j])
      if (length(select_var) > 0) {
        cat(sprintf("  - \u03C6%d : %s\n", j, paste(select_var,
          collapse = ", "
        )))
      } else {
        cat(sprintf("  - \u03C6%d : (none)\n", j))
      }

      all_select_var <- c(all_select_var, select_var)
    }

    all_select_var_unique <- unique(all_select_var)

    beta_est <- res@est_mle[[selected_model]]$beta
    nb_phi <- dim(beta_est)[2]
    gamma_est <- res@est_mle[[selected_model]]$gamma

    colnames(beta_est) <- paste0("\u03C6", 1:nb_phi)
    n <- nrow(beta_est)
    n_v <- length(all_select_var_unique)
    if (is.null(rownames(beta_est))) {
      rownames(beta_est) <- paste0("Row", 1:n)
    }
    last_idx <- (n - n_v + 1):n
    rownames(beta_est)[last_idx] <- paste0("v", all_select_var_unique)
    first_idx <- setdiff(1:n, last_idx)
    if (length(first_idx) == 1) {
      rownames(beta_est)[first_idx] <- "\u03BC"
    } else if (length(first_idx) > 1) {
      rownames(beta_est)[first_idx] <- c(
        "\u03BC",
        paste0("w", seq_along(first_idx[-1]))
      )
    }
    colnames(gamma_est) <- paste0("\u03C6", 1:nb_phi)
    rownames(gamma_est) <- paste0("\u03C6", 1:nb_phi)

    cat("\n")
    cat("\n")
    cat("---- Estimated parameters ----\n")
    cat("------------------------------\n")

    cat("\n")
    cat("---- Coefficients ----\n")
    print(round(beta_est, dec))
    cat("\n")
    cat("---- Covariance matrix ----\n")
    cat("\n")
    print(round(gamma_est, dec))
    cat("\n")
  }
)
