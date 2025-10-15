#' Compute threshold for spike-and-slab prior
#'
#' Computes the threshold value for variable selection in a spike-and-slab prior,
#'              based on prior variances and inclusion probability.
#' @param nu1 Numeric. Variance of the "slab" component (large variance).
#' @param nu0 Numeric. Variance of the "spike" component (small variance).
#' @param alpha Numeric. Prior inclusion probability of the variable.
#' @return Numeric. Threshold value used in spike-and-slab selection.
#' @keywords internal
threshold <- function(nu1, nu0, alpha) {
  value <- 2 * nu0 * nu1 * log(sqrt(nu1 / nu0) * (1 - alpha) / alpha) / (nu1 - nu0)
  return(sqrt(value))
}

#' Compute posterior inclusion probability p*
#'
#' Computes the posterior probability that a coefficient belongs to the "slab" component
#'              of a spike-and-slab prior. This is used for fast variable selection updates.
#' @param beta Numeric matrix. Current coefficient estimates (rows: covariates, columns: responses/parameters).
#' @param alpha Numeric vector. Prior inclusion probabilities for each coefficient.
#' @param nu0 Numeric. Variance of the spike (small variance) component.
#' @param nu1 Numeric. Variance of the slab (large variance) component.
#' @return Numeric matrix of the same dimension as \code{beta}. Each entry is the posterior probability p*.
#' @keywords internal
p_star <- function(beta, alpha, nu0, nu1) {
  norm1 <- stats::dnorm(beta, mean = 0, sd = sqrt(nu1))
  norm0 <- stats::dnorm(beta, mean = 0, sd = sqrt(nu0))
  num <- norm1 * rep(alpha, each = nrow(beta))
  denom <- num + norm0 * rep(1 - alpha, each = nrow(beta))
  p_star <- num / denom
  return(p_star)
}

#' Determine algorithmic case based on model configuration
#'
#' Determines the operational case for the algorithm (e.g., MAP vs MLE, full vs partial selection,
#' presence of fixed parameters) based on the provided configuration.
#' @param config List. Configuration object containing the following elements:
#'   \itemize{
#'     \item \code{method_type} Character. Either "map" or "mle".
#'     \item \code{num_parameters_not_to_select} Integer. Number of parameters not subject to selection.
#'     \item \code{fixed_parameters_indices} Integer vector. Indices of parameters fixed a priori.
#'   }
#' @return Character string indicating the operational case. Possible values:
#'   \itemize{
#'     \item "map_full_select": MAP with all parameters subject to selection.
#'     \item "map_part_select_nofixed": MAP with partial selection, no fixed parameters.
#'     \item "map_part_select_fixed": MAP with partial selection and some fixed parameters.
#'     \item "mle_nofixed": MLE with no fixed parameters.
#'     \item "mle_fixed": MLE with some fixed parameters.
#'   }
#' @keywords internal
get_case <- function(config) {
  if (config$method_type == "map") {
    if (config$num_parameters_not_to_select == 0) {
      "map_full_select"
    } else if (length(config$fixed_parameters_indices) == 0) {
      "map_part_select_nofixed"
    } else {
      "map_part_select_fixed"
    }
  } else {
    if (length(config$fixed_parameters_indices) > 0) "mle_fixed" else "mle_nofixed"
  }
}



estimate_phi_individuals <- function(data, model, init, maxit = 1000) {
  g <- model@model_func
  n_phi <- model@phi_dim
  n <- length(data@y_series)
  
  est_mat <- matrix(NA_real_, nrow = n, ncol = n_phi)
  
  sse_phi <- function(phi, y, t) {
    pred <- g(phi, t)
    sum((y - pred)^2)
  }
  
  # bornes simples basées sur init@intercept
  lower <- rep(0, n_phi)           # tous positifs pour ce modèle
  upper <- init@intercept * 10     # bornes larges
  
  for (i in seq_len(n)) {
    y <- data@y_series[[i]]
    t <- data@t_series[[i]]
    p0 <- init@intercept
    
    opt <- optim(
      par = p0,
      fn = sse_phi,
      y = y, t = t,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = maxit)
    )
    
    est_mat[i, ] <- opt$par
  }
  
  colnames(est_mat) <- paste0("phi", seq_len(n_phi))
  est_mat
}


build_init_from_phi_lasso <- function(est_indiv,
                                      data_processed,
                                      model,
                                      init,
                                      lasso_lambda = NULL) {
  phi_est_matrix <- est_indiv #est_indiv$phi_est

  beta_candidates <- NULL
  beta_forced <- NULL


  if (!is.matrix(phi_est_matrix)) stop("phi_est_matrix must be a matrix N x n_phi")
  n_phi <- ncol(phi_est_matrix)

  # Indices
  phi_to_select_idx <- model@phi_to_select_idx
  phi_not_to_select_idx <- setdiff(seq_len(n_phi), phi_to_select_idx)

  beta_sel <- NULL

  if (length(phi_to_select_idx) > 0 &&
    !is.null(data_processed@x_phi_to_select)) {
    x_sel <- data_processed@x_phi_to_select[, -1, drop = FALSE]
    # remove intercept column
    beta_sel <- matrix(0, ncol = length(phi_to_select_idx), nrow = ncol(x_sel))
    for (i in seq_along(phi_to_select_idx)) {
      phi_sel <- phi_est_matrix[, phi_to_select_idx[i]]
      # Lasso fit
      fit <- glmnet::glmnet(
        x = x_sel, y = phi_sel, alpha = 1,
        standardize = TRUE
      )
      s_val <- if (is.null(lasso_lambda)) min(fit$lambda) else lasso_lambda
      coefs_list <- coef(fit, s = s_val)
      beta_sel[, i] <- as.matrix(coefs_list[-1])
    }
  }


  # === 5. Beta and covariance for phi_not_to_select_idx ===

  beta_not_sel <- NULL
  x_not_sel <- data_processed@x_phi_not_to_select[, -1, drop = FALSE] # remove intercept column∏
  if (length(phi_not_to_select_idx) > 0) {
    if (!is_empty_matrix(x_not_sel)) {
      beta_not_sel <- matrix(0, ncol = length(phi_not_to_select_idx), nrow = ncol(x_not_sel))
      for (i in seq_along(phi_not_to_select_idx)) {
        phi_not_sel <- phi_est_matrix[, phi_not_to_select_idx[i]]
        fit <- lm(phi_not_sel ~ x_not_sel)
        coefs <- coef(fit)
        beta_not_sel[, i] <- coefs[-1]
      }
    }
  }


  # === Build beta_candidates and beta_forced ===

  supp <- model@x_forced_support
  if (is_empty_support(supp)) {
    beta_forced <- NULL
    beta_candidates <- matrix(0, nrow = ncol(x_sel), ncol = n_phi)
    if (!is.null(beta_sel)) beta_candidates[, phi_to_select_idx] <- beta_sel
  } else {
    x_forced_idx <- seq(1, nrow(supp))
    # beta_forced <- matrix(0, nrow = ncol(x_not_sel), ncol = n_phi)
    beta_forced <- matrix(0, nrow = nrow(supp), ncol = n_phi)
    if (!all(supp[, phi_not_to_select_idx] == 0)) {
      beta_forced[, phi_not_to_select_idx] <- beta_not_sel
    }
    if (!all(supp[, phi_to_select_idx] == 0)) {
      beta_forced[, phi_to_select_idx] <- beta_sel[nrow(supp), ]
    }
    beta_candidates <- matrix(0, nrow = ncol(x_sel[, -x_forced_idx]), ncol = n_phi)
    beta_candidates[, phi_to_select_idx] <- beta_sel[-seq(1, nrow(supp)), ]
  }


  # === 6. Return saemvsInit ===
  methods::new("saemvsInit",
    intercept = init@intercept,
    beta_candidates = beta_candidates,
    beta_forced = beta_forced,
    cov_re = init@cov_re,
    sigma2 = init@sigma2
  )
}
