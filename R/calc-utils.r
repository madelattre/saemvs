#' Compute threshold for spike-and-slab prior
#'
#' Computes the threshold value for variable selection with the Gaussian
#'  spike-and-slab prior, based on prior variances and inclusion probability.
#' @param nu1 Numeric. Variance of the "slab" component (large variance).
#' @param nu0 Numeric. Variance of the "spike" component (small variance).
#' @param alpha Numeric. Prior inclusion probability of the variable.
#' @return Numeric. Threshold value used in spike-and-slab selection.
#' @keywords internal
threshold <- function(nu1, nu0, alpha) {
  value <- 2 * nu0 * nu1 * log(
    sqrt(nu1 / nu0) * (1 - alpha) / alpha
  ) / (nu1 - nu0)
  sqrt(value)
}

#' Compute posterior inclusion probability p*
#'
#' Computes the posterior probability that a coefficient belongs to the
#' "slab" component of the spike-and-slab prior.
#' @param beta Numeric matrix. Current coefficient estimates
#' (rows: covariates, columns: responses/parameters).
#' @param alpha Numeric vector. Prior inclusion probabilities for each
#' coefficient.
#' @param nu0 Numeric. Variance of the spike (small variance) component.
#' @param nu1 Numeric. Variance of the slab (large variance) component.
#' @return Numeric matrix of the same dimension as \code{beta}. Each entry
#' is the posterior probability p*.
#' @keywords internal
p_star <- function(beta, alpha, nu0, nu1) {
  norm1 <- stats::dnorm(beta, mean = 0, sd = sqrt(nu1))
  norm0 <- stats::dnorm(beta, mean = 0, sd = sqrt(nu0))
  num <- norm1 * rep(alpha, each = nrow(beta))
  denom <- num + norm0 * rep(1 - alpha, each = nrow(beta))
  num / denom
}

#' Determine algorithmic case based on model configuration
#'
#' Determines the operational case for the algorithm (e.g., MAP vs MLE,
#' full vs partial selection, presence of fixed parameters) based on the
#' provided configuration.
#' @param config List. Configuration object containing in particular the
#'  following elements:
#'   \itemize{
#'     \item \code{method_type} Character. Either "map" or "mle".
#'     \item \code{num_parameters_not_to_select} Integer. Number of parameters
#'  not subject to selection.
#'     \item \code{fixed_parameters_indices} Integer vector. Indices of
#'  parameters without random variability.
#'   }
#' @return Character string indicating the operational case. Possible values:
#'   \itemize{
#'     \item "map_full_select": MAP with all parameters subject to selection.
#'     \item "map_part_select_nofixed": MAP with some (but not all) of the
#'  parameters subject to selection, no fixed parameters.
#'     \item "map_part_select_fixed": MAP with some (but not all) of the
#'  parameters subject to selection and some fixed parameters.
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
    if (length(config$fixed_parameters_indices) > 0) {
      "mle_fixed"
    } else {
      "mle_nofixed"
    }
  }
}


#' Estimate individual-level phi parameters
#'
#' This internal function estimates the individual-specific \eqn{\varphi}
#' parameters for each series in the provided dataset, given a model structure
#' and initial parameter values. Optimization is performed independently for
#' each individual using the L-BFGS-B method.
#'
#' @param data An object of class \code{saemvsData} containing the observed
#' data. It is expected to provide the time series via the slots
#' \code{y_series} (list of response vectors) and \code{t_series} (list of
#' corresponding time vectors).
#' @param model An object of class \code{saemvsModel} representing the model
#' structure. It should provide at least the following slots:
#'   \describe{
#'     \item{\code{model_func}}{A function \eqn{g(\varphi, t)} that returns the
#'     model prediction given parameters and time.}
#'     \item{\code{phi_dim}}{An integer specifying the dimension of the
#'     parameter vector \eqn{\varphi}.}
#'   }
#' @param init An object of class \code{saemvsInit} providing initial values
#' for optimization. It is expected to include an \code{intercept} slot
#' used to initialize and bound the search.
#' @param maxit Integer. Maximum number of iterations allowed for the optimizer.
#' Defaults to 1000.
#'
#' @details
#' For each individual \eqn{i}, the function minimizes the sum of squared
#' errors:
#' \deqn{SSE(\varphi_i) = \sum_t (y_{i,t} - g(\varphi_i, t))^2}
#' using the \code{optim()} function with box constraints. Lower bounds are set
#' to zero, and upper bounds are ten times the initial intercept value.
#'
#' The optimization results are stored in a matrix with one row per individual
#' and one column per parameter.
#'
#' @return A numeric matrix of dimension \eqn{n \times q}, where
#'   \eqn{n} is the number of individuals and \eqn{q} is the number of
#'   estimated \eqn{\varphi} parameters. Column names are of the form
#'   \code{"phi1"}, \code{"phi2"}, etc.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' est <- estimate_phi_individuals(
#'   data = my_data, model = my_model,
#'   init = my_init
#' )
#' }
#'
estimate_phi_individuals <- function(data, model, init, maxit = 1000) {
  g <- model@model_func
  n_phi <- model@phi_dim
  n <- length(data@y_series)

  est_mat <- matrix(NA_real_, nrow = n, ncol = n_phi)

  sse_phi <- function(phi, y, t) {
    pred <- g(phi, t)
    sum((y - pred)^2)
  }

  lower <- rep(0, n_phi)
  upper <- init@intercept * 10

  for (i in seq_len(n)) {
    y <- data@y_series[[i]]
    t <- data@t_series[[i]]
    phi0 <- init@intercept

    opt <- stats::optim(
      par = phi0,
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

#' Build initial values from individual phi estimates using Lasso
#'
#' This internal function constructs an initial \code{saemvsInit} object for
#' the SAEMVS algorithm based on individual-level \eqn{\phi} estimates. Lasso
#' regression is applied to initialize coefficients corresponding to the
#' parameters that are subject to selection, while ordinary least squares is
#' used for the parameters that are not subject to selection.
#'
#' @param est_indiv A numeric matrix of dimension \eqn{n \times q}, where
#' \eqn{n} is the number of individuals and \eqn{q} the number of
#' \eqn{\varphi} parameters. Typically obtained from
#' \code{estimate_phi_individuals}.
#' @param data_processed An object of class \code{saemvsDataProcessed} (or
#' similar) containing pre-processed covariates for \eqn{\varphi} parameters
#' subject and not suject to selection. Expected to provide at least the slots
#' \code{x_phi_to_select} and \code{x_phi_not_to_select}.
#' @param model An object of class \code{saemvsModel} representing the model
#'   structure, with slots:
#'   \describe{
#'     \item{\code{phi_to_select_idx}}{Indices of \eqn{\varphi} parameters on
#' which to apply the Lasso selection.}
#'     \item{\code{x_forced_support}}{Support (binary) matrix indicating the
#' position of forced covariates.}
#'   }
#' @param init An object of class \code{saemvsInit} providing initial values
#'   for the intercept, covariance matrix, and residual variance.
#' @param lasso_lambda Optional numeric value specifying the lambda to use for
#'   Lasso regression. If \code{NULL}, the minimum lambda from the Lasso path
#'   is used.
#'
#' @details
#' For parameters designated for selection, Lasso regression is applied to
#' the corresponding covariates. For non-selected parameters, ordinary least
#' squares regression is used. Forced covariates indicated by the model are
#' incorporated into the \code{beta_forced} matrix, and the remaining
#' coefficients are stored in \code{beta_candidates}.
#'
#' The resulting \code{saemvsInit} object contains:
#' \itemize{
#'   \item \code{intercept}: initial intercept values.
#'   \item \code{beta_candidates}: candidate coefficients for selection.
#'   \item \code{beta_forced}: coefficients for forced covariates.
#'   \item \code{cov_re}: covariance matrix of random effects.
#'   \item \code{sigma2}: residual variance.
#' }
#'
#' @return An object of class \code{saemvsInit} containing initial values
#'   for the SAEM-VS algorithm.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' init_obj <- build_init_from_phi_lasso(
#'   est_indiv = phi_est_matrix,
#'   data_processed = processed_data,
#'   model = my_model,
#'   init = init_values
#' )
#' }
build_init_from_phi_lasso <- function(est_indiv,
                                      data_processed,
                                      model,
                                      init,
                                      lasso_lambda = NULL) {
  phi_est_matrix <- est_indiv

  beta_candidates <- NULL
  beta_forced <- NULL

  n_phi <- ncol(phi_est_matrix)

  phi_to_select_idx <- model@phi_to_select_idx
  phi_not_to_select_idx <- setdiff(seq_len(n_phi), phi_to_select_idx)

  beta_sel <- NULL

  if (
    length(phi_to_select_idx) > 0 && !is.null(data_processed@x_phi_to_select)
  ) {
    x_sel <- data_processed@x_phi_to_select[, -1, drop = FALSE]
    # remove intercept column
    beta_sel <- matrix(0, ncol = length(phi_to_select_idx), nrow = ncol(x_sel))
    for (i in seq_along(phi_to_select_idx)) {
      phi_sel <- phi_est_matrix[, phi_to_select_idx[i]]
      fit <- glmnet::glmnet(
        x = x_sel, y = phi_sel, alpha = 1,
        standardize = TRUE
      )
      s_val <- if (is.null(lasso_lambda)) min(fit$lambda) else lasso_lambda
      coefs_list <- stats::coef(fit, s = s_val)
      beta_sel[, i] <- as.matrix(coefs_list[-1])
    }
  }


  beta_not_sel <- NULL

  x_not_sel <- data_processed@x_phi_not_to_select[, -1, drop = FALSE]
  # remove intercept column
  if (length(phi_not_to_select_idx) > 0) {
    if (!is_empty_matrix(x_not_sel)) {
      beta_not_sel <- matrix(
        0,
        ncol = length(phi_not_to_select_idx), nrow = ncol(x_not_sel)
      )
      for (i in seq_along(phi_not_to_select_idx)) {
        fit <- stats::lm(
          phi_est_matrix[, phi_not_to_select_idx[i]] ~ x_not_sel
        )
        coefs <- stats::coef(fit)
        beta_not_sel[, i] <- coefs[-1]
      }
    }
  }

  supp <- model@x_forced_support

  if (is_empty_support(supp)) {
    beta_forced <- NULL
    beta_candidates <- matrix(0, nrow = ncol(x_sel), ncol = n_phi)
    if (!is.null(beta_sel)) beta_candidates[, phi_to_select_idx] <- beta_sel
  } else {
    x_forced_idx <- seq(1, nrow(supp))
    beta_forced <- matrix(0, nrow = nrow(supp), ncol = n_phi)
    if (!all(supp[, phi_not_to_select_idx] == 0)) {
      beta_forced[, phi_not_to_select_idx] <- beta_not_sel
    }
    if (!all(supp[, phi_to_select_idx] == 0)) {
      beta_forced[, phi_to_select_idx] <- beta_sel[nrow(supp), ]
    }
    beta_candidates <- matrix(
      0,
      nrow = ncol(x_sel[, -x_forced_idx]), ncol = n_phi
    )
    beta_candidates[, phi_to_select_idx] <- beta_sel[-seq(1, nrow(supp)), ]
  }

  methods::new("saemvsInit",
    intercept = init@intercept,
    beta_candidates = beta_candidates,
    beta_forced = beta_forced,
    cov_re = init@cov_re,
    sigma2 = init@sigma2
  )
}
