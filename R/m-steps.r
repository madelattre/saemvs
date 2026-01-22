#' Internal M-step: update parameters not subject to selection
#'
#' Performs the M-step updates for parameters that are not subject to variable
#' selection. Updates regression coefficients (`beta_not_to_select`), covariance
#' matrices (`gamma_not_to_select`), and related sufficient statistics for the
#' "not_to_select" parameters.
#'
#' @param config List containing model configuration, including:
#'   \itemize{
#'     \item \code{x_phi_not_to_select_list}: list of design matrices for each
#'  series,
#'     \item \code{num_series}: number of independent series,
#'     \item \code{num_parameters_not_to_select}: number of parameters not
#'  subject to selection,
#'     \item \code{num_covariates_not_to_select}: number of covariates for
#'  not-to-select parameters,
#'     \item \code{forced_covariates_indices}: indices of covariates that are
#'  forced into the model,
#'     \item \code{covariance_decay}: scaling factor used for shrinkage of
#'  gamma matrices.
#'   }
#' @param k Integer, current iteration index (1-based), used to index and update
#'   \code{state$beta_not_to_select} and \code{state$gamma_not_to_select}.
#' @param state List containing the current MCMC state, including:
#'   \itemize{
#'     \item \code{phi}: list of latent parameter matrices per iteration,
#'     \item \code{s2_not_to_select}, \code{s3_not_to_select}: stochastic
#' approximation statistics,
#'     \item \code{beta_mat_0}: current vector or matrix of regression
#'  coefficients,
#'     \item \code{beta_not_to_select}: list of beta matrices for not-to-select
#'  parameters,
#'     \item \code{gamma_not_to_select}: list of covariance matrices for
#'  not-to-select parameters.
#'   }
#'   These components are read and updated in the M-step for parameters not
#'  subject to selection.
#' @param backend A list of compiled model functions,
#'   passed for interface consistency with the SAEM iteration framework.
#'
#' @return Updated \code{state} list with:
#'   \itemize{
#'     \item \code{beta_not_to_select[[k + 1]]}: updated regression
#'  coefficients,
#'     \item \code{gamma_not_to_select[[k + 1]]}: updated covariance matrix,
#'     \item other components of \code{state} remain unchanged.
#'   }
#'
#' @keywords internal
#' @noRd
#' @note This function performs the M-step for parameters not subject to
#' selection, enforcing symmetry of covariance matrices and applying shrinkage
#' if necessary.
#' It is intended for internal use within the SAEM-MCMC algorithm.

m_step_not_to_select <- function(config, k, state, backend) {
  old_gamma <- state$gamma_not_to_select[[k]]
  s3 <- state$s3_not_to_select[[k + 1]]
  s2 <- state$s2_not_to_select[[k + 1]]

  gamma_inv <- mat_inv(old_gamma)

  sum_xgx <- Reduce(`+`, lapply(
    config$x_phi_not_to_select_list, function(xi) t(xi) %*% gamma_inv %*% xi
  ))

  sum_xgs3 <- Reduce(`+`, lapply(1:config$num_series, function(i) {
    t(config$x_phi_not_to_select_list[[i]]) %*% gamma_inv %*% s3[i, ]
  }))

  beta_vector <- solve_linear_syst(sum_xgx, sum_xgs3)
  beta_outer <- beta_vector %*% t(beta_vector)
  state$beta_mat_0[config$forced_covariates_indices] <- beta_vector

  state$beta_not_to_select[[k + 1]] <- matrix(
    state$beta_mat_0,
    ncol = config$num_parameters_not_to_select,
    nrow = config$num_covariates_not_to_select + 1,
    byrow = FALSE
  )

  sum_bx <- Reduce(`+`, lapply(config$x_phi_not_to_select_list, function(xi) {
    xi %*% beta_outer %*% t(xi)
  }))

  sum_xbs3 <- Reduce(`+`, lapply(1:config$num_series, function(i) {
    config$x_phi_not_to_select_list[[i]] %*% beta_vector %*% t(s3[i, ])
  }))

  gamma_proposed <- (s2 + sum_bx - sum_xbs3 - t(sum_xbs3)) / config$num_series
  gamma_proposed <- (gamma_proposed + t(gamma_proposed)) / 2 # enforce symmetry
  gamma_scaled <- config$covariance_decay * old_gamma

  use_scaled <- sum(gamma_scaled^2) > sum(gamma_proposed^2)

  state$gamma_not_to_select[[k + 1]] <- if (use_scaled) {
    gamma_scaled
  } else {
    gamma_proposed
  }

  return(state) # nolint : return-linter
}

#' Internal M-step: update parameters subject to selection
#'
#' Performs the M-step updates for parameters that are subject to variable
#' selection. Updates regression coefficients (`beta_to_select`), covariance
#' matrices (`gamma_to_select`), and inclusion probabilities (`alpha`) for the
#' "to_select" parameters.
#'
#' @param config List containing model configuration, including:
#'   \itemize{
#'     \item \code{x_phi_to_select}: design matrix for parameters subject to
#'  selection,
#'     \item \code{tx_x_phi_to_select}, \code{kron_tx_x_phi_to_select}:
#' precomputed matrices for M-step,
#'     \item \code{num_series}: number of independent series,
#'     \item \code{num_covariates_to_select}: number of covariates for
#'  to-select parameters,
#'     \item \code{num_parameters_to_select}: number of parameters subject to
#'  selection,
#'     \item \code{phi_intercept_prior_variance}: prior variance for intercept
#'  terms,
#'     \item \code{spike_parameter}, \code{slab_parameter}: spike-and-slab
#'  prior hyperparameters,
#'     \item \code{inclusion_prob_prior_a}, \code{inclusion_prob_prior_b}:
#' hyperparameters for inclusion probabilities,
#'     \item \code{cov_re_prior_scale}, \code{cov_re_prior_df}: prior parameters
#'  for gamma,
#'     \item \code{covariance_decay}: scaling factor used for shrinkage of
#'  gamma matrices.
#'   }
#'
#' @param k Integer, current iteration index (1-based), used to index and update
#'   \code{state$beta_to_select}, \code{state$gamma_to_select}, and
#'  \code{state$alpha}.
#'
#' @param state List containing the current MCMC state, including:
#'   \itemize{
#'     \item \code{phi}: list of latent parameter matrices per iteration,
#'     \item \code{s2_to_select}, \code{s3_to_select}: stochastic approximation
#'  statistics for to-select parameters,
#'     \item \code{beta_to_select}: list of regression coefficient matrices
#'  for parameters to select,
#'     \item \code{gamma_to_select}: list of covariance matrices for parameters
#'  to select,
#'     \item \code{alpha}: vector of inclusion probabilities per iteration.
#'   }
#'   These components are read and updated in the M-step for parameters subject
#'  to selection.
#'
#' @param backend A list of compiled model functions and helper routines,
#'   passed for interface consistency with the SAEM iteration framework.
#'
#' @return Updated \code{state} list with:
#'   \itemize{
#'     \item \code{beta_to_select[[k + 1]]}: updated regression coefficients,
#'     \item \code{gamma_to_select[[k + 1]]}: updated covariance matrix
#'  (symmetric, possibly shrinked),
#'     \item \code{alpha[[k + 1]]}: updated inclusion probabilities,
#'     \item other components of \code{state} remain unchanged.
#'   }
#'
#' @keywords internal
#' @noRd
#' @note This function performs the M-step for parameters subject to variable
#'  selection,
#' enforcing symmetry of covariance matrices and applying shrinkage if
#'  necessary.
#' Inclusion probabilities are updated according to a Beta prior and
#'  spike-and-slab posterior.
#' It is intended for internal use within the SAEM-MCMC algorithm.
m_step_to_select <- function(config, k, state, backend) {
  old_beta <- state$beta_to_select[[k]]
  old_gamma <- state$gamma_to_select[[k]]
  old_alpha <- state$alpha[[k]]
  s2 <- state$s2_to_select[[k + 1]]
  s3 <- state$s3_to_select[[k + 1]]

  p_star_k <- p_star(
    old_beta[-1, , drop = FALSE],
    old_alpha,
    config$spike_parameter,
    config$slab_parameter
  )

  d_tilde_star <- rbind(
    rep(
      1 / config$phi_intercept_prior_variance,
      config$num_parameters_to_select
    ),
    (1 - p_star_k) / config$spike_parameter + p_star_k / config$slab_parameter
  )

  sum_pstar_k <- apply(p_star_k, 2, sum)

  state$alpha[[k + 1]] <- (sum_pstar_k + config$inclusion_prob_prior_a - 1) /
    (config$num_covariates_to_select + config$inclusion_prob_prior_b +
       config$inclusion_prob_prior_a - 2)

  s_xgx <- config$kron_tx_x_phi_to_select +
    kronecker_gamma_diag_mult(
      old_gamma, c(d_tilde_star),
      config$num_covariates_to_select + 1
    )

  s_xgs3 <- matrix(t(config$x_phi_to_select) %*% s3, ncol = 1)

  beta_vector <- solve_linear_syst(s_xgx, s_xgs3)

  state$beta_to_select[[k + 1]] <- matrix(
    beta_vector,
    nrow = config$num_covariates_to_select + 1,
    ncol = config$num_parameters_to_select,
    byrow = FALSE
  )

  sum_bx <- t(state$beta_to_select[[k + 1]]) %*% config$tx_x_phi_to_select %*%
    state$beta_to_select[[k + 1]]
  sum_xbs3 <- t(config$x_phi_to_select %*% state$beta_to_select[[k + 1]]) %*% s3

  gamma_proposed <- (config$cov_re_prior_scale + s2 + sum_bx - sum_xbs3 -
                       t(sum_xbs3)) /
    (config$num_series + config$cov_re_prior_df +
       config$num_parameters_to_select + 1)
  gamma_proposed <- (gamma_proposed + t(gamma_proposed)) / 2 # enforce symmetry
  gamma_scaled <- config$covariance_decay * old_gamma

  use_scaled <- sum(gamma_scaled^2) > sum(gamma_proposed^2)

  state$gamma_to_select[[k + 1]] <- if (use_scaled) {
    gamma_scaled
  } else {
    gamma_proposed
  }
  return(state) # nolint: return-linter
}

#' Internal M-step: update MLE parameters
#'
#' Performs the M-step updates for Maximum Likelihood Estimation (MLE) of
#'  parameters.
#' Updates regression coefficients, covariance matrices, and residual
#'  variance (\code{sigma2})
#' for parameters not subject to selection.
#'
#' @param config List containing model configuration and hyperparameters,
#'  including:
#'   \itemize{
#'     \item \code{covariance_decay}: scaling factor for shrinkage of
#'  covariance matrices,
#'     \item \code{total_observations}: total number of observations for
#'  computing residual variance.
#'   }
#'   Other fields required by \code{m_step_not_to_select} are also expected.
#' @param k Integer, current iteration index (1-based), used to index and update
#'   \code{state$beta_not_to_select}, \code{state$gamma_not_to_select},
#'  and \code{state$sigma2}.
#' @param state List containing the current MCMC state, including:
#'   \itemize{
#'     \item \code{beta_not_to_select}, \code{gamma_not_to_select}:
#'  parameters not subject to selection,
#'     \item \code{s1}: stochastic approximation statistic used to compute
#'  residual variance,
#'     \item \code{sigma2}: residual variance per iteration.
#'   }
#' @param backend A list of compiled model functions passed for interface
#'  consistency.
#'
#' @return Updated \code{state} list with:
#'   \itemize{
#'     \item \code{beta_not_to_select[[k + 1]]},
#'  \code{gamma_not_to_select[[k + 1]]}: updated MLE parameters,
#'     \item \code{sigma2[k + 1]}: updated residual variance,
#'     \item other components of \code{state} remain unchanged.
#'   }
#'
#' @keywords internal
#' @noRd
#' @note Intended for internal use in the SAEM-MCMC framework when computing
#'  MLE estimates.
m_step_mle <- function(config, k, state, backend) {
  state <- m_step_not_to_select(config, k, state, backend)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$covariance_decay * old_sigma2,
    state$s1[k + 1] / config$total_observations
  )
  return(state) # nolint: return-linter
}


#' Internal M-step: update MAP estimates for parameters subject to selection
#'
#' Performs the M-step updates for Maximum A Posteriori (MAP) estimation of
#'  parameters
#' subject to variable selection. Updates regression coefficients
#'  (\code{beta_to_select}),
#' covariance matrices (\code{gamma_to_select}), inclusion probabilities
#'  (\code{alpha}),
#' and residual variance (\code{sigma2}).
#'
#' @param config List containing model configuration and hyperparameters,
#'  including:
#'   \itemize{
#'     \item \code{covariance_decay}: scaling factor for shrinkage of
#'  covariance matrices,
#'     \item \code{residual_variance_prior_shape},
#' \code{residual_variance_prior_rate}: prior parameters for residual variance,
#'     \item \code{total_observations}: total number of observations for
#'  computing residual variance.
#'   }
#'   Other fields required by \code{m_step_to_select} are also expected.
#' @param k Integer, current iteration index (1-based), used to index and update
#'   \code{state$beta_to_select}, \code{state$gamma_to_select},
#'  \code{state$alpha}, and \code{state$sigma2}.
#' @param state List containing the current MCMC state, including:
#'   \itemize{
#'     \item \code{beta_to_select}, \code{gamma_to_select}: parameters subject
#'  to selection,
#'     \item \code{alpha}: inclusion probabilities,
#'     \item \code{s1}: stochastic approximation statistic used to compute
#'  residual variance,
#'     \item \code{sigma2}: residual variance per iteration.
#'   }
#' @param backend A list of compiled model functions passed
#'  for interface consistency.
#'
#' @return Updated \code{state} list with:
#'   \itemize{
#'     \item \code{beta_to_select[[k + 1]]}, \code{gamma_to_select[[k + 1]]},
#'  \code{alpha[[k + 1]]}: updated MAP parameters,
#'     \item \code{sigma2[k + 1]}: updated residual variance,
#'     \item other components of \code{state} remain unchanged.
#'   }
#'
#' @keywords internal
#' @noRd
#' @note Intended for internal use in the SAEM-MCMC framework when computing
#'  MAP estimates.

m_step_map_to_select <- function(config, k, state, backend) {
  state <- m_step_to_select(config, k, state, backend)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$covariance_decay * old_sigma2,
    (config$residual_variance_prior_shape *
       config$residual_variance_prior_rate +
       state$s1[k + 1]) /
      (config$total_observations + config$residual_variance_prior_shape + 2)
  )
  return(state) # nolint: return-linter
}

#' Internal M-step: update MAP for all parameters
#'
#' Performs a full Maximum A Posteriori (MAP) update for both parameters
#' subject to selection (`to_select`) and not subject to selection
#'  (`not_to_select`).
#' Updates regression coefficients (`beta`), covariance matrices (`gamma`),
#'  inclusion probabilities (`alpha`), and residual variance (`sigma2`) for
#'  all parameters.
#'
#' @param config List containing model configuration and hyperparameters,
#'  including:
#'   \itemize{
#'     \item \code{covariance_decay}: scaling factor used for shrinkage of
#'  covariance matrices,
#'     \item \code{residual_variance_prior_shape},
#'  \code{residual_variance_prior_rate}: prior parameters for residual variance,
#'     \item \code{total_observations}: total number of observations,
#'     \item other fields required by \code{m_step_to_select} and
#'  \code{m_step_not_to_select}.
#'   }
#'
#' @param k Integer, current iteration index (1-based), used to index and update
#'   all relevant \code{state} components.
#'
#' @param state List containing the current MCMC state, including:
#'   \itemize{
#'     \item \code{beta_to_select}, \code{gamma_to_select}, \code{alpha}:
#'  parameters subject to selection,
#'     \item \code{beta_not_to_select}, \code{gamma_not_to_select}:
#'  parameters not subject to selection,
#'     \item \code{s1}: stochastic approximation statistic used to compute
#'  residual variance,
#'     \item \code{sigma2}: residual variance per iteration.
#'   }
#'
#' @param backend A list of compiled model functions passed
#'  for interface consistency.
#'
#' @return Updated \code{state} list with:
#'   \itemize{
#'     \item \code{beta_to_select[[k + 1]]}, \code{gamma_to_select[[k + 1]]},
#' \code{alpha[[k + 1]]}: updated MAP parameters for to-select,
#'     \item \code{beta_not_to_select[[k + 1]]},
#'  \code{gamma_not_to_select[[k + 1]]}: updated MAP parameters for
#'  not-to-select,
#'     \item \code{sigma2[k + 1]}: updated residual variance,
#'     \item other components of \code{state} remain unchanged.
#'   }
#'
#' @keywords internal
#' @noRd
#' @note Intended for internal use in the SAEM-MCMC framework when computing MAP
#'  estimates for all parameters simultaneously.
m_step_map_all <- function(config, k, state, backend) {
  state <- m_step_to_select(config, k, state, backend)
  state <- m_step_not_to_select(config, k, state, backend)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$covariance_decay * old_sigma2,
    (config$residual_variance_prior_shape *
       config$residual_variance_prior_rate +
       state$s1[k + 1]) /
      (config$total_observations + config$residual_variance_prior_shape + 2)
  )
  return(state) # nolint: return-linter
}
