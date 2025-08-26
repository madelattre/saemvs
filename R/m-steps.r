#' Internal M-step: update parameters not subject to selection
#' 
#' Performs the M-step updates for parameters that are not subject to variable selection.
#'              Updates beta, gamma, and related sufficient statistics for the "not_to_select" parameters.
#' @param config List containing model configuration, including design matrices and hyperparameters.
#' @param k Current iteration index (integer).
#' @param state List containing the current state of the algorithm (phi, beta, gamma, s1, s2, s3, etc.).
#' @return Updated \code{state} list with new beta_not_to_select, gamma_not_to_select, and s1/s2/s3 entries.
#' @keywords internal
m_step_not_to_select <- function(config, k, state) {
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

  state$gamma_not_to_select[[k + 1]] <- if (use_scaled) gamma_scaled else gamma_proposed

  return(state)
}

#' Internal M-step: update parameters subject to selection
#' 
#' Performs the M-step updates for parameters that are subject to variable selection.
#'              Updates beta, gamma, and inclusion probabilities (alpha) for the "to_select" parameters.
#' @param config List containing model configuration, including design matrices and hyperparameters.
#' @param k Current iteration index (integer).
#' @param state List containing the current state of the algorithm (phi, beta, gamma, alpha, s1, s2, s3, etc.).
#' @return Updated \code{state} list with new beta_to_select, gamma_to_select, alpha, and related statistics.
#' @keywords internal
m_step_to_select <- function(config, k, state) {
  old_beta <- state$beta_to_select[[k]]
  old_gamma <- state$gamma_to_select[[k]]
  old_alpha <- state$alpha[[k]]
  s2 <- state$s2_to_select[[k + 1]]
  s3 <- state$s3_to_select[[k + 1]]

  p_star_k <- p_star(
    old_beta[-1, ],
    old_alpha,
    config$spike_parameter,
    config$slab_parameter
  )

  d_tilde_star <- rbind(
    rep(1 / config$phi_intercept_prior_variance, config$num_parameters_to_select),
    (1 - p_star_k) / config$spike_parameter + p_star_k / config$slab_parameter
  )

  sum_pstar_k <- apply(p_star_k, 2, sum)

  state$alpha[[k + 1]] <- (sum_pstar_k + config$inclusion_prob_prior_a - 1) /
    (config$num_covariates_to_select + config$inclusion_prob_prior_b + config$inclusion_prob_prior_a - 2)

  s_xgx <- config$kron_tx_x_phi_to_select +
    kronecker_gamma_diag_mult(old_gamma, c(d_tilde_star), config$num_covariates_to_select + 1)

  s_xgs3 <- matrix(t(config$x_phi_to_select) %*% s3, ncol = 1)

  beta_vector <- solve_linear_syst(s_xgx, s_xgs3)

  state$beta_to_select[[k + 1]] <- matrix(
    beta_vector,
    nrow = config$num_covariates_to_select + 1,
    ncol = config$num_parameters_to_select,
    byrow = FALSE
  )

  sum_bx <- t(state$beta_to_select[[k + 1]]) %*% config$tx_x_phi_to_select %*% state$beta_to_select[[k + 1]]
  sum_xbs3 <- t(config$x_phi_to_select %*% state$beta_to_select[[k + 1]]) %*% s3

  gamma_proposed <- (config$cov_re_prior_scale + s2 + sum_bx - sum_xbs3 - t(sum_xbs3)) /
    (config$num_series + config$cov_re_prior_df + config$num_parameters_to_select + 1)
  gamma_proposed <- (gamma_proposed + t(gamma_proposed)) / 2 # enforce symmetry
  gamma_scaled <- config$covariance_decay * old_gamma
  use_scaled <- sum(gamma_scaled^2) > sum(gamma_proposed^2)

  state$gamma_to_select[[k + 1]] <- if (use_scaled) gamma_scaled else gamma_proposed

  return(state)
}

#' Internal M-step: update MLE parameters
#' 
#' Updates the residual variance (sigma2) and parameters not subject to selection (MLE case).
#' @param config List containing model configuration and hyperparameters.
#' @param k Current iteration index (integer).
#' @param state List containing the current state of the algorithm.
#' @return Updated \code{state} list with new sigma2 and beta/gamma for non-selected parameters.
#' @keywords internal
m_step_mle <- function(config, k, state) {
  state <- m_step_not_to_select(config, k, state)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$covariance_decay * old_sigma2,
    state$s1[k + 1] / config$total_observations
  )
  return(state)
}

#' Internal M-step: update MAP estimates for parameters subject to selection
#' 
#' Performs MAP update for parameters subject to selection, including beta, gamma, alpha, and sigma2.
#' @param config List containing model configuration and hyperparameters.
#' @param k Current iteration index (integer).
#' @param state List containing the current state of the algorithm.
#' @return Updated \code{state} list with new beta_to_select, gamma_to_select, alpha, and sigma2.
#' @keywords internal
m_step_map_to_select <- function(config, k, state) {
  state <- m_step_to_select(config, k, state)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$covariance_decay * old_sigma2,
    (config$residual_variance_prior_shape * config$residual_variance_prior_rate + state$s1[k + 1]) /
      (config$total_observations + config$residual_variance_prior_rate + 2)
  )
  return(state)
}

#' Internal M-step: update MAP for all parameters
#' 
#' @description Performs a full MAP update for both parameters subject and not subject to selection,
#'              combining the updates for beta, gamma, alpha, and sigma2.
#' @param config List containing model configuration and hyperparameters.
#' @param k Current iteration index (integer).
#' @param state List containing the current state of the algorithm.
#' @return Updated \code{state} list with new beta, gamma, alpha, and sigma2 for all parameters.
#' @keywords internal
m_step_map_all <- function(config, k, state) {
  state <- m_step_to_select(config, k, state)
  state <- m_step_not_to_select(config, k, state)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$covariance_decay * old_sigma2,
    (config$residual_variance_prior_shape * config$residual_variance_prior_rate + state$s1[k + 1]) /
      (config$total_observations + config$residual_variance_prior_rate + 2)
  )
  return(state)
}