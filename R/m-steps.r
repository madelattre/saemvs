m_step_ldim <- function(config, k, state) {
  old_gamma <- state$gamma_not_to_select[[k]]
  s3 <- state$s3_not_to_select[[k + 1]]
  s2 <- state$s2_not_to_select[[k + 1]]


  gamma_inv <- mat_inv(old_gamma)

  s_xgx <- Reduce(`+`, lapply(
    config$x_phi_not_to_select_list, function(xi) t(xi) %*% gamma_inv %*% xi
  ))

  s_xgs3 <- Reduce(`+`, lapply(1:config$num_series, function(i) {
    t(config$x_phi_not_to_select_list[[i]]) %*% gamma_inv %*% s3[i, ]
  }))

  beta_vect <- solve_linear_syst(s_xgx, s_xgs3)
  cbeta_outer <- beta_vect %*% t(beta_vect)
  state$beta_mat_0[config$forced_covariates_indices] <- beta_vect

  state$beta_not_to_select[[k + 1]] <-
    matrix(state$beta_mat_0,
      ncol = config$num_parameters_not_to_select,
      nrow = config$num_covariates_not_to_select + 1, byrow = FALSE
    )


  sum_bx <- Reduce(`+`, lapply(config$x_phi_not_to_select_list, function(xi) {
    xi %*% cbeta_outer %*% t(xi)
  }))

  sum_xbs3 <- Reduce(`+`, lapply(1:config$num_series, function(i) {
    config$x_phi_not_to_select_list[[i]] %*% beta_vect %*% t(s3[i, ])
  }))

  gamma_prop <- (s2 + sum_bx - sum_xbs3 - t(sum_xbs3)) / config$num_series
  gamma_prop <- (gamma_prop + t(gamma_prop)) / 2 # Forcing symmetry
  gamma_scaled <- config$covariance_decay * old_gamma
  cond_gamma <- (sum(gamma_scaled^2) > sum(gamma_prop^2))
  gamma_scaled <- config$covariance_decay * old_gamma
  state$gamma_not_to_select[[k + 1]] <- if (cond_gamma) gamma_scaled else gamma_prop

  return(state)
}

# Update coefficients and covariance matrix for the high-dimensional part of
# the model
m_step_hdim <- function(config, k, state) {
  old_beta <- state$beta_to_select[[k]]
  old_gamma <- state$gamma_to_select[[k]]
  old_alpha <- state$alpha[[k]]
  s2 <- state$s2_to_select[[k + 1]]
  s3 <- state$s3_to_select[[k + 1]]

  p_star_k <- p_star_fast(
    old_beta[-1, ],
    old_alpha, config$spike_parameter, config$slab_parameter
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


  beta_vect <- solve_linear_syst(s_xgx, s_xgs3)

  state$beta_to_select[[k + 1]] <-
    matrix(beta_vect, nrow = config$num_covariates_to_select + 1, ncol = config$num_parameters_to_select, byrow = FALSE)

  sum_bx <- t(state$beta_to_select[[k + 1]]) %*% config$tx_x_phi_to_select %*%
    state$beta_to_select[[k + 1]]

  sum_xbs3 <- t(config$x_phi_to_select %*% state$beta_to_select[[k + 1]]) %*% s3

  gamma_prop <- (config$cov_re_prior_scale + s2 + sum_bx - sum_xbs3 -
    t(sum_xbs3)) / (config$num_series + config$cov_re_prior_df + config$num_parameters_to_select + 1)
  gamma_prop <- (gamma_prop + t(gamma_prop)) / 2 # Forcing symmetry
  gamma_scaled <- config$covariance_decay * old_gamma
  cond_gamma <- (sum(gamma_scaled^2) > sum(gamma_prop^2))

  state$gamma_to_select[[k + 1]] <- if (cond_gamma) gamma_scaled else gamma_prop

  return(state)
}

m_step_mle <- function(config, k, state) {
  state <- m_step_ldim(config, k, state)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$covariance_decay * old_sigma2,
    state$s1[k + 1] / (config$total_observations)
  )
  return(state)
}

m_step_full_map <- function(config, k, state) {
  state <- m_step_hdim(config, k, state)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$covariance_decay * old_sigma2,
    (config$residual_variance_prior_shape * config$residual_variance_prior_rate + state$s1[k + 1]) /
      (config$total_observations + config$residual_variance_prior_rate + 2)
  )

  return(state)
}

m_step_mix <- function(config, k, state) {
  state <- m_step_hdim(config, k, state)
  state <- m_step_ldim(config, k, state)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$covariance_decay * old_sigma2,
    (config$residual_variance_prior_shape * config$residual_variance_prior_rate + state$s1[k + 1]) /
      (config$total_observations + config$residual_variance_prior_rate + 2)
  )

  return(state)
}
