m_step_ldim <- function(config, k, state) {
  old_gamma <- state$gamma_ldim[[k]]
  s3 <- state$s3_ldim[[k + 1]]
  s2 <- state$s2_ldim[[k + 1]]


  gamma_inv <- mat_inv(old_gamma)

  s_xgx <- Reduce(`+`, lapply(
    config$x, function(xi) t(xi) %*% gamma_inv %*% xi
  ))

  s_xgs3 <- Reduce(`+`, lapply(1:config$n, function(i) {
    t(config$x[[i]]) %*% gamma_inv %*% s3[i, ]
  }))

  beta_vect <- solve_linear_syst(s_xgx, s_xgs3)
  cbeta_outer <- beta_vect %*% t(beta_vect)
  state$beta_mat_0[config$supp_index] <- beta_vect

  state$beta_ldim[[k + 1]] <-
    matrix(state$beta_mat_0,
      ncol = config$q_ldim,
      nrow = config$pw + 1, byrow = FALSE
    )


  sum_bx <- Reduce(`+`, lapply(config$x, function(xi) {
    xi %*% cbeta_outer %*% t(xi)
  }))

  sum_xbs3 <- Reduce(`+`, lapply(1:config$n, function(i) {
    config$x[[i]] %*% beta_vect %*% t(s3[i, ])
  }))

  gamma_prop <- (s2 + sum_bx - sum_xbs3 - t(sum_xbs3)) / config$n
  gamma_prop <- (gamma_prop + t(gamma_prop)) / 2 # Forcing symmetry
  gamma_scaled <- config$tau * old_gamma
  cond_gamma <- (sum(gamma_scaled^2) > sum(gamma_prop^2))
  gamma_scaled <- config$tau * old_gamma
  state$gamma_ldim[[k + 1]] <- if (cond_gamma) gamma_scaled else gamma_prop

  return(state)
}

# Update coefficients and covariance matrix for the high-dimensional part of
# the model
m_step_hdim <- function(config, k, state) {
  old_beta <- state$beta_hdim[[k]]
  old_gamma <- state$gamma_hdim[[k]]
  old_alpha <- state$alpha[[k]]
  s2 <- state$s2_hdim[[k + 1]]
  s3 <- state$s3_hdim[[k + 1]]

  p_star_k <- p_star_fast(
    old_beta[-1, ],
    old_alpha, config$nu0, config$nu1
  )

  d_tilde_star <- rbind(
    rep(1 / config$sigma2_mu, config$q_hdim),
    (1 - p_star_k) / config$nu0 + p_star_k / config$nu1
  )

  sum_pstar_k <- apply(p_star_k, 2, sum)

  state$alpha[[k + 1]] <- (sum_pstar_k + config$a - 1) /
    (config$pv + config$b + config$a - 2)

  s_xgx <- config$kron_tv_v +
    kronecker_gamma_diag_mult(old_gamma, c(d_tilde_star), config$pv + 1)

  s_xgs3 <- matrix(t(config$v) %*% s3, ncol = 1)


  beta_vect <- solve_linear_syst(s_xgx, s_xgs3)

  state$beta_hdim[[k + 1]] <-
    matrix(beta_vect, nrow = config$pv + 1, ncol = config$q_hdim, byrow = FALSE)

  sum_bx <- t(state$beta_hdim[[k + 1]]) %*% config$tv_v %*%
    state$beta_hdim[[k + 1]]

  sum_xbs3 <- t(config$v %*% state$beta_hdim[[k + 1]]) %*% s3

  gamma_prop <- (config$sgam + s2 + sum_bx - sum_xbs3 - t(sum_xbs3)) /
    (config$n + config$d + config$q_hdim + 1)
  gamma_prop <- (gamma_prop + t(gamma_prop)) / 2 # Forcing symmetry
  gamma_scaled <- config$tau * old_gamma
  cond_gamma <- (sum(gamma_scaled^2) > sum(gamma_prop^2))
  gamma_scaled <- config$tau * old_gamma
  state$gamma_hdim[[k + 1]] <- if (cond_gamma) gamma_scaled else gamma_prop

  return(state)
}

m_step_mle <- function(config, k, state) {
  state <- m_step_ldim(config, k, state)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$tau * old_sigma2,
    state$s1[k + 1] / (config$ntot)
  )
  return(state)
}

m_step_full_map <- function(config, k, state) {
  state <- m_step_hdim(config, k, state)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$tau * old_sigma2,
    (config$nsig * config$lsig + state$s1[k + 1]) /
      (config$ntot + config$nsig + 2)
  )

  return(state)
}

m_step_mix <- function(config, k, state) {
  state <- m_step_hdim(config, k, state)
  state <- m_step_ldim(config, k, state)
  old_sigma2 <- state$sigma2[k]
  state$sigma2[k + 1] <- max(
    config$tau * old_sigma2,
    (config$nsig * config$lsig + state$s1[k + 1]) /
      (config$ntot + config$nsig + 2)
  )

  return(state)
}
