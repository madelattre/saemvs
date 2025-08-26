single_iteration <- function(
    k, config, state, sa_func,
    m_func, mh_update_func) {
  state <- s_step(config, k, state)
  state <- sa_func(config, k, state)
  state <- m_func(config, k, state)

  state <- mh_update_func(config, k, state)


  return(state)
}

run_generic_saem <- function(config, state, sa_func, m_func, mh_update_func) {
  for (k in 1:config$num_iterations) {
    state <- single_iteration(k, config, state, sa_func, m_func, mh_update_func)
  }
  return(state)
}

run_fixed_saem <- function(config, state, sa_func, m_func, mh_update_func) {
  for (k in 1:config$num_burnin) {
    state <- single_iteration(k, config, state, sa_func, m_func, mh_update_func)
  }

  index_fixed <- which(config$parameters_not_to_select_indices %in% config$fixed_parameters_indices)

  alphas <- exp((-2 * log(10) -
    log(diag(state$gamma_ldim[[config$num_burnin]])[index_fixed])) /
    (config$num_iterations - config$num_burnin))

  for (k in (config$num_burnin + 1):config$num_iterations) {
    state <- single_iteration(k, config, state, sa_func, m_func, mh_update_func)
    state$gamma_ldim[[k + 1]] <- shrink_covariance_matrix(
      state$gamma_ldim[[k]], state$gamma_ldim[[k + 1]],
      index_fixed, alphas
    )
  }

  shrinked_gamma <- lapply(seq_len(config$num_iterations + 1), function(k) {
    zero_out_shrinked(state$gamma_ldim[[k]], index_fixed)
  })

  return(state)
}
