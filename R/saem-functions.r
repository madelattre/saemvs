#' Perform a single SAEM iteration
#'
#' Executes one iteration of the SAEM algorithm, including the Metropolis step,
#'              stochastic approximation (SA) step, M-step, and MH update.
#' @param k Integer. Current iteration index.
#' @param config List. Model configuration including design matrices,
#' hyperparameters, and algorithm settings.
#' @param state List. Current state of the algorithm (phi, beta, gamma, alpha,
#' etc.).
#' @param sa_func Function. Stochastic approximation step function (SA-step) to
#' update sufficient statistics.
#' @param m_func Function. M-step function to update model parameters.
#' @param mh_update_func Function. Optional function to update proposal
#' distributions for the Metropolis step.
#' @return Updated \code{state} list after one SAEM iteration.
#' @keywords internal
perform_saem_iteration <- function(k, config, state, sa_func, m_func,
                                   mh_update_func) {
  state <- metropolis_s_step(config, k, state)
  state <- sa_func(config, k, state)
  state <- m_func(config, k, state)
  state <- mh_update_func(config, k, state)
  return(state)
}

#' Run full SAEM algorithm
#'
#' Executes the full SAEM algorithm for all iterations, applying the specified
#'              SA-step, M-step, and Metropolis-Hastings update functions.
#' @param config List. Model configuration including number of iterations and
#' algorithm settings.
#' @param state List. Initial state of the algorithm.
#' @param sa_func Function. Stochastic approximation step function (SA-step) to
#' update sufficient statistics.
#' @param m_func Function. M-step function to update model parameters.
#' @param mh_update_func Function. Metropolis-Hastings update function.
#' @return Final \code{state} list after completing all SAEM iterations.
#' @keywords internal
run_saem_full <- function(config, state, sa_func, m_func, mh_update_func) {
  for (k in 1:config$num_iterations) {
    state <- perform_saem_iteration(
      k, config, state, sa_func, m_func,
      mh_update_func
    )
  }
  return(state)
}

#' Run SAEM algorithm with fixed parameters
#'
#' Executes the SAEM algorithm for models with both random and fixed effects.
#' Shrinkage is applied to the covariance of fixed parameters after burn-in
#' iterations.
#' @param config List. Model configuration including number of iterations,
#' burn-in period, and fixed parameter indices.
#' @param state List. Initial state of the algorithm.
#' @param sa_func Function. Stochastic approximation step function (SA-step) to
#' update sufficient statistics.
#' @param m_func Function. M-step function to update model parameters.
#' @param mh_update_func Function. Metropolis-Hastings update function.
#' @return Final \code{state} list after completing all SAEM iterations with
#' shrinkage for fixed effects.
#' @keywords internal
run_saem_with_fixed_effects <- function(
    config, state, sa_func, m_func, mh_update_func) {
  for (k in 1:config$num_burnin) {
    state <- perform_saem_iteration(
      k, config, state, sa_func, m_func,
      mh_update_func
    )
  }

  index_fixed <- which(
    config$parameters_not_to_select_indices %in%
      config$fixed_parameters_indices
  )
  alphas <-
    exp(
      (-2 * log(10) -
        log(diag(state$gamma_not_to_select[[config$num_burnin]])[index_fixed])
      ) / (config$num_iterations - config$num_burnin)
    )

  for (k in (config$num_burnin + 1):config$num_iterations) {
    state <- perform_saem_iteration(
      k, config, state, sa_func, m_func,
      mh_update_func
    )
    state$gamma_not_to_select[[k + 1]] <- shrink_covariance_matrix(
      state$gamma_not_to_select[[k]],
      state$gamma_not_to_select[[k + 1]],
      index_fixed,
      alphas
    )
  }

#  shrinked_gamma <- lapply(seq_len(config$num_iterations + 1), function(k) {
#    zero_out_shrinked(state$gamma_not_to_select[[k]], index_fixed)
#  })

  return(state)
}
