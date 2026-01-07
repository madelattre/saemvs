#' Perform a single SAEM iteration
#'
#' Executes one iteration of the SAEM algorithm, sequentially performing
#' the Metropolis simulation step, the stochastic approximation (SA) step,
#' the M-step, and a Metropolis–Hastings proposal update.
#' This function acts as a lightweight orchestrator and does not implement
#' model-specific logic itself; all computations are delegated to the
#' functions passed as arguments.
#' @param k Integer. Current iteration index.
#' @param config List. Model configuration including design matrices,
#' hyperparameters, and algorithm settings.
#' @param state List. Current state of the SAEM algorithm, containing latent
#'   variables, sufficient statistics, and model parameters (e.g. \code{phi},
#'   \code{beta}, \code{gamma}, \code{alpha}). This object is successively
#'   updated and returned by each step of the iteration.
#' @param sa_func Function. Stochastic approximation step function (SA-step) to
#' update sufficient statistics.
#' @param m_func Function. M-step function to update model parameters.
#' @param mh_update_func Function. Function used to update the proposal
#'   distributions for the Metropolis–Hastings step. This function is called
#'   at each iteration and must be provided.
#' @param backend A list containing the compiled model backend used to evaluate
#'   the likelihood and perform model-specific computations during the
#'   Metropolis, SA, and M-steps.
#' @return Updated \code{state} list after one SAEM iteration.
#' @keywords internal
perform_saem_iteration <- function(k, config, state, sa_func, m_func,
                                   mh_update_func, backend) {
  state <- metropolis_s_step(config, k, state, backend)
  state <- sa_func(config, k, state, backend)
  state <- m_func(config, k, state, backend)
  state <- mh_update_func(config, k, state)
  return(state) # nolint: return-linter
}

#' Run full SAEM algorithm
#'
#' At each iteration, this function delegates the actual computation to
#' \code{\link{perform_saem_iteration}}, which sequentially applies the
#' Metropolis, SA, M-step, and Metropolis–Hastings update functions.
#' The number of iterations is fixed in advance, and no convergence or early
#' stopping criterion is evaluated within this function.
#' @param config List. Model configuration object containing algorithm settings,
#'   including the total number of iterations (\code{num_iterations}).
#' @param state List. Initial state of the SAEM algorithm, containing latent
#'   variables, sufficient statistics, and model parameters. This state is
#'   updated in place (functionally) at each iteration.
#' @param sa_func Function. Stochastic approximation step function (SA-step) to
#' update sufficient statistics.
#' @param m_func Function. M-step function to update model parameters.
#' @param mh_update_func Function. Metropolis-Hastings update function.
#' @param backend A list containing the compiled model backend used throughout
#'   the SAEM iterations to evaluate the likelihood and perform model-specific
#'   computations.
#' @return Final \code{state} list after completing all SAEM iterations.
#' @keywords internal
run_saem_full <- function(config, state, sa_func, m_func,
                          mh_update_func, backend) {
  for (k in 1:config$num_iterations) {
    state <- perform_saem_iteration(
      k, config, state, sa_func, m_func,
      mh_update_func, backend
    )
  }
  return(state) # nolint: return-linter
}

#' Run SAEM algorithm with fixed parameters
#'
#' Executes the SAEM algorithm for models with both random and fixed effects.
#' Shrinkage is applied to the covariance of fixed parameters after burn-in
#' iterations.
#' @param config List. Model configuration object containing algorithm settings,
#'   including:
#'   \itemize{
#'     \item \code{num_iterations}: total number of SAEM iterations,
#'     \item \code{num_burnin}: number of burn-in iterations,
#'     \item \code{parameters_not_to_select_indices}: indices of parameters not
#'       subject to variable selection,
#'     \item \code{fixed_parameters_indices}: indices of parameters treated as
#'       fixed effects.
#'   }
#' @param state List. Initial state of the algorithm. This object is
#'   iteratively updated at each iteration
#' @param sa_func Function. Stochastic approximation step function (SA-step) to
#' update sufficient statistics.
#' @param m_func Function. M-step function to update model parameters.
#' @param mh_update_func Function. Metropolis-Hastings update function.
#' @return Final \code{state} list after completing all SAEM iterations with
#' shrinkage for fixed effects.
#' @param backend A list containing the compiled model backend used to evaluate
#' the likelihood and perform model-specific computations during the SAEM
#' iterations.
#' @details
#' In this context, "fixed effects" refers to parameters that are not subject
#' to variable selection and whose random-effect covariance matrix is
#' progressively shrunk after the burn-in phase. These parameters are not
#' fixed to constant values, but their variability is reduced over iterations.
#' After the burn-in phase, shrinkage is progressively applied to the
#' covariance matrix of the fixed-effect parameters. Shrinkage coefficients
#' are computed from the diagonal elements of the covariance matrix at the
#' end of burn-in and are applied multiplicatively at each subsequent
#' iteration via \code{shrink_covariance_matrix}.

#' @keywords internal
run_saem_with_fixed_effects <- function(
    config, state, sa_func, m_func, mh_update_func, backend) {
  for (k in 1:config$num_burnin) {
    state <- perform_saem_iteration(
      k, config, state, sa_func, m_func,
      mh_update_func, backend
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
      mh_update_func, backend
    )
    state$gamma_not_to_select[[k + 1]] <- shrink_covariance_matrix(
      state$gamma_not_to_select[[k]],
      state$gamma_not_to_select[[k + 1]],
      index_fixed,
      alphas
    )
  }

  return(state) # nolint: return-linter
}
