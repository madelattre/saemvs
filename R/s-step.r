#' @title Metropolis-Hastings Update and Proposal Management for Latent
#' Parameters
#' @description Internal functions to update latent parameters (\code{phi}) and
#'   the corresponding proposal distributions during the SAEM-MCMC algorithm.
#'   These functions are used internally for the Metropolis-Hastings (MH) steps.
#' @details
#' - `metropolis_s_step()` performs one MH update for all series at a given
#' iteration.
#' - `update_proposal_mh_all()` updates the mean and variance of the proposal
#' distribution for both parameters subject to selection (to-select) and
#' parameters not subject to selection (not-to-select).
#' - `update_proposal_mh_to_select()` updates the proposals only for parameters
#' subject to selection.
#' - `update_proposal_mh_not_to_select()` updates the proposals only for
#' parameters not subject to selection.
#' @param config Internal configuration object containing data, model matrices,
#' and MH parameters.
#' The \code{config} object must include:
#'   \itemize{
#'     \item \code{y_series}, \code{t_series}: observations and time points,
#'     \item \code{x_phi_to_select}, \code{x_phi_not_to_select}:
#' design matrices,
#'     \item \code{parameters_to_select_indices},
#'  \code{parameters_not_to_select_indices}: indices of parameters
#' for selection,
#'     \item \code{num_mh_iterations}: number of MH iterations per SAEM
#' iteration,
#'     \item \code{mh_proposal_scale}: scaling factor for proposal distribution,
#'     \item \code{mh_kernel_type}: kernel type for MH updates,
#'     \item \code{num_series}: number of independent series.
#'   }
#' @param iteration Integer specifying the current iteration of the SAEM
#' algorithm.
#' @param state List containing the current latent state and model parameters,
#'   including:
#'   \itemize{
#'     \item \code{phi}: list of latent parameter matrices per iteration,
#'     \item \code{mprop_mh}: list of proposal mean matrices per iteration,
#'     \item \code{vprop_mh}: list of proposal covariance matrices per
#' iteration,
#'     \item \code{beta_to_select}, \code{beta_not_to_select}: regression
#' coefficients,
#'     \item \code{gamma_to_select}, \code{gamma_not_to_select}: covariance
#' matrices.
#'   }
#'   This list is updated in place by the MH steps and proposal update
#' functions.
#' @param backend A list containing compiled model functions and helper
#'   routines used to perform Metropolis-Hastings updates for latent
#'   parameters.
#' @return Updated \code{state} list:
#'   \itemize{
#'     \item \code{metropolis_s_step()}: updates \code{phi} for the next
#' iteration,
#'     \item \code{update_proposal_mh_all()}: updates \code{mprop_mh} and
#' \code{vprop_mh} for all parameters,
#'     \item \code{update_proposal_mh_to_select()}: updates only the to-select
#'  parameters,
#'     \item \code{update_proposal_mh_not_to_select()}: updates only the
#' not-to-select parameters.
#'   }

#' @keywords internal
#' @note These functions are intended for internal use within the SAEM-MCMC
#' algorithm and should not be called directly in user code.
metropolis_s_step <- function(config, iteration, state, backend) {
  sampled_phi <- backend$metropolis_vector(
    y = config$y_series,
    t = config$t_series,
    phi_current = split(state$phi[[iteration]], row(state$phi[[iteration]])),
    mean_prop = split(
      state$mprop_mh[[iteration]],
      row(state$mprop_mh[[iteration]])
    ),
    var_prop = state$vprop_mh[[iteration]],
    sigma2 = state$sigma2[iteration],
    niter_mh = config$num_mh_iterations,
    kappa = config$mh_proposal_scale,
    kernel = config$mh_kernel_type
  )

  state$phi[[iteration + 1]] <- matrix(unlist(sampled_phi),
    nrow = config$num_series,
    byrow = TRUE
  )
  return(state) # nolint: return-linter
}

#' @rdname metropolis_s_step
update_proposal_mh_all <- function(config, iteration, state) {
  state$mprop_mh[[iteration + 1]][, config$parameters_to_select_indices] <-
    config$x_phi_to_select %*% state$beta_to_select[[iteration + 1]]

  state$mprop_mh[[iteration + 1]][, config$parameters_not_to_select_indices] <-
    config$x_phi_not_to_select %*% state$beta_not_to_select[[iteration + 1]]

  state$vprop_mh[[iteration + 1]][
    config$parameters_to_select_indices,
    config$parameters_to_select_indices
  ] <-
    state$gamma_to_select[[iteration + 1]]

  state$vprop_mh[[iteration + 1]][
    config$parameters_not_to_select_indices,
    config$parameters_not_to_select_indices
  ] <-
    state$gamma_not_to_select[[iteration + 1]]

  return(state) # nolint: return-linter
}

#' @rdname metropolis_s_step
update_proposal_mh_to_select <- function(config, iteration, state) {
  state$mprop_mh[[iteration + 1]][, config$parameters_to_select_indices] <-
    config$x_phi_to_select %*% state$beta_to_select[[iteration + 1]]

  state$vprop_mh[[iteration + 1]][
    config$parameters_to_select_indices,
    config$parameters_to_select_indices
  ] <-
    state$gamma_to_select[[iteration + 1]]

  return(state) # nolint: return-linter
}

#' @rdname metropolis_s_step
update_proposal_mh_not_to_select <- function(config, iteration, state) { # nolint: object_length_linter
  state$mprop_mh[[iteration + 1]][, config$parameters_not_to_select_indices] <-
    config$x_phi_not_to_select %*% state$beta_not_to_select[[iteration + 1]]

  state$vprop_mh[[iteration + 1]][
    config$parameters_not_to_select_indices,
    config$parameters_not_to_select_indices
  ] <-
    state$gamma_not_to_select[[iteration + 1]]

  return(state) # nolint: return-linter
}
