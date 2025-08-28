#' @title Metropolis-Hastings Update and Proposal Management for Latent Parameters
#' @description Internal functions to update latent parameters (\code{phi}) and
#'   the corresponding proposal distributions during the SAEM-MCMC algorithm.
#'   These functions are used internally for the Metropolis-Hastings (MH) steps.
#' @details
#' - `metropolis_s_step()` performs one MH update for all series at a given iteration.
#' - `update_proposal_mh_all()` updates the mean and variance of the proposal distribution
#'   for both parameters subject to selection (to-select) and parameters not subject to selection (not-to-select).
#' - `update_proposal_mh_to_select()` updates the proposals only for parameters subject to selection.
#' - `update_proposal_mh_not_to_select()` updates the proposals only for parameters not subject to selection.
#' @param config Internal configuration object containing data, model matrices, and MH parameters.
#' @param iteration Integer specifying the current iteration of the SAEM algorithm.
#' @param state List containing the current latent state, including:
#'   \itemize{
#'     \item \code{phi}: list of latent parameter matrices per iteration
#'     \item \code{mprop_mh}: list of proposal mean matrices per iteration
#'     \item \code{vprop_mh}: list of proposal covariance matrices per iteration
#'     \item \code{beta_to_select}: list of regression coefficient matrices for parameters to select
#'     \item \code{beta_not_to_select}: list of regression coefficient matrices for parameters not to select
#'     \item \code{gamma_to_select}: list of covariance matrices for parameters to select
#'     \item \code{gamma_not_to_select}: list of covariance matrices for parameters not to select
#'   }
#' @return Updated \code{state} list with new latent parameters or proposal matrices.
#' @keywords internal
metropolis_s_step <- function(config, iteration, state) {
  sampled_phi <- metropolis_vector_cpp(
    y = config$y_series,
    t = config$t_series,
    phi_current = split(state$phi[[iteration]], row(state$phi[[iteration]])),
    mean_prop = split(state$mprop_mh[[iteration]], row(state$mprop_mh[[iteration]])),
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
  return(state)
}

#' @rdname metropolis_s_step
update_proposal_mh_all <- function(config, iteration, state) {
  # Update mean proposals for to-select and not-to-select parameters
  state$mprop_mh[[iteration + 1]][, config$parameters_to_select_indices] <-
    config$x_phi_to_select %*% state$beta_to_select[[iteration + 1]]

  state$mprop_mh[[iteration + 1]][, config$parameters_not_to_select_indices] <-
    config$x_phi_not_to_select %*% state$beta_not_to_select[[iteration + 1]]

  # Update covariance proposals for both groups
  state$vprop_mh[[iteration + 1]][config$parameters_to_select_indices, config$parameters_to_select_indices] <-
    state$gamma_to_select[[iteration + 1]]

  state$vprop_mh[[iteration + 1]][config$parameters_not_to_select_indices, config$parameters_not_to_select_indices] <-
    state$gamma_not_to_select[[iteration + 1]]

  return(state)
}

#' @rdname metropolis_s_step
update_proposal_mh_to_select <- function(config, iteration, state) {
  # Update mean and covariance proposals for parameters to select only
  state$mprop_mh[[iteration + 1]][, config$parameters_to_select_indices] <-
    config$x_phi_to_select %*% state$beta_to_select[[iteration + 1]]

  state$vprop_mh[[iteration + 1]][config$parameters_to_select_indices, config$parameters_to_select_indices] <-
    state$gamma_to_select[[iteration + 1]]

  return(state)
}

#' @rdname metropolis_s_step
update_proposal_mh_not_to_select <- function(config, iteration, state) {
  # Update mean and covariance proposals for parameters not to select only
  state$mprop_mh[[iteration + 1]][, config$parameters_not_to_select_indices] <-
    config$x_phi_not_to_select %*% state$beta_not_to_select[[iteration + 1]]

  state$vprop_mh[[iteration + 1]][config$parameters_not_to_select_indices, config$parameters_not_to_select_indices] <-
    state$gamma_not_to_select[[iteration + 1]]

  return(state)
}
