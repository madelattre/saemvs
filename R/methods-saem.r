#' Run the SAEM algorithm
#'
#' Internal generic method for executing the Stochastic Approximation
#' Expectation-Maximization (SAEM) algorithm in the context of variable
#' selection with spike-and-slab priors or maximum likelihood estimation.
#'
#' This method orchestrates the full SAEM procedure, including:
#' \itemize{
#'   \item Validation and preprocessing of input data and model.
#'   \item Initialization of algorithmic state.
#'   \item Configuration of hyperparameters and tuning parameters.
#'   \item Selection of appropriate stochastic approximation (SA) and
#'         maximization (M) steps depending on the model type and presence of
#'         fixed effects.
#'   \item Iterative execution of the SAEM algorithm until convergence.
#' }
#'
#' @param data An object of class \code{saemvsProcessedData} containing the
#'   observed data, design matrices, and related metadata.
#' @param model An object of class \code{saemvsProcessedModel} specifying the
#'   model structure, covariates, random effects, and indices of parameters
#'   subject or not subject to selection.
#' @param init An object of class \code{saemvsProcessedInit} providing
#'   initialization for the algorithm (parameters, variances, latent states,
#'  etc.).
#' @param tuning_algo An object of class \code{saemvsTuning} defining
#'  algorithmic
#'   parameters such as number of iterations, burn-in period, step-size
#'  sequence, and parallelization settings.
#' @param hyperparam An object of class \code{saemvsHyperSpikeAndSlab}
#' containing hyperparameters for the spike-and-slab prior or other priors used
#'  for variable selection.
#' @param backend A list of compiled model functions and helper routines, used
#'   internally for fast computation during SAEM-MCMC iterations.
#'
#' @return A \code{state} list containing the results of the SAEM iterations,
#'  including:
#'   \itemize{
#'     \item \code{phi}: latent parameter matrices per iteration.
#'     \item \code{beta_to_select}, \code{beta_not_to_select}: regression
#'  coefficients.
#'     \item \code{gamma_to_select}, \code{gamma_not_to_select}: covariance
#'  matrices.
#'     \item \code{alpha}: posterior inclusion probabilities for parameters to
#'  select.
#'     \item \code{sigma2}: residual variance per iteration.
#'     \item Stochastic approximation statistics (s1, s2, s3) used internally.
#'   }
#'
#' @details
#' Depending on the configuration of the model, hyperparameters, and presence of
#'  fixed effects,
#' the function selects the appropriate SA and M-step implementations:
#' \itemize{
#'   \item \code{"map_full_select"}: MAP estimation with all parameters
#'  subject to selection.
#'   \item \code{"map_part_select_nofixed"}: MAP with partial selection,
#'  no fixed parameters.
#'   \item \code{"map_part_select_fixed"}: MAP with partial selection and
#'  fixed parameters.
#'   \item \code{"mle_nofixed"}: Maximum Likelihood Estimation with no fixed
#'  parameters.
#'   \item \code{"mle_fixed"}: Maximum Likelihood Estimation with fixed
#'  parameters.
#' }
#'
#' @note This method is internal to the package and should not be called
#'  directly by end-users.
#' High-level functions (e.g., \code{saemvs}) provide the proper interface
#'  for running the algorithm.
#'
#' @keywords internal
#' @noRd
setGeneric(
  "run_saem",
  function(data, model, init, tuning_algo, hyperparam, backend) {
    standardGeneric("run_saem")
  }
)

setMethod(
  "run_saem",
  signature(
    data = "saemvsProcessedData",
    model = "saemvsProcessedModel",
    init = "saemvsProcessedInit",
    tuning_algo = "saemvsTuning",
    hyperparam = "saemvsHyperSpikeAndSlab",
    backend = "list"
  ),
  function(data, model, init, tuning_algo, hyperparam, backend) {
    check_hyper(hyperparam, model, tuning_algo)
    hyperparam <- prepare_hyper(hyperparam, data, model)

    config <- make_config(
      data, model, tuning_algo, init, hyperparam
    )
    state <- init_state(config)
    case <- get_case(config)

    switch(case,
      map_full_select = {
        state <- run_saem_full(
          config, state, sa_step_to_select, m_step_map_to_select,
          update_proposal_mh_to_select, backend
        )
      },
      map_part_select_nofixed = {
        state <- run_saem_full(
          config, state, sa_step_all, m_step_map_all, update_proposal_mh_all,
          backend
        )
      },
      map_part_select_fixed = {
        state <- run_saem_with_fixed_effects(
          config, state, sa_step_all, m_step_map_all, update_proposal_mh_all,
          backend
        )
      },
      mle_nofixed = {
        state <- run_saem_full(
          config, state, sa_step_not_to_select, m_step_mle,
          update_proposal_mh_not_to_select, backend
        )
      },
      mle_fixed = {
        state <- run_saem_with_fixed_effects(
          config, state, sa_step_not_to_select, m_step_mle,
          update_proposal_mh_not_to_select, backend
        )
      }
    )

    return(state) # nolint: return-linter
  }
)
