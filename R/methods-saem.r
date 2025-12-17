#' Run the SAEM algorithm
#'
#' Internal generic method for running the Stochastic Approximation
#' Expectation-Maximization (SAEM) algorithm in the context of variable
#' selection with spike-and-slab priors or maximum likelihood estimation.
#'
#' This method orchestrates the full SAEM procedure, including:
#' \itemize{
#'   \item Data validation and preprocessing
#'   \item Initialization of algorithmic state
#'   \item Configuration of hyperparameters and tuning parameters
#'   \item Selection of appropriate stochastic approximation (SA) and
#'   maximization (M) steps depending on the model type and presence of fixed
#'   effects
#'   \item Iterative execution of the SAEM algorithm until convergence
#' }
#'
#' The function is dispatched on specific classes of data, model,
#' initialization, tuning, and hyperparameter objects, and should not be called
#' directly by the end-user.
#'
#' @param data An object of class \code{saemvsData}.
#'   Contains the observed data and related metadata.
#' @param model An object of class \code{saemvsProcessedModel}.
#'   Describes the statistical model, including covariates and random effect
#'   structure.
#' @param init An object of class \code{saemvsInit}.
#'   Contains initialization values for the algorithm (parameters, variances,
#'   etc.).
#' @param tuning_algo An object of class \code{saemvsTuning}.
#'   Defines algorithmic tuning parameters (number of iterations, burn-in,
#'   decay rates, etc.).
#' @param hyperparam An object of class \code{saemvsHyperSpikeAndSlab}.
#'   Hyperparameters for the spike-and-slab prior used in variable selection.
#'
#' @return An updated \code{state} object containing the results of the SAEM
#' iterations, including parameter estimates, variances, and posterior inclusion
#' probabilities.
#'
#' @details
#' Depending on the configuration of the model and hyperparameters, the function
#' automatically selects the appropriate SA and M step implementations:
#' \itemize{
#'   \item \code{"map_full_select"}: MAP estimation with all parameters subject
#'   to selection.
#'   \item \code{"map_part_select_nofixed"}: MAP with partial selection and no
#'   fixed parameters.
#'   \item \code{"map_part_select_fixed"}: MAP with partial selection and some
#'   fixed parameters.
#'   \item \code{"mle_nofixed"}: MLE with no fixed parameters.
#'   \item \code{"mle_fixed"}: MLE with some fixed parameters.
#' }
#'
#' @note This method is internal to the package and not intended to be called
#' directly by end-users. Instead, high-level user functions should be used to
#' set up and run the algorithm.
#'
#' @keywords internal
setGeneric(
  "run_saem",
  function(data, model, init, tuning_algo, hyperparam) {
    standardGeneric("run_saem")
  }
)

setGeneric(
  "run_saem",
  function(data, model, init, tuning_algo, hyperparam) {
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
    hyperparam = "saemvsHyperSpikeAndSlab"
  ),
  function(data, model, init, tuning_algo, hyperparam) {
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
          update_proposal_mh_to_select
        )
      },
      map_part_select_nofixed = {
        state <- run_saem_full(
          config, state, sa_step_all, m_step_map_all, update_proposal_mh_all
        )
      },
      map_part_select_fixed = {
        state <- run_saem_with_fixed_effects(
          config, state, sa_step_all, m_step_map_all, update_proposal_mh_all
        )
      },
      mle_nofixed = {
        state <- run_saem_full(
          config, state, sa_step_not_to_select, m_step_mle,
          update_proposal_mh_not_to_select
        )
      },
      mle_fixed = {
        state <- run_saem_with_fixed_effects(
          config, state, sa_step_not_to_select, m_step_mle,
          update_proposal_mh_not_to_select
        )
      }
    )

    return(state) # nolint: return-linter
  }
)
