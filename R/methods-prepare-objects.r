#' Prepare processed design matrices for SAEMVS model fitting
#'
#' Internal method that constructs design matrices and related structures
#' required for fitting a \code{saemvsModel} to a \code{saemvsData} object.
#'
#' Specifically, it:
#' \itemize{
#'   \item separates covariates for parameters subject to selection from those for parameters not subject to selection,
#'   \item builds corresponding design matrices (with intercept),
#'   \item computes Gram matrices and block-diagonal expansions for covariates related to parameters subject to selection,
#'   \item constructs a list representation of design matrices for parameters not subject to selection for efficient use in SAEM updates,
#'   \item and returns a \code{saemvsProcessedData} object encapsulating all results.
#' }
#' @param data A \code{\linkS4class{saemvsData}} object.
#'   Contains observed time series \code{y_series}, covariates
#'   to be selected (\code{x_candidates}), and forced covariates
#'   (\code{x_forced}).
#'
#' @param model A \code{\linkS4class{saemvsModel}} object.
#'   Provides model dimension (\code{phi_dim}), indices of parameters
#'   subject to selection (\code{phi_to_select_idx}), and the forced
#'   support structure (\code{x_forced_support}).
#'
#' @return A \code{\linkS4class{saemvsProcessedData}} object containing:
#' \describe{
#'   \item{\code{y_series}}{Original response time series.}
#'   \item{\code{t_series}}{Associated time points.}
#'   \item{\code{x_candidates}}{Design matrix of candidate covariates.}
#'   \item{\code{x_forced}}{Design matrix of forced covariates.}
#'   \item{\code{x_phi_to_select}}{Final design matrix for parameters subject to selection, including intercept.}
#'   \item{\code{x_phi_not_to_select}}{Final design matrix for parameters not subject to selection, including intercept.}
#'   \item{\code{tx_x_phi_to_select}}{Gram matrix \eqn{X'X} for the parameters subject to selection.}
#'   \item{\code{kron_tx_x_phi_to_select}}{Block-diagonal expansion of the Gram matrix, one block per parameter subject to selection.}
#'   \item{\code{x_phi_not_to_select_list}}{List of matrices representing the design for parameters not subject to selection, expanded by support.}
#' }
#'
#' @note
#' This method is internal to the package and not intended for direct user calls.
#' It assumes that \code{data} and \code{model} are consistent, and may error
#' if the forced support structure does not match model dimensions.
#'
#' @seealso \code{\linkS4class{saemvsData}}, \code{\linkS4class{saemvsModel}}, \code{\linkS4class{saemvsProcessedData}}
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' data_obj <- methods::new("saemvsData", ...) # create data
#' model_obj <- methods::new("saemvsModel", ...) # create model
#' processed <- prepare_data(data_obj, model_obj)
#' }
setGeneric(
  "prepare_data",
  function(data, model) {
    standardGeneric("prepare_data")
  }
)

setMethod(
  "prepare_data",
  signature(data = "saemvsData", model = "saemvsModel"),
  function(data, model) {
    n_obs <- length(data@y_series)

    ## === Extract model indices and supports ===
    phi_to_select_indices <- model@phi_to_select_idx
    phi_not_to_select_indices <- setdiff(seq_len(model@phi_dim), phi_to_select_indices)
    forced_support <- model@x_forced_support

    ## === 1. Build design matrix for parameters subject to selection ===
    x_to_select_design <- NULL

    if (length(phi_to_select_indices) > 0) {
      if (is_empty_support(forced_support)) {
        ## No forced covariates at all → intercept + candidates
        x_to_select_design <- cbind(1, data@x_candidates)
      } else if (is_empty_support(forced_support[, phi_to_select_indices])) {
        ## No forced covariates among the candidate subset → intercept + candidates
        x_to_select_design <- cbind(1, data@x_candidates)
      } else {
        ## Some forced covariates correspond to params subject to selection
        forced_support_params_to_select <- matrix(
          forced_support[, phi_to_select_indices],
          ncol = length(phi_to_select_indices)
        )
        forced_cov_params_to_select_idx <- extract_rows_with_ones(forced_support_params_to_select)
        x_to_select_design <- cbind(
          1,
          data@x_forced[, forced_cov_params_to_select_idx, drop = FALSE],
          data@x_candidates
        )
      }
    }

    ## === 2. Compute Gram and block-diagonal matrices ===
    if (is.null(x_to_select_design)) {
      tx_x_to_select <- NULL
      kron_tx_x_to_select <- NULL
    } else {
      tx_x_to_select <- t(x_to_select_design) %*% x_to_select_design
      kron_tx_x_to_select <- as.matrix(Matrix::bdiag(replicate(
        length(phi_to_select_indices),
        tx_x_to_select,
        simplify = FALSE
      )))
    }

    ## === 3. Build design matrix for unselected parameters ===
    x_not_to_select_design <- NULL
    forced_support_params_not_to_select <- NULL

    if (is_empty_support(forced_support)) {
      x_not_to_select_design <- matrix(1, nrow = n_obs)
    } else if (length(phi_not_to_select_indices) == 0 ||
      is_empty_support(forced_support[, phi_not_to_select_indices])) {
      ## No unselected parameters or no forced covariates among them
      x_not_to_select_design <- matrix(1, nrow = n_obs)
    } else {
      forced_support_params_not_to_select <- matrix(
        forced_support[, phi_not_to_select_indices],
        ncol = length(phi_not_to_select_indices)
      )
      forced_cov_params_not_to_select_idx <- extract_rows_with_ones(forced_support_params_not_to_select)
      x_not_to_select_design <- cbind(
        1,
        data@x_forced[, forced_cov_params_not_to_select_idx, drop = FALSE]
      )
    }

    ## === 4. Build list representation for unselected design ===
    if (ncol(x_not_to_select_design) == 1) {
      x_not_to_select_list <- lapply(
        seq_len(length(data@y_series)),
        function(i) diag(1, length(phi_not_to_select_indices))
      )
    } else {
      if (is.null(forced_support_params_not_to_select)) {
        stop("Internal error: 'support' is NULL while x_design has >1 column.")
      }
      x_not_to_select_list <- expand_to_list(
        x_not_to_select_design,
        rbind(1, forced_support_params_not_to_select), # include intercept
        length(phi_not_to_select_indices)
      )
    }

    ## === 5. Construct processed data object ===
    data_processed <- methods::new("saemvsProcessedData",
      y_series                 = data@y_series,
      t_series                 = data@t_series,
      x_candidates             = data@x_candidates,
      x_forced                 = data@x_forced,
      x_phi_to_select          = x_to_select_design,
      x_phi_not_to_select      = x_not_to_select_design,
      tx_x_phi_to_select       = tx_x_to_select,
      kron_tx_x_phi_to_select  = kron_tx_x_to_select,
      x_phi_not_to_select_list = x_not_to_select_list
    )

    return(data_processed)
  }
)

#' Prepare initial values for the SAEMVS algorithm
#'
#' This internal method constructs processed initial values for the SAEMVS
#' algorithm, splitting parameters into those subject to variable selection
#' and those not. It combines intercepts, forced and candidate
#' betas into matrices ready for algorithmic use.
#'
#' @param init An object of class \code{saemvsInit} containing initial values
#'   for intercepts, candidate betas, forced betas, covariance of random effects,
#'   and residual variance.
#' @param model An object of class \code{saemvsModel} defining the model structure,
#'   including indices of parameters subject to selection (\code{phi_to_select_idx}) and
#'   forced covariate support (\code{x_forced_support}).
#'
#' @return An object of class \code{saemvsProcessedInit} containing:
#' \describe{
#'   \item{\code{inclusion_prob}}{Vector of prior inclusion probabilities for
#'         parameters subject to selection (default 0.5).}
#'   \item{\code{sigma2}}{Residual variance, copied from \code{init}.}
#'   \item{\code{gamma_to_select}}{Covariance matrix of random effects for
#'         parameters subject to selection.}
#'   \item{\code{gamma_not_to_select}}{Covariance matrix of random effects for
#'         parameters not subject to selection.}
#'   \item{\code{beta_to_select}}{Matrix of coefficients for parameters subject
#'         to selection, combining intercept, forced and candidate betas.}
#'   \item{\code{beta_not_to_select}}{Matrix of coefficients for parameters not
#'         subject to selection, combining intercept and forced covariates.}
#' }
#'
#' @details
#' The function handles the following cases explicitly:
#' \itemize{
#'   \item No parameters are subject to selection (\code{phi_to_select_idx} is empty),
#'         resulting in \code{beta_to_select}, \code{gamma_to_select}, and
#'         \code{inclusion_prob} being \code{NULL}.
#'   \item No forced covariates for either selected or non-selected parameters,
#'         in which case only intercepts and candidate betas are used.
#'   \item Forced covariates present, in which case only rows corresponding to
#'         active support (rows with ones) are included in the beta matrices.
#' }
#'
#' The method ensures all matrices maintain proper dimensions using \code{drop = FALSE}.
#'
#' @keywords internal
setGeneric(
  "prepare_init",
  function(init, model) {
    standardGeneric("prepare_init")
  }
)

setMethod(
  "prepare_init",
  signature(init = "saemvsInit", model = "saemvsModel"),
  function(init, model) {
    ## === 0. Extract indices ===
    phi_to_select_idx <- model@phi_to_select_idx
    phi_not_to_select_idx <- setdiff(seq_len(model@phi_dim), phi_to_select_idx)
    forced_support <- model@x_forced_support

    ## === 1. Prepare beta and gamma for parameters subject to selection ===
    if (length(phi_to_select_idx) == 0) {
      beta_to_select <- NULL
      gamma_to_select <- NULL
      inclusion_prob <- NULL
    } else {
      # Extract forced support for parameters subject to selection
      if (is_empty_support(forced_support)) {
        xf_sel <- NULL
      } else {
        xf_sel <- matrix(forced_support[, phi_to_select_idx],
          ncol = length(phi_to_select_idx)
        )
      }

      # Compose beta_to_select
      if (is_empty_support(xf_sel)) {
        beta_to_select <- rbind(
          init@intercept[phi_to_select_idx],
          init@beta_candidates[, phi_to_select_idx, drop = FALSE]
        )
      } else {
        rows_forced_sel <- extract_rows_with_ones(xf_sel)
        beta_forced_sel <- init@beta_forced[rows_forced_sel, phi_to_select_idx, drop = FALSE]
        beta_to_select <- rbind(
          init@intercept[phi_to_select_idx],
          beta_forced_sel,
          init@beta_candidates[, phi_to_select_idx, drop = FALSE]
        )
      }

      # gamma_to_select
      gamma_to_select <- matrix(
        init@cov_re[phi_to_select_idx, phi_to_select_idx, drop = FALSE],
        ncol = length(phi_to_select_idx)
      )

      # gamma_to_select <- matrix(
      #   cov_re[phi_to_select_idx, phi_to_select_idx, drop = FALSE],
      #   ncol = length(phi_to_select_idx)
      # )

      # inclusion probabilities (default 0.5)
      inclusion_prob <- rep(0.5, length(phi_to_select_idx))
    }

    ## === 2. Prepare beta and gamma for parameters NOT subject to selection ===
    if (length(phi_not_to_select_idx) == 0) {
      beta_not_to_select <- NULL
      gamma_not_to_select <- NULL
    } else {
      if (is_empty_support(forced_support)) {
        xf_not_sel <- NULL
      } else {
        xf_not_sel <- matrix(forced_support[, phi_not_to_select_idx],
          ncol = length(phi_not_to_select_idx)
        )
      }

      # Compose beta_not_to_select
      if (is_empty_support(xf_not_sel)) {
        beta_not_to_select <- matrix(init@intercept[phi_not_to_select_idx],
          ncol = length(phi_not_to_select_idx)
        )
      } else {
        rows_forced_not_sel <- extract_rows_with_ones(xf_not_sel)
        beta_forced_not_sel <- matrix(
          init@beta_forced[rows_forced_not_sel, phi_not_to_select_idx, drop = FALSE],
          ncol = length(phi_not_to_select_idx)
        )
        beta_not_to_select <- rbind(
          init@intercept[phi_not_to_select_idx],
          beta_forced_not_sel
        )
      }

      # gamma_not_to_select

      # replace fixed effect variances with strictly positive values for sampling (exponentialization trick)
      cov_re <- init@cov_re
      phi_idx <- model@phi_fixed_idx

      if (length(phi_idx) > 0) {
        cov_re[cbind(phi_idx, phi_idx)] <-
          pmax(0.5 * (init@intercept[phi_idx])**2, 1)
      }



      # gamma_not_to_select <- matrix(
      #   init@cov_re[phi_not_to_select_idx, phi_not_to_select_idx, drop = FALSE],
      #   ncol = length(phi_not_to_select_idx)
      # )

      gamma_not_to_select <- matrix(
        cov_re[phi_not_to_select_idx, phi_not_to_select_idx, drop = FALSE],
        ncol = length(phi_not_to_select_idx)
      )

    }


    ## === 3. Construct processed init object ===
    init_processed <- methods::new("saemvsProcessedInit",
      inclusion_prob      = inclusion_prob,
      sigma2              = init@sigma2,
      gamma_to_select     = gamma_to_select,
      gamma_not_to_select = gamma_not_to_select,
      beta_to_select      = beta_to_select,
      beta_not_to_select  = beta_not_to_select
    )

    return(init_processed)
  }
)


#' Prepare and Complete Hyperparameters for SAEM Variable Selection
#'
#' This internal function completes and validates the hyperparameters of a
#' \code{saemvsHyperSlab} object according to the processed data and model.
#' For parameters subject to selection (MAP case), it ensures that the prior
#' inclusion probabilities \code{a} and \code{b} are set correctly.
#' In the case where no parameter is selected (MLE case), it returns
#' a default hyperparameter object with NULL fields.
#'
#' @param hyper A \code{saemvsHyperSlab} object containing the current hyperparameters.
#' @param data A \code{saemvsProcessedData} object containing the processed data for the algorithm.
#' @param model A \code{saemvsModel} object representing the model structure.
#'
#' @return A \code{saemvsHyperSlab} object with completed and validated hyperparameters.
#'
#' @details
#' - If no parameters are subject to selection (\code{length(model@phi_to_select_idx) == 0}),
#'   the function returns \code{saemvsHyperSlab(NULL, NULL, NULL)}.
#' - If some parameters are subject to selection, it fills missing \code{inclusion_prob_prior_a}
#'   and \code{inclusion_prob_prior_b} with default values (\code{1} and number of candidates, respectively),
#'   and checks their length consistency.
#' - Stops with a descriptive error if the provided prior vectors have inconsistent lengths.
#'
#' @examples
#' \dontrun{
#' hyper <- saemvsHyperSlab(NULL, NULL, NULL)
#' data_processed <- prepare_data(data, model)
#' hyper_completed <- prepare_hyper(hyper, data_processed, model)
#' }
#'
#' @keywords internal
setGeneric(
  "prepare_hyper",
  function(hyper, data, model) {
    standardGeneric("prepare_hyper")
  }
)

setMethod(
  "prepare_hyper",
  signature(
    hyper = "saemvsHyperSlab",
    data = "saemvsProcessedData",
    model = "saemvsModel"
  ),
  function(hyper, data, model) {
    # Number of parameters subject to selection
    n_selected <- length(model@phi_to_select_idx)

    # MLE case: no parameter is selected
    if (n_selected == 0) {
      return(saemvsHyperSlab(NULL, NULL, NULL))
    }


    # MAP case: parameters are selected
    n_candidates <- ncol(data@x_candidates)

    # Complete inclusion probability prior 'a' if missing
    if (is.null(hyper@inclusion_prob_prior_a)) {
      hyper@inclusion_prob_prior_a <- rep(1, n_selected)
    } else if (length(hyper@inclusion_prob_prior_a) != n_selected) {
      stop(sprintf(
        "Length of 'inclusion_prob_prior_a' (%d) does not match the number of selected parameters (%d).",
        length(hyper@inclusion_prob_prior_a), n_selected
      ))
    }

    # Complete inclusion probability prior 'b' if missing
    if (is.null(hyper@inclusion_prob_prior_b)) {
      hyper@inclusion_prob_prior_b <- rep(n_candidates, n_selected)
    } else if (length(hyper@inclusion_prob_prior_b) != n_selected) {
      stop(sprintf(
        "Length of 'inclusion_prob_prior_b' (%d) does not match the number of selected parameters (%d).",
        length(hyper@inclusion_prob_prior_b), n_selected
      ))
    }

    return(hyper)
  }
)

#' Create a unified configuration list for the SAEMVS algorithm
#'
#' This function constructs a configuration list containing processed data, model information,
#' initialization, hyperparameters, and tuning parameters, with fields aligned with the
#' original object slots.
#'
#' @param data A \code{saemvsProcessedData} object.
#' @param model A \code{saemvsModel} object.
#' @param tuning_algo A \code{saemvsTuning} object.
#' @param init A \code{saemvsProcessedInit} object.
#' @param hyperparam A \code{saemvsHyperSlab} object.
#'
#' @return A named list containing all information necessary to run the SAEMVS algorithm.
#' @keywords internal
setGeneric(
  "make_config",
  function(data, model, tuning_algo, init, hyperparam) {
    standardGeneric("make_config")
  }
)

setMethod(
  "make_config",
  signature(
    data = "saemvsProcessedData", model = "saemvsModel", tuning_algo = "saemvsTuning",
    init = "saemvsProcessedInit", hyperparam = "saemvsHyperSlab"
  ),
  function(data, model, tuning_algo, init, hyperparam) {
    # Dimensions
    total_phi <- model@phi_dim
    phi_to_select_idx <- model@phi_to_select_idx
    phi_not_to_select_idx <- setdiff(seq_len(total_phi), phi_to_select_idx)
    n_phi_to_select <- length(phi_to_select_idx)
    n_phi_not_to_select <- total_phi - n_phi_to_select

    # Method type
    method_type <- if (n_phi_to_select == 0) "mle" else "map"

    # MAP hyperparameters
    spike_parameter <- slab_parameter <- residual_variance_prior_shape <-
      residual_variance_prior_rate <- inclusion_prob_prior_a <-
      inclusion_prob_prior_b <- phi_intercept_prior_variance <-
      cov_re_prior_scale <- cov_re_prior_df <- NULL

    if (method_type == "map") {
      spike_parameter <- hyperparam@spike_parameter
      slab_parameter <- hyperparam@slab_parameter
      residual_variance_prior_shape <- hyperparam@residual_variance_prior_shape
      residual_variance_prior_rate <- hyperparam@residual_variance_prior_rate
      inclusion_prob_prior_a <- hyperparam@inclusion_prob_prior_a
      inclusion_prob_prior_b <- hyperparam@inclusion_prob_prior_b
      phi_intercept_prior_variance <- hyperparam@phi_intercept_prior_variance
      cov_re_prior_scale <- hyperparam@cov_re_prior_scale
      cov_re_prior_df <- hyperparam@cov_re_prior_df
    }

    # Support for unselected parameters
    phi_not_to_select_forced_support <- extract_sub_support(
      model@x_forced_support,
      phi_not_to_select_idx
    )
    x_support_phi_not_to_select <- if (!is_empty_support(phi_not_to_select_forced_support)) {
      matrix(model@x_forced_support[, phi_not_to_select_idx], ncol = n_phi_not_to_select)
    } else {
      NULL
    }

    # Build configuration list with explicit names
    config <- list(
      # === Data ===
      y_series = data@y_series,
      t_series = data@t_series,
      x_candidates = data@x_candidates, # Covariates available for selection
      x_forced = data@x_forced, # Forced covariates
      x_phi_to_select = data@x_phi_to_select, # Covariates for parameters to select
      x_phi_not_to_select = data@x_phi_not_to_select, # Covariates for parameters not to select
      tx_x_phi_to_select = data@tx_x_phi_to_select, # t(X) %*% X for selected parameters
      kron_tx_x_phi_to_select = data@kron_tx_x_phi_to_select, # Block-diagonal matrix
      x_phi_not_to_select_list = data@x_phi_not_to_select_list, # List for sequences
      num_series = length(data@y_series), # <- n
      series_lengths = lengths(data@y_series), # <- ni
      total_observations = sum(lengths(data@y_series)), # <- ntot
      num_covariates_to_select = if (!is.null(data@x_phi_to_select)) dim(data@x_phi_to_select)[2] - 1 else 0,
      num_covariates_not_to_select = if (!is.null(data@x_phi_not_to_select)) dim(data@x_phi_not_to_select)[2] - 1 else 0,

      # === Model info ===
      total_parameters = total_phi,
      parameters_to_select_indices = phi_to_select_idx,
      parameters_not_to_select_indices = phi_not_to_select_idx,
      fixed_parameter_indices = model@phi_fixed_idx,
      model_function = model@model_func,
      x_support_phi_not_to_select = x_support_phi_not_to_select,
      forced_covariates_indices = if (!is.null(x_support_phi_not_to_select)) which(rbind(1, x_support_phi_not_to_select) == 1) else seq(1, length(phi_not_to_select_idx)),


      # === Initialization ===
      init_parameters = init,

      # === Tuning ===
      num_iterations = tuning_algo@niter,
      num_burnin = tuning_algo@nburnin,
      num_mh_iterations = tuning_algo@niter_mh,
      mh_kernel_type = tuning_algo@kernel_mh,
      mh_proposal_scale = tuning_algo@mh_proposal_scale,
      step_size = tuning_algo@step,
      covariance_decay = tuning_algo@covariance_decay,

      # === Hyperparameters (MAP only) ===
      method_type = method_type,
      spike_parameter = spike_parameter, # Spike prior value
      slab_parameter = slab_parameter, # Slab prior value
      residual_variance_prior_shape = residual_variance_prior_shape,
      residual_variance_prior_rate = residual_variance_prior_rate,
      inclusion_prob_prior_a = inclusion_prob_prior_a,
      inclusion_prob_prior_b = inclusion_prob_prior_b,
      phi_intercept_prior_variance = phi_intercept_prior_variance,
      cov_re_prior_scale = cov_re_prior_scale,
      cov_re_prior_df = cov_re_prior_df,

      # === Dimension info ===
      num_parameters_to_select = n_phi_to_select,
      num_parameters_not_to_select = n_phi_not_to_select
    )

    return(config)
  }
)
