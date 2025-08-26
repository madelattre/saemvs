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
#'
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
#' data_obj <- new("saemvsData", ...) # create data
#' model_obj <- new("saemvsModel", ...) # create model
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
        function(i) matrix(1, nrow = 1, ncol = 1)
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
    data_processed <- new("saemvsProcessedData",
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
    phi_to_select_idx   <- model@phi_to_select_idx
    phi_not_to_select_idx <- setdiff(seq_len(model@phi_dim), phi_to_select_idx)
    forced_support       <- model@x_forced_support

    ## === 1. Prepare beta and gamma for parameters subject to selection ===
    if (length(phi_to_select_idx) == 0) {
      beta_to_select  <- NULL
      gamma_to_select <- NULL
      inclusion_prob  <- NULL
    } else {
      # Extract forced support for parameters subject to selection
      if (is_empty_support(forced_support)) {
        xf_sel <- NULL
      } else {
        xf_sel <- matrix(forced_support[, phi_to_select_idx],
                         ncol = length(phi_to_select_idx))
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

      # inclusion probabilities (default 0.5)
      inclusion_prob <- rep(0.5, length(phi_to_select_idx))
    }

    ## === 2. Prepare beta and gamma for parameters NOT subject to selection ===
    if (length(phi_not_to_select_idx) == 0) {
      beta_not_to_select  <- NULL
      gamma_not_to_select <- NULL
    } else {
      if (is_empty_support(forced_support)) {
        xf_not_sel <- NULL
      } else {
        xf_not_sel <- matrix(forced_support[, phi_not_to_select_idx],
                             ncol = length(phi_not_to_select_idx))
      }

      # Compose beta_not_to_select
      if (is_empty_support(xf_not_sel)) {
        beta_not_to_select <- matrix(init@intercept[phi_not_to_select_idx],
                                     ncol = length(phi_not_to_select_idx))
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
      gamma_not_to_select <- matrix(
        init@cov_re[phi_not_to_select_idx, phi_not_to_select_idx, drop = FALSE],
        ncol = length(phi_not_to_select_idx)
      )
    }

    ## === 3. Construct processed init object ===
    init_processed <- new("saemvsProcessedInit",
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

## -- Complete the hyperparameters according to the model

setGeneric(
  "prepare_hyper",
  function(hyper, data, model) {
    standardGeneric("prepare_hyper")
  }
)

setMethod(
  "prepare_hyper",
  signature(hyper = "saemvsHyperSlab", data = "saemvsProcessedData", model = "saemvsModel"),
  function(hyper, data, model) {
    nbs <- length(model@phi_to_select_idx)

    if (nbs == 0) { # mle
      hyper <- saemvsHyperSlab(NULL, NULL, NULL)
    } else { # map

      p <- dim(data@x_candidates)[2] #- 1

      if (is.null(hyper@inclusion_prob_prior_a)) {
        hyper@inclusion_prob_prior_a <- rep(1, nbs)
      }

      if (is.null(hyper@inclusion_prob_prior_b)) {
        hyper@inclusion_prob_prior_b <- rep(p, nbs)
      }
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
      forced_covariates_indices = if (!is.null(x_support_phi_not_to_select)) which(rbind(1, x_support_phi_not_to_select) == 1) else matrix(1, ncol = length(phi_not_to_select_idx), nrow = 1),


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


## -- Prepare data, model and initialization from map to mle

setGeneric(
  "from_map_to_mle_data",
  function(data, cand_support) {
    standardGeneric("from_map_to_mle_data")
  }
)

setMethod(
  "from_map_to_mle_data",
  signature(
    data = "saemvsData", cand_support = "matrix"
  ),
  function(data, cand_support) {
    n <- length(data@y_series)

    lines_with_ones <- which(rowSums(cand_support) > 0)

    v_restricted <- matrix(
      cbind(1, data@x_candidates)[, lines_with_ones],
      nrow = n
    )[, -1] # -1 pour ne pas contenir l'intercept

    # Si on ne force pas l'inclusion de variables dans phi_sel
    new_x <- matrix(
      cbind(data@x_forced, v_restricted),
      nrow = n
    )

    new_data <- saemvsData(
      y = data@y_series,
      t = data@t_series,
      x_forced = new_x
    )

    return(new_data)
  }
)

setGeneric(
  "from_map_to_mle_model",
  function(model, cand_support) {
    standardGeneric("from_map_to_mle_model")
  }
)

setMethod(
  "from_map_to_mle_model",
  signature(
    model = "saemvsModel", cand_support = "matrix"
  ),
  function(model, cand_support) {
    # Ici il faut tenir compte des indices des phi_s et des phi_ns

    all_phi_index <- seq(1, model@phi_dim)
    index_unselect <- setdiff(all_phi_index, model@phi_to_select_idx)
    perm <- c(model@phi_to_select_idx, index_unselect)
    inv_perm <- match(seq_along(perm), perm)

    nb_phi_s <- length(model@phi_to_select_idx)

    lines_with_ones <- which(rowSums(cand_support) > 0)

    covariate_restricted_support <- cand_support[lines_with_ones, ]

    if (is.matrix(covariate_restricted_support)) {
      selected_support_matrix <- matrix(
        covariate_restricted_support[-1, ],
        ncol = nb_phi_s,
        nrow = length(covariate_restricted_support[-1, ]) / nb_phi_s
      )
    } else {
      selected_support_matrix <- NULL
    }

    new_covariate_support <- merge_support(
      model@x_forced_support,
      selected_support_matrix,
      # data@w,
      nb_phi_s,
      model@phi_dim - nb_phi_s,
      inv_perm
    )

    new_model <- saemvsModel(
      g = model@model_func,
      phi_dim = model@phi_dim,
      phi_fixed_idx = model@phi_fixed_idx,
      x_forced_support = new_covariate_support
    )

    return(new_model)
  }
)


setGeneric(
  "from_map_to_mle_init",
  function(init, model, cand_support) {
    standardGeneric("from_map_to_mle_init")
  }
)

setMethod(
  "from_map_to_mle_init",
  signature(
    init = "saemvsInit", model = "saemvsModel", cand_support = "matrix"
  ),
  # Ici il faut tenir compte des indices des phi_s et des phi_ns
  function(init, model, cand_support) {
    lines_with_ones <- which(rowSums(cand_support[-1, ]) > 0)

    new_beta_init <- rbind(
      init@beta_forced,
      init@beta_candidates[lines_with_ones, ]
    )

    # if (is.null(init@gamma_ldim)) {
    #   new_gamma_init <- init@gamma_hdim
    # } else {
    #   new_gamma_init <- as.matrix(Matrix::bdiag(
    #     init@gamma_hdim,
    #     init@gamma_ldim
    #   ))
    # }

    new_init <- saemvsInit(
      intercept = init@intercept,
      beta_forced = new_beta_init,
      cov_re = init@cov_re,
      sigma2 = init@sigma2
    )

    return(new_init)
  }
)
