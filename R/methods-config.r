## -- Prepare data according to the model

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
    n <- length(data@y_series)

    phi_sel_idx <- model@phi_to_select_idx
    x_forced_support <- model@x_forced_support

    # On remplit x_phi_sel
    if (length(phi_sel_idx) > 0) {
      if (is_empty_support(x_forced_support) == TRUE) {
        x_phi_sel <- cbind(1, data@x_candidates)
      } else {
        if (is_empty_support(x_forced_support[, phi_sel_idx]) == TRUE) {
          x_phi_sel <- cbind(1, data@x_candidates)
        } else {
          xf_supp_phi_sel <-
            matrix(x_forced_support[, phi_sel_idx], ncol = length(phi_sel_idx))
          x_sel_forced_idx <-
            extract_rows_with_ones(xf_supp_phi_sel)
          x_phi_sel <- cbind(
            1, data@x_forced[, x_sel_forced_idx],
            data@x_candidates
          )
        }
      }
    } else {
      x_phi_sel <- NULL
    }

    # On en déduit tx_x_phi_sel et kron_tx_x_phi_sel

    if (is.null(x_phi_sel)) {
      tx_x_phi_sel <- kron_tx_x_phi_sel <- NULL
    } else {
      tx_x_phi_sel <- t(x_phi_sel) %*% x_phi_sel
      kron_tx_x_phi_sel <- as.matrix(Matrix::bdiag(replicate(
        length(phi_sel_idx),
        tx_x_phi_sel,
        simplify = FALSE
      )))
    }

    # On remplit x_phi_insel
    phi_insel_idx <- setdiff(seq(1, model@phi_dim), phi_sel_idx)
    if (is_empty_support(x_forced_support) == TRUE) {
      x_phi_insel <- matrix(1, nrow = n)
    } else if (is_empty_support(x_forced_support[, phi_insel_idx]) == TRUE) {
      x_phi_insel <- matrix(1, nrow = n)
    } else {
      xf_supp_phi_insel <- matrix(x_forced_support[, phi_insel_idx],
        ncol = model@phi_dim - length(phi_sel_idx)
      )
      x_insel_forced_idx <-
        extract_rows_with_ones(xf_supp_phi_insel)
      x_phi_insel <- cbind(1, data@x_forced[, x_insel_forced_idx])
    }

    # On en déduit x_phi_insel_list
    if (ncol(x_phi_insel) == 1) {
      x_phi_insel_list <- lapply(
        seq_len(length(data@y_series)),
        function(i) {
          matrix(1, nrow = 1, ncol = 1)
        }
      )
    } else {
      x_phi_insel_list <- expand_to_list(
        x_phi_insel,
        rbind(1, xf_supp_phi_insel), # besoin de l'intercept
        model@phi_dim - length(phi_sel_idx)
      )
    }

    data_alg <- new("saemvsProcessedData",
      y_series = data@y_series,
      t_series = data@t_series,
      x_candidates = data@x_candidates,
      x_forced = data@x_forced,
      x_phi_to_select = x_phi_sel,
      x_phi_not_to_select = x_phi_insel,
      tx_x_phi_to_select = tx_x_phi_sel,
      kron_tx_x_phi_to_select = kron_tx_x_phi_sel,
      x_phi_not_to_select_list = x_phi_insel_list
    )
    return(data_alg)
  }
)


## -- Prepare the initial values according to the model

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
    # A compléter, initialisation de beta_sel par défaut?

    # On remplit beta_hdim, gamma_hdim et alpha


    # S'assurer qu'on ait bien des matrices partout

    if (length(model@phi_to_select_idx) == 0) {
      beta_hdim <- NULL
      gamma_hdim <- NULL
      alpha <- NULL
    } else {
      if (is_empty_support(model@x_forced_support) == TRUE) {
        xf_supp_phi_sel <- NULL
      } else {
        xf_supp_phi_sel <- matrix(
          model@x_forced_support[, model@phi_to_select_idx],
          ncol = length(model@phi_to_select_idx)
        )
      }
      if (is_empty_support(xf_supp_phi_sel) == TRUE) {
        beta_hdim <- rbind(
          init@intercept[model@phi_to_select_idx],
          init@beta_candidates[, model@phi_to_select_idx]
        )
      } else {
        raws_bf_phi_sel <- extract_rows_with_ones(xf_supp_phi_sel)
        bf_phi_sel <- init@beta_forced[raws_bf_phi_sel, model@phi_to_select_idx]
        beta_hdim <- rbind(
          init@intercept[model@phi_to_select_idx],
          bf_phi_sel,
          init@beta_candidates[, model@phi_to_select_idx]
        )
      }
      gamma_hdim <- matrix(
        init@cov_re[model@phi_to_select_idx, model@phi_to_select_idx],
        ncol = length(model@phi_to_select_idx)
      )

      alpha <- rep(0.5, length(model@phi_to_select_idx))
    }

    phi_insel_idx <- setdiff(seq(1, model@phi_dim), model@phi_to_select_idx)

    if (length(phi_insel_idx) == 0) {
      beta_ldim <- NULL
      gamma_ldim <- NULL
    } else {
      if (is_empty_support(model@x_forced_support) == TRUE) {
        xf_supp_phi_insel <- NULL
      } else {
        xf_supp_phi_insel <- matrix(
          model@x_forced_support[, phi_insel_idx],
          ncol = length(phi_insel_idx)
        )
      }
      if (is_empty_support(xf_supp_phi_insel) == TRUE) {
        beta_ldim <- matrix(init@intercept[phi_insel_idx],
          ncol = length(phi_insel_idx)
        )
      } else {
        raws_bf_phi_insel <- extract_rows_with_ones(xf_supp_phi_insel)
        bf_phi_insel <- matrix(
          init@beta_forced[raws_bf_phi_insel, phi_insel_idx],
          ncol = length(phi_insel_idx)
        )
        beta_ldim <- rbind(
          init@intercept[phi_insel_idx],
          bf_phi_insel
        )
      }

      gamma_ldim <- matrix(
        init@cov_re[phi_insel_idx, phi_insel_idx],
        ncol = length(phi_insel_idx)
      )
    }


    init_alg <- new("saemvsProcessedInit",
      inclusion_prob = alpha,
      sigma2 = init@sigma2,
      gamma_to_select = gamma_hdim,
      gamma_not_to_select = gamma_ldim,
      beta_to_select = beta_hdim,
      beta_not_to_select = beta_ldim
    )

    return(init_alg)
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
    phi_not_to_select_forced_support <- extract_sub_support(model@x_forced_support, 
    phi_not_to_select_idx)
    x_support_phi_not_to_select <- if (!is_empty_support(phi_not_to_select_forced_support)) {
      matrix(model@x_forced_support[, phi_not_to_select_idx], ncol = n_phi_not_to_select)
    } else NULL

    # Build configuration list with explicit names
    config <- list(
      # === Data ===
      y_series = data@y_series,
      t_series = data@t_series,
      x_candidates = data@x_candidates,           # Covariates available for selection
      x_forced = data@x_forced,                   # Forced covariates
      x_phi_to_select = data@x_phi_to_select,     # Covariates for parameters to select
      x_phi_not_to_select = data@x_phi_not_to_select, # Covariates for parameters not to select
      tx_x_phi_to_select = data@tx_x_phi_to_select,   # t(X) %*% X for selected parameters
      kron_tx_x_phi_to_select = data@kron_tx_x_phi_to_select, # Block-diagonal matrix
      x_phi_not_to_select_list = data@x_phi_not_to_select_list, # List for sequences
      num_series = length(data@y_series),            # <- n
      series_lengths = lengths(data@y_series),      # <- ni
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
      spike_parameter = spike_parameter,                     # Spike prior value
      slab_parameter = slab_parameter,                       # Slab prior value
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


# ## -- Config

# setGeneric(
#   "make_config",
#   function(data, model, tuning_algo, init, hyperparam) {
#     standardGeneric("make_config")
#   }
# )

# setMethod(
#   "make_config",
#   signature(
#     data = "saemvsProcessedData", model = "saemvsModel", tuning_algo = "saemvsTuning",
#     init = "saemvsProcessedInit", hyperparam = "saemvsHyperSlab"
#   ),
#   function(data, model, tuning_algo, init, hyperparam) {
#     q_phi <- model@phi_dim
#     index_select <- model@phi_to_select_idx


#     if (length(index_select) == 0) {
#       method <- "mle"
#       q_hdim <- 0
#       q_ldim <- q_phi
#       index_unselect <- seq(1, q_phi)
#       nu0 <- nu1 <- nsig <- lsig <- a <- b <- sigma2_mu <- sgam <- d <- NULL
#     } else {
#       method <- "map"
#       q_hdim <- length(index_select)
#       q_ldim <- q_phi - length(index_select)
#       index_unselect <- setdiff(seq(1, q_phi), index_select)
#       nu0 <- hyperparam@spike_parameter
#       nu1 <- hyperparam@slab_parameter
#       nsig <- hyperparam@residual_variance_prior_shape
#       lsig <- hyperparam@residual_variance_prior_rate
#       a <- hyperparam@inclusion_prob_prior_a
#       b <- hyperparam@inclusion_prob_prior_b
#       sigma2_mu <- hyperparam@phi_intercept_prior_variance
#       sgam <- hyperparam@cov_re_prior_scale
#       d <- hyperparam@cov_re_prior_df
#     }

#     unselect_support <-
#       extract_sub_support(model@x_forced_support, index_unselect)

#     if (!is_empty_support(unselect_support)) {
#       x_support_insel <- matrix(
#         model@x_forced_support[, index_unselect],
#         ncol = model@phi_dim - length(model@phi_to_select_idx)
#       )
#     } else {
#       x_support_insel <- NULL
#     }

#     return(list(
#       method = method,
#       index_select = index_select,
#       index_unselect = index_unselect,
#       q_phi = q_phi,
#       q_hdim = q_hdim,
#       q_ldim = q_ldim,
#       yi = data@y_series,
#       n = length(data@y_series),
#       ni = lengths(data@y_series),
#       ntot = sum(lengths(data@y_series)),
#       ti = data@t_series,
#       v = data@x_phi_to_select,
#       tv_v = data@tx_x_phi_to_select,
#       kron_tv_v = data@kron_tx_x_phi_to_select,
#       w = data@x_phi_not_to_select,
#       x = data@x_phi_not_to_select_list,
#       pv = dim(data@x_phi_to_select)[2] - 1,
#       pw = dim(data@x_phi_not_to_select)[2] - 1,
#       g = model@model_func,
#       support = x_support_insel,
#       # est-ce qu'on a besoin de support? OUI...
#       supp_index = which(rbind(
#         1,
#         x_support_insel
#       ) == 1),
#       index_fixed = model@phi_fixed_idx,
#       q_fixed = length(model@phi_fixed_idx),
#       niter = tuning_algo@niter,
#       nburnin = tuning_algo@nburnin,
#       niter_mh = tuning_algo@niter_mh,
#       kernel_mh = tuning_algo@kernel_mh,
#       kappa = tuning_algo@mh_proposal_scale,
#       step = tuning_algo@step,
#       param_init = init,
#       tau = tuning_algo@covariance_decay,
#       nu0 = nu0,
#       nu1 = nu1,
#       nsig = nsig,
#       lsig = lsig,
#       a = a,
#       b = b,
#       sigma2_mu = sigma2_mu,
#       sgam = sgam,
#       d = d
#     ))
#   }
# )

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
