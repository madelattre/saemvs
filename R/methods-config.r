## -- Prepare data according to the model

setGeneric(
  "prepare_data",
  function(data, model) {
    standardGeneric("prepare_data")
  }
)


setMethod(
  "prepare_data",
  signature(data = "dataC", model = "modelC"),
  function(data, model) {
    n <- length(data@y_list)

    phi_sel_idx <- model@phi_sel_idx
    x_forced_support <- model@x_forced_support

    # On remplit x_phi_sel
    if (length(phi_sel_idx) > 0) {
      if (empty_support(x_forced_support) == TRUE) {
        x_phi_sel <- cbind(1, data@x_sel)
      } else {
        if (empty_support(x_forced_support[, phi_sel_idx]) == TRUE) {
          x_phi_sel <- cbind(1, data@x_sel)
        } else {
          xf_supp_phi_sel <-
            matrix(x_forced_support[, phi_sel_idx], ncol = length(phi_sel_idx))
          x_sel_forced_idx <-
            extract_raws_with_ones(xf_supp_phi_sel)
          x_phi_sel <- cbind(1, data@x_forced[, x_sel_forced_idx], data@x_sel)
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
    if (empty_support(x_forced_support) == TRUE) {
      x_phi_insel <- matrix(1, nrow = n)
    } else if (empty_support(x_forced_support[, phi_insel_idx]) == TRUE) {
      x_phi_insel <- matrix(1, nrow = n)
    } else {
      xf_supp_phi_insel <- matrix(x_forced_support[, phi_insel_idx],
        ncol = model@phi_dim - length(phi_sel_idx)
      )
      x_insel_forced_idx <-
        extract_raws_with_ones(xf_supp_phi_insel)
      x_phi_insel <- cbind(1, data@x_forced[, x_insel_forced_idx])
    }

    # On en déduit x_phi_insel_list
    if (ncol(x_phi_insel) == 1) {
      x_phi_insel_list <- lapply(
        seq_len(length(data@y_list)),
        function(i) {
          matrix(1, nrow = 1, ncol = 1)
        }
      )
    } else {
      x_phi_insel_list <- from_v_to_x(
        x_phi_insel,
        rbind(1, xf_supp_phi_insel), # besoin de l'intercept
        model@phi_dim - length(phi_sel_idx)
      )
    }

    data_alg <- new("dataAlgo",
      y_list = data@y_list,
      t_list = data@t_list,
      x_sel = data@x_sel,
      x_forced = data@x_forced,
      x_phi_sel = x_phi_sel,
      x_phi_insel = x_phi_insel,
      tx_x_phi_sel = tx_x_phi_sel,
      kron_tx_x_phi_sel = kron_tx_x_phi_sel,
      x_phi_insel_list = x_phi_insel_list
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
  signature(init = "initC", model = "modelC"),
  function(init, model) {
    # A compléter, initialisation de beta_sel par défaut?

    # On remplit beta_hdim, gamma_hdim et alpha


    # S'assurer qu'on ait bien des matrices partout

    if (length(model@phi_sel_idx) == 0) {
      beta_hdim <- NULL
      gamma_hdim <- NULL
      alpha <- NULL
    } else {
      if (empty_support(model@x_forced_support) == TRUE) {
        xf_supp_phi_sel <- NULL
      } else {
        xf_supp_phi_sel <- matrix(
          model@x_forced_support[, model@phi_sel_idx],
          ncol = length(model@phi_sel_idx)
        )
      }
      if (empty_support(xf_supp_phi_sel) == TRUE) {
        beta_hdim <- rbind(
          init@intercept[model@phi_sel_idx],
          init@beta_sel[, model@phi_sel_idx]
        )
      } else {
        raws_bf_phi_sel <- extract_raws_with_ones(xf_supp_phi_sel)
        bf_phi_sel <- init@beta_forced[raws_bf_phi_sel, model@phi_sel_idx]
        beta_hdim <- rbind(
          init@intercept[model@phi_sel_idx],
          bf_phi_sel,
          init@beta_sel[, model@phi_sel_idx]
        )
      }
      gamma_hdim <- matrix(
        init@gamma[model@phi_sel_idx, model@phi_sel_idx],
        ncol = length(model@phi_sel_idx)
      )

      alpha <- rep(0.5, length(model@phi_sel_idx))
    }

    phi_insel_idx <- setdiff(seq(1, model@phi_dim), model@phi_sel_idx)

    if (length(phi_insel_idx) == 0) {
      beta_ldim <- NULL
      gamma_ldim <- NULL
    } else {
      if (empty_support(model@x_forced_support) == TRUE) {
        xf_supp_phi_insel <- NULL
      } else {
        xf_supp_phi_insel <- matrix(
          model@x_forced_support[, phi_insel_idx],
          ncol = length(phi_insel_idx)
        )
      }
      if (empty_support(xf_supp_phi_insel) == TRUE) {
        beta_ldim <- matrix(init@intercept[phi_insel_idx],
          ncol = length(phi_insel_idx)
        )
      } else {
        raws_bf_phi_insel <- extract_raws_with_ones(xf_supp_phi_insel)
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
        init@gamma[phi_insel_idx, phi_insel_idx],
        ncol = length(phi_insel_idx)
      )
    }


    init_alg <- new("initAlgo",
      alpha = alpha,
      sigma2 = init@sigma2,
      gamma_hdim = gamma_hdim,
      gamma_ldim = gamma_ldim,
      beta_hdim = beta_hdim,
      beta_ldim = beta_ldim
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
  signature(hyper = "hyperC", data = "dataAlgo", model = "modelC"),
  function(hyper, data, model) {
    nbs <- length(model@phi_sel_idx)

    if (nbs == 0) { # mle
      hyper <- hyperC(NULL, NULL, NULL)
    } else { # map

      p <- dim(data@x_sel)[2] #- 1

      if (is.null(hyper@a)) {
        hyper@a <- rep(1, nbs)
      }

      if (is.null(hyper@b)) {
        hyper@b <- rep(p, nbs)
      }
    }
    return(hyper)
  }
)


## -- Config

setGeneric(
  "make_config",
  function(data, model, tuning_algo, init, hyperparam) {
    standardGeneric("make_config")
  }
)

setMethod(
  "make_config",
  signature(
    data = "dataAlgo", model = "modelC", tuning_algo = "tuningC",
    init = "initAlgo", hyperparam = "hyperC"
  ),
  function(data, model, tuning_algo, init, hyperparam) {
    q_phi <- model@phi_dim
    index_select <- model@phi_sel_idx
    

    if (length(index_select) == 0) {
      method <- "mle"
      q_hdim <- 0
      q_ldim <- q_phi
      index_unselect <- seq(1, q_phi)
      nu0 <- nu1 <- nsig <- lsig <- a <- b <- sigma2_mu <- sgam <- d <- NULL
    } else {
      method <- "map"
      q_hdim <- length(index_select)
      q_ldim <- q_phi - length(index_select)
      index_unselect <- setdiff(seq(1,q_phi),index_select)
      nu0 <- hyperparam@nu0
      nu1 <- hyperparam@nu1
      nsig <- hyperparam@nsig # nu_sigma
      lsig <- hyperparam@lsig # lb_sigma
      a <- hyperparam@a
      b <- hyperparam@b
      sigma2_mu <- hyperparam@sigma2_mu
      sgam <- hyperparam@sgam # Sigma_Gamma
      d <- hyperparam@d # qu'est-ce que d?
    }

    unselect_support <- 
    extract_sub_support(model@x_forced_support, index_unselect)

    if (!empty_support(unselect_support)) {
      x_support_insel <- matrix(
        model@x_forced_support[, index_unselect],
        ncol = model@phi_dim - length(model@phi_sel_idx)
      )
    } else {
      x_support_insel <- NULL
    }

    return(list(
      method = method,
      index_select = index_select,
      index_unselect = index_unselect,
      q_phi = q_phi,
      q_hdim = q_hdim,
      q_ldim = q_ldim,
      yi = data@y_list,
      n = length(data@y_list),
      ni = lengths(data@y_list),
      ntot = sum(lengths(data@y_list)),
      ti = data@t_list,
      v = data@x_phi_sel,
      tv_v = data@tx_x_phi_sel,
      kron_tv_v = data@kron_tx_x_phi_sel,
      w = data@x_phi_insel,
      x = data@x_phi_insel_list,
      pv = dim(data@x_phi_sel)[2] - 1,
      pw = dim(data@x_phi_insel)[2] - 1,
      g = model@model_func,
      support = x_support_insel,
      # est-ce qu'on a besoin de support? OUI...
      supp_index = which(rbind(
        1,
        x_support_insel
      ) == 1),
      index_fixed = model@phi_fixed_idx,
      q_fixed = length(model@phi_fixed_idx),
      niter = tuning_algo@niter,
      nburnin = tuning_algo@nburnin,
      niter_mh = tuning_algo@niter_mh,
      kernel_mh = tuning_algo@kernel_mh,
      kappa = tuning_algo@kappa_mh,
      step = tuning_algo@step,
      param_init = init,
      tau = tuning_algo@tau,
      nu0 = nu0,
      nu1 = nu1,
      nsig = nsig,
      lsig = lsig,
      a = a,
      b = b,
      sigma2_mu = sigma2_mu,
      sgam = sgam,
      d = d
    ))
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
    data = "dataC", cand_support = "matrix"
  ),
  function(data, cand_support) {
    n <- length(data@y_list)

    lines_with_ones <- which(rowSums(cand_support) > 0)

    v_restricted <- matrix(
      cbind(1, data@x_sel)[, lines_with_ones],
      nrow = n
    )[, -1] # -1 pour ne pas contenir l'intercept

    # Si on ne force pas l'inclusion de variables dans phi_sel
    new_x <- matrix(
      cbind(data@x_forced, v_restricted),
      nrow = n
    )

    new_data <- dataC(
      y = data@y_list,
      t = data@t_list,
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
    model = "modelC", cand_support = "matrix"
  ),
  function(model, cand_support) {
    # Ici il faut tenir compte des indices des phi_s et des phi_ns

    all_phi_index <- seq(1, model@phi_dim)
    index_unselect <- setdiff(all_phi_index, model@phi_sel_idx)
    perm <- c(model@phi_sel_idx, index_unselect)
    inv_perm <- match(seq_along(perm), perm)

    nb_phi_s <- length(model@phi_sel_idx)

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

    new_model <- modelC(
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
    init = "initC", model = "modelC", cand_support = "matrix"
  ),
  # Ici il faut tenir compte des indices des phi_s et des phi_ns
  function(init, model, cand_support) {
    # all_phi_index <- seq(1, model@phi_dim)
    # index_unselect <- setdiff(all_phi_index, model@phi_sel_idx)
    # perm <- c(model@phi_sel_idx, index_unselect)
    # inv_perm <- match(seq_along(perm), perm)

    lines_with_ones <- which(rowSums(cand_support[-1, ]) > 0)

    new_beta_init <- rbind(
      init@beta_forced,
      init@beta_sel[lines_with_ones, ]
    )

    # if (is.null(init@gamma_ldim)) {
    #   new_gamma_init <- init@gamma_hdim
    # } else {
    #   new_gamma_init <- as.matrix(Matrix::bdiag(
    #     init@gamma_hdim,
    #     init@gamma_ldim
    #   ))
    # }

    new_init <- initC(
      intercept = init@intercept,
      beta_forced = new_beta_init,
      gamma = init@gamma,
      # beta_ns = new_beta_init[, inv_perm],
      # gamma_ns = new_gamma_init[inv_perm, inv_perm],
      sigma2 = init@sigma2
    )

    return(new_init)
  }
)
