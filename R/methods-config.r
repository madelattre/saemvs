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
    if (length(model@index_select) == 0) {
      # mle

      v <- tv_v <- kron_tv_v <- NULL

      if (is.null(data@w)) { # No fixed covariate, intercept only
        data@w <- matrix(1, nrow = length(data@y_list), ncol = 1)
        data@x <- lapply(
          seq_len(length(data@y_list)),
          function(i) {
            matrix(1, nrow = 1, ncol = 1)
          }
        )
      } else {
        data@w <- cbind(1, data@w)
        if (is.null(model@covariate_support)) {
          supp <- matrix(1, ncol = model@q_phi, nrow = 1)
        } else {
          supp <- rbind(1, model@covariate_support)
        }
        data@x <- from_v_to_x(
          data@w, supp,
          model@q_phi
        )
      }
    } else if (length(model@index_select) == model@q_phi) {
      # Covariate matrix for the selection part of the model
      data@v <- cbind(1, data@v)
      data@tv_v <- t(data@v) %*% data@v
      data@kron_tv_v <- as.matrix(Matrix::bdiag(replicate(
        length(model@index_select),
        data@tv_v,
        simplify = FALSE
      )))

      data@w <- NULL
      data@x <- NULL
    } else {
      # map

      # Covariate matrix for the selection part of the model
      data@v <- cbind(1, data@v)
      data@tv_v <- t(data@v) %*% data@v
      data@kron_tv_v <- as.matrix(Matrix::bdiag(replicate(
        length(model@index_select),
        data@tv_v,
        simplify = FALSE
      )))

      # Covariate matrix for the effects with no covariate selection
      if (is.null(data@w)) { # No fixed covariate, intercept only
        data@w <- matrix(1, nrow = length(data@y_list), ncol = 1)
        data@x <- lapply(
          seq_len(length(data@y_list)),
          function(i) {
            matrix(1, nrow = 1, ncol = 1)
          }
        )
      } else {
        data@w <- cbind(1, data@w)
        if (is.null(model@covariate_support)) {
          supp <- matrix(1, nrow = model@q_phi -
            length(model@index_select), ncol = 1)
        } else {
          supp <- rbind(1, model@covariate_support)
        }
        data@x <- from_v_to_x(
          data@w, supp,
          model@q_phi -
            length(model@index_select)
        )
      }
    }

    return(data)
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
    if (length(model@index_select) > 0) {
      init@alpha <- rep(0.5, length(model@index_select))
    }

    if ((length(model@index_select) != model@q_phi) &&
      (!is.null(model@covariate_support))) {
      index <- which(model@covariate_support == 0)
      # if (init@beta_ldim[index] != 0) {
      #   init@beta_ldim[index] <- 0
      #   # Afficher un warning?
      # }
    }

    return(init)
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
  signature(hyper = "hyperC", data = "dataC", model = "modelC"),
  function(hyper, data, model) {
    nbs <- length(model@index_select)

    if (nbs == 0) { # mle
      hyper <- hyperC(NULL, NULL, NULL)
    } else { # map

      p <- dim(data@v)[2] - 1

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
    data = "dataC", model = "modelC", tuning_algo = "tuningC",
    init = "initC", hyperparam = "hyperC"
  ),
  function(data, model, tuning_algo, init, hyperparam) {
    q_phi <- model@q_phi
    index_select <- model@index_select

    if (length(index_select) == 0) {
      method <- "mle"
      q_hdim <- 0
      q_ldim <- q_phi
      index_unselect <- seq(1, q_phi)
      nu0 <- nu1 <- nsig <- lsig <- a <- b <- sigma2_mu <- sgam <- d <- NULL
    } else {
      method <- "map"
      q_hdim <- length(index_select)
      q_ldim <- model@q_phi - length(index_select)
      index_unselect <- seq(1, q_phi)[-index_select]
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
      v = data@v,
      tv_v = data@tv_v,
      kron_tv_v = data@kron_tv_v,
      w = data@w,
      x = data@x,
      pv = dim(data@v)[2] - 1,
      pw = dim(data@w)[2] - 1,
      g = model@model_func,
      support = model@covariate_support, # est-ce qu'on a besoin de support?
      supp_index = which(rbind(1, model@covariate_support) == 1),
      index_fixed = model@index_fixed,
      q_fixed = length(model@index_fixed),
      # covariance_structure = model$covariance_structure, # to be removed?,
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
