## -- Check that the data fits the model

setGeneric(
  "check_data",
  function(data, model) {
    standardGeneric("check_data")
  }
)

setMethod(
  "check_data",
  signature(data = "dataC", model = "modelC"),
  function(data, model) {
    if (is.null(data@v) && (length(model@index_select) > 0)) {
      stop(
        paste0(
          "Parameter selection is requested (index_select is not empty), ",
          "but the covariate matrix 'v' is missing. Please provide 'v' with ",
          "one row per sequence in 'y'."
        )
      )
    }

    if (!is.null(data@w) && !is.null(model@covariate_support)) {
      if (ncol(data@w) != nrow(model@covariate_support)) {
        stop(
          paste0(
            "The number of columns in matrix 'w' must be equal ",
            "to the number of lines in 'covariate_support'."
          )
        )
      }
    }
  }
)

## -- Check that the initial values fit the model

setGeneric(
  "check_init",
  function(init, data, model) {
    standardGeneric("check_init")
  }
)


setMethod(
  "check_init",
  signature(init = "initC", data = "dataC", model = "modelC"),
  function(init, data, model) {
    if (length(model@index_select) == 0) { # mle
      if (!is.null(init@beta_hdim) || !is.null(init@gamma_hdim)) {
        #|| !is.null(init@alpha)
        warning(
          paste0(
            "As 'index_select' is empty, the initial ",
            "values for 'beta_s' and 'gamma_s' ",
            "will not be used."
          )
        )
      }

      if (is.null(init@beta_ldim) || is.null(init@gamma_ldim)) {
        stop(
          paste0(
            "As 'index_select' is empty, initial values for 'beta_ns' ",
            "and 'gamma_ns' must be provided."
          )
        )
      }

      check_beta_gamma(init@beta_ldim, init@gamma_ldim, model@q_phi)
      check_beta_support(init@beta_ldim, model@covariate_support)
    } else { # map

      if (is.null(init@beta_hdim) || is.null(init@gamma_hdim)) {
        #||is.null(init@alpha)
        stop(
          paste0(
            "As 'index_select' is not empty, initial values for 'beta_s'",
            " and 'gamma_s' must be provided."
          )
        )
      }

      if (length(model@index_select) == model@q_phi) {
        if (!is.null(init@beta_ldim) || !is.null(init@gamma_ldim)) {
          warning(
            paste0(
              "As all the parameters are subject to selection",
              "('index_select' = ", seq(1, model@q_phi), "),",
              "the initial values for 'beta_ns' and 'gamma_ns'",
              "will not be used."
            )
          )
        }

        check_beta_gamma(init@beta_hdim, init@gamma_hdim, model@q_phi,
          case = "s"
        )
        check_beta_hdim(init@beta_hdim, ncol(data@v))
      } else {
        check_beta_gamma(
          init@beta_hdim, init@gamma_hdim,
          length(model@index_select),
          case = "s"
        )
        check_beta_gamma(
          init@beta_ldim, init@gamma_ldim,
          model@q_phi - length(model@index_select)
        )
        check_beta_hdim(init@beta_hdim, ncol(data@v))
        check_beta_support(init@beta_ldim, model@covariate_support)
      }

      if (length(init@alpha) != length(model@index_select)) {
        stop(
          paste0(
            "The initial value for 'alpha' must contain ",
            length(model@index_select), " elements (",
            length(model@index_select), " elements in 'index_select'"
          )
        )
      }
    }
  }
)

## -- Check that the hyperparameters fit the model

setGeneric(
  "check_hyper",
  function(hyper, model, tuning) {
    standardGeneric("check_hyper")
  }
)


setMethod(
  "check_hyper",
  signature(hyper = "hyperC", model = "modelC", tuning = "tuningC"),
  function(hyper, model, tuning) {
    nbs <- length(model@index_select)

    # if (nbs == 0) { # mle
    #   if (!is.null(hyper@nu0) || !is.null(hyper@nu1) ||
    #     !is.null(hyper@nsig) || !is.null(hyper@lsig) ||
    #     !is.null(hyper@a) || !is.null(hyper@b) ||
    #     !is.null(hyper@sigma2_mu) || !is.null(hyper@sgam) ||
    #     !is.null(hyper@d)) {
    #     warning(paste0(
    #       "As no parameter is subject to selection,",
    #       " no hyperparameter is required, so the ",
    #       "hyperparameter slots will be set to NULL."
    #     ))
    #   }
    # }
    if (nbs != 0) { # map
      if (!is.null(hyper@a) && (length(hyper@a) != nbs)) {
        stop(paste0(
          "As ", nbs, " parameters are subject to selection,",
          "hyperparameter 'a' must contain ", nbs, " components."
        ))
      }

      if (!is.null(hyper@b) && (length(hyper@a) != nbs)) {
        stop(paste0(
          "As ", nbs, " parameters are subject to selection,",
          "hyperparameter 'b' must contain ", nbs, " components."
        ))
      }

      if ((dim(hyper@sgam)[1] != nbs) || (dim(hyper@sgam)[2] != nbs)) {
        stop(paste0(
          "As ", nbs, " parameters are subject to selection,",
          "hyperparameter 'sgam' must be a squared matrix with ",
          nbs, "columns and ", nbs, " rows."
        ))
      }
    }

    # Check that all nu0 values are smaller than nu1

    if ((!is.null(tuning@nu0_grid)) && (!all(tuning@nu0_grid < hyper@nu1))) {
      stop(
        paste0(
          "All spike parameter values in 'nu0_grid' must be smaller than ",
          "slab hyperparameter 'nu1'."
        )
      )
    }
  }
)
