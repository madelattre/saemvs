## -- Check that the data fits the model

setGeneric(
  "check_data",
  function(data, model) {
    standardGeneric("check_data")
  }
)

setMethod(
  "check_data",
  signature(data = "saemvsData", model = "saemvsModel"),
  function(data, model) {
    if (is.null(data@x_candidates) && (length(model@phi_to_select_idx) > 0)) {
      stop(
        paste0(
          "Parameter selection is requested ('phi_to_select_idx' is not empty), ",
          "but the covariate matrix 'x_candidates' is missing. Please provide ",
          "'x_candidates' with one row per sequence in 'y'."
        )
      )
    }

    if (!is.null(data@x_forced) && !is.null(model@x_forced_support)) {
      if (ncol(data@x_forced) != nrow(model@x_forced_support)) {
        stop(
          paste0(
            "The number of columns in matrix 'x_forced' must be equal ",
            "to the number of lines in 'x_forced_support'."
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
  signature(init = "initC", data = "saemvsData", model = "saemvsModel"),
  function(init, data, model) {
    # A réécrire
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
  signature(hyper = "hyperC", model = "saemvsModel", tuning = "tuningC"),
  function(hyper, model, tuning) {
    nbs <- length(model@phi_to_select_idx)

    if (nbs != 0) { # map
      if (!is.null(hyper@a) && (length(hyper@a) != nbs)) {
        stop(paste0(
          "As ", nbs, " parameters are subject to selection,",
          "hyperparameter 'a' must contain ", nbs, " components."
        ))
      }

      if (!is.null(hyper@b) && (length(hyper@b) != nbs)) {
        stop(paste0(
          "As ", nbs, " parameters are subject to selection,",
          "hyperparameter 'b' must contain ", nbs, " components."
        ))
      }

      # Check !is.null(hyper@sgam) pour éviter les erreurs si vide?
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
      # Vérifier que tuning@nu0_grid n'est pas vide? (length(tuning@nu_grid)>0)
      stop(
        paste0(
          "All spike parameter values in 'nu0_grid' must be smaller than ",
          "slab hyperparameter 'nu1'."
        )
      )
    }
  }
)
