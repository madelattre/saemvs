#' Sparse Variable Selection with SAEM
#'
#' The `saemvs` function performs stochastic approximation EM (SAEM)
#' combined with variable selection using a spike-and-slab prior.
#' It supports both BIC and extended BIC (e-BIC) as model selection criteria.
#'
#' This method explores a grid of spike variances (`nu0`) in parallel,
#' fits MAP estimates for each candidate support, and then evaluates unique
#' supports with the chosen information criterion to select the best model.
#'
#' @param data An object of class \link[=saemvsData-class]{saemvsData}
#'   containing the dataset (observations and design matrices).
#' @param model An object of class \link[=saemvsModel-class]{saemvsModel}
#'   specifying the model structure and the indices of parameters
#'   that are candidates for selection (`phi_to_select_idx`).
#' @param init An object of class \link[=saemvsInit-class]{saemvsInit}
#'   providing the initialization for the SAEM algorithm.
#' @param tuning_algo An object of class \link[=saemvsTuning-class]{saemvsTuning}
#'   containing algorithmic tuning parameters such as
#'   the spike grid (`spike_values_grid`), number of workers, and random seed.
#' @param hyperparam An object of class \link[=saemvsHyperSlab-class]{saemvsHyperSlab}
#'   containing the hyperparameters for the slab prior.
#' @param pen Character string indicating the model selection criterion.
#'   Must be either `"BIC"` or `"e-BIC"`.
#'
#' @details
#' The procedure is carried out in two steps:
#' \enumerate{
#'   \item For each spike variance candidate (`nu0`), a MAP estimation
#'   is performed in parallel, producing supports, thresholds, and MAP estimates.
#'   \item Unique supports are identified, and each support is re-evaluated
#'   with the chosen information criterion (`pen`) to select the most
#'   appropriate model.
#' }
#'
#' Parallel execution is handled using the `future` and `furrr` packages.
#' The model function (C++ code) must be compiled on each worker before execution.
#'
#' @return An object of class \code{\linkS4class{saemvsResults}} containing:
#' \itemize{
#'   \item `criterion`: The selection criterion used (`"BIC"` or `"e-BIC"`).
#'   \item `criterion_values`: Values of the criterion for each unique support.
#'   \item `thresholds`: Thresholds computed for each spike value.
#'   \item `beta_map`: MAP estimates of regression parameters for each spike value.
#'   \item `mle_estimates`: MLE estimates for each unique support.
#'   \item `support`: List of supports obtained across spike values.
#'   \item `unique_support`: List of unique supports identified.
#'   \item `support_mapping`: Mapping from each run to the corresponding unique support.
#'   \item `spike_values_grid`: The grid of spike variances used.
#'   \item `phi_fixed_idx`: Indices of fixed parameters (not subject to selection).
#'   \item `forced_variables_idx`: For each unique support, indices of covariates
#'     that were forced into the model (always included).
#'   \item `selected_variables_idx`: For each unique support, indices of covariates
#'     that were actively selected by the algorithm (excluding forced ones).
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming data_obj, model_obj, init_obj, tuning_obj, hyper_obj are created
#' result <- saemvs(
#'   data = data_obj,
#'   model = model_obj,
#'   init = init_obj,
#'   tuning_algo = tuning_obj,
#'   hyperparam = hyper_obj,
#'   pen = "e-BIC"
#' )
#' }
#'
#' @seealso \code{\linkS4class{saemvsResults}}, \code{\link{run_saem}}
#' @export
setGeneric(
  "saemvs",
  function(data, model, init, tuning_algo, hyperparam, pen) {
    standardGeneric("saemvs")
  }
)

#' @rdname saemvs
#' @exportMethod saemvs
setMethod(
  "saemvs",
  signature(
    data = "saemvsData", model = "saemvsModel", init = "saemvsInit",
    tuning_algo = "saemvsTuning", hyperparam = "saemvsHyperSlab",
    pen = "character"
  ),
  function(data, model, init, tuning_algo, hyperparam, pen) {
    # --- Step 1: Argument checks ---
    if (!pen %in% c("e-BIC", "BIC")) {
      stop(
        paste0(
          "Criterion ", pen, " is not supported. ",
          "'pen' must be either 'e-BIC' or 'BIC'."
        )
      )
    }

    if (is.null(model@phi_to_select_idx) || (length(model@phi_to_select_idx) == 0)) {
      stop(
        paste0(
          "'phi_to_select_idx' must contain at least one parameter index",
          " for variable selection."
        )
      )
    }

    # --- Step 2: Parallelization setup ---
    future::plan(future::multisession, workers = tuning_algo@nb_workers)

    spike_values_grid <- tuning_algo@spike_values_grid

    # --- Step 3: First SAEMVS phase (MAP estimation) ---
    message("SAEMVS: step 1/2")
    progressr::with_progress({
      p <- progressr::progressor(along = spike_values_grid)

      map_runs <- furrr::future_map(
        spike_values_grid,
        function(nu0) {
          # C++ model must be compiled on each worker
          compile_model(model@model_func)
          fh <- saemvsHyperSpikeAndSlab(nu0, hyperparam)
          res <- saemvs_one_map_run(data, model, init, tuning_algo, fh)
          p()
          res
        },
        .options = furrr::furrr_options(seed = tuning_algo@seed)
      )
    })

    # Extract MAP results
    support <- lapply(
      seq_along(map_runs),
      function(x) map_runs[[x]]$support
    )

    thresholds <- lapply(
      seq_along(map_runs),
      function(x) map_runs[[x]]$threshold
    )

    beta_map <- lapply(
      seq_along(map_runs),
      function(x) map_runs[[x]]$beta
    )

    # Deduplicate supports
    support_hashes <- sapply(support, digest::digest)
    unique_support_indices <- which(!duplicated(support_hashes) == TRUE)
    hash_to_compact_index <- stats::setNames(
      seq_along(unique_support_indices),
      support_hashes[unique_support_indices]
    )


    support_index_mapping <- unname(sapply(
      support_hashes,
      function(h) hash_to_compact_index[[h]]
    ))


    # --- Step 4: Second SAEMVS phase (criterion evaluation: BIC/e-BIC) ---
    message("SAEMVS: step 2/2")
    progressr::with_progress({
      p <- progressr::progressor(along = length(unique_support_indices))
      criterion_results <- furrr::future_map(
        unique_support_indices,
        function(k) {
          compile_model(model@model_func)
          res <- saemvs_one_ic_run(
            k, support, data, model, init, tuning_algo, pen
          )
          p()
          res
        },
        .options = furrr::furrr_options(seed = tuning_algo@seed)
      )
    })

    # Extract evaluation results
    criterion_values <- unlist(lapply(
      seq_along(criterion_results),
      function(x) criterion_results[[x]]$ll
    ))

    mle <- lapply(
      seq_along(criterion_results),
      function(x) criterion_results[[x]]$mle_param
    )

    forced_variables_idx <- lapply(
      seq_along(criterion_results),
      function(x) criterion_results[[x]]$forced_variables_idx
    )

    selected_variables_idx <- lapply(
      seq_along(criterion_results),
      function(x) criterion_results[[x]]$selected_variables_idx
    )

    # --- Step 5: Assemble results object including forced and selected variable indices ---

    res <- methods::new(
      "saemvsResults",
      criterion = pen,
      criterion_values = criterion_values,
      thresholds = thresholds,
      beta_map = beta_map,
      mle_estimates = mle,
      support = support,
      unique_support = support[unique_support_indices],
      support_mapping = support_index_mapping,
      spike_values_grid = tuning_algo@spike_values_grid,
      phi_fixed_idx = model@phi_fixed_idx,
      forced_variables_idx = forced_variables_idx,
      selected_variables_idx = selected_variables_idx
    )

    return(res)
  }
)



#' Internal: Single MAP Run of SAEMVS
#'
#' Performs one MAP estimation run of the SAEM algorithm for a given
#' spike-and-slab prior configuration.
#' This function is used internally by \code{\link{saemvs}} to explore
#' different spike variances and extract candidate supports.
#'
#' @param data An object of class \link[=saemvsData-class]{saemvsData}
#'   containing the dataset (observations and design matrices).
#' @param model An object of class \link[=saemvsModel-class]{saemvsModel}
#'   specifying the model structure and indices of parameters subject
#'   to selection (\code{phi_to_select_idx}).
#' @param init An object of class \link[=saemvsInit-class]{saemvsInit}
#'   providing initial values for the SAEM algorithm.
#' @param tuning_algo An object of class \link[=saemvsTuning-class]{saemvsTuning}
#'   defining algorithmic parameters such as the number of iterations (\code{niter}).
#' @param hyperparam An object of class \link[=saemvsHyperSpikeAndSlab-class]{saemvsHyperSpikeAndSlab}
#'   containing hyperparameters of the spike-and-slab prior
#'   (spike variance, slab variance, and mixing proportion).
#'
#' @details
#' The function runs the SAEM algorithm via \code{\link{run_saem}},
#' then constructs the MAP estimates of coefficients and derives
#' the support set by thresholding \eqn{\beta} coefficients.
#'
#' Special handling is performed for parameters that are forced
#' to be in the model (\code{x_forced_support}).
#' These are always marked as selected in the support matrix.
#'
#' The returned support matrix includes the intercept term as the first row.
#'
#' @return A \code{list} with elements:
#' \itemize{
#'   \item \code{threshold} Vector of threshold values applied to coefficients.
#'   \item \code{support} Logical matrix indicating selected variables
#'         across iterations (including intercept).
#'   \item \code{beta} Matrix of MAP coefficient estimates.
#' }
#'
#' @keywords internal
#' @seealso \code{\link{saemvs}}, \code{\link{run_saem}}
setGeneric(
  "saemvs_one_map_run",
  function(data, model, init, tuning_algo, hyperparam) {
    standardGeneric("saemvs_one_map_run")
  }
)

setMethod(
  "saemvs_one_map_run",
  signature(
    data = "saemvsData", model = "saemvsModel", init = "saemvsInit",
    tuning_algo = "saemvsTuning", hyperparam = "saemvsHyperSpikeAndSlab"
  ),
  function(data, model, init, tuning_algo, hyperparam) {
    map <- run_saem(data, model, init, tuning_algo, hyperparam)
    forced_support <- extract_sub_support(
      model@x_forced_support,
      model@phi_to_select_idx
    )

    n_phi_select <- length(model@phi_to_select_idx)
    niter <- tuning_algo@niter
    if (is_empty_support(forced_support) == TRUE) {
      p <- dim(data@x_candidates)[2]
    } else {
      p <- dim(data@x_candidates)[2] + dim(forced_support)[1]
    }


    beta_map <- map$beta_to_select[[niter + 1]][-1, ]
    alpha_map <- map$alpha[[niter + 1]]
    threshold_matrix <- matrix(
      rep(threshold(
        hyperparam@slab_parameter,
        hyperparam@spike_parameter, alpha_map
      ), p),
      nrow = p, byrow = TRUE
    )
    support <- rbind(
      rep(TRUE, n_phi_select),
      abs(beta_map) >= threshold_matrix
    )

    if (!is_empty_support(forced_support)) {
      support[2:(1 + dim(forced_support)[1]), ] <- TRUE
    }

    res <- list(
      threshold = threshold_matrix[1, ],
      support = support,
      beta = beta_map
    )
    return(res)
  }
)

#' Internal: EBIC/BIC Evaluation for a Candidate Support
#'
#' Performs a maximum likelihood estimation (MLE) of the model restricted to
#' a given support set of parameters, and computes the corresponding
#' information criterion (BIC or extended BIC).
#'
#' This function is used internally by \code{\link{saemvs}} in the second step
#' of the SAEMVS algorithm to rank candidate supports obtained from MAP runs.
#'
#' @param k Integer index of the candidate support to be evaluated.
#' @param support A list of logical or numeric matrices, each encoding a
#'   candidate support obtained from MAP estimation (\code{saemvs_one_map_run}).
#' @param data An object of class \link[=saemvsData-class]{saemvsData}, the dataset.
#' @param model An object of class \link[=saemvsModel-class]{saemvsModel}, specifying
#'   model structure and indices of parameters subject to selection.
#' @param init An object of class \link[=saemvsInit-class]{saemvsInit}, providing
#'   initialization for SAEM.
#' @param tuning_algo An object of class \link[=saemvsTuning-class]{saemvsTuning},
#'   containing algorithmic tuning parameters (e.g. number of iterations).
#' @param pen Character string, the information criterion to use.
#'   Must be either \code{"BIC"} or \code{"e-BIC"}.
#'
#' @details
#' The function restricts the model to the variables specified in
#' \code{support[[k]]}, taking into account variables that are always forced
#' into the model (\code{x_forced_support}).
#' The SAEM algorithm is then rerun under this restricted model to obtain
#' MLE estimates of the parameters.
#' Finally, the log-likelihood penalized by the chosen criterion (BIC/e-BIC)
#' is computed via \code{\link{loglik}}.
#'
#' @return A \code{list} with elements:
#' \itemize{
#'   \item \code{ll} Numeric, the penalized log-likelihood value (BIC/e-BIC).
#'   \item \code{mle_param} List of estimated parameters under the restricted model:
#'     \itemize{
#'       \item \code{beta} Estimated regression coefficients.
#'       \item \code{gamma} Estimated variances for the random effects.
#'       \item \code{sigma2} Estimated residual variance.
#'     }
#' }
#'
#' @keywords internal
#' @seealso \code{\link{saemvs}}, \code{\link{saemvs_one_map_run}}, \code{\link{loglik}}
setGeneric(
  "saemvs_one_ic_run",
  function(k, support, data, model, init, tuning_algo, pen) {
    standardGeneric("saemvs_one_ic_run")
  }
)

setMethod(
  "saemvs_one_ic_run",
  signature(
    k = "integer", support = "list",
    data = "saemvsData", model = "saemvsModel", init = "saemvsInit",
    tuning_algo = "saemvsTuning", pen = "character"
  ),
  function(k, support, data, model, init, tuning_algo, pen) {
    p <- dim(data@x_candidates)[2]

    forced_support <- extract_sub_support(
      model@x_forced_support,
      model@phi_to_select_idx
    )

    if (is_empty_support(forced_support)) {
      cand_support <- matrix(as.numeric(support[[k]]),
        nrow = nrow(support[[k]])
      )
      forced_rows <- integer(0) # no forced covariates
      # selected_rows <- seq_len(nrow(cand_support))
      active_candidate_idx <- which(rowSums(cand_support[-1, , drop = FALSE]) > 0)
      selected_rows <- active_candidate_idx
    } else {
      cand_support <- matrix(as.numeric(support[[k]]),
        nrow = nrow(support[[k]])
      )
      idx_forced_phi_sel <- seq(1, dim(forced_support)[1])
      cand_support <- cand_support[-idx_forced_phi_sel, ]
      forced_rows <- idx_forced_phi_sel
      active_candidate_idx <- which(rowSums(cand_support[-1, , drop = FALSE]) > 0)
      selected_rows <- active_candidate_idx
    }


    new_data <- map_to_mle_data(data, cand_support)
    new_model <- map_to_mle_model(model, cand_support)
    new_init <- map_to_mle_init(init, model, cand_support)
    new_hyperparam <- saemvsHyperSlab(NULL, NULL, NULL)
    new_full_hyperparam <- saemvsHyperSpikeAndSlab(NULL, new_hyperparam)

    mle <- run_saem(
      new_data, new_model, new_init, tuning_algo, new_full_hyperparam
    )

    mle_param <- list(
      beta   = mle$beta_not_to_select[[tuning_algo@niter + 1]],
      gamma  = mle$gamma_not_to_select[[tuning_algo@niter + 1]],
      sigma2 = mle$sigma2[tuning_algo@niter + 1]
    )

    ll <- loglik(new_data, new_model, tuning_algo, mle_param, pen, p)


    return(list(ll = ll, mle_param = mle_param, forced_variables_idx = forced_rows, selected_variables_idx = selected_rows))
  }
)


#' Run a Single SAEMVS Test Fit
#'
#' Executes a simplified SAEMVS procedure for testing purposes.
#'
#' @param data A \link[=saemvsData-class]{saemvsData} object containing the observed response series.
#' @param model A \link[=saemvsModel-class]{saemvsModel} object defining the model structure.
#' @param init A \link[=saemvsInit-class]{saemvsInit} object providing initial parameter values.
#' @param tuning_algo A \link[=saemvsTuning-class]{saemvsTuning} object specifying algorithm hyperparameters.
#' @param hyperparam A \link[=saemvsHyperSlab-class]{saemvsHyperSlab} object defining the slab prior.
#'
#' @return A \link[=saemResults-class]{saemResults} object containing estimated parameters.
#'
#' @details
#' This function allows a quick fit of SAEMVS to inspect the evolution of parameter estimates.
#' The function performs the following steps:
#' \enumerate{
#'   \item Compile the R model function to C++ using `compile_model`.
#'   \item Construct a spike-and-slab hyperparameter object.
#'   \item Run the SAEM algorithm using `run_saem`.
#'   \item Return parameter estimates packaged in a `saemResults` object.
#' }
#'
#' @examples
#' \dontrun{
#' test_results <- test_saemvs(data_obj, model_obj, init_obj, tuning_obj, hyper_obj)
#' summary(test_results)
#' }
#'
#' @export
setGeneric(
  "test_saemvs",
  function(data, model, init, tuning_algo, hyperparam) {
    standardGeneric("test_saemvs")
  }
)

#' @rdname test_saemvs
#' @exportMethod test_saemvs
setMethod(
  "test_saemvs",
  signature(
    data = "saemvsData", model = "saemvsModel", init = "saemvsInit",
    tuning_algo = "saemvsTuning", hyperparam = "saemvsHyperSlab"
  ),
  function(data, model, init, tuning_algo, hyperparam) {
    # Compile the model function to C++ for efficiency
    compile_model(model@model_func)

    # Construct full spike-and-slab hyperparameter object for SAEM
    full_hyperparam <- saemvsHyperSpikeAndSlab(
      tuning_algo@spike_values_grid[1],
      hyperparam
    )

    # Run the SAEM algorithm to estimate parameters
    saem_state <- run_saem(data, model, init, tuning_algo, full_hyperparam)

    # Package results into saemResults object
    res <- methods::new(
      "saemResults",
      beta_to_select = saem_state$beta_to_select,
      beta_not_to_select = saem_state$beta_not_to_select,
      gamma_to_select = saem_state$gamma_to_select,
      gamma_not_to_select = saem_state$gamma_not_to_select,
      sigma2 = saem_state$sigma2
    )

    return(res)
  }
)
