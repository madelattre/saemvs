#' Sparse Variable Selection with SAEM
#'
#' The `saemvs` function performs stochastic approximation EM (SAEM)
#' combined with variable selection using a spike-and-slab prior.
#' It supports both BIC and extended BIC (e-BIC) as model selection criteria.
#'
#' This method explores a grid of spike variances (`nu0`) in parallel,
#' fits MAP estimates for each spike value, derives candidate supports,
#' and then evaluates unique supports with the chosen information
#' criterion to select the best model.
#'
#' @param data An object of class \link[=saemvsData-class]{saemvsData}
#'   containing the dataset (observations and design matrices).
#' @param model An object of class \link[=saemvsModel-class]{saemvsModel}
#'   specifying the model structure and the indices of parameters
#'   that are candidates for selection (`phi_to_select_idx`).
#' @param init An object of class \link[=saemvsInit-class]{saemvsInit}
#'   providing the initialization for the SAEM algorithm.
#' @param tuning_algo An object of class
#' \link[=saemvsTuning-class]{saemvsTuning} containing algorithmic tuning
#' parameters such as the spike grid (`spike_values_grid`), number of workers,
#' and random seed.
#' @param hyperparam An object of class
#' \link[=saemvsHyperSlab-class]{saemvsHyperSlab} containing the
#' hyperparameters for the slab prior.
#' @param pen Character string indicating the model selection criterion. Must be
#' either `"BIC"` or `"e-BIC"`.
#' @param use_cpp Logical. If `TRUE` (default), the C++ backend of the model
#' function is compiled and used. If `FALSE`, a non-C++ (R-based) backend
#' may be used. The model function is compiled on each worker when running
#' in parallel.

#'
#' @details
#' The procedure is carried out in two steps:
#' \enumerate{
#'   \item For each spike variance candidate (`nu0`), a MAP estimation is
#'   performed in parallel, producing supports, thresholds, and MAP estimates.
#'   \item Unique supports are identified, and each support is re-evaluated
#'   with the chosen information criterion (`pen`) to select the most
#'   appropriate model.
#' }
#' A support corresponds to a specific subset of parameters (among
#' `phi_to_select_idx`) that are included in the model, typically represented
#' as a binary inclusion pattern obtained after thresholding the MAP estimates.
#' For a given support, covariates are partitioned into:
#' \itemize{
#'   \item \emph{forced variables}, which are always included in the model
#'   (e.g., fixed effects or parameters not subject to spike-and-slab
#'  selection),
#'   \item \emph{selected variables}, which are actively selected by the
#'   spike-and-slab procedure among the candidate parameters.
#' }

#'
#' Parallel execution is handled using the `future` and `furrr` packages.
#' When `use_cpp = TRUE`, the model C++ code is compiled multiple times:
#' once before the first step, and again within each parallel worker during
#' both the MAP estimation step and the information-criterion evaluation step.
#' This ensures that each worker has access to a valid compiled backend.

#' @return An object of class \code{\linkS4class{saemvsResults}} containing:
#' \itemize{
#'   \item `criterion`: The selection criterion used (`"BIC"` or `"e-BIC"`).
#'   \item `criterion_values`: Values of the selected information criterion
#'   (BIC or e-BIC), computed from the log-likelihood of the MLE associated
#'   with each unique support.
#'   \item `thresholds`: Thresholds computed for each spike value.
#'   \item `beta_map`: MAP estimates of regression parameters for each spike
#'   value.
#'   \item `mle_estimates`: MLE estimates for each unique support.
#'   \item `support`: List of supports obtained across spike values.
#'   \item `unique_support`: List of unique supports identified.
#'   \item `support_mapping`: Mapping from each run to the corresponding unique
#'   support.
#'   \item `spike_values_grid`: The grid of spike variances used.
#'   \item `phi_fixed_idx`: Indices of fixed parameters.
#'   \item `phi_to_select_idx`: Indices of parameters subject to selection.
#'   \item `forced_variables_idx`: For each unique support, indices of
#'   covariates that were forced into the model (always included).
#'   \item `selected_variables_idx`: For each unique support, indices of
#'   covariates that were actively selected by the algorithm (excluding forced
#'   ones).
#' }

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
  function(
    data, model, init, tuning_algo, hyperparam, pen = "e-BIC", use_cpp = TRUE
  ) {
    standardGeneric("saemvs")
  }
)


#' @rdname saemvs
#' @exportMethod saemvs
setMethod(
  "saemvs",
  signature(
    data = "saemvsData", model = "saemvsModel", init = "saemvsInit",
    tuning_algo = "saemvsTuning", hyperparam = "saemvsHyperSlab"
  ),
  function(
    data, model, init, tuning_algo, hyperparam, pen = "e-BIC", use_cpp = TRUE
  ) {
    # --- Step 1: Argument checks ---
    if (!pen %in% c("e-BIC", "BIC")) {
      stop(
        paste0(
          "Criterion ", pen, " is not supported. ",
          "'pen' must be either 'e-BIC' or 'BIC'."
        )
      )
    }

    if (is.null(model@phi_to_select) ||
          (length(model@phi_to_select) == 0)) {
      stop(
        paste0(
          "'phi_to_select' must contain at least one parameter",
          " for variable selection."
        )
      )
    }

    spike_values_grid <- tuning_algo@spike_values_grid

    nb_workers <- min(tuning_algo@nb_workers, length(spike_values_grid))
    tuning_algo@nb_workers <- nb_workers

    check_data_and_model(data, model)
    model_processed <- prepare_model(data, model)
    data_processed <- prepare_data(data, model_processed)
    check_init(init, data_processed, model_processed)
    init_alg <- prepare_init(init, model_processed, data_processed)

    compile_model(model_processed@model_func, use_cpp)

    message("SAEMVS: step 1/2")
    progressr::with_progress({
      p <- progressr::progressor(along = spike_values_grid)
      map_runs <- safe_future_map(
        spike_values_grid,
        function(nu0) {
          backend <- compile_model(
            model_processed@model_func, use_cpp,
            silent = TRUE
          )
          fh <- saemvsHyperSpikeAndSlab(nu0, hyperparam)
          res <- saemvs_one_map_run(
            data_processed, model_processed, init_alg, tuning_algo, fh, backend
          )
          p()
          res
        },
        workers = tuning_algo@nb_workers,
        seed = tuning_algo@seed
      )
    })

    support <- lapply(map_runs, `[[`, "support")
    thresholds <- lapply(map_runs, `[[`, "threshold")
    beta_map <- lapply(map_runs, `[[`, "beta")

    support_hashes <- vapply(support, digest::digest, character(1))
    unique_support_indices <- which(!duplicated(support_hashes))
    hash_to_compact_index <- stats::setNames(
      seq_along(unique_support_indices),
      support_hashes[unique_support_indices]
    )
    support_index_mapping <- unname(vapply(
      support_hashes,
      function(h) hash_to_compact_index[[h]],
      integer(1)
    ))

    nb_workers <- min(tuning_algo@nb_workers, length(unique_support_indices))
    tuning_algo@nb_workers <- nb_workers

    message("SAEMVS: step 2/2")
    progressr::with_progress({
      criterion_results <- safe_future_map(
        unique_support_indices,
        function(k) {
          backend <- compile_model(
            model_processed@model_func, use_cpp,
            silent = TRUE
          )
          saemvs_one_ic_run(
            k, support, data_processed, model_processed, init_alg,
            tuning_algo, pen, backend
          )
        },
        workers = tuning_algo@nb_workers,
        seed = tuning_algo@seed
      )
    })

    criterion_values <- vapply(criterion_results, `[[`, numeric(1), "ll")
    mle <- lapply(criterion_results, `[[`, "mle_param")
    forced_variables_idx <- lapply(
      criterion_results, `[[`,
      "forced_variables_idx"
    )

    selected_variables_idx <- lapply(
      criterion_results, `[[`,
      "selected_variables_idx"
    )

    phi <- lapply(
      criterion_results, `[[`,
      "phi"
    )

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
      phi_fixed_idx = model_processed@phi_fixed_idx,
      phi_to_select_idx = model_processed@phi_to_select_idx,
      forced_variables_idx = forced_variables_idx,
      selected_variables_idx = selected_variables_idx,
      phi_names = model@phi_names,
      x_candidates_names = data@x_candidates_names,
      x_forced_names = data@x_forced_names,
      phi = phi
    )

    return(res) # nolint : return_linter
  }
)

#' Safe wrapper around `furrr::future_map` with controlled future plan
#'
#' This function provides a safe and reproducible way to map a function over
#' a list or vector in parallel using `furrr::future_map`. It restores the
#' previous `future` plan on exit and captures errors to provide a clean,
#' informative message. Optionally, a random seed can be set for
#' reproducibility.
#'
#' @param .x A vector or list to iterate over.
#' @param .f A function to apply to each element of `.x`.
#' @param ... Additional arguments passed to `.f`.
#' @param workers Integer. Number of parallel workers. Default is 1
#' (sequential).
#' @param seed Optional integer. Seed for reproducibility of parallel
#' operations.
#'
#' @return A list of results returned by applying `.f` over `.x`.
#'
#' @keywords internal
#' @noRd

safe_future_map <- function(.x, .f, ..., workers = 1, seed = NULL) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (workers > 1) {
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::sequential)
  }

  tryCatch(
    {
      furrr::future_map(
        .x,
        .f,
        ...,
        .options = furrr::furrr_options(seed = seed)
      )
    },
    error = function(e) {
      msg <- conditionMessage(e)
      msg <- sub("^Progress interrupted by [^:]+: ?", "", msg)

      stop(paste0("Parallelization interrupted: ", msg), call. = FALSE)
    }
  )
}


#' Internal: Single MAP Run of SAEMVS
#'
#' Performs one MAP estimation run of the SAEM algorithm for a given
#' spike-and-slab prior configuration.
#' This function is used internally by \code{\link{saemvs}} to explore
#' different spike variances and extract candidate supports.
#'
#' @param data An object of class
#' \link[=saemvsProcessedData-class]{saemvsProcessedData}
#' containing the processed dataset (observations and design matrices),
#' as returned by internal preprocessing steps.
#' @param model An object of class
#' \link[=saemvsProcessedModel-class]{saemvsProcessedModel}
#'   specifying the processed model structure and indices of parameters subject
#' @param init An object of class
#' \link[=saemvsProcessedInit-class]{saemvsProcessedInit}
#' providing processed initial values for the SAEM algorithm.
#' @param tuning_algo An object of class
#' \link[=saemvsTuning-class]{saemvsTuning}
#'   defining algorithmic parameters such as the number of iterations
#'   (\code{niter}).
#' @param hyperparam An object of class
#' \link[=saemvsHyperSpikeAndSlab-class]{saemvsHyperSpikeAndSlab}
#'   containing hyperparameters of the spike-and-slab prior (spike variance,
#'   slab variance, and mixing proportion).
#' @param backend A list containing the compiled model backend (typically
#'   generated by \code{\link{compile_model}}), used by \code{\link{run_saem}}
#'   to evaluate the model likelihood and its derivatives.

#'
#' @details
#' The function runs the SAEM algorithm via \code{\link{run_saem}}, then
#' constructs the MAP estimates of coefficients and derives the support set by
#' thresholding \eqn{\beta} coefficients.
#' Threshold values are computed from the spike-and-slab hyperparameters
#' (spike variance, slab variance, and mixing proportions) and are applied
#' uniformly across parameters to determine inclusion based on the magnitude
#' of the MAP estimates.
#' Forced covariates (specified in \code{x_forced_support}) are inserted
#' at the top of the support matrix (after the intercept) and are always
#' marked as selected, regardless of the thresholding step.
#'
#' The returned support matrix includes the intercept term as the first row.
#'
#' @return A \code{list} with elements:
#' \itemize{
#'   \item \code{threshold} Vector of threshold values applied to coefficients.
#'   \item \code{support} Logical matrix indicating which variables are
#'   selected in the final MAP solution. Rows correspond to variables
#'   subject to selection (with the intercept as the first row), and
#'   columns correspond to model parameters. Variables that are forced
#'   into the model are always marked as \code{TRUE}.
#'  \item \code{beta} Matrix of MAP coefficient estimates.
#' }
#'
#' @keywords internal
#' @seealso \code{\link{saemvs}}, \code{\link{run_saem}}
setGeneric(
  "saemvs_one_map_run",
  function(data, model, init, tuning_algo, hyperparam, backend) {
    standardGeneric("saemvs_one_map_run")
  }
)

setMethod(
  "saemvs_one_map_run",
  signature(
    data = "saemvsProcessedData",
    model = "saemvsProcessedModel",
    init = "saemvsProcessedInit",
    tuning_algo = "saemvsTuning",
    hyperparam = "saemvsHyperSpikeAndSlab",
    backend = "list"
  ),
  function(data, model, init, tuning_algo, hyperparam, backend) {
    map <- run_saem(data, model, init, tuning_algo, hyperparam, backend)
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


    beta_map <- map$beta_to_select[[niter + 1]][-1, , drop = FALSE]
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
    return(res) # nolint : return-linter
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
#' @param k Integer index indicating which element of the \code{support} list
#'   is evaluated (i.e., \code{support[[k]]}).
#' @param support A list of logical or numeric matrices, each representing
#'   a candidate support. Rows correspond to covariates subject to selection
#'   (with the intercept as the first row), and columns correspond to model
#'   parameters. Forced covariates may be included as leading rows and are
#'   handled separately.
#' @param data An object of class
#' \link[=saemvsProcessedData-class]{saemvsProcessedData},
#' containing the processed dataset restricted to candidate variables.
#' @param model An object of class
#'  \link[=saemvsProcessedModel-class]{saemvsProcessedModel},
#' specifying model structure and indices of parameters subject to selection.
#' @param init An object of class
#' \link[=saemvsProcessedInit-class]{saemvsProcessedInit},
#' providing processed initialization values for SAEM under the restricted
#' model.
#' @param tuning_algo An object of class
#' \link[=saemvsTuning-class]{saemvsTuning}, containing algorithmic tuning
#' parameters (e.g. number of iterations).
#' @param pen Character string, the information criterion to use. Must be either
#' \code{"BIC"} or \code{"e-BIC"}.
#' @param backend A list containing the compiled model backend (typically
#'   produced by \code{\link{compile_model}}), used for likelihood evaluation
#'   and SAEM-based MLE estimation.
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
#'   \item \code{mle_param} List of estimated parameters under the restricted
#'   model:
#'     \itemize{
#'       \item \code{beta} Estimated regression coefficients.
#'       \item \code{gamma} Estimated variances for the random effects.
#'       \item \code{sigma2} Estimated residual variance.
#'   \item \code{forced_variables_idx} Integer vector giving the row indices
#'   (within the support matrix) of covariates that were forced into the model.
#'   \item \code{selected_variables_idx} Integer vector giving the row indices
#'   of covariates that were actively selected among the candidate covariates
#'   (excluding forced ones).
#'     }
#' }
#'
#' @keywords internal
#' @seealso \code{\link{saemvs}}, \code{\link{saemvs_one_map_run}},
#' \code{\link{loglik}}
setGeneric(
  "saemvs_one_ic_run",
  function(k, support, data, model, init, tuning_algo, pen, backend) {
    standardGeneric("saemvs_one_ic_run")
  }
)

setMethod(
  "saemvs_one_ic_run",
  signature(
    k = "integer",
    support = "list",
    data = "saemvsProcessedData",
    model = "saemvsProcessedModel",
    init = "saemvsProcessedInit",
    tuning_algo = "saemvsTuning",
    pen = "character",
    backend = "list"
  ),
  function(k, support, data, model, init, tuning_algo, pen, backend) {
    p <- dim(data@x_candidates)[2]

    if (!is_empty_support(model@x_forced_support)) {
      nb_forced_beta <- sum(model@x_forced_support)
    } else {
      nb_forced_beta <- 0
    }


    forced_support <- extract_sub_support(
      model@x_forced_support,
      model@phi_to_select_idx
    )

    if (is_empty_support(forced_support)) {
      cand_support <- matrix(as.numeric(support[[k]]),
        nrow = nrow(support[[k]])
      )
      forced_rows <- integer(0)
      active_candidate_idx <- which(rowSums(cand_support[-1, ,
        drop = FALSE # nolint: indent_linter
      ]) > 0)
      selected_rows <- active_candidate_idx
    } else {
      cand_support <- matrix(as.numeric(support[[k]]),
        nrow = nrow(support[[k]])
      )
      idx_forced_phi_sel <- seq(1, dim(forced_support)[1])
      cand_support <- cand_support[-idx_forced_phi_sel, , drop = FALSE]
      forced_rows <- idx_forced_phi_sel
      active_candidate_idx <- which(rowSums(cand_support[-1, ,
        drop = FALSE # nolint: indent_linter
      ]) > 0)
      selected_rows <- active_candidate_idx
    }
    new_data <- map_to_mle_data(data, cand_support)
    new_model <- map_to_mle_model(model, cand_support)
    new_data_processed <- prepare_data(new_data, new_model)
    new_init <- map_to_mle_init(init, model, cand_support)
    new_processed_init <- prepare_init(new_init, new_model, new_data_processed)
    new_hyperparam <- saemvsHyperSlab(NULL, NULL, NULL)
    new_full_hyperparam <- saemvsHyperSpikeAndSlab(NULL, new_hyperparam)


    mle <- run_saem(
      new_data_processed, new_model, new_processed_init, tuning_algo,
      new_full_hyperparam, backend
    )

    mle_param <- list(
      beta   = mle$beta_not_to_select[[tuning_algo@niter + 1]],
      gamma  = mle$gamma_not_to_select[[tuning_algo@niter + 1]],
      sigma2 = mle$sigma2[tuning_algo@niter + 1]
    )

    ll <- loglik(
      new_data, new_model, tuning_algo, mle_param, pen, p,
      model@phi_to_select_idx, nb_forced_beta, backend
    )

    return(list( # nolint : return_linter
      ll = ll,
      mle_param = mle_param,
      forced_variables_idx = forced_rows,
      selected_variables_idx = selected_rows,
      phi = mle$phi
    ))
  }
)


#' Run a Single SAEMVS Test Fit
#'
#' Executes a simplified SAEMVS procedure for testing purposes.
#'
#' @param data A \link[=saemvsData-class]{saemvsData} object containing the
#' observed response series.
#' @param model A \link[=saemvsModel-class]{saemvsModel} object defining the
#' model structure.
#' @param init A \link[=saemvsInit-class]{saemvsInit} object providing initial
#' parameter values.
#' @param tuning_algo A \link[=saemvsTuning-class]{saemvsTuning} object
#' specifying algorithm hyperparameters.
#' @param hyperparam A \link[=saemvsHyperSlab-class]{saemvsHyperSlab} object
#' defining the slab prior. Internally combined with a spike value to form
#' a \code{saemvsHyperSpikeAndSlab} object.
#' @param use_cpp Logical. If TRUE (default), the model function is compiled
#'   to C++ for faster execution. If FALSE, a pure R backend is used.

#'
#' @return A \link[=saemResults-class]{saemResults} object containing estimated
#' parameters.
#'
#' @details
#' This function allows a quick fit of SAEMVS to inspect the evolution of
#' parameter estimates. The function performs the following steps:
#' \enumerate{
#'   \item Compile the R model function to C++ using `compile_model`.
#'   \item Construct a spike-and-slab hyperparameter object.
#'   \item Run the SAEM algorithm using `run_saem`.
#'   \item Return parameter estimates packaged in a `saemResults` object.
#' }
#'
#' @note
#' This function is intended for testing and debugging purposes only.
#' It does not perform variable selection across a grid of spike values
#' and should not be used for final model comparison.
#'
#' @examples
#' \dontrun{
#' test_results <- test_saemvs(
#'   data_obj, model_obj, init_obj, tuning_obj,
#'   hyper_obj
#' )
#' summary(test_results)
#' }
#'
#' @export
setGeneric(
  "test_saemvs",
  function(data, model, init, tuning_algo, hyperparam, use_cpp = TRUE) {
    standardGeneric("test_saemvs")
  }
)

#' @rdname test_saemvs
#' @exportMethod test_saemvs
setMethod(
  "test_saemvs",
  signature(
    data = "saemvsData",
    model = "saemvsModel",
    init = "saemvsInit",
    tuning_algo = "saemvsTuning",
    hyperparam = "saemvsHyperSlab"
  ),
  function(data, model, init, tuning_algo, hyperparam, use_cpp = TRUE) {
    check_data_and_model(data, model)
    model_processed <- prepare_model(data, model)
    backend <- compile_model(
      model_processed@model_func, use_cpp,
      silent = TRUE
    )

    full_hyperparam <- saemvsHyperSpikeAndSlab(
      tuning_algo@spike_values_grid[1],
      hyperparam
    )

    data_processed <- prepare_data(data, model_processed)
    check_init(init, data_processed, model_processed)
    init_alg <- prepare_init(init, model_processed, data_processed)

    saem_state <- run_saem(
      data_processed, model_processed, init_alg, tuning_algo, full_hyperparam,
      backend
    )

    res <- methods::new(
      "saemResults",
      beta_to_select = saem_state$beta_to_select,
      beta_not_to_select = saem_state$beta_not_to_select,
      gamma_to_select = saem_state$gamma_to_select,
      gamma_not_to_select = saem_state$gamma_not_to_select,
      sigma2 = saem_state$sigma2,
      phi_to_select_idx = model_processed@phi_to_select_idx,
      phi_not_to_select_idx = setdiff(
        seq_len(model_processed@phi_dim),
        model_processed@phi_to_select_idx
      ),
      phi_names = model@phi_names,
      x_candidates_names = data@x_candidates_names,
      x_forced_names = data@x_forced_names,
      phi = saem_state$phi
    )

    return(res) # nolint : return_linter
  }
)
