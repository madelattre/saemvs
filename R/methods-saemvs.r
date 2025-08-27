#' Sparse Variable Selection with SAEM
#'
#' The \code{saemvs} function performs stochastic approximation EM (SAEM) 
#' combined with variable selection using a spike-and-slab prior. 
#' It supports both BIC and extended BIC (e-BIC) as model selection criteria.
#'
#' This method explores a grid of spike variances (\eqn{\nu_0}) in parallel,
#' fits MAP estimates for each candidate support, and then evaluates unique 
#' supports with the chosen information criterion to select the best model.
#'
#' @param data An object of class \code{\link{saemvsData}} 
#'   containing the dataset (observations and design matrices).
#' @param model An object of class \code{\link{saemvsModel}} 
#'   specifying the model structure and the indices of parameters 
#'   that are candidates for selection (\code{phi_to_select_idx}).
#' @param init An object of class \code{\link{saemvsInit}} 
#'   providing the initialization for the SAEM algorithm.
#' @param tuning_algo An object of class \code{\link{saemvsTuning}} 
#'   containing algorithmic tuning parameters such as 
#'   the spike grid (\code{spike_values_grid}), number of workers, and random seed.
#' @param hyperparam An object of class \code{\link{saemvsHyperSlab}} 
#'   containing the hyperparameters for the slab prior.
#' @param pen Character string indicating the model selection criterion.
#'   Must be either \code{"BIC"} or \code{"e-BIC"}.
#'
#' @details
#' The procedure is carried out in two steps:
#' \enumerate{
#'   \item For each spike variance candidate (\eqn{\nu_0}), a MAP estimation
#'   is performed in parallel, producing supports, thresholds, and MAP estimates.
#'   \item Unique supports are identified, and each support is re-evaluated
#'   with the chosen information criterion (\code{pen}) to select the most
#'   appropriate model.
#' }
#'
#' Parallel execution is handled using the \pkg{future} and \pkg{furrr} packages.
#' The model function (C++ code) must be compiled on each worker before execution.
#'
#' @return An object of class \code{\link{saemvsResults}} containing:
#' \itemize{
#'   \item \code{criterion} The selection criterion used (\code{"BIC"} or \code{"e-BIC"}).
#'   \item \code{criterion_values} Values of the criterion for each unique support.
#'   \item \code{thresholds} Thresholds computed for each spike value.
#'   \item \code{beta_map} MAP estimates of regression parameters for each spike value.
#'   \item \code{mle_estimates} MLE estimates for each unique support.
#'   \item \code{support} List of supports obtained across spike values.
#'   \item \code{unique_support} List of unique supports identified.
#'   \item \code{support_mapping} Mapping from each run to the corresponding unique support.
#'   \item \code{spike_values_grid} The grid of spike variances used.
#'   \item \code{phi_fixed_idx} Indices of fixed parameters (not subject to selection).
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
#' @seealso \code{\link{saemvsResults}}, \code{\link{run_saem}}
#'
 
#' @export
setGeneric(
  "saemvs",
  function(data, model, init, tuning_algo, hyperparam, pen) {
    standardGeneric("saemvs")
  }
)

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
    hash_to_compact_index <- setNames(
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
          res <- saemvs_one_ebic_run(
            k, support, data, model, init, tuning_algo, hyperparam, pen
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
      function(x) criterion_results[[x]]$param
    )

    # --- Step 5: Assemble results object ---
    res <- new(
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
      phi_fixed_idx = model@phi_fixed_idx
    )

    return(res)
  }
)


#' @export
setGeneric(
  "test_saemvs",
  function(data, model, init, tuning_algo, hyperparam) {
    standardGeneric("test_saemvs")
  }
)

#' @exportMethod test_saemvs
setMethod(
  "test_saemvs",
  signature(
    data = "saemvsData", model = "saemvsModel", init = "saemvsInit",
    tuning_algo = "saemvsTuning", hyperparam = "saemvsHyperSlab"
  ),
  function(data, model, init, tuning_algo, hyperparam) {
    compile_model(model@model_func)
    full_hyperparam <- saemvsHyperSpikeAndSlab(
      tuning_algo@spike_values_grid[1],
      hyperparam
    )
    state <- run_saem(data, model, init, tuning_algo, full_hyperparam)

    res <- new(
      "saemResults",
      beta_to_select = state$beta_to_select,
      beta_not_to_select = state$beta_not_to_select,
      gamma_to_select = state$gamma_to_select,
      gamma_not_to_select = state$gamma_not_to_select,
      sigma2 = state$sigma2
    )
    return(res)
  }
)

## -- Sub-methods called in the main saemvs method

## ---- Computation of the map at each point on the grid

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
    supp_forced_phi_sel <- extract_sub_support(
      model@x_forced_support,
      model@phi_to_select_idx
    )

    q <- length(model@phi_to_select_idx)
    niter <- tuning_algo@niter
    if (is_empty_support(supp_forced_phi_sel) == TRUE) {
      p <- dim(data@x_candidates)[2]
    } else {
      p <- dim(data@x_candidates)[2] + dim(supp_forced_phi_sel)[1]
    }

    # et non -1, car les data sont transformées dans run_saem
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
      rep(TRUE, q),
      abs(beta_map) >= threshold_matrix
    ) # Attention, ici, le support contient l'intercept

    if (!is_empty_support(supp_forced_phi_sel)) {
      support[2:(1 + dim(supp_forced_phi_sel)[1]), ] <- TRUE
    }

    res <- list(
      threshold = threshold_matrix[1, ],
      support = support,
      beta = beta_map
    )
    return(res)
  }
)

## ---- Computation of the mle and of the log-likelihood for each support

setGeneric(
  "saemvs_one_ebic_run",
  function(k, support, data, model, init, tuning_algo, hyperparam, pen) {
    standardGeneric("saemvs_one_ebic_run")
  }
)

setMethod(
  "saemvs_one_ebic_run",
  signature(
    k = "numeric", support = "numeric",
    data = "saemvsData", model = "saemvsModel", init = "saemvsInit",
    tuning_algo = "saemvsTuning", hyperparam = "saemvsHyperSpikeAndSlab",
    pen = "character"
  ),
  saemvs_one_ebic_run <- function(k, support, data, model, init, tuning_algo,
                                  hyperparam, pen) {
    p <- dim(data@x_candidates)[2]

    supp_forced_phi_sel <- extract_sub_support(
      model@x_forced_support,
      model@phi_to_select_idx
    )

    if (is_empty_support(supp_forced_phi_sel)) {
      cand_support <- matrix(as.numeric(support[[k]]),
        nrow = nrow(support[[k]])
      )
      # contient l'intercept et ne concerne que les phi_sel de la première étape
    } else {
      cand_support <- matrix(as.numeric(support[[k]]),
        nrow = nrow(support[[k]])
      )
      idx_forced_phi_sel <- seq(1, dim(supp_forced_phi_sel)[1]) + 1
      cand_support <- cand_support[-idx_forced_phi_sel, ]
    }


    new_data <- map_to_mle_data(data, cand_support)
    new_model <- map_to_mle_model(model, cand_support)
    new_init <- map_to_mle_init(init, model, cand_support)
    new_hyperparam <- saemvsHyperSlab(NULL, NULL, NULL)
    new_full_hyperparam <- saemvsHyperSpikeAndSlab(NULL, new_hyperparam)

    mle <- run_saem(
      new_data, new_model, new_init, tuning_algo, new_full_hyperparam
    )

    param <- list(
      beta   = mle$beta_not_to_select[[tuning_algo@niter + 1]],
      gamma  = mle$gamma_not_to_select[[tuning_algo@niter + 1]],
      sigma2 = mle$sigma2[tuning_algo@niter + 1]
    )

    ll <- loglik(new_data, new_model, tuning_algo, param, pen, p)


    return(list(ll = ll, param = param))
  }
)
