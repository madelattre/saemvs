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

    # Plan de parallélisation
    future::plan(future::multisession, workers = tuning_algo@nb_workers)

    nu0_list <- tuning_algo@spike_values_grid

    # Lancer les calculs en parallèle

    message("SAEMVS: step 1/2")
    progressr::with_progress({
      p <- progressr::progressor(along = nu0_list)

      map_results <- furrr::future_map(
        nu0_list,
        function(nu0) {
          # Important : le code C++ doit être dispo sur chaque worker
          compile_model(model@model_func)
          fh <- saemvsHyperSpikeAndSlab(nu0, hyperparam)
          res <- saemvs_one_map_run(data, model, init, tuning_algo, fh)
          p()
          res
        },
        .options = furrr::furrr_options(seed = tuning_algo@seed)
      )
    })


    support <- lapply(
      seq_along(map_results),
      function(x) map_results[[x]]$support
    )

    thresholds <- lapply(
      seq_along(map_results),
      function(x) map_results[[x]]$threshold
    )

    beta_map <- lapply(
      seq_along(map_results),
      function(x) map_results[[x]]$beta
    )

    # recherche des supports uniques
    support_hashes <- sapply(support, digest::digest)
    unique_support <- which(!duplicated(support_hashes) == TRUE)
    hash_to_compact_index <- setNames(
      seq_along(unique_support),
      support_hashes[unique_support]
    )


    map_to_unique_support <- unname(sapply(
      support_hashes,
      function(h) hash_to_compact_index[[h]]
    ))


    # Lancer les calculs en parallèle
    message("SAEMVS: step 2/2")
    progressr::with_progress({
      p <- progressr::progressor(along = nu0_list)
      ebic_res <- furrr::future_map(
        unique_support,
        function(k) {
          # Important : le code C++ doit être dispo sur chaque worker
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

    ebic <- unlist(lapply(
      seq_along(ebic_res),
      function(x) ebic_res[[x]]$ll
    ))

    est_mle <- lapply(
      seq_along(ebic_res),
      function(x) ebic_res[[x]]$param
    )

    # ebic <- 0

    # est_mle <- list()

    res <- new(
      "saemvsResults",
      criterion = pen,
      criterion_values = ebic,
      thresholds = thresholds,
      beta_map = beta_map,
      mle_estimates = est_mle,
      support = support,
      unique_support = support[unique_support],
      support_mapping = map_to_unique_support,
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
      "resSAEM",
      beta_s = state$beta_hdim,
      beta_ns = state$beta_ldim,
      gamma_s = state$gamma_hdim,
      gamma_ns = state$gamma_ldim,
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
    beta_map <- map$beta_hdim[[niter + 1]][-1, ]
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


    new_data <- from_map_to_mle_data(data, cand_support)
    new_model <- from_map_to_mle_model(model, cand_support)
    new_init <- from_map_to_mle_init(init, model, cand_support)
    new_hyperparam <- saemvsHyperSlab(NULL, NULL, NULL)
    new_full_hyperparam <- saemvsHyperSpikeAndSlab(NULL, new_hyperparam)

    mle <- run_saem(
      new_data, new_model, new_init, tuning_algo, new_full_hyperparam
    )

    param <- list(
      beta   = mle$beta_ldim[[tuning_algo@niter + 1]],
      gamma  = mle$gamma_ldim[[tuning_algo@niter + 1]],
      sigma2 = mle$sigma2[tuning_algo@niter + 1]
    )

    ll <- loglik(new_data, new_model, tuning_algo, param, pen, p)


    return(list(ll = ll, param = param))
  }
)
