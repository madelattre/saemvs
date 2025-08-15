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
    data = "dataC", model = "modelC", init = "initC",
    tuning_algo = "tuningC", hyperparam = "hyperC",
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

    # Plan de parallélisation
    future::plan(future::multisession, workers = tuning_algo@nb_workers)

    nu0_list <- tuning_algo@nu0_grid

    # Lancer les calculs en parallèle
    map_results <- furrr::future_map(
      nu0_list,
      function(nu0) {
        # Important : le code C++ doit être dispo sur chaque worker
        compile_model(model@model_func)

        fh <- fullHyperC(nu0, hyperparam)
        saemvs_one_map_run(data, model, init, tuning_algo, fh)
      },
      .options = furrr::furrr_options(seed = tuning_algo@seed)
    )

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
    ebic_res <- furrr::future_map(
      unique_support,
      function(k) {
        # Important : le code C++ doit être dispo sur chaque worker
        compile_model(model@model_func)
        saemvs_one_ebic_run(
          k, support, data, model, init, tuning_algo, hyperparam, pen
        )
      },
      .options = furrr::furrr_options(seed = tuning_algo@seed)
    )

    ebic <- unlist(ebic_res)

    list(
      ebic = ebic,
      thresholds = thresholds,
      beta = beta_map,
      support = support,
      unique_support = unique_support,
      nu0_grid = tuning_algo@nu0_grid,
      index_select = model@index_select,
      map_to_unique_support = map_to_unique_support,
      pen = pen
    )

    res <- new(
      "resSAEMVS",
      pen = pen,
      crit_values = ebic,
      thresholds = thresholds,
      beta = beta_map,
      support = support,
      map_to_unique_support = map_to_unique_support,
      nu0_grid = tuning_algo@nu0_grid
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
    data = "dataC", model = "modelC", init = "initC",
    tuning_algo = "tuningC", hyperparam = "fullHyperC"
  ),
  function(data, model, init, tuning_algo, hyperparam) {
    map <- run_saem(data, model, init, tuning_algo, hyperparam)
    q <- length(model@index_select)
    niter <- tuning_algo@niter
    p <- dim(data@v)[2]
    # et non -1, car les data sont transformées dans run_saem
    beta_map <- map$beta_hdim[[niter + 1]][-1, ]
    alpha_map <- map$alpha[[niter + 1]]
    threshold_matrix <- matrix(
      rep(threshold(hyperparam@nu1, hyperparam@nu0, alpha_map), p),
      nrow = p, byrow = TRUE
    )
    support <- rbind(
      rep(TRUE, q),
      abs(beta_map) >= threshold_matrix
    ) # Attention, ici, le support contient l'intercept

    res <- list(
      threshold = threshold_matrix[1, ],
      support = support,
      beta = beta_map
    )
    # On garde beta pour le plot mais les autres paramètres sont-ils utiles?
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
    data = "dataC", model = "modelC", init = "initC",
    tuning_algo = "tuningC", hyperparam = "fullHyperC",
    pen = "character"
  ),
  saemvs_one_ebic_run <- function(k, support, data, model, init, tuning_algo, hyperparam, pen) {
    p <- dim(data@v)[2]

    nb_phi_s <- length(model@index_select)
    n <- length(data@y_list)

    covariate_support <- matrix(as.numeric(support[[k]]),
      nrow = nrow(support[[k]])
    )

    lines_with_ones <- which(rowSums(covariate_support) > 0)

    covariate_restricted_support <- covariate_support[lines_with_ones, ]
    # contient l'intercept

    v_restricted <- matrix(
      cbind(1, data@v)[, lines_with_ones],
      nrow = n
    )

    new_w <- matrix(
      cbind(data@w, v_restricted)[, -1],
      nrow = n
    )
    # -1 pour ne pas contenir l'intercept (cf v_restricted)


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
      model@covariate_support,
      selected_support_matrix,
      data@w,
      nb_phi_s,
      model@q_phi - nb_phi_s
    )

    new_model <- modelC(
      g = model@model_func,
      nphi = model@q_phi,
      phi_fixed = model@index_fixed,
      support = new_covariate_support
    )

    new_data <- dataC(y = data@y_list, t = data@t_list, w = new_w)

    new_beta_init <- merge_init_beta(
      init@beta_hdim, init@beta_ldim, lines_with_ones, data@w
    )

    new_gamma_init <- as.matrix(Matrix::bdiag(init@gamma_hdim, init@gamma_ldim))

    new_init <- initC(
      beta_ns = new_beta_init,
      gamma_ns = new_gamma_init,
      sigma2 = init@sigma2
    )

    new_hyperparam <- hyperC(NULL, NULL, NULL)
    new_full_hyperparam <- fullHyperC(NULL, new_hyperparam)

    mle <- run_saem(
      new_data, new_model, new_init, tuning_algo, new_full_hyperparam
    )

    param <- list(
      beta   = mle$beta_ldim[[tuning_algo@niter + 1]],
      gamma  = mle$gamma_ldim[[tuning_algo@niter + 1]],
      sigma2 = mle$sigma2[tuning_algo@niter + 1]
    )

    ll <- loglik(new_data, new_model, tuning_algo, param, pen, p)

    return(ll)
  }
)
