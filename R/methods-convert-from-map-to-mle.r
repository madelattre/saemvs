
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
