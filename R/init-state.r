init_state <- function(config) {
  sigma2 <- rep(0, config$num_iterations + 1)

  phi <- lapply(
    seq_len(config$num_iterations + 1),
    function(k) matrix(0, config$num_series, config$total_parameters)
  )

  s1 <- rep(0, config$num_series + 1)


  # For the estimation of parameters subject to selection (hdim)

  if (config$num_parameters_to_select != 0) {
    beta_hdim <- lapply(
      seq_len(config$num_iterations + 1),
      function(k) matrix(0, config$num_covariates_to_select + 1, config$num_parameters_to_select)
    )

    gamma_hdim <- lapply(
      seq_len(config$num_iterations + 1),
      function(k) matrix(0, config$num_parameters_to_select, config$num_parameters_to_select)
    )

    alpha <- lapply(
      seq_len(config$num_iterations + 1),
      function(k) matrix(0, config$num_parameters_to_select, 1)
    )

    s2_hdim <- lapply(
      seq_len(config$num_iterations + 1),
      function(k) matrix(0, config$num_parameters_to_select, config$num_parameters_to_select)
    )
    s3_hdim <- lapply(
      seq_len(config$num_iterations + 1),
      function(k) matrix(0, config$num_series, config$num_parameters_to_select)
    )
  } else {
    beta_hdim <- NULL
    gamma_hdim <- NULL
    alpha <- NULL
    s2_hdim <- NULL
    s3_hdim <- NULL
  }


  # For the estimation of parameters non subject to selection (ldim)

  if (config$num_parameters_not_to_select != 0) {
    beta_ldim <- lapply(
      seq_len(config$num_iterations + 1),
      function(k) matrix(0, config$num_covariates_not_to_select + 1, config$num_parameters_not_to_select)
    )

    gamma_ldim <- lapply(
      seq_len(config$num_iterations + 1),
      function(k) matrix(0, config$num_parameters_not_to_select, config$num_parameters_not_to_select)
    )

    s2_ldim <- lapply(
      seq_len(config$num_iterations + 1),
      function(k) matrix(0, config$num_parameters_not_to_select, config$num_parameters_not_to_select)
    )
    s3_ldim <- lapply(
      seq_len(config$num_iterations + 1),
      function(k) matrix(0, config$num_series, config$num_parameters_not_to_select)
    )

    beta_mat_0 <- rep(0, config$num_parameters_not_to_select * (config$num_covariates_not_to_select + 1))
  } else {
    beta_ldim <- NULL
    gamma_ldim <- NULL
    s2_ldim <- lapply(
      seq_len(config$num_iterations + 1),
      function(k) NULL
    )
    s3_ldim <- lapply(
      seq_len(config$num_iterations + 1),
      function(k) NULL
    )
    beta_mat_0 <- NULL
  }




  sigma2[1] <- config$init_parameters@sigma2


  if (length(config$parameters_to_select_indices) > 0) {
    beta_hdim[[1]] <- config$init_parameters@beta_to_select
    gamma_hdim[[1]] <- config$init_parameters@gamma_to_select
    phi[[1]][, config$parameters_to_select_indices] <- config$x_phi_to_select %*% beta_hdim[[1]]
    alpha[[1]] <- config$init_parameters@inclusion_prob
  }

  if (length(config$parameters_not_to_select_indices) > 0) {
    beta_ldim[[1]] <- config$init_parameters@beta_not_to_select
    gamma_ldim[[1]] <- config$init_parameters@gamma_not_to_select
    phi[[1]][, config$parameters_not_to_select_indices] <- config$x_phi_not_to_select %*% beta_ldim[[1]]
  }

  # -- For the proposal distribution in the MH step
  mprop_mh <- phi
  vprop_mh <- lapply(
    seq_len(config$num_iterations + 1),
    function(k) matrix(0, config$total_parameters, config$total_parameters)
  )
  vprop_mh[[1]][config$parameters_to_select_indices, config$parameters_to_select_indices] <- gamma_hdim[[1]]
  if (length(config$parameters_not_to_select_indices) > 0) {
    vprop_mh[[1]][config$parameters_not_to_select_indices, config$parameters_not_to_select_indices] <-
      gamma_ldim[[1]]
  }

  list(
    beta_hdim = beta_hdim,
    gamma_hdim = gamma_hdim,
    beta_ldim = beta_ldim,
    gamma_ldim = gamma_ldim,
    beta_mat_0 = beta_mat_0,
    sigma2 = sigma2,
    alpha = alpha,
    phi = phi,
    s1 = s1,
    s2_hdim = s2_hdim,
    s3_hdim = s3_hdim,
    s2_ldim = s2_ldim,
    s3_ldim = s3_ldim,
    mprop_mh = mprop_mh,
    vprop_mh = vprop_mh
  )
}
