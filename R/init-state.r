init_state <- function(config) {
  sigma2 <- rep(0, config$niter + 1)

  phi <- lapply(
    seq_len(config$niter + 1),
    function(k) matrix(0, config$n, config$q_phi)
  )

  s1 <- rep(0, config$niter + 1)


  # For the estimation of parameters subject to selection (hdim)

  if (config$q_hdim != 0) {
    beta_hdim <- lapply(
      seq_len(config$niter + 1),
      function(k) matrix(0, config$pv + 1, config$q_hdim)
    )

    gamma_hdim <- lapply(
      seq_len(config$niter + 1),
      function(k) matrix(0, config$q_hdim, config$q_hdim)
    )

    alpha <- lapply(
      seq_len(config$niter + 1),
      function(k) matrix(0, config$q_hdim, 1)
    )

    s2_hdim <- lapply(
      seq_len(config$niter + 1),
      function(k) matrix(0, config$q_hdim, config$q_hdim)
    )
    s3_hdim <- lapply(
      seq_len(config$niter + 1),
      function(k) matrix(0, config$n, config$q_hdim)
    )
  } else {
    beta_hdim <- NULL
    gamma_hdim <- NULL
    alpha <- NULL
    s2_hdim <- NULL
    s3_hdim <- NULL
  }


  # For the estimation of parameters non subject to selection (ldim)

  if (config$q_ldim != 0) {
    beta_ldim <- lapply(
      seq_len(config$niter + 1),
      function(k) matrix(0, config$pw + 1, config$q_ldim)
    )

    gamma_ldim <- lapply(
      seq_len(config$niter + 1),
      function(k) matrix(0, config$q_ldim, config$q_ldim)
    )

    s2_ldim <- lapply(
      seq_len(config$niter + 1),
      function(k) matrix(0, config$q_ldim, config$q_ldim)
    )
    s3_ldim <- lapply(
      seq_len(config$niter + 1),
      function(k) matrix(0, config$n, config$q_ldim)
    )

    beta_mat_0 <- rep(0, config$q_ldim * (config$pw + 1))
  } else {
    beta_ldim <- NULL
    gamma_ldim <- NULL
    s2_ldim <- lapply(
      seq_len(config$niter + 1),
      function(k) NULL
    )
    s3_ldim <- lapply(
      seq_len(config$niter + 1),
      function(k) NULL
    )
    beta_mat_0 <- NULL
  }




  sigma2[1] <- config$param_init@sigma2


  if (length(config$index_select) > 0) {
    beta_hdim[[1]] <- config$param_init@beta_to_select
    gamma_hdim[[1]] <- config$param_init@gamma_to_select
    phi[[1]][, config$index_select] <- config$v %*% beta_hdim[[1]]
    alpha[[1]] <- config$param_init@inclusion_prob
  }

  if (length(config$index_unselect) > 0) {
    beta_ldim[[1]] <- config$param_init@beta_not_to_select
    gamma_ldim[[1]] <- config$param_init@gamma_not_to_select
    phi[[1]][, config$index_unselect] <- config$w %*% beta_ldim[[1]]
  }

  # -- For the proposal distribution in the MH step
  mprop_mh <- phi
  vprop_mh <- lapply(
    seq_len(config$niter + 1),
    function(k) matrix(0, config$q_phi, config$q_phi)
  )
  vprop_mh[[1]][config$index_select, config$index_select] <- gamma_hdim[[1]]
  if (length(config$index_unselect) > 0) {
    vprop_mh[[1]][config$index_unselect, config$index_unselect] <-
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
