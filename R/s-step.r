s_step <- function(config, k, state) {
  new_phi <- metropolis_vector_cpp(
    y = config$y_series,
    t = config$t_series,
    phi_current = split(state$phi[[k]], row(state$phi[[k]])),
    mean_prop = split(
      state$mprop_mh[[k]],
      row(state$mprop_mh[[k]])
    ),
    var_prop = state$vprop_mh[[k]],
    sigma2 = state$sigma2[k],
    niter_mh = config$num_mh_iterations,
    kappa = config$mh_proposal_scale,
    kernel = config$mh_kernel_type
  )

  state$phi[[k + 1]] <- matrix(unlist(new_phi), nrow = config$num_series, byrow = TRUE)

  return(state)
}


update_prop_mh_mix <- function(config, k, state) {
  state$mprop_mh[[k + 1]][, config$parameters_to_select_indices] <-
    config$x_phi_to_select %*% state$beta_hdim[[k + 1]]
  state$mprop_mh[[k + 1]][, -config$parameters_to_select_indices] <-
    config$x_phi_not_to_select %*% state$beta_ldim[[k + 1]]

  state$vprop_mh[[k + 1]][config$parameters_to_select_indices, config$parameters_to_select_indices] <-
    state$gamma_hdim[[k + 1]]
  state$vprop_mh[[k + 1]][-config$parameters_to_select_indices, -config$parameters_to_select_indices] <-
    state$gamma_ldim[[k + 1]]

  return(state)
}

update_prop_mh_hdim <- function(config, k, state) {
  state$mprop_mh[[k + 1]][, config$parameters_to_select_indices] <-
    config$x_phi_to_select %*% state$beta_hdim[[k + 1]]

  state$vprop_mh[[k + 1]][config$parameters_to_select_indices, config$parameters_to_select_indices] <- state$gamma_hdim[[k + 1]]

  return(state)
}

update_prop_mh_ldim <- function(config, k, state) {
  state$mprop_mh[[k + 1]][, config$parameters_not_to_select_indices] <-
    config$x_phi_not_to_select %*% state$beta_ldim[[k + 1]]

  state$vprop_mh[[k + 1]][config$parameters_not_to_select_indices, config$parameters_not_to_select_indices] <-
    state$gamma_ldim[[k + 1]]

  return(state)
}
