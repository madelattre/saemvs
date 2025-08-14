s_step <- function(config, k, state) {
  new_phi <- metropolis_vector_cpp(
    y = config$yi,
    t = config$ti,
    phi_current = split(state$phi[[k]], row(state$phi[[k]])),
    mean_prop = split(
      state$mprop_mh[[k]],
      row(state$mprop_mh[[k]])
    ),
    var_prop = state$vprop_mh[[k]],
    sigma2 = state$sigma2[k],
    niter_mh = config$niter_mh,
    kappa = config$kappa,
    kernel = config$kernel_mh
  )

  state$phi[[k + 1]] <- matrix(unlist(new_phi), nrow = config$n, byrow = TRUE)

  return(state)
}


update_prop_mh_mix <- function(config, k, state) {
  state$mprop_mh[[k + 1]][, config$index_select] <-
    config$v %*% state$beta_hdim[[k + 1]]
  state$mprop_mh[[k + 1]][, -config$index_select] <-
    config$w %*% state$beta_ldim[[k + 1]]

  state$vprop_mh[[k + 1]][config$index_select, config$index_select] <-
    state$gamma_hdim[[k + 1]]
  state$vprop_mh[[k + 1]][-config$index_select, -config$index_select] <-
    state$gamma_ldim[[k + 1]]

  return(state)
}

update_prop_mh_hdim <- function(config, k, state) {
  state$mprop_mh[[k + 1]][, config$index_select] <-
    config$v %*% state$beta_hdim[[k + 1]]

  state$vprop_mh[[k + 1]][config$index_select, config$index_select] <-
    state$gamma_hdim[[k + 1]]

  return(state)
}

update_prop_mh_ldim <- function(config, k, state) {
  state$mprop_mh[[k + 1]][, config$index_unselect] <-
    config$w %*% state$beta_ldim[[k + 1]]

  state$vprop_mh[[k + 1]][config$index_unselect, config$index_unselect] <-
    state$gamma_ldim[[k + 1]]

  return(state)
}
