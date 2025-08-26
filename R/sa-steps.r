sa_step_ldim <- function(config, k, state) {
  step <- config$step_size[k]

  phi <- state$phi[[k + 1]]

  errs <- vapply(seq_along(config$y_series), function(i) {
    sum((config$y_series[[i]] - config$model_function(phi[i, ], config$t_series[[i]]))^2)
  }, numeric(1))

  mco <- sum(errs)

  state$s1[k + 1] <- state$s1[k] + step * (mco - state$s1[k])
  state$s2_not_to_select[[k + 1]] <- state$s2_not_to_select[[k]] +
    step * (t(phi) %*% phi - state$s2_not_to_select[[k]])
  state$s3_not_to_select[[k + 1]] <- state$s3_not_to_select[[k]] +
    step * (phi - state$s3_not_to_select[[k]])

  return(state)
}

sa_step_hdim <- function(config, k, state) {
  step <- config$step_size[k]

  phi <- state$phi[[k + 1]]

  errs <- vapply(seq_along(config$y_series), function(i) {
    sum((config$y_series[[i]] - config$model_function(phi[i, ], config$t_series[[i]]))^2)
  }, numeric(1))

  mco <- sum(errs)

  state$s1[k + 1] <- state$s1[k] + step * (mco - state$s1[k])
  state$s2_to_select[[k + 1]] <- state$s2_to_select[[k]] +
    step * (t(phi) %*% phi - state$s2_to_select[[k]])
  state$s3_to_select[[k + 1]] <- state$s3_to_select[[k]] +
    step * (phi - state$s3_to_select[[k]])

  return(state)
}

sa_step_split <- function(config, k, state) {
  step <- config$step_size[k]

  phi <- state$phi[[k + 1]]


  errs <- vapply(seq_along(config$y_series), function(i) {
    sum((config$y_series[[i]] - g_vector_cpp(phi[i, ], config$t_series[[i]]))^2)
  }, numeric(1))

  mco <- sum(errs)

  phih <- as.matrix(phi[, config$parameters_to_select_indices])
  phil <- as.matrix(phi[, -config$parameters_to_select_indices])

  state$s1[k + 1] <- state$s1[k] + step * (mco - state$s1[k])
  state$s2_to_select[[k + 1]] <- state$s2_to_select[[k]] +
    step * (t(phih) %*% phih - state$s2_to_select[[k]])
  state$s3_to_select[[k + 1]] <- state$s3_to_select[[k]] +
    step * (phih - state$s3_to_select[[k]])
  state$s2_not_to_select[[k + 1]] <- state$s2_not_to_select[[k]] +
    step * (t(phil) %*% phil - state$s2_not_to_select[[k]])
  state$s3_not_to_select[[k + 1]] <- state$s3_not_to_select[[k]] +
    step * (phil - state$s3_not_to_select[[k]])

  return(state)
}
