sa_step_ldim <- function(config, k, state) {
  step <- config$step[k]

  phi <- state$phi[[k + 1]]

  errs <- vapply(seq_along(config$yi), function(i) {
    sum((config$yi[[i]] - config$g(phi[i, ], config$ti[[i]]))^2)
  }, numeric(1))

  mco <- sum(errs)

  state$s1[k + 1] <- state$s1[k] + step * (mco - state$s1[k])
  state$s2_ldim[[k + 1]] <- state$s2_ldim[[k]] +
    step * (t(phi) %*% phi - state$s2_ldim[[k]])
  state$s3_ldim[[k + 1]] <- state$s3_ldim[[k]] +
    step * (phi - state$s3_ldim[[k]])

  return(state)
}

sa_step_hdim <- function(config, k, state) {
  step <- config$step[k]

  phi <- state$phi[[k + 1]]

  errs <- vapply(seq_along(config$yi), function(i) {
    sum((config$yi[[i]] - config$g(phi[i, ], config$ti[[i]]))^2)
  }, numeric(1))

  mco <- sum(errs)

  state$s1[k + 1] <- state$s1[k] + step * (mco - state$s1[k])
  state$s2_hdim[[k + 1]] <- state$s2_hdim[[k]] +
    step * (t(phi) %*% phi - state$s2_hdim[[k]])
  state$s3_hdim[[k + 1]] <- state$s3_hdim[[k]] +
    step * (phi - state$s3_hdim[[k]])

  return(state)
}

sa_step_split <- function(config, k, state) {
  step <- config$step[k]

  phi <- state$phi[[k + 1]]


  # errs <- vapply(seq_along(config$yi), function(i) {
  #   sum((config$yi[[i]] - config$g(phi[i, ], config$ti[[i]]))^2)
  # }, numeric(1))


  errs <- vapply(seq_along(config$yi), function(i) {
    sum((config$yi[[i]] - g_vector_cpp(phi[i, ], config$ti[[i]]))^2)
  }, numeric(1))

  mco <- sum(errs)

  phih <- as.matrix(phi[, config$index_select])
  phil <- as.matrix(phi[, -config$index_select])

  state$s1[k + 1] <- state$s1[k] + step * (mco - state$s1[k])
  state$s2_hdim[[k + 1]] <- state$s2_hdim[[k]] +
    step * (t(phih) %*% phih - state$s2_hdim[[k]])
  state$s3_hdim[[k + 1]] <- state$s3_hdim[[k]] +
    step * (phih - state$s3_hdim[[k]])
  state$s2_ldim[[k + 1]] <- state$s2_ldim[[k]] +
    step * (t(phil) %*% phil - state$s2_ldim[[k]])
  state$s3_ldim[[k + 1]] <- state$s3_ldim[[k]] +
    step * (phil - state$s3_ldim[[k]])

  return(state)
}
