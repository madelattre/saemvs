#' Stochastic Approximation Step for parameters subject to selection
#' 
#' Update the stochastic approximation statistics for parameters that are subject to selection.
#' This function updates s1, s2_to_select, and s3_to_select in the MCMC state for a given iteration.
#' @param config List containing model configuration, including step sizes, series data, parameter indices, and the model function.
#' @param iteration Integer, the current iteration index (1-based).
#' @param state List representing the current MCMC state, including phi, s1, s2_to_select, and s3_to_select.
#' @return Updated \code{state} list with statistics updated for the parameters to select.
#' @keywords internal
sa_step_to_select <- function(config, iteration, state) {
  step_size <- config$step_size[iteration]
  phi_matrix <- state$phi[[iteration + 1]]

  squared_errors <- vapply(
    seq_along(config$y_series),
    function(i) {
      sum((config$y_series[[i]] - config$model_function(phi_matrix[i, ], config$t_series[[i]]))^2)
    },
    numeric(1)
  )
  total_error <- sum(squared_errors)

  state$s1[iteration + 1] <- state$s1[iteration] + step_size * (total_error - state$s1[iteration])
  state$s2_to_select[[iteration + 1]] <- state$s2_to_select[[iteration]] +
    step_size * (t(phi_matrix) %*% phi_matrix - state$s2_to_select[[iteration]])
  state$s3_to_select[[iteration + 1]] <- state$s3_to_select[[iteration]] +
    step_size * (phi_matrix - state$s3_to_select[[iteration]])

  return(state)
}

#' Stochastic Approximation Step for parameters not subject to selection
#' 
#' Update the stochastic approximation statistics for parameters that are not subject to selection.
#' This function updates s1, s2_not_to_select, and s3_not_to_select in the MCMC state for a given iteration.
#' @param config List containing model configuration, including step sizes, series data, parameter indices, and the model function.
#' @param iteration Integer, the current iteration index (1-based).
#' @param state List representing the current MCMC state, including phi, s1, s2_not_to_select, and s3_not_to_select.
#' @return Updated \code{state} list with statistics updated for the parameters not to select.
#' @keywords internal
sa_step_not_to_select <- function(config, iteration, state) {
  step_size <- config$step_size[iteration]
  phi_matrix <- state$phi[[iteration + 1]]

  squared_errors <- vapply(
    seq_along(config$y_series),
    function(i) {
      sum((config$y_series[[i]] - config$model_function(phi_matrix[i, ], config$t_series[[i]]))^2)
    },
    numeric(1)
  )
  total_error <- sum(squared_errors)

  state$s1[iteration + 1] <- state$s1[iteration] + step_size * (total_error - state$s1[iteration])
  state$s2_not_to_select[[iteration + 1]] <- state$s2_not_to_select[[iteration]] +
    step_size * (t(phi_matrix) %*% phi_matrix - state$s2_not_to_select[[iteration]])
  state$s3_not_to_select[[iteration + 1]] <- state$s3_not_to_select[[iteration]] +
    step_size * (phi_matrix - state$s3_not_to_select[[iteration]])

  return(state)
}

#' Stochastic Approximation Step for All Parameters
#' 
#' Update the stochastic approximation statistics for all parameters simultaneously.
#' This function updates s1, s2_to_select, s3_to_select, s2_not_to_select, and s3_not_to_select in the MCMC state for a given iteration.
#' @param config List containing model configuration, including step sizes, series data, parameter indices, and the model function.
#' @param iteration Integer, the current iteration index (1-based).
#' @param state List representing the current MCMC state, including phi, s1, s2_to_select, s3_to_select, s2_not_to_select, and s3_not_to_select.
#' @return Updated \code{state} list with statistics updated for all parameters.
#' @keywords internal
sa_step_all <- function(config, iteration, state) {
  step_size <- config$step_size[iteration]
  phi_matrix <- state$phi[[iteration + 1]]

  squared_errors <- vapply(
    seq_along(config$y_series),
    function(i) {
      sum((config$y_series[[i]] - g_vector_cpp(phi_matrix[i, ], config$t_series[[i]]))^2)
    },
    numeric(1)
  )
  total_error <- sum(squared_errors)

  phi_to_select <- as.matrix(phi_matrix[, config$parameters_to_select_indices])
  phi_not_to_select <- as.matrix(phi_matrix[, -config$parameters_to_select_indices])

  state$s1[iteration + 1] <- state$s1[iteration] + step_size * (total_error - state$s1[iteration])
  state$s2_to_select[[iteration + 1]] <- state$s2_to_select[[iteration]] +
    step_size * (t(phi_to_select) %*% phi_to_select - state$s2_to_select[[iteration]])
  state$s3_to_select[[iteration + 1]] <- state$s3_to_select[[iteration]] +
    step_size * (phi_to_select - state$s3_to_select[[iteration]])
  state$s2_not_to_select[[iteration + 1]] <- state$s2_not_to_select[[iteration]] +
    step_size * (t(phi_not_to_select) %*% phi_not_to_select - state$s2_not_to_select[[iteration]])
  state$s3_not_to_select[[iteration + 1]] <- state$s3_not_to_select[[iteration]] +
    step_size * (phi_not_to_select - state$s3_not_to_select[[iteration]])

  return(state)
}