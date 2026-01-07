#' Stochastic Approximation Step for parameters subject to selection
#'
#' Update the stochastic approximation statistics for parameters that are
#' subject to selection.
#' This function updates s1, s2_to_select, and s3_to_select in the MCMC state
#' for a given iteration.
#' @param config List containing the model configuration, including:
#'   \itemize{
#'     \item \code{step_size}: vector of SA step sizes per iteration,
#'     \item \code{y_series}: list of observed data series,
#'     \item \code{t_series}: list of corresponding time points per series,
#'     \item \code{model_function}: function mapping a latent parameter vector
#'           \code{phi_i} and time points \code{t_i} to predicted observations.
#'   }
#' @param iteration Integer, the current iteration index (1-based).
#' @param state List representing the current MCMC state, including:
#'   \itemize{
#'     \item \code{phi}: list of latent parameter matrices per iteration,
#'     \item \code{s1}: vector of scalar statistics per iteration,
#'     \item \code{s2_to_select}: list of covariance-like matrices for
#' parameters subject to selection,
#'     \item \code{s3_to_select}: list of mean-like vectors for parameters
#'           subject to selection.
#'   }
#'   These statistics are updated in place according to the stochastic
#' approximation step.
#' @param backend A list containing compiled model functions, passed for
#' consistency with the SAEM iteration framework.

#' @return Updated \code{state} list with:
#'   \itemize{
#'     \item \code{s1[iteration + 1]} updated scalar statistic,
#'     \item \code{s2_to_select[[iteration + 1]]} updated covariance statistic,
#'     \item \code{s3_to_select[[iteration + 1]]} updated mean statistic.
#'   }

#' @keywords internal
#' @note This function is intended for internal use within the SAEM-MCMC
#' algorithm and should not be called directly in user code.
sa_step_to_select <- function(config, iteration, state, backend) {
  step_size <- config$step_size[iteration]
  phi_matrix <- state$phi[[iteration + 1]]

  squared_errors <- vapply(
    seq_along(config$y_series),
    function(i) {
      sum((config$y_series[[i]] -
             config$model_function(phi_matrix[i, ], config$t_series[[i]]))^2)
    },
    numeric(1)
  )
  total_error <- sum(squared_errors)

  state$s1[iteration + 1] <- state$s1[iteration] +
    step_size * (total_error - state$s1[iteration])
  state$s2_to_select[[iteration + 1]] <- state$s2_to_select[[iteration]] +
    step_size * (t(phi_matrix) %*% phi_matrix - state$s2_to_select[[iteration]])
  state$s3_to_select[[iteration + 1]] <- state$s3_to_select[[iteration]] +
    step_size * (phi_matrix - state$s3_to_select[[iteration]])

  return(state) # nolint: return-linter
}

#' Stochastic Approximation Step for parameters not subject to selection
#'
#' Update the stochastic approximation statistics for parameters that are not
#' subject to selection.
#' This function updates s1, s2_not_to_select, and s3_not_to_select in the MCMC
#' state for a given iteration.
#' @param config List containing the model configuration, including:
#'   \itemize{
#'     \item \code{step_size}: vector of SA step sizes per iteration,
#'     \item \code{y_series}: list of observed data series,
#'     \item \code{t_series}: list of corresponding time points per series,
#'     \item \code{model_function}: function mapping a latent parameter vector
#'           \code{phi_i} and time points \code{t_i} to predicted observations.
#'   }
#' @param iteration Integer, the current iteration index (1-based).
#' @param state List representing the current MCMC state, including:
#'   \itemize{
#'     \item \code{phi}: list of latent parameter matrices per iteration,
#'     \item \code{s1}: vector of scalar statistics per iteration,
#'     \item \code{s2_not_to_select}: list of covariance-like matrices for
#' parameters not subject to selection,
#'     \item \code{s3_not_to_select}: list of mean-like vectors for parameters
#'           not subject to selection.
#'   }
#'   These statistics are updated in place according to the stochastic
#'  approximation step.
#' @param backend A list containing compiled model functions, passed for
#' consistency with the SAEM iteration framework.
#' @return Updated \code{state} list with:
#'   \itemize{
#'     \item \code{s1[iteration + 1]} updated scalar statistic,
#'     \item \code{s2_not_to_select[[iteration + 1]]} updated covariance
#' statistic,
#'     \item \code{s3_not_to_select[[iteration + 1]]} updated mean statistic.
#'   }
#' @keywords internal
#' @note This function is intended for internal use within the SAEM-MCMC
#' algorithm and should not be called directly in user code.
sa_step_not_to_select <- function(config, iteration, state, backend) {
  step_size <- config$step_size[iteration]
  phi_matrix <- state$phi[[iteration + 1]]

  squared_errors <- vapply(
    seq_along(config$y_series),
    function(i) {
      sum((config$y_series[[i]] -
             config$model_function(phi_matrix[i, ], config$t_series[[i]]))^2)
    },
    numeric(1)
  )
  total_error <- sum(squared_errors)

  state$s1[iteration + 1] <- state$s1[iteration] +
    step_size * (total_error - state$s1[iteration])
  state$s2_not_to_select[[iteration + 1]] <-
    state$s2_not_to_select[[iteration]] +
    step_size * (t(phi_matrix) %*% phi_matrix -
                 state$s2_not_to_select[[iteration]])
  state$s3_not_to_select[[iteration + 1]] <-
    state$s3_not_to_select[[iteration]] +
    step_size * (phi_matrix - state$s3_not_to_select[[iteration]])

  return(state) # nolint: return-linter
}

#' Stochastic Approximation Step for All Parameters
#'
#' Update the stochastic approximation statistics for all parameters
#' simultaneously.
#' This function updates s1, s2_to_select, s3_to_select, s2_not_to_select, and
#' s3_not_to_select in the MCMC state for a given iteration.
#' @param config List containing the model configuration, including:
#'   \itemize{
#'     \item \code{step_size}: vector of SA step sizes per iteration,
#'     \item \code{y_series}: list of observed data series,
#'     \item \code{t_series}: list of corresponding time points per series,
#'     \item \code{parameters_to_select_indices}: indices of parameters
#' subject to selection.
#'   }
#' @param iteration Integer, the current iteration index (1-based).
#' @param state List representing the current MCMC state, including:
#'   \itemize{
#'     \item \code{phi}: list of latent parameter matrices per iteration,
#'     \item \code{s1}: vector of scalar statistics per iteration,
#'     \item \code{s2_to_select}, \code{s3_to_select}: lists of covariance-like
#'           and mean-like statistics for parameters subject to selection,
#'     \item \code{s2_not_to_select}, \code{s3_not_to_select}: lists of
#' covariance-like and mean-like statistics for parameters not subject to
#'  selection.
#'   }
#'   All these statistics are updated according to the stochastic approximation
#'  step.
#' @param backend A list containing compiled model functions,
#'   including \code{g_vector}, which computes predicted observations given
#'   a latent parameter vector. Used to compute squared errors for SA updates.
#' @return Updated \code{state} list with:
#'   \itemize{
#'     \item \code{s1[iteration + 1]} updated scalar statistic,
#'     \item \code{s2_to_select[[iteration + 1]]},
#'  \code{s3_to_select[[iteration + 1]]} updated statistics for
#' parameters to select,
#'     \item \code{s2_not_to_select[[iteration + 1]]},
#'  \code{s3_not_to_select[[iteration + 1]]} updated statistics
#'           for parameters not to select.
#'   }
#' @keywords internal
#' @note This function is intended for internal use within the SAEM-MCMC
#' algorithm and should not be called directly in user code.
sa_step_all <- function(config, iteration, state, backend) {
  step_size <- config$step_size[iteration]
  phi_matrix <- state$phi[[iteration + 1]]

  squared_errors <- vapply(
    seq_along(config$y_series),
    function(i) {
      sum((config$y_series[[i]] -
             backend$g_vector(phi_matrix[i, ], config$t_series[[i]]))^2)
    },
    numeric(1)
  )
  total_error <- sum(squared_errors)

  phi_to_select <- as.matrix(phi_matrix[
    ,
    config$parameters_to_select_indices
  ])
  phi_not_to_select <-
    as.matrix(phi_matrix[, -config$parameters_to_select_indices])

  state$s1[iteration + 1] <- state$s1[iteration] +
    step_size * (total_error - state$s1[iteration])
  state$s2_to_select[[iteration + 1]] <- state$s2_to_select[[iteration]] +
    step_size * (t(phi_to_select) %*% phi_to_select -
                   state$s2_to_select[[iteration]])
  state$s3_to_select[[iteration + 1]] <- state$s3_to_select[[iteration]] +
    step_size * (phi_to_select - state$s3_to_select[[iteration]])
  state$s2_not_to_select[[iteration + 1]] <-
    state$s2_not_to_select[[iteration]] +
    step_size * (t(phi_not_to_select) %*% phi_not_to_select -
                 state$s2_not_to_select[[iteration]])
  state$s3_not_to_select[[iteration + 1]] <-
    state$s3_not_to_select[[iteration]] +
    step_size * (phi_not_to_select - state$s3_not_to_select[[iteration]])

  return(state) # nolint: return-linter
}
