#' Initialize algorithm state for SAEMVS
#'
#' Internal function to initialize all necessary matrices and lists for the
#' SAEMVS algorithm according to the provided configuration.  
#' Supports parameters subject to selection (`_to_select`) and
#' parameters not subject to selection (`_not_to_select`).  
#' Includes initial values for `beta`, `gamma`, `phi`, `sigma2`, inclusion probabilities (`alpha`),
#' sufficient statistics, and proposal distributions for the Metropolis-Hastings step.
#'
#' @param config A `saemvsConfig` object containing:
#'   \describe{
#'     \item{num_iterations}{Total number of iterations.}
#'     \item{num_series}{Number of series/individuals.}
#'     \item{total_parameters}{Total number of parameters in the model.}
#'     \item{num_parameters_to_select}{Number of parameters subject to selection.}
#'     \item{num_covariates_to_select}{Number of covariates for parameters subject to selection.}
#'     \item{num_parameters_not_to_select}{Number of parameters not subject to selection.}
#'     \item{num_covariates_not_to_select}{Number of covariates for parameters not subject to selection.}
#'     \item{init_parameters}{Initial values as a `saemvsProcessedInit` object.}
#'     \item{parameters_to_select_indices}{Indices of parameters subject to select in the full vector.}
#'     \item{parameters_not_to_select_indices}{Indices of parameters not subject to selection.}
#'     \item{x_phi_to_select}{Design matrix for parameters subject to selection.}
#'     \item{x_phi_not_to_select}{Design matrix for parameters not subject to selection.}
#'   }
#'
#' @return A named list containing:
#'   \describe{
#'     \item{beta_to_select}{List of beta matrices for parameters subject to selection.}
#'     \item{gamma_to_select}{List of gamma covariance matrices for parameters subject to selection.}
#'     \item{beta_not_to_select}{List of beta matrices for parameters not subject to select.}
#'     \item{gamma_not_to_select}{List of gamma covariance matrices for parameters not subject to selection.}
#'     \item{beta_mat_0}{Initial vectorized beta for active covariates for parameters not subject to selection.}
#'     \item{sigma2}{Vector of residual variance across iterations.}
#'     \item{alpha}{List of inclusion probabilities.}
#'     \item{phi}{List of parameter values for all series.}
#'     \item{s1}{Vector of sufficient statistics for residuals.}
#'     \item{s2_to_select, s3_to_select}{Sufficient statistics for parameters subject to selection.}
#'     \item{s2_not_to_select, s3_not_to_select}{Sufficient statistics for parameters not subject to selection.}
#'     \item{mprop_mh}{List of mean proposal matrices for MH step.}
#'     \item{vprop_mh}{List of variance proposal matrices for MH step.}
#'   }
#'
#' @details
#' Performs consistent initialization for all elements required in the SAEMVS
#' algorithm. Ensures that all lists have length `num_iterations + 1` and
#' that initial values are properly assigned. Handles the case when either
#' parameters subject to selection or parameters not subject to selection are absent.
#'
#' @keywords internal
init_state <- function(config) {
  n_iter <- config$num_iterations + 1
  n_series <- config$num_series
  n_total_params <- config$total_parameters

  # Initialize residual variance and phi values
  sigma2 <- rep(0, n_iter)
  phi <- lapply(seq_len(n_iter), function(k) matrix(0, n_series, n_total_params))
  s1 <- rep(0, n_series + 1)

  # --- Parameters to select
  if (config$num_parameters_to_select > 0) {
    beta_to_select <- lapply(seq_len(n_iter),
                             function(k) matrix(0, config$num_covariates_to_select + 1, config$num_parameters_to_select))
    gamma_to_select <- lapply(seq_len(n_iter),
                              function(k) matrix(0, config$num_parameters_to_select, config$num_parameters_to_select))
    alpha <- lapply(seq_len(n_iter),
                    function(k) matrix(0, config$num_parameters_to_select, 1))
    s2_to_select <- lapply(seq_len(n_iter),
                           function(k) matrix(0, config$num_parameters_to_select, config$num_parameters_to_select))
    s3_to_select <- lapply(seq_len(n_iter),
                           function(k) matrix(0, n_series, config$num_parameters_to_select))
  } else {
    beta_to_select <- gamma_to_select <- alpha <- s2_to_select <- s3_to_select <- NULL
  }

  # --- Parameters not to select
  if (config$num_parameters_not_to_select > 0) {
    beta_not_to_select <- lapply(seq_len(n_iter),
                                 function(k) matrix(0, config$num_covariates_not_to_select + 1, config$num_parameters_not_to_select))
    gamma_not_to_select <- lapply(seq_len(n_iter),
                                  function(k) matrix(0, config$num_parameters_not_to_select, config$num_parameters_not_to_select))
    s2_not_to_select <- lapply(seq_len(n_iter),
                               function(k) matrix(0, config$num_parameters_not_to_select, config$num_parameters_not_to_select))
    s3_not_to_select <- lapply(seq_len(n_iter),
                               function(k) matrix(0, n_series, config$num_parameters_not_to_select))
    beta_mat_0 <- rep(0, config$num_parameters_not_to_select * (config$num_covariates_not_to_select + 1))
  } else {
    beta_not_to_select <- gamma_not_to_select <- NULL
    s2_not_to_select <- s3_not_to_select <- lapply(seq_len(n_iter), function(k) NULL)
    beta_mat_0 <- NULL
  }

  # --- Set initial values
  sigma2[1] <- config$init_parameters@sigma2

  if (length(config$parameters_to_select_indices) > 0) {
    beta_to_select[[1]] <- config$init_parameters@beta_to_select
    gamma_to_select[[1]] <- config$init_parameters@gamma_to_select
    phi[[1]][, config$parameters_to_select_indices] <- config$x_phi_to_select %*% beta_to_select[[1]]
    alpha[[1]] <- config$init_parameters@inclusion_prob
  }

  if (length(config$parameters_not_to_select_indices) > 0) {
    beta_not_to_select[[1]] <- config$init_parameters@beta_not_to_select
    gamma_not_to_select[[1]] <- config$init_parameters@gamma_not_to_select
    phi[[1]][, config$parameters_not_to_select_indices] <- config$x_phi_not_to_select %*% beta_not_to_select[[1]]
  }

  # --- Initialize proposal distributions for Metropolis-Hastings
  mprop_mh <- phi
  vprop_mh <- lapply(seq_len(n_iter), function(k) matrix(0, n_total_params, n_total_params))
  if (!is.null(beta_to_select)) {
    vprop_mh[[1]][config$parameters_to_select_indices, config$parameters_to_select_indices] <- gamma_to_select[[1]]
  }
  if (!is.null(beta_not_to_select)) {
    vprop_mh[[1]][config$parameters_not_to_select_indices, config$parameters_not_to_select_indices] <- gamma_not_to_select[[1]]
  }

  # --- Return the initialized state
  list(
    beta_to_select = beta_to_select,
    gamma_to_select = gamma_to_select,
    beta_not_to_select = beta_not_to_select,
    gamma_not_to_select = gamma_not_to_select,
    beta_mat_0 = beta_mat_0,
    sigma2 = sigma2,
    alpha = alpha,
    phi = phi,
    s1 = s1,
    s2_to_select = s2_to_select,
    s3_to_select = s3_to_select,
    s2_not_to_select = s2_not_to_select,
    s3_not_to_select = s3_not_to_select,
    mprop_mh = mprop_mh,
    vprop_mh = vprop_mh
  )
}




# init_state <- function(config) {
#   sigma2 <- rep(0, config$num_iterations + 1)

#   phi <- lapply(
#     seq_len(config$num_iterations + 1),
#     function(k) matrix(0, config$num_series, config$total_parameters)
#   )

#   s1 <- rep(0, config$num_series + 1)


#   # For the estimation of parameters subject to selection (hdim)

#   if (config$num_parameters_to_select != 0) {
#     beta_hdim <- lapply(
#       seq_len(config$num_iterations + 1),
#       function(k) matrix(0, config$num_covariates_to_select + 1, config$num_parameters_to_select)
#     )

#     gamma_hdim <- lapply(
#       seq_len(config$num_iterations + 1),
#       function(k) matrix(0, config$num_parameters_to_select, config$num_parameters_to_select)
#     )

#     alpha <- lapply(
#       seq_len(config$num_iterations + 1),
#       function(k) matrix(0, config$num_parameters_to_select, 1)
#     )

#     s2_hdim <- lapply(
#       seq_len(config$num_iterations + 1),
#       function(k) matrix(0, config$num_parameters_to_select, config$num_parameters_to_select)
#     )
#     s3_hdim <- lapply(
#       seq_len(config$num_iterations + 1),
#       function(k) matrix(0, config$num_series, config$num_parameters_to_select)
#     )
#   } else {
#     beta_hdim <- NULL
#     gamma_hdim <- NULL
#     alpha <- NULL
#     s2_hdim <- NULL
#     s3_hdim <- NULL
#   }


#   # For the estimation of parameters non subject to selection (ldim)

#   if (config$num_parameters_not_to_select != 0) {
#     beta_ldim <- lapply(
#       seq_len(config$num_iterations + 1),
#       function(k) matrix(0, config$num_covariates_not_to_select + 1, config$num_parameters_not_to_select)
#     )

#     gamma_ldim <- lapply(
#       seq_len(config$num_iterations + 1),
#       function(k) matrix(0, config$num_parameters_not_to_select, config$num_parameters_not_to_select)
#     )

#     s2_ldim <- lapply(
#       seq_len(config$num_iterations + 1),
#       function(k) matrix(0, config$num_parameters_not_to_select, config$num_parameters_not_to_select)
#     )
#     s3_ldim <- lapply(
#       seq_len(config$num_iterations + 1),
#       function(k) matrix(0, config$num_series, config$num_parameters_not_to_select)
#     )

#     beta_mat_0 <- rep(0, config$num_parameters_not_to_select * (config$num_covariates_not_to_select + 1))
#   } else {
#     beta_ldim <- NULL
#     gamma_ldim <- NULL
#     s2_ldim <- lapply(
#       seq_len(config$num_iterations + 1),
#       function(k) NULL
#     )
#     s3_ldim <- lapply(
#       seq_len(config$num_iterations + 1),
#       function(k) NULL
#     )
#     beta_mat_0 <- NULL
#   }




#   sigma2[1] <- config$init_parameters@sigma2


#   if (length(config$parameters_to_select_indices) > 0) {
#     beta_hdim[[1]] <- config$init_parameters@beta_to_select
#     gamma_hdim[[1]] <- config$init_parameters@gamma_to_select
#     phi[[1]][, config$parameters_to_select_indices] <- config$x_phi_to_select %*% beta_hdim[[1]]
#     alpha[[1]] <- config$init_parameters@inclusion_prob
#   }

#   if (length(config$parameters_not_to_select_indices) > 0) {
#     beta_ldim[[1]] <- config$init_parameters@beta_not_to_select
#     gamma_ldim[[1]] <- config$init_parameters@gamma_not_to_select
#     phi[[1]][, config$parameters_not_to_select_indices] <- config$x_phi_not_to_select %*% beta_ldim[[1]]
#   }

#   # -- For the proposal distribution in the MH step
#   mprop_mh <- phi
#   vprop_mh <- lapply(
#     seq_len(config$num_iterations + 1),
#     function(k) matrix(0, config$total_parameters, config$total_parameters)
#   )
#   vprop_mh[[1]][config$parameters_to_select_indices, config$parameters_to_select_indices] <- gamma_hdim[[1]]
#   if (length(config$parameters_not_to_select_indices) > 0) {
#     vprop_mh[[1]][config$parameters_not_to_select_indices, config$parameters_not_to_select_indices] <-
#       gamma_ldim[[1]]
#   }

#   list(
#     beta_hdim = beta_hdim,
#     gamma_hdim = gamma_hdim,
#     beta_ldim = beta_ldim,
#     gamma_ldim = gamma_ldim,
#     beta_mat_0 = beta_mat_0,
#     sigma2 = sigma2,
#     alpha = alpha,
#     phi = phi,
#     s1 = s1,
#     s2_hdim = s2_hdim,
#     s3_hdim = s3_hdim,
#     s2_ldim = s2_ldim,
#     s3_ldim = s3_ldim,
#     mprop_mh = mprop_mh,
#     vprop_mh = vprop_mh
#   )
# }
