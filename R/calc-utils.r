#' Compute threshold for spike-and-slab prior
#'
#' Computes the threshold value for variable selection in a spike-and-slab prior,
#'              based on prior variances and inclusion probability.
#' @param nu1 Numeric. Variance of the "slab" component (large variance).
#' @param nu0 Numeric. Variance of the "spike" component (small variance).
#' @param alpha Numeric. Prior inclusion probability of the variable.
#' @return Numeric. Threshold value used in spike-and-slab selection.
#' @keywords internal
threshold <- function(nu1, nu0, alpha) {
  value <- 2 * nu0 * nu1 * log(sqrt(nu1 / nu0) * (1 - alpha) / alpha) / (nu1 - nu0)
  return(sqrt(value))
}

#' Compute posterior inclusion probability p*
#'
#' Computes the posterior probability that a coefficient belongs to the "slab" component
#'              of a spike-and-slab prior. This is used for fast variable selection updates.
#' @param beta Numeric matrix. Current coefficient estimates (rows: covariates, columns: responses/parameters).
#' @param alpha Numeric vector. Prior inclusion probabilities for each coefficient.
#' @param nu0 Numeric. Variance of the spike (small variance) component.
#' @param nu1 Numeric. Variance of the slab (large variance) component.
#' @return Numeric matrix of the same dimension as \code{beta}. Each entry is the posterior probability p*.
#' @keywords internal
p_star <- function(beta, alpha, nu0, nu1) {
  norm1 <- stats::dnorm(beta, mean = 0, sd = sqrt(nu1))
  norm0 <- stats::dnorm(beta, mean = 0, sd = sqrt(nu0))
  num <- norm1 * rep(alpha, each = nrow(beta))
  denom <- num + norm0 * rep(1 - alpha, each = nrow(beta))
  p_star <- num / denom
  return(p_star)
}

#' Determine algorithmic case based on model configuration
#'
#' Determines the operational case for the algorithm (e.g., MAP vs MLE, full vs partial selection,
#' presence of fixed parameters) based on the provided configuration.
#' @param config List. Configuration object containing the following elements:
#'   \itemize{
#'     \item \code{method_type} Character. Either "map" or "mle".
#'     \item \code{num_parameters_not_to_select} Integer. Number of parameters not subject to selection.
#'     \item \code{fixed_parameters_indices} Integer vector. Indices of parameters fixed a priori.
#'   }
#' @return Character string indicating the operational case. Possible values:
#'   \itemize{
#'     \item "map_full_select": MAP with all parameters subject to selection.
#'     \item "map_part_select_nofixed": MAP with partial selection, no fixed parameters.
#'     \item "map_part_select_fixed": MAP with partial selection and some fixed parameters.
#'     \item "mle_nofixed": MLE with no fixed parameters.
#'     \item "mle_fixed": MLE with some fixed parameters.
#'   }
#' @keywords internal
get_case <- function(config) {
  if (config$method_type == "map") {
    if (config$num_parameters_not_to_select == 0) {
      "map_full_select"
    } else if (length(config$fixed_parameters_indices) == 0) {
      "map_part_select_nofixed"
    } else {
      "map_part_select_fixed"
    }
  } else {
    if (length(config$fixed_parameters_indices) > 0) "mle_fixed" else "mle_nofixed"
  }
}
