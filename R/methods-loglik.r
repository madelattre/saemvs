#' Compute log-likelihood (or penalized log-likelihood) for a SAEMVS model
#'
#' This internal function computes the (penalized) log-likelihood of the data
#' given the model and parameter estimates. The method supports different
#' penalty types including "BIC" and "e-BIC".
#'
#' @param data A \code{saemvsData} object containing the response and design matrices.
#' @param model A \code{saemvsModel} object defining the model structure.
#' @param tuning_algo A \code{saemvsTuning} object containing tuning settings (e.g., number of Monte Carlo samples).
#' @param param A list with parameter estimates. Must contain:
#'   \describe{
#'     \item{beta}{Matrix of fixed-effect coefficients.}
#'     \item{gamma}{Covariance matrix of random effects.}
#'     \item{sigma2}{Residual variance.}
#'   }
#' @param pen Character string specifying the penalty. One of "BIC", "e-BIC", or no penalty.
#' @param p Numeric, number of candidate covariates (for e-BIC computation). Consider renaming to \code{num_covariates} for clarity.
#' @return Numeric value of the penalized log-likelihood (BIC or e-BIC).  
#' @keywords internal
loglik <- function(data, model, tuning_algo, param, pen, p) {
  
  # Prepare the data
  data <- prepare_data(data, model)
  
  # Number of Monte Carlo samples for marginal likelihood approximation
  num_samples <- tuning_algo@n_is_samples
  
  # Forced support for phi parameters
  support_matrix <- model@x_forced_support
  supp_index <- which(c(t(support_matrix)) == 1)
  
  yi <- data@y_series  # List of response vectors for each individual
  ti <- data@t_series  # List of time vectors for each individual
  
  n <- length(yi)      # Number of individuals
  ni <- lengths(yi)    # Observations per individual
  
  # Compute X %*% beta for non-selected parameters
  beta_x <- data@x_phi_not_to_select %*% param$beta
  beta_x_list <- split(beta_x, row(beta_x))
  
  gamma <- param$gamma
  sigma2 <- param$sigma2
  
  # Simulate phi samples from multivariate normal
  phi_samples <- lapply(beta_x_list, function(mu_i) {
    mvnfast::rmvn(num_samples, mu = mu_i, sigma = gamma)
  })
  
  # Compute marginal likelihood for individual i
  log_lik_i <- function(i) {
    # Sum over simulated phi
    contribution <- sum(apply(phi_samples[[i]], 1, function(phi_i) {
      exp(-sum((yi[[i]] - g_vector_cpp(phi_i, ti[[i]]))^2) / (2 * sigma2))
    }))
    
    # Contribution to marginal likelihood
    log((2 * pi * sigma2)^(-ni[i] / 2) * contribution / num_samples)
  }
  
  # Sum over individuals
  loglike <- sum(sapply(seq_along(yi), log_lik_i))
  
  nb_beta <- length(supp_index)
  
  # Apply penalty if requested
  biclike <- switch(
    pen,
    "e-BIC" = -2 * loglike + (nb_beta - model@phi_dim) * log(n) +
      2 * log(choose(p * model@phi_dim, nb_beta - model@phi_dim)),
    "BIC" = -2 * loglike + nb_beta * log(n),
    -2 * loglike
  )
  
  return(biclike)
}
