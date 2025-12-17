#' Compute log-likelihood (or penalized log-likelihood) for a SAEMVS model
#'
#' This internal function computes the (penalized) log-likelihood of the data
#' given the model and parameter estimates. The method supports different
#' penalty types including "BIC" and "e-BIC".
#'
#' @param data A \code{saemvsData} object containing the response and design
#' matrices.
#' @param model A \code{saemvsProcessedModel} object defining the model
#'  structure.
#' @param tuning_algo A \code{saemvsTuning} object containing tuning settings
#' (e.g., number of Monte Carlo samples).
#' @param param A list with parameter estimates. Must contain:
#'   \describe{
#'     \item{beta}{Matrix of fixed-effect coefficients.}
#'     \item{gamma}{Covariance matrix of random effects.}
#'     \item{sigma2}{Residual variance.}
#'   }
#' @param pen Character string specifying the penalty. One of "BIC", "e-BIC",
#' or no penalty.
#' @param p Numeric, number of candidate covariates (for e-BIC computation).
#' Consider renaming to \code{num_covariates} for clarity.
#' @param phi_to_select_idx Integer vector indicating which covariates
#' (columns of the design matrix) are subject to selection. Used for computing
#' the penalization term in BIC/e-BIC.
#' @param nb_forced_beta Integer specifying the number of covariates that are
#' forced to be included in the model. This is used to adjust the penalization
#' term in the e-BIC computation.
#' @return Numeric value of the penalized log-likelihood (BIC or e-BIC).
#' @keywords internal
#'
setGeneric(
  "loglik",
  function(data,
           model,
           tuning_algo,
           param,
           pen,
           p,
           phi_to_select_idx,
           nb_forced_beta) {
    standardGeneric("loglik")
  }
)

setMethod(
  "loglik",
  signature(
    data = "saemvsData",
    model = "saemvsProcessedModel",
    tuning_algo = "saemvsTuning",
    param = "list",
    pen = "character",
    p = "numeric",
    phi_to_select_idx = "numeric",
    nb_forced_beta = "numeric"
  ),
  function(data,
           model,
           tuning_algo,
           param,
           pen,
           p,
           phi_to_select_idx,
           nb_forced_beta) {
    data <- prepare_data(data, model)

    num_samples <- tuning_algo@n_is_samples

    yi <- data@y_series
    ti <- data@t_series

    n <- length(yi)
    ni <- lengths(yi)


    beta_x <- data@x_phi_not_to_select %*% param$beta
    beta_x_list <- split(beta_x, row(beta_x))

    gamma <- param$gamma
    sigma2 <- param$sigma2

    phi_samples <- lapply(beta_x_list, function(mu_i) {
      mvnfast::rmvn(num_samples, mu = mu_i, sigma = gamma)
    })

    log_lik_i <- function(i) {
      contribution <- sum(apply(phi_samples[[i]], 1, function(phi_i) {
        exp(-sum((yi[[i]] - g_vector_cpp(phi_i, ti[[i]]))^2) / (2 * sigma2))
      }))

      log((2 * pi * sigma2)^(-ni[i] / 2) * contribution / num_samples)
    }

    loglike <- sum(sapply(seq_along(yi), log_lik_i))

    nb_selected_beta <- sum(model@x_forced_support[, phi_to_select_idx])

    biclike <- switch(pen,
      "e-BIC" = -2 * loglike + (nb_selected_beta) * log(n) +
        2 * log(choose(
          p * length(phi_to_select_idx),
          (nb_selected_beta - nb_forced_beta)
        )),
      "BIC" = -2 * loglike + nb_selected_beta * log(n),
      -2 * loglike
    )

    return(biclike) # nolint: return-linter
  }
)
