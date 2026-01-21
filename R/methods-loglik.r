#' Compute log-likelihood (or penalized log-likelihood) for a SAEMVS model
#'
#' Internal function to compute the (penalized) log-likelihood of the data
#' given the model and parameter estimates. Supports different penalties
#' such as "BIC" and "e-BIC".
#'
#' @param data A \code{saemvsData} object containing the response variables
#'   and design matrices.
#' @param model A \code{saemvsProcessedModel} object defining the model
#'  structure,
#'   including covariates, random effects, and indices of parameters to select.
#' @param tuning_algo A \code{saemvsTuning} object specifying algorithmic
#'   tuning parameters such as number of Monte Carlo samples.
#' @param param A list of parameter estimates with components:
#'   \describe{
#'     \item{beta}{Matrix of fixed-effect coefficients.}
#'     \item{gamma}{Covariance matrix of random effects.}
#'     \item{sigma2}{Residual variance.}
#'   }
#' @param pen Character string specifying the penalty type. One of
#'   \code{"BIC"}, \code{"e-BIC"}, or \code{""} for no penalty.
#' @param p Numeric. Total number of candidate covariates (used for e-BIC).
#' @param phi_to_select_idx Integer vector of indices of covariates subject
#'   to selection (used to adjust the penalization term).
#' @param nb_forced_beta Integer. Number of covariates forced to be included
#'   in the model (affects penalization in e-BIC).
#' @param backend A list of compiled functions (e.g., \code{g_vector}) for
#'   evaluating model predictions efficiently.
#'
#' @return Numeric. The penalized log-likelihood (BIC or e-BIC) value.
#'
#' @details
#' The function proceeds as follows:
#' \enumerate{
#'   \item Preprocess the data using \code{prepare_data()}.
#'   \item For each subject/series, generate Monte Carlo samples of the latent
#'         parameters (\code{phi}) from a multivariate normal distribution
#'         with mean \eqn{\beta \times X} and covariance \eqn{\gamma}.
#'   \item Compute the likelihood contribution for each series using
#'         the model function in \code{backend}.
#'   \item Sum over all series to obtain the total log-likelihood.
#'   \item Apply the specified penalty ("BIC" or "e-BIC") using the number of
#'         selected covariates and forced covariates.
#' }
#'
#' @keywords internal
setGeneric(
  "loglik",
  function(data,
           model,
           tuning_algo,
           param,
           pen,
           p,
           phi_to_select_idx,
           nb_forced_beta,
           backend) {
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
    nb_forced_beta = "numeric",
    backend = "list"
  ),
  function(data,
           model,
           tuning_algo,
           param,
           pen,
           p,
           phi_to_select_idx,
           nb_forced_beta, backend) {
    data <- prepare_data(data, model)

    num_samples <- tuning_algo@n_is_samples

    yi <- data@y_series
    ti <- data@t_series

    n <- length(yi)

    beta_x <- data@x_phi_not_to_select %*% param$beta
    beta_x_list <- split(beta_x, row(beta_x))

    gamma <- param$gamma
    sigma2 <- param$sigma2

    phi_samples <- lapply(beta_x_list, function(mu_i) {
      backend$rmvnorm_mat(mu_i, gamma, num_samples)
    })

    loglike <- backend$ll(yi, ti, phi_samples, sigma2)

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
