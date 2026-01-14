#' Compute individual-level predictions for one support
#'
#' Internal utility function used by \code{predict()} and
#' \code{predict_support()}. For a given support, it averages the last
#' \eqn{k} SAEM iterations of the individual parameters \eqn{\phi}.
#'
#' @param object A \code{saemvsResults} object.
#' @param support_idx Integer. Index of the support to use.
#' @param k Integer. Number of last SAEM iterations to average.
#'
#' @return A numeric matrix with:
#' \itemize{
#'   \item rows corresponding to individuals
#'   \item columns corresponding to \eqn{\phi} parameters
#' }
#'
#' @details
#' Let \eqn{\phi^{(t)}_{i,p}} denote the value of parameter \eqn{p} for
#' individual \eqn{i} at SAEM iteration \eqn{t}. This function returns:
#' \deqn{
#' \hat{\phi}_{i,p} = \frac{1}{k} \sum_{t=T-k+1}^{T} \phi^{(t)}_{i,p},
#' }
#' where \eqn{T} is the total number of iterations for the selected support.
#'
#' @keywords internal
.predict_one_support <- function(object, support_idx, k = 20) {
  phi_sup <- object@phi[[support_idx]]

  if (length(phi_sup) < k) {
    stop("Not enough iterations: k = ", k,
         ", available = ", length(phi_sup))
  }

  iter_idx <- seq.int(length(phi_sup) - k + 1, length(phi_sup))

  phi_array <- simplify2array(phi_sup[iter_idx])

  phi_mean <- apply(phi_array, c(1, 2), mean)

  colnames(phi_mean) <- object@phi_names
  rownames(phi_mean) <- rownames(phi_sup[[1]])

  phi_mean
}


#' Prediction of individual parameters from a SAEMVS fit
#'
#' Computes individual-level predictions of the random effects
#' \eqn{\phi} by averaging the last \eqn{k} SAEM iterations.
#'
#' \describe{
#'   \item{\code{predict()}}{Returns predictions for the best support,
#'   defined as the one minimizing the selection criterion.}
#'   \item{\code{predict_support()}}{Returns predictions for a user-specified
#'   support.}
#' }
#'
#' @param object A \code{saemvsResults} object.
#' @param support_idx Integer. Index of the support to use (only for
#'   \code{predict_support()}).
#' @param k Integer (default = 20). Number of last SAEM iterations to average.
#' @param ... Not used. Included for method consistency.
#'
#' @return
#' A numeric matrix with one row per individual and one column per
#' \eqn{\phi} parameter.
#'
#' @details
#' Predictions are obtained by averaging the individual parameters
#' \eqn{\phi} over the last \eqn{k} SAEM iterations of the selected support.
#' This provides a Monte Carlo approximation of the posterior mean of
#' \eqn{\phi}.
#'
#' @seealso
#' \code{\link{summary}}, \code{\link{summary_support}}
#'
#' @examples
#' \dontrun{
#' # Best support
#' phi_hat <- predict(res)
#'
#' # More smoothing
#' phi_hat <- predict(res, k = 50)
#'
#' # Specific support
#' phi_sup2 <- predict_support(res, support_idx = 2, k = 30)
#' }
#'
#' @rdname predict
#' @export
predict_support <- function(object, support_idx, k = 20) {
  stopifnot(inherits(object, "saemvsResults"))

  n_sup <- length(object@phi)
  if (support_idx < 1 || support_idx > n_sup) {
    stop("support_idx must be between 1 and ", n_sup)
  }

  .predict_one_support(object, support_idx, k)
}


#' @rdname predict
#' @export
setMethod(
  "predict",
  signature(object = "saemvsResults"),
  function(object, ..., k = 20) {
    best_idx <- which.min(object@criterion_values)
    .predict_one_support(object, best_idx, k)
  }
)