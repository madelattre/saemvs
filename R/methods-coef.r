#' Extract labeled regression coefficients for one support
#'
#' Internal utility function used by \code{coef()} and
#' \code{coef_support()} to extract the fixed-effects regression
#' coefficients \eqn{\beta} associated with a given support.
#'
#' The function returns the coefficient matrix with informative row
#' and column names:
#' \itemize{
#'   \item rows correspond to the intercept, forced covariates, and
#'         selected covariates
#'   \item columns correspond to the \eqn{\phi} parameters
#' }
#'
#' @param object A \code{saemvsResults} object.
#' @param support_idx Integer. Index of the support to use.
#'
#' @return
#' A numeric matrix of fixed-effects coefficients \eqn{\beta} with
#' labeled rows and columns.
#'
#' @keywords internal
.coef_one_support <- function(object, support_idx) {

  beta_est <- object@mle_estimates[[support_idx]]$beta
  if (is.null(beta_est)) {
    stop("No beta coefficients found for support ", support_idx)
  }

  phi_names <- object@phi_names
  colnames(beta_est) <- phi_names

  x_candidates_names <- object@x_candidates_names %||% character(0)
  x_forced_names <- object@x_forced_names %||% character(0)

  forced_variables_idx <- object@forced_variables_idx[[support_idx]]
  selected_variables_idx <- object@selected_variables_idx[[support_idx]]

  ## ---- Build row labels ----
  row_labels <- character(nrow(beta_est))
  row_labels[1] <- "\u03BC"  # intercept

  # Forced covariates
  if (length(forced_variables_idx) > 0) {
    forced_rows <- 1 + seq_along(forced_variables_idx)
    for (i in seq_along(forced_rows)) {
      idx <- forced_variables_idx[i]
      row_labels[forced_rows[i]] <-
        x_forced_names[idx] %||% paste0("forced_", idx)
    }
  } else {
    forced_rows <- integer(0)
  }

  # Selected covariates
  if (length(selected_variables_idx) > 0) {
    candidate_rows <- setdiff(2:nrow(beta_est), forced_rows)
    for (i in seq_along(selected_variables_idx)) {
      idx <- selected_variables_idx[i]
      row_labels[candidate_rows[i]] <-
        x_candidates_names[idx] %||% paste0("selected_", idx)
    }
  }

  rownames(beta_est) <- row_labels
  beta_est
}


#' Extract regression coefficients from a SAEMVS fit
#'
#' \describe{
#'   \item{\code{coef()}}{Returns the fixed-effects regression coefficients
#'   associated with the best support, defined as the one minimizing the
#'   selection criterion.}
#'   \item{\code{coef_support()}}{Returns the fixed-effects regression
#'   coefficients associated with a user-specified support.}
#' }
#'
#' @param object A \code{saemvsResults} object.
#' @param support_idx Integer. Index of the support to use (only for
#'   \code{coef_support()}).
#' @param ... Not used. Included for method consistency.
#'
#' @return
#' A numeric matrix of fixed-effects coefficients \eqn{\beta}, where:
#' \itemize{
#'   \item rows correspond to the intercept, forced covariates, and
#'         selected covariates
#'   \item columns correspond to the individual parameters \eqn{\phi}
#' }
#'
#' @details
#' The ordering and labeling of the rows matches the coefficient tables
#' displayed by \code{\link{summary}} and \code{\link{summary_support}},
#' ensuring consistency between printed summaries and programmatic
#' access to the estimates.
#'
#' @seealso
#' \code{\link{summary}}, \code{\link{predict}}, \code{\link{summary_support}}
#'
#' @examples
#' \dontrun{
#' # Coefficients for the best support
#' coef(res)
#'
#' # Coefficients for a specific support
#' coef_support(res, support_idx = 2)
#' }
#' @name coef
#' @rdname coef
#' @export
coef_support <- function(object, support_idx, ...) {
  stopifnot(inherits(object, "saemvsResults"))

  n_sup <- length(object@mle_estimates)
  if (support_idx < 1 || support_idx > n_sup) {
    stop("support_idx must be between 1 and ", n_sup)
  }

  .coef_one_support(object, support_idx)
}

#' @rdname coef
#' @export
setMethod(
  "coef",
  signature(object = "saemvsResults"),
  function(object, ...) {
    best_idx <- which.min(object@criterion_values)
    .coef_one_support(object, best_idx)
  }
)
