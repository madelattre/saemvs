############################################################
# Print methods for each exported class
############################################################

# --- saemvsData print method ---
#' Print method for objects of class \code{saemvsData}
#'
#' This S4 method displays a concise summary of a \code{saemvsData} object,
#' including:
#' \itemize{
#'   \item Number of individuals (length of \code{y_series})
#'   \item Observations (\code{y_series}), showing only the first values of
#'  each individual
#'   \item Within-subject repeated factor (\code{t_series}), truncated for
#'  readability
#'   \item Covariates to be selected (\code{x_candidates}), showing a
#'  truncated matrix
#'   \item Covariates forced into the model (\code{x_forced}), if any
#' }
#' Long vectors and matrices are truncated to improve readability. The full
#'  data can
#' be accessed directly via the object's slots.
#'
#' @param x An object of class \code{saemvsData}.
#'
#'
#' @examples
#' # Assume 'data' is a saemvsData object
#' \dontrun{
#' print(data)
#' }
#'
#' @aliases print,saemvsData-method
#' @export
setMethod(
  "print",
  signature(x = "saemvsData"),
  function(x) {
    object <- x
    cat("Object of class 'saemvsData'\n\n")
    cat("- Number of individuals (y_series):", length(object@y_series), "\n\n")
    cat("- Observations:\n\n")
    cat(.format_list_of_vectors(object@y_series), "\n\n")
    cat("- Within-suject repeated factor (t_series):\n\n")
    cat(.format_list_of_vectors(object@t_series), "\n\n")
    cat(
      "- Covariates to be selected (x_candidates):\n",
      .format_matrix(object@x_candidates, rownames = FALSE), "\n\n"
    )
    cat(
      "- Unselected covariates to include in the model (x_forced):\n",
      .format_matrix(object@x_forced, rownames = FALSE), "\n"
    )
    invisible(object)
  }
)

#' Print method for objects of class \code{saemvsModel}
#'
#' This S4 method displays a concise summary of a \code{saemvsModel} object.
#' It prints the model function, the parameter names, parameters subject to
#'  selection, non-random parameters, and any support for forced included
#'  covariates.
#' Vectors are printed on a single line, truncated when necessary for
#' readability.
#'
#' @param x An object of class \code{saemvsModel}.
#'
#'
#' @examples
#' # Assuming 'mod' is a saemvsModel object
#' \dontrun{
#' print(mod)
#' }
#' @aliases print,saemvsModel-method
#' @export
setMethod(
  "print",
  signature(x = "saemvsModel"),
  function(x) {
    object <- x
    cat("Object of class 'saemvsModel'\n\n")

    cat("- Model function:\n")
    f_body <- deparse(object@model_func)
    f_body <- paste(f_body, collapse = "\n")
    f_body <- gsub("^", "  ", f_body)
    cat(f_body, "\n\n")

    cat("- Parameter names (phi_names):\n")
    if (is.null(object@phi_names)) {
      cat("  NULL\n\n")
    } else {
      cat("  ", paste(object@phi_names, collapse = " "), "\n\n", sep = "")
    }

    cat("- Parameters subject to selection (phi_to_select):\n")
    if (is.null(object@phi_to_select)) {
      cat("  NULL\n\n")
    } else {
      cat("  ", paste(object@phi_to_select, collapse = " "), "\n\n", sep = "")
    }

    cat("- Non-random parameters (phi_fixed):\n")
    if (is.null(object@phi_fixed)) {
      cat("  NULL\n\n")
    } else {
      cat("  ", paste(object@phi_fixed, collapse = " "), "\n\n", sep = "")
    }

    cat("- Support for forced included covariates (x_forced_support):\n")
    if (is.null(object@x_forced_support) ||
          length(object@x_forced_support) == 0) {
      cat("  NULL\n")
    } else {
      cat("  ", paste(names(object@x_forced_support), collapse = " "), "\n")
    }

    invisible(object)
  }
)


# --- saemvsHyperSlab print method ---
#' Print method for objects of class \code{saemvsHyperSlab}
#'
#' Displays a concise summary of a \code{saemvsHyperSlab} object, including the
#'  prior distributions for model parameters. The output shows:
#' \itemize{
#'   \item The slab parameter controlling the prior variance of covariate
#'  coefficients
#'   \item Covariance of parameters subject to selection, assumed to follow an
#'  inverse Wishart prior
#'   \item Residual variance, assumed to follow an inverse Gamma prior (with
#'  parametrization
#'    shape = residual_variance_prior_shape / 2 and
#'    rate = residual_variance_prior_shape * residual_variance_prior_rate / 2)
#'   \item Prior variance for the intercept parameter
#' }
#' Matrices and vectors are formatted for readability.
#'
#' @param x An object of class \code{saemvsHyperSlab}.
#'
#'
#' @examples
#' # Assuming 'hyper_slab' is a saemvsHyperSlab object
#' \dontrun{
#' print(hyper_slab)
#' }
#' @aliases print,saemvsHyperSlab-method
#' @export
setMethod(
  "print",
  signature(x = "saemvsHyperSlab"),
  function(x) {
    object <- x
    res_shape <- object@residual_variance_prior_shape
    res_rate <- object@residual_variance_prior_rate
    cat("Object of class 'saemvsHyperSlab'\n\n")
    cat(
      "- Slab parameter for covariate coefficients (slab_parameter):",
      object@slab_parameter, "\n\n"
    )
    cat("- Covariance of parameters subject to selection (inverse Wishart prior):\n") # nolint : line-length-linter
    cat(
      "  - Degrees of freedom (cov_re_prior_df):",
      object@cov_re_prior_df, "\n"
    )
    cat("  - Scale matrix (cov_re_prior_scale):\n")
    cat("    ", .format_matrix(object@cov_re_prior_scale,
        rownames = FALSE, colnames = FALSE
      ),
      "\n\n",
      sep = ""
    )
    shape_param <- res_shape / 2
    rate_param <- res_shape * res_rate / 2
    cat("- Residual variance (inverse Gamma prior):\n")
    cat(
      "  - Shape parameter (shape = residual_variance_prior_shape / 2):",
      shape_param, "\n"
    )
    cat(
      "  - Rate parameter (rate = residual_variance_prior_shape * residual_variance_prior_rate / 2):", # nolint : line-length-linter
      rate_param, "\n\n"
    )
    cat(
      "- Intercept phi prior variance (phi_intercept_prior_variance):",
      object@phi_intercept_prior_variance, "\n\n"
    )

    invisible(object)
  }
)



# --- saemvsInit print method ---
#' Print method for objects of class \code{saemvsInit}
#'
#' Displays the initial values used by the SAEMVS algorithm. The output
#' provides a structured overview of the starting values for model parameters,
#'  including:
#' \itemize{
#'   \item The intercept parameter
#'   \item Regression coefficients for forced covariates
#'   \item Regression coefficients for covariates subject to selection
#'   \item Covariance matrix of the random effects
#'   \item Residual variance
#'   \item Indicator specifying whether default initialisation was used for
#'         covariate coefficients
#' }
#'
#' Matrices are formatted for readability and truncated when necessary.
#' Column names corresponding to model parameters are prefixed with \code{phi}
#' when displayed.
#'
#' @param x An object of class \code{saemvsInit}.
#'
#' @examples
#' # Assuming 'init' is a saemvsInit object
#' \dontrun{
#' print(init)
#' }
#'
#' @aliases print,saemvsInit-method
#' @export
setMethod(
  "print",
  signature(
    x = "saemvsInit"
  ),
  function(x) {
    object <- x
    cat("Object of class 'saemvsInit'\n\n")
    cat("- Intercept:", .format_vector(object@intercept), "\n\n")
    cat(
      "- Coefficients for forced covariates (beta_forced):\n",
      .format_matrix(object@beta_forced, phi_colnames = TRUE), "\n\n"
    )
    cat(
      "- Coefficients for covariates to be selected (beta_candidates):\n",
      .format_matrix(object@beta_candidates, phi_colnames = TRUE), "\n\n"
    )
    cat(
      "- Covariance of the random effects (cov_re):\n",
      .format_matrix(object@cov_re, phi_colnames = TRUE), "\n\n"
    )
    cat("- Residual variance (sigma2):", object@sigma2, "\n\n")
    cat(
      "default initialisation of covariate coefficients:",
      object@default, "\n\n"
    )
    invisible(object)
  }
)

# --- saemvsTuning print method ---
#' Print method for objects of class \code{saemvsTuning}
#'
#' Displays the tuning parameters used by the SAEMVS algorithm. The output
#' provides a structured overview of algorithmic settings, including:
#' \itemize{
#'   \item SAEM iteration and burn-in parameters
#'   \item Metropolis-Hastings settings used in the stochastic approximation
#'  step
#'   \item Spike-and-slab hyperparameter grid
#'   \item Importance sampling settings for log-likelihood computation
#'   \item Reproducibility and parallelisation options
#' }
#'
#'
#' @param x An object of class \code{saemvsTuning}.
#' @examples
#' # Assuming 'tuning' is a saemvsTuning object
#' \dontrun{
#' print(tuning)
#' }
#' @aliases print,saemvsTuning-method
#' @export
setMethod(
  "print",
  signature(x = "saemvsTuning"),
  function(x) {
    object <- x

    cat("Object of class 'saemvsTuning'\n\n")

    cat("- SAEM settings:\n")
    cat("  - Number of iterations (niter):", object@niter, "\n")
    cat("  - Burn-in period (nburnin):", object@nburnin, "\n\n")

    cat("- Metropolis-Hastings settings(S-step):\n")
    cat("  - Number of MH iterations (niter_mh):", object@niter_mh, "\n")
    cat("  - Proposal kernel (kernel_mh):", object@kernel_mh, "\n")
    cat(
      "  - Proposal scale (mh_proposal_scale):",
      object@mh_proposal_scale, "\n\n"
    )

    cat("- Spike-and-slab tuning:\n")
    cat(
      "  - Spike values grid (spike_values_grid):",
      .format_vector(object@spike_values_grid, full = TRUE), "\n\n"
    )

    cat("- Importance sampling (log-likelihood comutation):\n")
    cat("  - Number of IS samples (n_is_samples):", object@n_is_samples, "\n\n")

    cat("- Reproducibility and parallelisation:\n")
    cat("  - Seed:", object@seed, "\n")
    cat("  - Number of workers (nb_workers):", object@nb_workers, "\n\n")

    invisible(object)
  }
)


#' Print a Summary of SAEM Variable Selection Results
#'
#' Displays a concise overview of the results from a SAEM variable selection
#' procedure. The output provides high-level information about the model
#' selection process and identifies the best model according to the chosen
#' selection criterion.
#'
#' The printed summary includes:
#' \itemize{
#'   \item The selection criterion used (BIC or e-BIC).
#'   \item The number of spike parameter values explored.
#'   \item The number of candidate supports evaluated.
#'   \item The number of unique supports retained after removing duplicates.
#'   \item The criterion value of the best selected model.
#'   \item The number of covariates selected in the best model.
#' }
#'
#' For detailed parameter estimates and selected variables, use
#' \code{\link{summary}}.
#'
#' @param x An object of class \code{saemvsResults} containing the output of
#'   a SAEM variable selection run.
#'
#'
#' @examples
#' \dontrun{
#' res <- saemvs(...)
#' print(res)
#' }
#'
#' @aliases print,saemvsResults-method
#' @export
setMethod(
  "print",
  signature(x = "saemvsResults"),
  function(x) {
    object <- x

    best_idx <- which.min(object@criterion_values)

    n_spike <- length(object@spike_values_grid)
    n_unique <- length(object@unique_support)

    n_selected <- length(object@selected_variables_idx[[best_idx]])

    cat("Object of class 'saemvsResults'\n\n")

    cat(
      "- Selection criterion:",
      object@criterion, "\n"
    )
    cat(
      "- Number of spike values explored:",
      n_spike, "\n"
    )
    cat(
      "- Number of unique supports:",
      n_unique, "\n\n"
    )

    cat("Best model:\n")
    cat(
      "- Criterion value:",
      round(object@criterion_values[best_idx], 3), "\n"
    )
    cat(
      "- Number of selected covariates:",
      n_selected, "\n\n"
    )

    cat("Call summary() for detailed results.\n")

    invisible(object)
  }
)
