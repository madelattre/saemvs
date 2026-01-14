############################################################
# Print methods for SAEMVS package classes
############################################################

#' Print methods for SAEMVS package classes
#'
#' Concise display of the main slots for each SAEMVS class:
#' \code{saemvsData}, \code{saemvsModel}, \code{saemvsHyperSlab},
#' \code{saemvsHyperSpikeAndSlab}, \code{saemvsInit}, \code{saemvsTuning}.
#'
#' @param object An object of the corresponding SAEMVS class.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @name print
#' @rdname print
#' @export
#' @examples
#' # print(my_data)
#' # print(my_model)


############################################################
# Utility functions (internal, not exported)
############################################################

format_matrix <- function(mat, max_rows = 5, max_cols = 5, digits = 3,
                          rownames = TRUE, colnames = TRUE,
                          phi_colnames = FALSE, phi_prefix = "phi") {
  if (is.null(mat)) {
    return("NULL")
  }

  nr <- nrow(mat)
  nc <- ncol(mat)

  mat_disp <- mat
  if (nr > max_rows) mat_disp <- mat_disp[1:max_rows, , drop = FALSE]
  if (nc > max_cols) mat_disp <- mat_disp[, 1:max_cols, drop = FALSE]

  mat_disp <- as.data.frame(mat_disp)
  for (j in seq_along(mat_disp)) {
    x <- mat_disp[[j]]
    x_fmt <- formatC(x, digits = digits, format = "f")
    w <- max(nchar(x_fmt), na.rm = TRUE)
    mat_disp[[j]] <- formatC(x, digits = digits, format = "f", width = w)
  }

  if (colnames) {
    cn <- colnames(mat_disp)

    if (phi_colnames) {
      default_cn <- paste0("V", seq_len(ncol(mat_disp)))

      if (is.null(cn) || identical(cn, default_cn)) {
        cn <- paste0(phi_prefix, seq_len(ncol(mat_disp)))
      } else {
        cn <- paste0(phi_prefix, "_", cn)
      }
    }

    colnames(mat_disp) <- cn
  }
  mat_disp <- as.matrix(mat_disp)

  if (!rownames) rownames(mat_disp) <- NULL
  if (!colnames) colnames(mat_disp) <- NULL

  txt <- capture.output(print(mat_disp, quote = FALSE, right = TRUE))

  if (nr > max_rows) txt <- c(txt, "...")

  header <- paste0("(", nr, " rows x ", nc, " columns)")
  paste(c(header, txt), collapse = "\n")
}


format_vector <- function(
    x, max_length = 5, digits = 3, indent = "",
    full = FALSE) {
  n <- length(x)
  if (full || n <= max_length) {
    x_disp <- x
  } else {
    x_disp <- x[1:max_length]
  }

  x_fmt <- formatC(x_disp, digits = digits, format = "f")
  w <- max(nchar(x_fmt), na.rm = TRUE)
  x_fmt <- formatC(x_disp, digits = digits, format = "f", width = w)

  txt <- paste(x_fmt, collapse = ", ")
  if (!full && n > max_length) {
    txt <- paste0(txt, ", ... [length=", n, "]")
  }

  txt <- gsub("\n", paste0("\n", indent), txt)
  paste0("c(", txt, ")")
}


format_list_of_vectors <- function(lst,
                                   max_items = 3,
                                   max_vec_length = 5,
                                   digits = 3,
                                   indent = "  ") {
  if (length(lst) == 0) {
    return("list()")
  }

  out <- vapply(
    seq_len(min(length(lst), max_items)),
    function(i) {
      vec_str <- format_vector(
        lst[[i]],
        max_length = max_vec_length,
        digits = digits,
        indent = paste0(indent, "  ")
      )
      paste0(indent, "[[", i, "]] ", vec_str)
    },
    character(1)
  )

  if (length(lst) > max_items) {
    out <- c(
      out,
      paste0(indent, "... (", length(lst) - max_items, " more elements)")
    )
  }

  paste(out, collapse = "\n")
}


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
#' @param ... Additional arguments (currently unused) for compatibility.
#'
#' @return The function prints a summary to the console invisibly. It returns
#'  the input
#' object invisibly (typical behavior for \code{print()} methods).
#'
#' @examples
#' # Assume 'data' is a saemvsData object
#' print(data)
#'
#' @rdname print
#' @aliases print,saemvsData-method
#' @export
setMethod(
  "print",
  signature(x = "saemvsData"),
  function(x, ...) {
    object <- x
    cat("Object of class 'saemvsData'\n\n")
    cat("- Number of individuals (y_series):", length(object@y_series), "\n\n")
    cat("- Observations:\n\n")
    cat(format_list_of_vectors(object@y_series), "\n\n")
    cat("- Within-suject repeated factor (t_series):\n\n")
    cat(format_list_of_vectors(object@t_series), "\n\n")
    cat(
      "- Covariates to be selected (x_candidates):\n",
      format_matrix(object@x_candidates, rownames = FALSE), "\n\n"
    )
    cat(
      "- Unselected covariates to include in the model (x_forced):\n",
      format_matrix(object@x_forced, rownames = FALSE), "\n"
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
#' The print method is intended for interactive display and provides a quick
#'  overview of the model structure and associated parameters. It does not
#'  modify the object.
#'
#' @param x An object of class \code{saemvsModel}.
#' @param ... Additional arguments (currently ignored), for compatibility with
#'  the S4 generic \code{print}.
#'
#' @return Invisibly returns the input \code{x} object. The main purpose is
#'  side-effect printing.
#'
#' @examples
#' # Assuming 'mod' is a saemvsModel object
#' print(mod)
#'
#' @rdname print
#' @aliases print,saemvsModel-method
#' @export
setMethod(
  "print",
  signature(x = "saemvsModel"),
  function(x, ...) {
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
#' @param ... Additional arguments (currently ignored), for compatibility with
#'  the S4 generic \code{print}.
#'
#' @return Invisibly returns the input object. The main purpose is side-effect
#'  printing.
#'
#' @examples
#' # Assuming 'hyper_slab' is a saemvsHyperSlab object
#' print(hyper_slab)
#'
#' @rdname print
#' @aliases print,saemvsHyperSlab-method
#' @export
setMethod(
  "print",
  signature(x = "saemvsHyperSlab"),
  function(x, ...) {
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
    cat("    ", format_matrix(object@cov_re_prior_scale,
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
#' Displays the initial values used by the SAEM-VS algorithm. The output
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
#' @param ... Additional arguments (currently ignored), included for
#' compatibility with the S4 generic \code{print}.
#'
#' @return Invisibly returns the input object. The function is called for
#' its side effect of printing to the console.
#'
#' @examples
#' # Assuming 'init' is a saemvsInit object
#' print(init)
#'
#' @rdname print
#' @aliases print,saemvsInit-method
#' @export
setMethod(
  "print",
  signature(
    x = "saemvsInit"
  ),
  function(x, ...) {
    object <- x
    cat("Object of class 'saemvsInit'\n\n")
    cat("- Intercept:", format_vector(object@intercept), "\n\n")
    cat(
      "- Coefficients for forced covariates (beta_forced):\n",
      format_matrix(object@beta_forced, phi_colnames = TRUE), "\n\n"
    )
    cat(
      "- Coefficients for covariates to be selected (beta_candidates):\n",
      format_matrix(object@beta_candidates, phi_colnames = TRUE), "\n\n"
    )
    cat(
      "- Covariance of the random effects (cov_re):\n",
      format_matrix(object@cov_re, phi_colnames = TRUE), "\n\n"
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
#' Displays the tuning parameters used by the SAEM-VS algorithm. The output
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
#' The spike values grid is fully displayed, while other parameters are shown
#' as scalar values for clarity.
#'
#' @param x An object of class \code{saemvsTuning}.
#' @param ... Additional arguments (currently ignored), included for
#' compatibility with the S4 generic \code{print}.
#'
#' @return Invisibly returns the input object. This method is called for
#' its side effect of printing to the console.
#'
#' @examples
#' # Assuming 'tuning' is a saemvsTuning object
#' print(tuning)
#'
#' @rdname print
#' @aliases print,saemvsTuning-method
#' @export
setMethod(
  "print",
  signature(x = "saemvsTuning"),
  function(x, ...) {
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
      format_vector(object@spike_values_grid, full = TRUE), "\n\n"
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
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' The object \code{x}, invisibly. This method is called for its side effect
#' of printing a summary to the console.
#'
#' @seealso
#' \code{\link{summary}}, \code{\link{saemvsResults-class}}
#'
#' @examples
#' \dontrun{
#' res <- saemvs(...)
#' print(res)
#' }
#'
#' @aliases print,saemvsResults-method
#' @rdname print
#' @export
setMethod(
  "print",
  signature(x = "saemvsResults"),
  function(x, ...) {
    object <- x

    best_idx <- which.min(object@criterion_values)

    n_spike <- length(object@spike_values_grid)
    n_support <- length(object@support)
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
      "- Number of candidate supports:",
      n_support, "\n"
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
