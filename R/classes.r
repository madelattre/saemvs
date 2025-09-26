#' @title Class Union: matrix or NULL
#' @description Internal union class used in slots that may contain
#' a numeric matrix or be \code{NULL}.
#' @keywords internal
setClassUnion("matrixORNULL", c("matrix", "NULL"))

#' @title Class Union: numeric or NULL
#' @description Internal union class used in slots that may contain
#' a numeric vector or be \code{NULL}.
#' @keywords internal
setClassUnion("numericORNULL", c("numeric", "NULL"))

#' @title Class Union: list or NULL
#' @description Internal union class used in slots that may contain
#' a list or be \code{NULL}.
#' @keywords internal
setClassUnion("listORNULL", c("list", "NULL"))


#' @title saemvsData class
#' @description A class to store user-provided data for the SAEMVS algorithm.
#'
#' This is the main entry point for users: responses, time indices,
#' candidate covariates for selection, and forced covariates.
#'
#' @slot y_series list of numeric vectors. Each element contains observed
#' responses for one individual.
#' @slot t_series list of numeric vectors. Each element contains corresponding
#' time indices.
#' @slot x_candidates matrix or NULL. Optional matrix of candidate covariates
#' (rows = individuals, columns = covariates).
#' @slot x_forced matrix or NULL. Optional matrix of forced covariates
#' (rows = individuals, columns = covariates).
#'
#' @section Validity:
#' \itemize{
#'   \item \code{y_series} and \code{t_series} must have the same length.
#'   \item Each element of \code{y_series} and \code{t_series} must be numeric
#'  vectors of the same length.
#'   \item If provided, \code{x_candidates} and \code{x_forced} must each have
#'  as many rows as the length of \code{y_series}.
#'   \item All matrices must be numeric.
#' }
#'
#' @return An object of class \code{saemvsData}.
#' @exportClass saemvsData
#'
#' @examples
#' \dontrun{
#' y <- list(rnorm(10), rnorm(8))
#' t <- list(1:10, 1:8)
#' d <- saemvsData(y, t)
#' }
setClass(
  "saemvsData",
  slots = list(
    y_series     = "list",
    t_series     = "list",
    x_candidates = "matrixORNULL",
    x_forced     = "matrixORNULL"
  ),
  prototype = list(
    y_series     = list(),
    t_series     = list(),
    x_candidates = NULL,
    x_forced     = NULL
  ),
  validity = function(object) {
    n <- length(object@y_series)

    # Check lengths of y_series and t_series
    if (n != length(object@t_series)) {
      return("y_series and t_series must have the same number of elements.")
    }

    # Check each element is numeric and lengths match
    for (i in seq_len(n)) {
      if (!is.numeric(object@y_series[[i]])) {
        return(sprintf("y_series[[%d]] must be numeric.", i))
      }
      if (!is.numeric(object@t_series[[i]])) {
        return(sprintf("t_series[[%d]] must be numeric.", i))
      }
      if (length(object@y_series[[i]]) != length(object@t_series[[i]])) {
        return(sprintf("y_series[[%d]] and t_series[[%d]] do not have the same
        length.", i, i))
      }
    }

    # Check candidate and forced covariates
    for (mat_name in c("x_candidates", "x_forced")) {
      mat <- slot(object, mat_name)
      if (!is.null(mat)) {
        if (!is.matrix(mat)) {
          return(sprintf("'%s' must be a matrix or NULL.", mat_name))
        }
        if (!is.numeric(mat)) {
          return(sprintf("All entries in '%s' must be numeric.", mat_name))
        }
        if (nrow(mat) != n) {
          return(sprintf(
            "'%s' must have %d rows (matching y_series); currently %d.",
            mat_name, n, nrow(mat)
          ))
        }
      }
    }

    # Check that there are no identical columns between x_candidates and
    # x_forced
    if (!is.null(object@x_candidates) && !is.null(object@x_forced)) {
      candidates <- object@x_candidates
      forced <- object@x_forced

      # Loop over each column in x_forced
      for (j in seq_len(ncol(forced))) {
        col_forced <- forced[, j, drop = FALSE]
        for (i in seq_len(ncol(candidates))) {
          col_candidate <- candidates[, i, drop = FALSE]
          if (all(col_forced == col_candidate)) {
            return(paste0(
              "Column ", j,
              " of 'x_forced' is identical to column ", i,
              " of 'x_candidates'."
            ))
          }
        }
      }
    }


    TRUE
  }
)


#' @rdname saemvsData
#' @title Constructor for saemvsData
#'
#' @param y list of numeric vectors (responses).
#' @param t list of numeric vectors (time indices).
#' @param x_candidates matrix or NULL, candidate covariates for selection.
#' @param x_forced matrix or NULL, forced covariates.
#'
#' @return An object of class \code{saemvsData}.
#' @export
saemvsData <- function(y, t, x_candidates = NULL, x_forced = NULL) {
  methods::new("saemvsData",
    y_series     = y,
    t_series     = t,
    x_candidates = x_candidates,
    x_forced     = x_forced
  )
}

#' @title Constructor from data.frames for saemvsData
#'
#' @description
#' Creates a \code{saemvsData} object from a longitudinal dataset and a
#' covariate dataset.
#' The function splits the longitudinal data by individual to build
#' \code{y_series} and \code{t_series}, and constructs \code{x_candidates}
#' and \code{x_forced} matrices from the covariate data.
#' @param long_df data.frame containing the longitudinal measurements.
#' Must include columns for ID, response, and time.
#' @param covar_df data.frame containing fixed covariates (one row per
#' individual).
#' @param id_col name of the column identifying individuals (must exist in
#' both dataframes).
#' @param y_col name of the response column in \code{long_df}.
#' @param t_col name of the time column in \code{long_df}.
#' @param x_candidates_cols vector of column names in \code{covar_df} to use as
#' candidate covariates (optional). If NULL, all columns except \code{id_col}
#' and forced covariates are used.
#' @param x_forced_cols vector of column names in \code{covar_df} to use as
#' forced covariates (optional).
#'
#' @return An object of class \code{saemvsData}.
#'
#' @details
#' The function ensures that the order of individuals in \code{y_series}
#' matches the order of rows in \code{x_candidates} and \code{x_forced}.
#' Candidate covariates are assumed to be fixed per individual.
#'
#' @examples
#' # Example longitudinal data
#' long_df <- data.frame(
#'   id = rep(1:2, each = 5),
#'   t = rep(1:5, 2),
#'   y = rnorm(10)
#' )
#'
#' # Example covariate data
#' covar_df <- data.frame(
#'   id = 1:2,
#'   x1 = rnorm(2),
#'   x2 = rnorm(2)
#' )
#'
#' # Create saemvsData object
#' d <- saemvsDataFromDFs(
#'   long_df, covar_df,
#'   id_col = "id", y_col = "y", t_col = "t",
#'   x_candidates_cols = c("x1"),
#'   x_forced_cols = c("x2")
#' )
#'
#' @export
saemvsDataFromDFs <- function(long_df,
                              covar_df,
                              id_col,
                              y_col,
                              t_col,
                              x_candidates_cols = NULL,
                              x_forced_cols = NULL) {
  # Check that required columns exist
  for (nm in c(id_col, y_col, t_col)) {
    if (!nm %in% names(long_df)) {
      stop(sprintf("Column '%s' is missing in long_df.", nm))
    }
  }
  if (!id_col %in% names(covar_df)) {
    stop(sprintf("Column '%s' is missing in covar_df.", id_col))
  }

  # Split longitudinal data by individual
  split_data <- split(long_df, long_df[[id_col]])
  ids_long <- names(split_data)

  # Create y_series and t_series as lists
  y_series <- lapply(split_data, function(df) df[[y_col]])
  t_series <- lapply(split_data, function(df) df[[t_col]])

  # Reorder covariate dataframe to match the order of individuals in long_df
  covar_df_ordered <- covar_df[match(ids_long, covar_df[[id_col]]), ]
  if (any(is.na(covar_df_ordered[[id_col]]))) {
    stop("Some individuals in long_df are not present in covar_df.")
  }

  # Determine forced covariates
  if (!is.null(x_forced_cols)) {
    x_forced <- as.matrix(covar_df_ordered[, x_forced_cols, drop = FALSE])
    if (!is.numeric(x_forced)) stop("x_forced must be numeric.")
  } else {
    x_forced <- NULL
    x_forced_cols <- character(0)
  }

  # Determine candidate covariates
  if (is.null(x_candidates_cols)) {
    # Use all columns except id and forced covariates
    x_candidates_cols <- setdiff(
      colnames(covar_df_ordered),
      c(id_col, x_forced_cols)
    )
  }
  if (length(x_candidates_cols) > 0) {
    x_candidates <-
      as.matrix(covar_df_ordered[, x_candidates_cols, drop = FALSE])
    if (!is.numeric(x_candidates)) stop("x_candidates must be numeric.")
  } else {
    x_candidates <- NULL
  }

  # Build the final saemvsData object
  saemvsData(
    y = y_series,
    t = t_series,
    x_candidates = x_candidates,
    x_forced = x_forced
  )
}




#' Internal class: saemvsProcessedData
#'
#' A class to store processed design matrices ready for SAEMVS algorithms.
#' This class is internal to the package and is not intended to be
#' directly created or manipulated by the user.
#' It inherits from \code{saemvsData} and contains additional slots
#' representing design matrices of covariates associated with the model
#' parameters \code{phi}, split into components subject to selection and
#' components that are not subject to selection.
#'
#' @slot x_phi_to_select matrix or NULL. Design matrix of covariates for
#' \code{phi} parameters on which variable selection will be performed.
#' Rows correspond to individuals; columns to covariates.
#'
#' @slot x_phi_not_to_select matrix or NULL.
#' Design matrix of covariates for \code{phi} parameters that are not
#' subject to selection. Rows correspond to individuals; columns to
#' covariates.
#'
#' @slot tx_x_phi_to_select matrix or NULL.
#' Transformation of \code{x_phi_to_select} (e.g., multiplied by its transpose)
#' required for algorithm computations.
#'
#' @slot kron_tx_x_phi_to_select matrix or NULL. Kronecker product of
#' \code{tx_x_phi_to_select} used for vectorized calculations.
#'
#' @slot x_phi_not_to_select_list list or NULL.
#'   List of \code{x_phi_not_to_select} entries per individual. Each element
#'   corresponds to one individual and contains a numeric vector.
#'
#' @section Validity:
#'   \itemize{
#'     \item All matrices must be numeric if not NULL.
#'     \item \code{x_phi_not_to_select_list}, if not NULL, must have length
#'  equal to the number of individuals.
#'     \item Columns in \code{x_phi_to_select} and \code{x_phi_not_to_select}
#'  should not overlap (warning issued if they do).
#'     \item Inherits all validity checks from \code{saemvsData} (e.g., length
#'  and type checks for y_series, t_series, x_candidates, x_forced).
#'   }
#'
#' @keywords internal
setClass(
  "saemvsProcessedData",
  slots = list(
    x_phi_to_select          = "matrixORNULL",
    x_phi_not_to_select      = "matrixORNULL",
    tx_x_phi_to_select       = "matrixORNULL",
    kron_tx_x_phi_to_select  = "matrixORNULL",
    x_phi_not_to_select_list = "listORNULL"
  ),
  prototype = list(
    x_phi_to_select          = NULL,
    x_phi_not_to_select      = NULL,
    tx_x_phi_to_select       = NULL,
    kron_tx_x_phi_to_select  = NULL,
    x_phi_not_to_select_list = NULL
  ),
  contains = "saemvsData",
  validity = function(object) {
    n_ind <- length(object@y_series)

    # Check matrices are numeric if present
    mats <- list(
      x_phi_to_select         = object@x_phi_to_select,
      x_phi_not_to_select     = object@x_phi_not_to_select
    )
    for (mat_name in names(mats)) {
      mat <- mats[[mat_name]]
      if (!is.null(mat)) {
        if (!is.matrix(mat)) {
          return(sprintf(
            "'%s' must be a matrix or NULL.", mat_name
          ))
        }
        if (!is.numeric(mat)) {
          return(sprintf(
            "All entries in '%s' must be numeric.", mat_name
          ))
        }
        if (nrow(mat) != n_ind) {
          return(sprintf(
            "'%s' must have %d rows (matching y_series); currently %d.",
            mat_name, n_ind, nrow(mat)
          ))
        }
      }
    }

    # Check list length
    if (!is.null(object@x_phi_not_to_select_list) &&
      length(object@x_phi_not_to_select_list) != n_ind) {
      return(sprintf(
        "'x_phi_not_to_select_list' must have %d elements (matching y_series).",
        n_ind
      ))
    }

    TRUE
  }
)

#' Internal constructor for saemvsProcessedData
#' @keywords internal
saemvsProcessedData <- function(x_phi_to_select = NULL,
                                x_phi_not_to_select = NULL,
                                tx_x_phi_to_select = NULL,
                                kron_tx_x_phi_to_select = NULL,
                                x_phi_not_to_select_list = NULL,
                                ...) {
  methods::new("saemvsProcessedData",
    x_phi_to_select          = x_phi_to_select,
    x_phi_not_to_select      = x_phi_not_to_select,
    tx_x_phi_to_select       = tx_x_phi_to_select,
    kron_tx_x_phi_to_select  = kron_tx_x_phi_to_select,
    x_phi_not_to_select_list = x_phi_not_to_select_list,
    ...
  )
}

#' @title saemvsModel class
#' @description Represents a model for SAEMVS, including the model function and
#' indexing information for variable selection on parameters phi.
#'
#' @slot model_func A function of form function(phi, t) returning predicted
#'  values.
#' @slot phi_dim Integer: total number of phi parameters.
#' @slot phi_to_select_idx Integer vector or NULL: indices of phi parameters
#'   on which variable selection will be performed.
#' @slot phi_fixed_idx Integer vector or NULL: indices of phi parameters
#'   that are fixed (i.e. without random variability).
#' @slot x_forced_support Numeric matrix or NULL: design matrix for covariates
#'   that are forced in the model. Must have phi_dim columns.
#'
#' @exportClass saemvsModel
setClass(
  "saemvsModel",
  slots = list(
    model_func = "function",
    phi_dim = "numeric",
    phi_to_select_idx = "numericORNULL",
    phi_fixed_idx = "numericORNULL",
    x_forced_support = "matrixORNULL"
  ),
  prototype = list(
    model_func = function(phi, t) numeric(0),
    phi_dim = as.integer(1),
    phi_to_select_idx = integer(0),
    phi_fixed_idx = integer(0),
    x_forced_support = matrix(numeric(0), nrow = 0, ncol = 0)
  ),
  validity = function(object) {
    # Check model function arguments
    args <- names(formals(object@model_func))
    if (!all(c("phi", "t") %in% args)) {
      return("The model function must have arguments 'phi' and 't'.")
    }

    # Check phi_dim
    if (!(is.numeric(object@phi_dim) && length(object@phi_dim) == 1 &&
      object@phi_dim >= 1)) {
      return("'phi_dim' must be a positive integer of length 1.")
    }

    # Check indices
    if (!(is.null(object@phi_to_select_idx) ||
      all(object@phi_to_select_idx >= 1 &
        object@phi_to_select_idx <= object@phi_dim))) {
      return("
      'phi_to_select_idx' must contain integers between 1 and 'phi_dim'.
      ")
    }

    if (!(is.null(object@phi_fixed_idx) ||
      all(object@phi_fixed_idx >= 1 &
        object@phi_fixed_idx <= object@phi_dim))) {
      return("
      'phi_not_to_select_idx' must contain integers between 1 and 'phi_dim'.
      ")
    }

    # No overlap
    if (length(intersect(object@phi_to_select_idx, object@phi_fixed_idx)) > 0) {
      return("'phi_to_select_idx' and 'phi_fixed_idx' must not overlap.")
    }

    # x_forced_support
    supp <- object@x_forced_support
    if (!is_empty_support(supp)) {
      if (!is.matrix(supp)) {
        return("'x_forced_support' must be a matrix or NULL.")
      }
      if (ncol(supp) != object@phi_dim) {
        return(paste0(
          "'x_forced_support' must have ", object@phi_dim,
          " columns (one per phi)."
        ))
      }
    }

    TRUE
  }
)

#' @rdname saemvsModel
#' @title Constructor for saemvsModel
#'
#' @param g Function of form `function(phi, t)` returning predicted values.
#' @param phi_dim Integer: total number of phi parameters.
#' @param phi_to_select_idx Integer vector (optional): indices of phi
#'  parameters for variable selection.
#' @param phi_fixed_idx Integer vector (optional): indices of phi parameters
#'  fixed (no random variability).
#' @param x_forced_support Numeric matrix or NULL: design matrix for forced
#'  covariates (phi_dim columns).
#'
#' @return An object of class \code{saemvsModel}.
#' @export
saemvsModel <- function(
    g, phi_dim, phi_to_select_idx = c(), phi_fixed_idx = c(),
    x_forced_support = matrix(numeric(0), nrow = 0, ncol = 0)) {
  methods::new("saemvsModel",
    model_func = g,
    phi_dim = as.integer(phi_dim),
    phi_to_select_idx = as.integer(phi_to_select_idx),
    phi_fixed_idx = as.integer(phi_fixed_idx),
    x_forced_support = x_forced_support
  )
}


#' @title saemvsHyperSlab
#' @description Represents the hyperparameters for the slab component in a
#'  spike-and-slab prior.
#'
#' Only the slab parameter, random effects covariance prior, and degrees of
#'  freedom need to be specified by the user. Other hyperparameters are fixed
#'  internally.
#'
#' @slot slab_parameter Numeric: positive slab parameter (controls the slab
#'  prior).
#' @slot cov_re_prior_scale Numeric matrix: scale matrix of the Inverse-Wishart
#'  prior for random effects covariance.
#' @slot cov_re_prior_df Numeric: degrees of freedom of the Inverse-Wishart
#'  prior for random effects covariance.
#' @slot residual_variance_prior_shape Numeric: shape parameter of the
#'  inverse-gamma prior on residual variance (fixed internally).
#' @slot residual_variance_prior_rate Numeric: rate parameter of the
#'  inverse-gamma prior on residual variance (fixed internally).
#' @slot phi_intercept_prior_variance Numeric: variance of the Gaussian prior
#'  on each intercept of phi (fixed internally).
#' @slot inclusion_prob_prior_a Numeric or NULL: alpha parameter of the Beta
#'  prior on covariate inclusion probabilities.
#' @slot inclusion_prob_prior_b Numeric or NULL: beta parameter of the Beta
#'  prior on covariate inclusion probabilities.
#'
#' @exportClass saemvsHyperSlab
setClass(
  "saemvsHyperSlab",
  slots = list(
    slab_parameter = "numericORNULL",
    cov_re_prior_scale = "matrixORNULL",
    cov_re_prior_df = "numericORNULL",
    residual_variance_prior_shape = "numericORNULL",
    residual_variance_prior_rate = "numericORNULL",
    phi_intercept_prior_variance = "numericORNULL",
    inclusion_prob_prior_a = "numericORNULL",
    inclusion_prob_prior_b = "numericORNULL"
  ),
  prototype = list(
    slab_parameter = 12000,
    cov_re_prior_scale = matrix(numeric(0), 0, 0),
    cov_re_prior_df = 1,
    residual_variance_prior_shape = 1,
    residual_variance_prior_rate = 1,
    phi_intercept_prior_variance = 3000^2,
    inclusion_prob_prior_a = NULL,
    inclusion_prob_prior_b = NULL
  ),
  validity = function(object) {
    errors <- c()

    errors <- c(
      errors,
      check_positive_or_null_slot(object@slab_parameter, "slab_parameter")
    )
    errors <- c(
      errors,
      check_positive_or_null_slot(object@cov_re_prior_df, "cov_re_prior_df")
    )
    errors <- c(
      errors,
      check_positive_or_null_slot(
        object@residual_variance_prior_shape,
        "residual_variance_prior_shape"
      )
    )
    errors <- c(
      errors,
      check_positive_or_null_slot(
        object@residual_variance_prior_rate,
        "residual_variance_prior_rate"
      )
    )
    errors <- c(
      errors,
      check_positive_or_null_slot(
        object@phi_intercept_prior_variance,
        "phi_intercept_prior_variance"
      )
    )

    if (!is.null(object@cov_re_prior_scale)) {
      errors <- c(
        errors,
        check_covariance(object@cov_re_prior_scale, "cov_re_prior_scale")
      )
    }

    if (length(errors) == 0) TRUE else errors
  }
)

#' @rdname saemvsHyperSlab
#' @title Constructor for saemvsHyperSlab
#'
#' @description Create a \code{saemvsHyperSlab} object specifying
#'  hyperparameters for the slab component in a spike-and-slab prior.
#'
#' Only the slab parameter, the scale matrix for the random effects covariance,
#' and the degrees of freedom need to be specified by the user. All other
#' hyperparameters are fixed internally.
#'
#' @param slab_parameter Numeric, positive. Controls the variance of the slab
#' prior. Default is 12000.
#' @param cov_re_prior_scale Numeric matrix. Scale matrix of the Inverse-Wishart
#'  prior for the random effects covariance. Must be square with dimension
#'  matching the number of phi parameters.
#' @param cov_re_prior_df Numeric, positive. Degrees of freedom of the
#'  Inverse-Wishart prior for the random effects covariance. Default is 1.
#'
#' @return An object of class \code{saemvsHyperSlab}.
#'
#' @examples
#' \dontrun{
#' cov_scale <- diag(2)
#' h <- saemvsHyperSlab(
#'   slab_parameter = 1000,
#'   cov_re_prior_scale = cov_scale, cov_re_prior_df = 3
#' )
#' }
#'
#' @export
saemvsHyperSlab <- function(slab_parameter = 12000,
                            cov_re_prior_scale,
                            cov_re_prior_df = 1) {
  methods::new("saemvsHyperSlab",
    slab_parameter = slab_parameter,
    cov_re_prior_scale = cov_re_prior_scale,
    cov_re_prior_df = cov_re_prior_df,
    residual_variance_prior_shape = 1,
    residual_variance_prior_rate = 1,
    phi_intercept_prior_variance = 3000^2,
    inclusion_prob_prior_a = NULL,
    inclusion_prob_prior_b = NULL
  )
}

#' Class "saemvsHyperSpikeAndSlab"
#'
#' An S4 class that extends \code{\linkS4class{saemvsHyperSlab}} by adding
#' a spike hyperparameter for the spike-and-slab prior.
#'
#' @slot spike_parameter \code{numeric}. The spike parameter, must be strictly
#'  positive.
#'
#' @section Prototype:
#' Defaults to \code{spike_parameter = 0.1}.
#'
#' @section Validity:
#' Ensures that \code{spike_parameter} is strictly positive.
#'
#' @seealso \code{\linkS4class{saemvsHyperSlab}}
#'
#' @keywords internal
setClass(
  "saemvsHyperSpikeAndSlab",
  slots = list(
    spike_parameter = "numericORNULL"
  ),
  prototype = list(
    spike_parameter = 0.1
  ),
  contains = "saemvsHyperSlab",
  validity = function(object) {
    errors <- c()

    errors <- c(
      errors,
      check_positive_or_null_slot(object@spike_parameter, "spike_parameter")
    )

    if (length(errors) == 0) TRUE else errors
  }
)

#' Internal constructor for saemvsHyperSpikeAndSlab
#'
#' @param spike_parameter numeric, strictly positive
#' @param hyper_slab an object of class \code{saemvsHyperSlab}
#' @keywords internal
#' @return An object of class \code{saemvsHyperSpikeAndSlab}
#' @name saemvsHyperSpikeAndSlab
saemvsHyperSpikeAndSlab <- function(spike_parameter,
                                    hyper_slab) {
  methods::new("saemvsHyperSpikeAndSlab",
    spike_parameter = spike_parameter,
    slab_parameter = hyper_slab@slab_parameter,
    cov_re_prior_scale = hyper_slab@cov_re_prior_scale,
    cov_re_prior_df = hyper_slab@cov_re_prior_df,
    residual_variance_prior_shape = hyper_slab@residual_variance_prior_shape,
    residual_variance_prior_rate = hyper_slab@residual_variance_prior_rate,
    phi_intercept_prior_variance = hyper_slab@phi_intercept_prior_variance,
    inclusion_prob_prior_a = hyper_slab@inclusion_prob_prior_a,
    inclusion_prob_prior_b = hyper_slab@inclusion_prob_prior_b
  )
}

#' @title saemvsInit class
#' @description Initialization of user-provided parameters
#'
#' Stores the initial values of parameters provided by the user for the SAEMVS
#'  algorithm.
#'
#' @slot intercept Numeric vector of intercepts for each component of phi.
#'   Must have length equal to \code{phi_dim} in the model.
#' @slot beta_forced Optional matrix of regression coefficients for forced
#'  covariates (from \code{x_forced}). Each column corresponds to a component
#'  of phi.
#' @slot beta_candidates Optional matrix of regression coefficients for
#'  candidate covariates (from \code{x_candidates}). Same structure as
#' \code{beta_forced}.
#' @slot cov_re Initial covariance matrix of the random effects.
#' @slot sigma2 Numeric. Initial residual variance (must be strictly positive).
#' @slot default Logical. If TRUE, all slots except intercept must be "empty".
#' @details
#' Validity checks ensure:
#' \itemize{
#'   \item The number of columns in \code{beta_forced} and
#' \code{beta_candidates} matches the length of \code{intercept}.
#'   \item \code{cov_re} is a valid covariance matrix.
#'   \item \code{sigma2} is strictly positive.
#' }
#'
#' @exportClass saemvsInit
setClass(
  "saemvsInit",
  slots = list(
    intercept = "numeric",
    beta_forced = "matrixORNULL",
    beta_candidates = "matrixORNULL",
    cov_re = "matrix",
    sigma2 = "numeric",
    default = "logical"
  ),
  prototype = list(
    intercept = numeric(0),
    beta_forced = NULL,
    beta_candidates = NULL,
    cov_re = matrix(numeric(0), nrow = 0, ncol = 0),
    sigma2 = 1,
    default = FALSE
  ),
  validity = function(object) {
    n_phi <- length(object@intercept)

    if (n_phi == 0) {
      return("'intercept' must be provided and non-empty.")
    }
    if (object@default) {
      # If default = TRUE, every slot except intercept must be empty or NULL
      if (!is.null(object@beta_forced)) {
        return("'beta_forced' must be NULL when default = TRUE.")
      }
      if (!is.null(object@beta_candidates)) {
        return("'beta_candidates' must be NULL when default = TRUE.")
      }
      # if (length(object@cov_re) > 0) {
      #   return("'cov_re' must be empty when default = TRUE.")
      # }
      # if (is.na(object@sigma2) && object@sigma2 != 1) {
      #   return("'sigma2' must remain at its default when default = TRUE.")
      #   ## A revoir
      # }
    } else {
      ## Check beta_candidates
      if (!is.null(object@beta_candidates)) {
        if (ncol(object@beta_candidates) != n_phi) {
          return("Number of columns in 'beta_candidates' must equal length of 'intercept'.")
        }
      }

      ## Check beta_forced
      if (!is.null(object@beta_forced)) {
        if (ncol(object@beta_forced) != n_phi) {
          return("Number of columns in 'beta_forced' must equal length of 'intercept'.")
        }
      }

      ## Check cov_re (possibly degenerated covariance matrix)
      if (length(object@cov_re) > 0) {
        if (!is.matrix(object@cov_re)) {
          return("'cov_re' must be a matrix.")
        }
        if (nrow(object@cov_re) != ncol(object@cov_re)) {
          return("'cov_re' must be square.")
        }

        diag_vals <- diag(object@cov_re)

        if (any(diag_vals < 0)) {
          return("Diagonal entries of 'cov_re' must be non-negative.")
        }

        ## Random components with positive variance
        idx_active <- which(diag_vals > 0)
        if (length(idx_active) > 0) {
          submat <- object@cov_re[idx_active, idx_active, drop = FALSE]

          msg <- check_covariance(submat, "cov_re (reduced)")
          if (!is.null(msg)) {
            return(msg)
          }
        }

        ## Check that if diag=0, raw and column are all zero too
        idx_zero <- which(diag_vals == 0)
        if (length(idx_zero) > 0) {
          for (i in idx_zero) {
            if (any(object@cov_re[i, ] != 0) || any(object@cov_re[, i] != 0)) {
              return(paste0("Row/column ", i, " of 'cov_re' must be all zeros if diagonal entry is zero."))
            }
          }
        }
      }

      ## Check sigma2
      if (length(object@sigma2) != 1 || object@sigma2 <= 0) {
        return("'sigma2' must be a single strictly positive value.")
      }
    }

    TRUE
  }
)


#' @rdname saemvsInit
#' @title Constructor for saemvsInit
#'
#' @description Create an object of class \code{saemvsInit} specifying initial
#'  values for the SAEMVS algorithm.
#'
#' @param intercept Numeric vector. Intercepts for each component of phi.
#' @param beta_forced Numeric matrix or NULL. Coefficients for forced
#'  covariates.
#' @param beta_candidates Numeric matrix or NULL. Coefficients for candidate
#'  covariates.
#' @param cov_re Numeric square matrix. Initial covariance of random effects.
#' @param sigma2 Numeric, positive. Initial residual variance (default = 1).
#' @param default Logical. If TRUE, ignore all slots except intercept.
#' @return An object of class \code{saemvsInit}.
#' @export
saemvsInit <- function(intercept,
                       beta_forced = NULL,
                       beta_candidates = NULL,
                       cov_re = matrix(numeric(0), nrow = 0, ncol = 0),
                       sigma2 = 1,
                       default = FALSE) {
  methods::new("saemvsInit",
    intercept = intercept,
    beta_forced = if (default) NULL else beta_forced,
    beta_candidates = if (default) NULL else beta_candidates,
    cov_re = cov_re,
    sigma2 = sigma2,
    default = default
  )
}

#' Internal classe saemvsProcessedInit
#'
#' Internal processed initialization of parameters
#'
#' Stores processed initial values of parameters transformed into a format
#'  compatible with the SAEMVS algorithm.
#'
#' @slot beta_to_select Matrix of regression coefficients for phi components
#'  subject to variable selection.
#' @slot beta_not_to_select Matrix of regression coefficients for phi
#'  components not subject to variable selection.
#' @slot gamma_to_select Covariance_matrix associated with phi_to_select.
#' @slot gamma_not_to_select Covariance_matrix associated with
#'  phi_not_to_select.
#' @slot sigma2 Numeric. Residual variance.
#' @slot inclusion_prob Numeric vector of inclusion probabilities, consistent
#'  with the processed dimension of phi.
#'
#' @note This class is created internally by the package. Users should not
#'  create or modify objects of this class directly.
#'
#' @keywords internal
setClass(
  "saemvsProcessedInit",
  slots = list(
    beta_to_select = "matrixORNULL",
    beta_not_to_select = "matrixORNULL",
    gamma_to_select = "matrixORNULL",
    gamma_not_to_select = "matrixORNULL",
    sigma2 = "numeric",
    inclusion_prob = "numericORNULL"
  ),
  prototype = list(
    beta_to_select = NULL,
    beta_not_to_select = NULL,
    gamma_to_select = NULL,
    gamma_not_to_select = NULL,
    sigma2 = 1,
    inclusion_prob = numeric(0)
  ),
  contains = "saemvsInit"
  # No validity: objects created automatically in the code
)

#' Internal constructor for saemvsProcessedInit
#' @keywords internal
#' @param beta_to_select Numeric matrix or NULL. Coefficients for phi
#'  components subject to selection.
#' @param beta_not_to_select Numeric matrix or NULL. Coefficients for phi
#'  components not subject to selection.
#' @param gamma_to_select Numeric matrix or NULL. Covariance matrix for phi
#'  components to select.
#' @param gamma_not_to_select Numeric matrix or NULL. Covariance matrix for phi
#'  components not to select.
#' @param sigma2 Numeric. Residual variance.
#' @param inclusion_prob Numeric vector or NULL. Inclusion probabilities for
#'  phi components.
#' @return An object of class \code{saemvsProcessedInit}.
saemvsProcessedInit <- function(beta_to_select = NULL,
                                beta_not_to_select = NULL,
                                gamma_to_select = NULL,
                                gamma_not_to_select = NULL,
                                sigma2 = 1,
                                inclusion_prob = numeric(0),
                                intercept,
                                beta_forced,
                                beta_candidates,
                                cov_re,
                                default) {
  methods::new("saemvsProcessedInit",
    beta_to_select = beta_to_select,
    beta_not_to_select = beta_not_to_select,
    gamma_to_select = gamma_to_select,
    gamma_not_to_select = gamma_not_to_select,
    sigma2 = sigma2,
    inclusion_prob = inclusion_prob,
    intercept = intercept,
    beta_forced = beta_forced,
    beta_candidates = beta_candidates,
    cov_re = cov_re,
    default = default
  )
}

#' @title saemvsTuning class
#' @description Algorithm tuning parameters for SAEMVS
#'
#' This class contains all parameters controlling the behavior of the SAEMVS
#'  algorithm, including iteration numbers, Metropolis-Hastings steps,
#'  importance sampling, step-size schedule, and the grid of spike values for
#'  variable selection.
#'
#' @slot niter Total number of SAEM iterations.
#' @slot nburnin Number of burn-in iterations.
#' @slot step Step-size schedule vector for stochastic approximation
#'  (automatically computed).
#' @slot niter_mh Number of Metropolis-Hastings iterations per S-step.
#' @slot kernel_mh Type of MH kernel: either "random_walk" or "pop".
#' @slot covariance_decay Decay factor controlling how quickly variance and
#'  covariance estimates are updated during the SAEMVS iterations. Values close
#'  to 1 lead to slow updates, ensuring stability of the Metropolis-Hastings
#'  kernel, while lower values increase the adaptation rate.
#' @slot mh_proposal_scale Scaling factor for the MH proposal variance
#'  (strictly positive).
#' @slot spike_values_grid Numeric vector of spike values (must be positive and
#'  non-empty).
#' @slot n_is_samples Number of importance sampling iterations for
#'  log-likelihood estimation.
#' @slot seed Integer random seed for reproducibility.
#' @slot nb_workers Number of parallel workers for SAEMVS.
#'
#' @exportClass saemvsTuning
setClass(
  "saemvsTuning",
  slots = list(
    niter = "numeric",
    nburnin = "numeric",
    step = "numeric",
    niter_mh = "numeric",
    kernel_mh = "character",
    covariance_decay = "numeric",
    mh_proposal_scale = "numeric",
    spike_values_grid = "numericORNULL",
    n_is_samples = "numeric",
    seed = "numeric",
    nb_workers = "numeric"
  ),
  prototype = list(
    niter = 500,
    nburnin = 350,
    step = numeric(0), # computed automatically in constructor
    niter_mh = 5,
    kernel_mh = "random_walk",
    covariance_decay = 0.98,
    mh_proposal_scale = 1.5,
    spike_values_grid = NULL,
    n_is_samples = 10000,
    seed = 220916,
    nb_workers = 4
  ),
  validity = function(object) {
    # niter
    if (!is.numeric(object@niter) || length(object@niter) != 1) {
      return("'niter' must be a single numeric value.")
    }
    if (object@niter <= 0 || object@niter != as.integer(object@niter)) {
      return("'niter' must be a positive integer.")
    }

    # nburnin
    if (!is.numeric(object@nburnin) || length(object@nburnin) != 1) {
      return("'nburnin' must be a single numeric value.")
    }
    if (object@nburnin < 0 || object@nburnin != as.integer(object@nburnin)) {
      return("'nburnin' must be a non-negative integer.")
    }
    if (object@nburnin > object@niter) {
      return("'nburnin' must be smaller than 'niter'.")
    }

    # step
    if (!is.numeric(object@step) || length(object@step) != object@niter) {
      return("'step' must be a numeric vector of length 'niter'.")
    }

    # niter_mh
    if (!is.numeric(object@niter_mh) || length(object@niter_mh) != 1) {
      return("'niter_mh' must be a single numeric value.")
    }
    if (object@niter_mh <= 0 ||
      object@niter_mh != as.integer(object@niter_mh)) {
      return("'niter_mh' must be a positive integer.")
    }

    # kernel_mh
    if (!is.character(object@kernel_mh) || length(object@kernel_mh) != 1) {
      return("'kernel_mh' must be a single character string.")
    }
    if (!object@kernel_mh %in% c("random_walk", "pop")) {
      return("'kernel_mh' must be either 'random_walk' or 'pop'.")
    }

    # covariance_decay
    if (!is.numeric(object@covariance_decay) ||
      length(object@covariance_decay) != 1) {
      return("'covariance_decay' must be a single numeric value.")
    }
    if (object@covariance_decay <= 0 || object@covariance_decay >= 1) {
      return("'covariance_decay' must be strictly between 0 and 1 (exclusive).")
    }

    # mh_proposal_scale
    if (!is.numeric(object@mh_proposal_scale) ||
      length(object@mh_proposal_scale) != 1) {
      return("'mh_proposal_scale' must be a single numeric value.")
    }
    if (object@mh_proposal_scale <= 0) {
      return("'mh_proposal_scale' must be strictly positive.")
    }

    # spike_values_grid
    if (!is.numeric(object@spike_values_grid) ||
      length(object@spike_values_grid) == 0) {
      return("'spike_values_grid' must be a non-empty numeric vector.")
    }
    if (any(object@spike_values_grid <= 0)) {
      return("All values in 'spike_values_grid' must be strictly positive.")
    }

    # n_is_samples
    if (!is.numeric(object@n_is_samples) || length(object@n_is_samples) != 1) {
      return("'n_is_samples' must be a single numeric value.")
    }
    if (object@n_is_samples <= 0 || object@n_is_samples != as.integer(object@n_is_samples)) {
      return("'n_is_samples' must be a positive integer.")
    }

    # seed
    if (!is.numeric(object@seed) || length(object@seed) != 1) {
      return("'seed' must be a single numeric value.")
    }
    if (object@seed != as.integer(object@seed)) {
      return("'seed' must be an integer.")
    }

    # nb_workers
    if (!is.numeric(object@nb_workers) || length(object@nb_workers) != 1) {
      return("'nb_workers' must be a single numeric value.")
    }
    if (object@nb_workers <= 0 ||
      object@nb_workers != as.integer(object@nb_workers)) {
      return("'nb_workers' must be a positive integer.")
    }

    TRUE
  }
)

#' @rdname saemvsTuning
#' @title Constructor for saemvsTuning
#'
#' @description Create a \code{saemvsTuning} object specifying algorithm tuning
#'  parameters for SAEMVS.
#'
#' @param niter Integer. Total number of SAEM iterations (positive).
#' @param nburnin Integer. Number of burn-in iterations (non-negative, less
#'  than niter).
#' @param niter_mh Integer. Number of Metropolis-Hastings iterations per S-step
#'  (positive).
#' @param kernel_mh Character. Type of MH kernel: "random_walk" or "pop".
#' @param covariance_decay Numeric between 0 and 1. Decay factor for covariance
#'  adaptation.
#' @param mh_proposal_scale Positive numeric. Scaling factor for MH proposal
#'  variance.
#' @param spike_values_grid Numeric vector, strictly positive and non-empty.
#' @param n_is_samples Integer. Number of importance sampling iterations
#'  (positive).
#' @param seed Integer. Random seed for reproducibility.
#' @param nb_workers Integer. Number of parallel workers (positive).
#'
#' @return An object of class \code{saemvsTuning}.
#' @export
saemvsTuning <- function(niter = 500,
                         nburnin = 350,
                         niter_mh = 5,
                         kernel_mh = "random_walk",
                         covariance_decay = 0.98,
                         mh_proposal_scale = 1.0,
                         spike_values_grid,
                         n_is_samples = 10000,
                         seed = 220916,
                         nb_workers = 4) {
  step <- c(
    rep(1, nburnin - 1),
    1 / ((1:(niter - nburnin + 1))^(2 / 3))
  )

  methods::new("saemvsTuning",
    niter = as.integer(niter),
    nburnin = as.integer(nburnin),
    step = step,
    niter_mh = as.integer(niter_mh),
    kernel_mh = kernel_mh,
    covariance_decay = covariance_decay,
    mh_proposal_scale = mh_proposal_scale,
    spike_values_grid = sort(spike_values_grid),
    n_is_samples = as.integer(n_is_samples),
    seed = as.integer(seed),
    nb_workers = as.integer(nb_workers)
  )
}


#' Class saemvsResults
#'
#' Results of SAEMVS model selection.
#'
#' Internal results container storing all outputs from a SAEMVS grid search and
#' model selection procedure. Objects of this class are created automatically
#' by the SAEMVS algorithm and are not intended to be constructed manually.
#'
#' @slot criterion Character scalar. The selection criterion used to identify
#'   the best model. Supported values: `"BIC"` or `"e-BIC"`.
#'
#' @slot criterion_values Numeric vector. Values of the selection criterion
#'   computed for each unique support in \code{unique_support}.
#'
#' @slot thresholds List. Thresholds applied to beta coefficients for each
#'   spike parameter value in \code{spike_values_grid}.
#'
#' @slot beta_map List. MAP (maximum a posteriori) estimates of beta
#'   coefficients for each spike value.
#'
#' @slot mle_estimates List. Maximum likelihood estimates of model parameters
#'   (beta and gamma) computed on each unique support.
#'
#' @slot support List. Raw supports selected across the spike grid
#'   (may contain duplicates).
#'
#' @slot unique_support List. Unique supports obtained from \code{support}
#'   (duplicates removed). Each entry corresponds to one element of
#'   \code{criterion_values}.
#'
#' @slot support_mapping Numeric vector. For each element of \code{support},
#'   gives the index of the corresponding entry in \code{unique_support}.
#'
#' @slot spike_values_grid Numeric vector. The grid of spike parameter values
#'   explored during the SAEMVS run.
#'
#' @slot phi_fixed_idx Numeric vector or \code{NULL}. Indices of fixed
#'   \eqn{\phi} components (non-estimated parameters), used for display in
#'   summaries.
#'
#' @slot forced_variables_idx List. For each unique support, indices of
#'   covariates that were forced into the model (always included).
#'
#' @slot selected_variables_idx List. For each unique support, indices of
#'   covariates that were actively selected by the variable selection
#'   procedure (excluding forced covariates).
#'
#' @details
#' The class is mainly used as an internal result structure and is returned by
#' functions implementing SAEMVS. Downstream methods such as
#' \code{\link{summary_saemvs}} make use of these slots for reporting.
#'
#' @keywords internal
setClass(
  "saemvsResults",
  slots = list(
    criterion = "character",
    criterion_values = "numeric",
    thresholds = "list",
    beta_map = "list",
    mle_estimates = "list",
    support = "list",
    unique_support = "list",
    support_mapping = "numeric",
    spike_values_grid = "numeric",
    phi_fixed_idx = "numericORNULL",
    forced_variables_idx = "list",
    selected_variables_idx = "list"
  ),
  prototype = list(
    criterion = character(0),
    criterion_values = numeric(0),
    thresholds = list(),
    beta_map = list(),
    mle_estimates = list(),
    support = list(),
    unique_support = list(),
    support_mapping = numeric(0),
    spike_values_grid = numeric(0),
    phi_fixed_idx = numeric(0),
    forced_variables_idx = list(),
    selected_variables_idx = list()
  ),
  validity = function(object) {
    if (!object@criterion %in% c("BIC", "e-BIC")) {
      return("Slot 'criterion' must be either 'BIC' or 'e-BIC'.")
    }

    if (length(object@criterion_values) != length(object@unique_support)) {
      return("The number of 'criterion_values' must match the number of
       'unique_support' entries.")
    }


    if (length(object@thresholds) != length(object@spike_values_grid)) {
      return("Slot 'thresholds' must have same length as 'spike_values_grid'.")
    }

    if (length(object@beta_map) != length(object@spike_values_grid)) {
      return("Slot 'beta_map' must have same length as 'spike_values_grid'.")
    }

    if (length(object@support_mapping) != length(object@support)) {
      return("Slot 'support_mapping' must map each support to a unique index.")
    }

    TRUE
  }
)




#' Classe saemResults
#'
#' Internal class storing intermediate estimation results produced by
#' the SAEM algorithm. Objects of this class are not meant
#' to be created by the user but are automatically produced by
#' \code{saemvs}.
#'
#' @slot beta_to_select A list of matrices containing the estimated regression
#'   coefficients for covariates that are related with the components of
#' `phi` that are subject to selection (\code{phi_to_select}) at each iteration.
#' @slot beta_not_to_select A list of matrices containing the estimated
#'  regression coefficients for covariates that are related with the components
#'  of `phi` that are not subject to selection (\code{phi_not_to_select}) at
#'  each iteration.
#' @slot gamma_to_select A list of covariance matrices estimated at each
#' iteration for the random effects associated with \code{phi_to_select}.
#' @slot gamma_not_to_select A list of covariance matrices estimated at each
#'  iteration for the random effects associated with \code{phi_not_to_select}.
#' @slot sigma2 A numeric vector of estimated residual variances at each
#'  iteration.
#'
#' @keywords internal
setClass(
  "saemResults",
  slots = list(
    beta_to_select = "listORNULL",
    beta_not_to_select = "listORNULL",
    gamma_to_select = "listORNULL",
    gamma_not_to_select = "listORNULL",
    sigma2 = "numericORNULL"
  ),
  prototype = list(
    beta_to_select = NULL,
    beta_not_to_select = NULL,
    gamma_to_select = NULL,
    gamma_not_to_select = NULL,
    sigma2 = NULL
  ),
  validity = function(object) {
    # Check beta_to_select
    if (!is.null(object@beta_to_select)) {
      if (!is.list(object@beta_to_select) ||
        !all(vapply(object@beta_to_select, is.matrix, logical(1)))) {
        return("'beta_to_select' must be a list of matrices.")
      }
    }

    # Check beta_not_to_select
    if (!is.null(object@beta_not_to_select)) {
      if (!is.list(object@beta_not_to_select) ||
        !all(vapply(object@beta_not_to_select, is.matrix, logical(1)))) {
        return("'beta_not_to_select' must be a list of matrices.")
      }
    }

    # Check gamma_to_select
    if (!is.null(object@gamma_to_select)) {
      if (!is.list(object@gamma_to_select) ||
        !all(vapply(object@gamma_to_select, is.matrix, logical(1)))) {
        return("'gamma_to_select' must be a list of matrices.")
      }
    }

    # Check gamma_not_to_select
    if (!is.null(object@gamma_not_to_select)) {
      if (!is.list(object@gamma_not_to_select) ||
        !all(vapply(object@gamma_not_to_select, is.matrix, logical(1)))) {
        return("'gamma_not_to_select' must be a list of matrices.")
      }
    }

    # Check sigma2
    if (!is.null(object@sigma2)) {
      if (!is.numeric(object@sigma2) || any(object@sigma2 <= 0)) {
        return("'sigma2' must be a positive numeric vector.")
      }
    }

    TRUE
  }
)
