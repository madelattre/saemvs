setClassUnion("matrixORNULL", c("matrix", "NULL"))
setClassUnion("numericORNULL", c("numeric", "NULL"))
setClassUnion("listORNULL", c("list", "NULL"))

#' Class saemvsData
#'
#' A class to store user-provided data for the SAEMVS algorithm.
#' This is the main entry point for users: responses, time indices,
#' candidate covariates for selection, and forced covariates.
#'
#' @slot y_series list of numeric vectors. Each element contains observed responses for one individual.
#' @slot t_series list of numeric vectors. Each element contains corresponding time indices.
#' @slot x_candidates matrix or NULL. Optional matrix of candidate covariates (rows = individuals, columns = covariates).
#' @slot x_forced matrix or NULL. Optional matrix of forced covariates (rows = individuals, columns = covariates).
#'
#' @section Validity:
#' \itemize{
#'   \item \code{y_series} and \code{t_series} must have the same length.
#'   \item Each element of \code{y_series} and \code{t_series} must be numeric vectors of the same length.
#'   \item If provided, \code{x_candidates} and \code{x_forced} must each have as many rows as the length of \code{y_series}.
#'   \item All matrices must be numeric.
#' }
#'
#' @return An object of class \code{saemvsData}.
#' @exportClass saemvsData
#'
#' @examples
#' y <- list(rnorm(10), rnorm(8))
#' t <- list(1:10, 1:8)
#' d <- saemvsData(y, t)
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
        return(sprintf("y_series[[%d]] and t_series[[%d]] do not have the same length.", i, i))
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

    # # Warning if overlap between to_select and not_to_select
    # if (!is.null(object@x_phi_to_select) &&
    # !is.null(object@x_phi_not_to_select)) {
    #   common_cols <- intersect(
    #     colnames(object@x_phi_to_select), colnames(object@x_phi_not_to_select)
    #     )
    #   if (length(common_cols) > 0) {
    #     warning("Some columns appear both in 'x_phi_to_select' and 'x_phi_not_to_select'.")
    #   }
    # }

    TRUE
  }
)



#' Constructor for saemvsData
#'
#' @param y list of numeric vectors (responses).
#' @param t list of numeric vectors (time indices).
#' @param x_candidates matrix or NULL, candidate covariates for selection.
#' @param x_forced matrix or NULL, forced covariates.
#'
#' @return An object of class \code{saemvsData}.
#' @export
saemvsData <- function(y, t, x_candidates = NULL, x_forced = NULL) {
  new("saemvsData",
    y_series     = y,
    t_series     = t,
    x_candidates = x_candidates,
    x_forced     = x_forced
  )
}


# Internal class: saemvsProcessedData
#' Internal class: saemvsProcessedData
#'
#' A class to store processed design matrices ready for SAEMVS algorithms.
#' This class is internal to the package and is not intended to be
#' directly created or manipulated by the user.
#' It inherits from \code{saemvsData} and contains additional slots
#' representing design matrices of covariates associated with the model parameters \code{phi},
#' split into components subject to selection and components that are not subject to selection.
#'
#' @slot x_phi_to_select matrix or NULL.
#'   Design matrix of covariates for \code{phi} parameters on which variable
#'   selection will be performed. Rows correspond to individuals; columns to covariates.
#'
#' @slot x_phi_not_to_select matrix or NULL.
#'   Design matrix of covariates for \code{phi} parameters that are not
#'   subject to selection. Rows correspond to individuals; columns to covariates.
#'
#' @slot tx_x_phi_to_select matrix or NULL.
#'   Transformation of \code{x_phi_to_select} (e.g., multiplied by its transpose) required for algorithm computations.
#'
#' @slot kron_tx_x_phi_to_select matrix or NULL.
#'   Kronecker product of \code{tx_x_phi_to_select} used for vectorized calculations.
#'
#' @slot x_phi_not_to_select_list list or NULL.
#'   List of \code{x_phi_not_to_select} entries per individual. Each element
#'   corresponds to one individual and contains a numeric vector.
#'
#' @section Validity:
#'   \itemize{
#'     \item All matrices must be numeric if not NULL.
#'     \item \code{x_phi_not_to_select_list}, if not NULL, must have length equal to the number of individuals.
#'     \item Columns in \code{x_phi_to_select} and \code{x_phi_not_to_select} should not overlap (warning issued if they do).
#'     \item Inherits all validity checks from \code{saemvsData} (e.g., length and type checks for y_series, t_series, x_candidates, x_forced).
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

# Internal constructor
saemvsProcessedData <- function(x_phi_to_select = NULL,
                                x_phi_not_to_select = NULL,
                                tx_x_phi_to_select = NULL,
                                kron_tx_x_phi_to_select = NULL,
                                x_phi_not_to_select_list = NULL,
                                ...) {
  new("saemvsProcessedData",
    x_phi_to_select          = x_phi_to_select,
    x_phi_not_to_select      = x_phi_not_to_select,
    tx_x_phi_to_select       = tx_x_phi_to_select,
    kron_tx_x_phi_to_select  = kron_tx_x_phi_to_select,
    x_phi_not_to_select_list = x_phi_not_to_select_list,
    ...
  )
}


#' Class saemvsModel
#'
#' Represents a model for SAEMVS, including the model function and
#' indexing information for variable selection on parameters phi.
#'
#' @slot model_func A function of form function(phi, t) returning predicted values.
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
      all(object@phi_to_select_idx >= 1 & object@phi_to_select_idx <= object@phi_dim))) {
      return("'phi_to_select_idx' must contain integers between 1 and 'phi_dim'.")
    }

    if (!(is.null(object@phi_fixed_idx) ||
      all(object@phi_fixed_idx >= 1 & object@phi_fixed_idx <= object@phi_dim))) {
      return("'phi_not_to_select_idx' must contain integers between 1 and 'phi_dim'.")
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
        return(paste0("'x_forced_support' must have ", object@phi_dim, " columns (one per phi)."))
      }
    }

    TRUE
  }
)

# Constructor
#' @export
saemvsModel <- function(
    g, phi_dim, phi_to_select_idx = c(), phi_fixed_idx = c(),
    x_forced_support = matrix(numeric(0), nrow = 0, ncol = 0)) {
  new("saemvsModel",
    model_func = g,
    phi_dim = as.integer(phi_dim),
    phi_to_select_idx = as.integer(phi_to_select_idx),
    phi_fixed_idx = as.integer(phi_fixed_idx),
    x_forced_support = x_forced_support
  )
}


#' Class saemvsHyperSlab
#'
#' Represents the hyperparameters for the slab component in a spike-and-slab prior.
#'
#' Only the slab parameter, random effects covariance prior, and degrees of freedom need
#' to be specified by the user. Other hyperparameters are fixed internally.
#'
#' @slot slab_parameter Numeric: positive slab parameter (controls the slab prior).
#' @slot cov_re_prior_scale Numeric matrix: scale matrix of the Inverse-Wishart prior for random effects covariance.
#' @slot cov_re_prior_df Numeric: degrees of freedom of the Inverse-Wishart prior for random effects covariance.
#' @slot residual_variance_prior_shape Numeric: shape parameter of the inverse-gamma prior on residual variance (fixed internally).
#' @slot residual_variance_prior_rate Numeric: rate parameter of the inverse-gamma prior on residual variance (fixed internally).
#' @slot phi_intercept_prior_variance Numeric: variance of the Gaussian prior on each intercept of phi (fixed internally).
#' @slot inclusion_prob_prior_a Numeric or NULL: alpha parameter of the Beta prior on covariate inclusion probabilities.
#' @slot inclusion_prob_prior_b Numeric or NULL: beta parameter of the Beta prior on covariate inclusion probabilities.
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
    check_positive_or_null_slot(object@slab_parameter, "slab_parameter")
    check_positive_or_null_slot(object@cov_re_prior_df, "cov_re_prior_df")
    check_positive_or_null_slot(
      object@residual_variance_prior_shape,
      "residual_variance_prior_shape"
    )
    check_positive_or_null_slot(
      object@residual_variance_prior_rate,
      "residual_variance_prior_rate"
    )
    check_positive_or_null_slot(
      object@phi_intercept_prior_variance,
      "phi_intercept_prior_variance"
    )
    if (!is.null(object@cov_re_prior_scale)) {
      check_covariance(object@cov_re_prior_scale, "cov_re_prior_scale")
    }
    TRUE
  }
)

# Constructor
#' @export
saemvsHyperSlab <- function(slab_parameter = 12000,
                            cov_re_prior_scale,
                            cov_re_prior_df = 1) {
  new("saemvsHyperSlab",
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
#' @slot spike_parameter \code{numeric}. The spike parameter, must be strictly positive.
#'
#' @section Prototype:
#' Defaults to \code{spike_parameter = 0.1}.
#'
#' @section Validity:
#' Ensures that \code{spike_parameter} is strictly positive.
#'
#' @seealso \code{\linkS4class{saemvsHyperSlab}}
#'
#' @exportClass saemvsHyperSpikeAndSlab
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
    check_positive_or_null_slot(object@spike_parameter, "spike_parameter")
    TRUE
  }
)

# Constructor
#' @export
saemvsHyperSpikeAndSlab <- function(spike_parameter,
                                    hyper_slab) {
  new("saemvsHyperSpikeAndSlab",
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



#' Class saemvsInit
#'
#' Initialization of user-provided parameters
#'
#' Stores the initial values of parameters provided by the user for the SAEMVS algorithm.
#'
#' @slot intercept Numeric vector of intercepts for each component of phi.
#'   Must have length equal to \code{phi_dim} in the model.
#' @slot beta_forced Optional matrix of regression coefficients for forced covariates
#'   (from \code{x_forced}). Each column corresponds to a component of phi.
#' @slot beta_candidates Optional matrix of regression coefficients for candidate covariates
#'   (from \code{x_candidates}). Same structure as \code{beta_forced}.
#' @slot cov_re Initial covariance matrix of the random effects.
#' @slot sigma2 Numeric. Initial residual variance (must be strictly positive).
#'
#' @details
#' Validity checks ensure:
#' \itemize{
#'   \item The number of columns in \code{beta_forced} and \code{beta_candidates}
#'         matches the length of \code{intercept}.
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
    sigma2 = "numeric"
  ),
  prototype = list(
    intercept = numeric(0),
    beta_forced = NULL,
    beta_candidates = NULL,
    cov_re = matrix(numeric(0), nrow = 0, ncol = 0),
    sigma2 = 1
  ),
  validity = function(object) {
    n_phi <- length(object@intercept)
    if (!is.null(object@beta_candidates)) {
      if (ncol(object@beta_candidates) != n_phi) {
        return("Number of columns in 'beta_candidates' must equal length of 'intercept'.")
      }
    }
    if (!is.null(object@beta_forced)) {
      if (ncol(object@beta_forced) != n_phi) {
        return("Number of columns in 'beta_forced' must equal length of 'intercept'.")
      }
    }
    check_covariance(object@cov_re, "cov_re")
    if (length(object@sigma2) != 1 || object@sigma2 <= 0) {
      return("'sigma2' must be a single strictly positive value.")
    }
    TRUE
  }
)

# Constructor
#' @export
saemvsInit <- function(intercept,
                       beta_forced = NULL,
                       beta_candidates = NULL,
                       cov_re,
                       sigma2 = 1) {
  new("saemvsInit",
    intercept = intercept,
    beta_forced = beta_forced,
    beta_candidates = beta_candidates,
    cov_re = cov_re,
    sigma2 = sigma2
  )
}


#' Internal processed initialization of parameters
#'
#' Stores processed initial values of parameters transformed into a format compatible with the SAEMVS algorithm.
#'
#' @slot beta_to_select Matrix of regression coefficients for phi components subject to variable selection.
#' @slot beta_not_to_select Matrix of regression coefficients for phi components not subject to variable selection.
#' @slot gamma_to_select Covariance_matrix associated with phi_to_select.
#' @slot gamma_not_to_select Covariance_matrix associated with phi_not_to_select.
#' @slot sigma2 Numeric. Residual variance.
#' @slot inclusion_prob Numeric vector of inclusion probabilities, consistent with the processed dimension of phi.
#'
#' @note This class is created internally by the package. Users should not create or modify objects of this class directly.
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
  )
  # No validity: objects created automatically in the code
)

# Constructor
saemvsProcessedInit <- function(beta_to_select = NULL,
                                beta_not_to_select = NULL,
                                gamma_to_select = NULL,
                                gamma_not_to_select = NULL,
                                sigma2 = 1,
                                inclusion_prob = numeric(0)) {
  new("saemvsProcessedInit",
      beta_to_select = beta_to_select,
      beta_not_to_select = beta_not_to_select,
      gamma_to_select = gamma_to_select,
      gamma_not_to_select = gamma_not_to_select,
      sigma2 = sigma2,
      inclusion_prob = inclusion_prob
  )
}


## -- Tuning parameters

#' @exportClass tuningC
setClass( # ne contient pas les paramètres initiaux ; nu0_grid est obligatoire
  "tuningC",
  slots = list(
    niter = "numeric",
    nburnin = "numeric",
    step = "numeric",
    niter_mh = "numeric",
    kernel_mh = "character",
    tau = "numeric",
    kappa_mh = "numeric",
    nu0_grid = "numericORNULL",
    nb_is = "numeric",
    seed = "numeric",
    nb_workers = "numeric"
  ),
  prototype = list(
    niter = integer(500),
    nburnin = integer(350),
    # step = c(rep(1, 350 - 1), 1 / ((1:(151))^(2 / 3))),
    niter_mh = integer(5),
    kernel_mh = "random_walk",
    # tau = 0.98,
    kappa_mh = 1.5,
    nu0_grid = NULL,
    nb_is = 10000,
    seed = 220916,
    nb_workers = 4
  ),
  validity = function(object) {
    if (object@kappa_mh <= 0) {
      return(
        paste0(
          "Parameter 'kappa_mh' used to tune the proposal ",
          "variance at S-step must of the SAEM algorithm be a ",
          "strictly positive value."
        )
      )
    }

    if (!object@kernel_mh %in% c("random_walk", "pop")) {
      return("Invalid kernel: should be 'random_walk' or 'pop'.")
    }

    check_positive_integer_slot(
      object@niter,
      "The total number of iterations of the SAEM algorithm 'niter' "
    )

    check_positive_integer_slot(
      object@nburnin,
      "The number of burn-in iterations of the SAEM algorithm 'nburnin' "
    )

    if (object@nburnin > object@niter) {
      return(
        paste0(
          "The number of burn-in iterations of the SAEM algorithm",
          " 'nburnin' must be smaller than the total number of iterations",
          " 'niter'."
        )
      )
    }

    check_positive_integer_slot(
      object@niter_mh,
      paste0(
        "The number of iterations of Metropolis-Hastings at each",
        " S-step 'niter_mh' "
      )
    )

    check_positive_integer_slot(
      object@nb_is,
      paste0(
        "The number of iterations of Importance Sampling 'nb_is'",
        " for log-likelihood estimation "
      )
    )

    if (!is.integer(object@seed)) {
      return(
        paste0("'seed' must be an integer value.")
      )
    }

    check_positive_integer_slot(
      object@nb_workers,
      paste0(
        "The number of workers 'nb_workers' for parallel processing",
        " of SAEMVS on 'nu0_grid' "
      )
    )


    if (!all(object@nu0_grid > 0)) {
      return(
        paste0(
          "All the spike parameter values in 'nu0_grid' must be positive."
        )
      )
    }

    # Ajouter que nu0_grid ne peut pas être vide

    TRUE
  }
)

#' @export
tuningC <- function(
    niter = 500,
    nburnin = 350,
    niter_mh = 5,
    kernel_mh = "random_walk",
    kappa_mh = 1.5,
    nu0_grid,
    nb_is = 10000,
    seed = 220916,
    nb_workers = 4) {
  new("tuningC",
    niter = as.integer(niter), nburnin = as.integer(nburnin),
    step = c(rep(1, nburnin - 1), 1 / ((1:(niter - nburnin + 1))^(2 / 3))),
    niter_mh = as.integer(niter_mh), kernel_mh = kernel_mh,
    tau = 0.98, kappa_mh = kappa_mh, nu0_grid = sort(nu0_grid),
    nb_is = as.integer(nb_is), seed = as.integer(seed),
    nb_workers = as.integer(nb_workers)
  )
}


## -- Results

#' @exportClass resSAEMVS
setClass(
  "resSAEMVS",
  slots = list(
    pen = "character",
    crit_values = "numeric",
    thresholds = "list",
    beta_map = "list",
    est_mle = "list",
    support = "list",
    unique_support = "list",
    map_to_unique_support = "numeric",
    nu0_grid = "numeric",
    index_fixed = "numericORNULL"
  ),
  prototype = list(
    pen = character(0),
    crit_values = numeric(0),
    thresholds = list(),
    beta_map = list(),
    est_mle = list(),
    support = list(),
    unique_support = list(),
    map_to_unique_support = numeric(0),
    nu0_grid = numeric(0),
    index_fixed = numeric(0)
  ),
  validity = function(object) {
    TRUE
  }
)


#' @exportClass resSAEM
setClass(
  "resSAEM",
  slots = list(
    beta_s = "listORNULL",
    beta_ns = "listORNULL",
    gamma_s = "listORNULL",
    gamma_ns = "listORNULL",
    sigma2 = "numericORNULL"
  ),
  prototype = list(
    beta_s = NULL,
    beta_ns = NULL,
    gamma_s = NULL,
    gamma_ns = NULL,
    sigma2 = NULL
  ),
  validity = function(object) {
    TRUE
  }
)
