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


## -- Initilialization of unknown parameters

#' @exportClass initC
setClass(
  "initC",
  slots = list(
    intercept = "numeric",
    beta_forced = "matrixORNULL",
    beta_sel = "matrixORNULL",
    gamma = "matrix",
    sigma2 = "numeric",
    alpha = "numericORNULL"
  ),
  prototype = list(
    intercept = numeric(0),
    beta_forced = NULL,
    beta_sel = NULL,
    gamma = matrix(numeric(0), nrow = 0, ncol = 0),
    sigma2 = numeric(0),
    alpha = NULL
  ),
  validity = function(object) {
    # Même nombre de composantes dans intercept que de colonnes dans beta_forced et beta_sel si ces deux derniers sont non nuls

    if (!is.null(object@beta_sel)) {
      if (ncol(object@beta_sel) != length(object@intercept)) {
        return(
          paste0(
            "The number of colums in 'beta_sel' must be equal to ",
            "the number of components in 'intercept'."
          )
        )
      }
    }


    if (!is.null(object@beta_forced)) {
      if (ncol(object@beta_forced) != length(object@intercept)) {
        return(
          paste0(
            "The number of colums in 'beta_forced' must be equal to ",
            "the number of components in 'intercept'."
          )
        )
      }
    }

    if (!is.null(object@gamma)) {
      check_covariance(object@gamma, "The initial value for 'gamma' ")
    }

    check_positive_slot(object@sigma2, "The initial value for 'sigma2' ")

    if ((!is.null(object@beta_sel)) && (!is.null(object@alpha))) {
      if (ncol(object@beta_sel) != length(object@alpha)) {
        return(
          paste0(
            "The number of elements in 'alpha' must be equal to the ",
            "number of columns in 'beta_sel'."
          )
        )
      }
    }

    TRUE
  }
)

#' @export
initC <- function(
    intercept, beta_forced = NULL, beta_sel = NULL,
    gamma, sigma2, alpha = NULL) {
  new("initC",
    intercept = intercept,
    beta_forced = beta_forced,
    beta_sel = beta_sel,
    gamma = gamma,
    sigma2 = sigma2,
    alpha = alpha
  )
}

#' @exportClass initC
setClass(
  "initAlgo",
  slots = list(
    beta_hdim = "matrixORNULL",
    beta_ldim = "matrixORNULL",
    gamma_hdim = "matrixORNULL",
    gamma_ldim = "matrixORNULL",
    sigma2 = "numeric",
    alpha = "numericORNULL"
  ),
  prototype = list(
    beta_hdim = NULL,
    beta_ldim = NULL,
    gamma_hdim = NULL,
    gamma_ldim = NULL,
    sigma2 = 10,
    alpha = NULL
  )
  # Pas besoin de validity, je crée ces objets dans le code
)

## -- Hyperparameters

#' @exportClass hyperC
setClass(
  "hyperC",
  slots = list(
    nu1 = "numericORNULL",
    nsig = "numericORNULL",
    lsig = "numericORNULL",
    a = "numericORNULL",
    b = "numericORNULL",
    sigma2_mu = "numericORNULL",
    sgam = "matrixORNULL",
    d = "numericORNULL"
  ),
  prototype = list(
    # nu0 = NULL,
    nu1 = 12000 # ,
    # nsig = 1,
    # lsig = 1,
    # a = NULL,
    # b = NULL,
    # sgam = NULL,
    # d = NULL,
    # sigma2_mu = 3000**2
  ),
  validity = function(object) {
    # if (!is.null(object@nu0)) {
    #   check_positive_slot(object@nu0, "Hyperparameter 'nu0' ")
    # }

    if (!is.null(object@nu1)) {
      check_positive_slot(object@nu1, "Hyperparameter 'nu1' ")
    }

    if (!is.null(object@nsig)) {
      check_positive_slot(object@nsig, "Hyperparameter 'nsig' ")
    }

    if (!is.null(object@lsig)) {
      check_positive_slot(object@lsig, "Hyperparameter 'lsig' ")
    }


    if (!is.null(object@ sigma2_mu)) {
      check_positive_slot(object@sigma2_mu, "Hyperparameter 'sigma2_mu' ")
    }

    if (!is.null(object@d)) {
      check_positive_slot(object@d, "Hyperparameter 'd' ")
    }

    if (!is.null(object@sgam)) {
      check_covariance(object@sgam, "Hyperparameter 'sgam' ")
    }
    TRUE
  }
)

#' @export
hyperC <- function(nu1 = 12000, sgam, d) {
  # L'utilisateur est obligé de rentrer des valeurs pour nu0, sgam, d
  # (cf prototype de la classe)
  new("hyperC",
    nu1 = nu1, nsig = 1, lsig = 1, a = NULL, b = NULL,
    sgam = sgam, d = d, sigma2_mu = 3000**2
  )
}

#' @exportClass fullHyperC
setClass(
  "fullHyperC",
  slots = list(
    nu0 = "numericORNULL"
  ),
  contains = "hyperC"
)

#' @export
fullHyperC <- function(nu0 = NULL, hc) {
  new(
    "fullHyperC",
    nu0 = nu0,
    nu1 = hc@nu1,
    nsig = hc@nsig,
    lsig = hc@lsig,
    a = hc@a,
    b = hc@b,
    sgam = hc@sgam,
    d = hc@d,
    sigma2_mu = hc@sigma2_mu
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
