setClassUnion("matrixORNULL", c("matrix", "NULL"))
setClassUnion("numericORNULL", c("numeric", "NULL"))
setClassUnion("listORNULL", c("list", "NULL"))
# setClassUnion("characterORNULL", c("character", "NULL"))

## -- Data

#' @exportClass dataC
setClass(
  "dataC",
  slots = list(
    y_list = "list",
    t_list = "list",
    v = "matrixORNULL",
    tv_v = "matrixORNULL",
    kron_tv_v = "matrixORNULL",
    w = "matrixORNULL",
    x = "listORNULL"
  ),
  prototype = list(
    y_list = list(),
    t_list = list(),
    v = NULL,
    tv_v = NULL,
    kron_tv_v = NULL,
    w = NULL,
    x = list()
  ),
  validity = function(object) {
    if (length(object@y_list) != length(object@t_list)) {
      return("Lists y et t must have the same number of elements.")
    }

    for (i in seq_along(object@y_list)) {
      if (length(object@y_list[[i]]) != length(object@t_list[[i]])) {
        return(
          paste0("y[[", i, "]] and t[[", i, "]] do not have the
          same number of elements.")
        )
      }
    }

    n <- length(object@y_list)

    if (!is.null(object@v) && nrow(object@v) != n) {
      return(
        paste0(
          "The matrix 'v' must have exactly", n,
          "rows, to match the", n, "individual sequences in 'y'.
          Currently, 'v' has,", nrow(v), "rows."
        )
      )
    }

    if (!is.null(object@w) && nrow(object@w) != n) {
      return(
        paste0(
          "The matrix 'w' must have exactly", n,
          "rows, to match the", n, "individual sequences in 'y'.
          Currently, 'w' has,", nrow(w), "rows."
        )
      )
    }

    TRUE
  }
)

#' @export
dataC <- function(
    y, t,
    v = NULL, w = NULL, tv_v = NULL, kron_tv_v = NULL, x = NULL) {
  new("dataC",
    y_list = y, t_list = t, v = v, w = w, tv_v = tv_v,
    kron_tv_v = kron_tv_v, x = x
  )
}

## -- Model

#' @exportClass modelC
setClass(
  "modelC",
  slots = list(
    model_func = "function",
    q_phi = "numeric",
    index_select = "numericORNULL",
    index_fixed = "numericORNULL",
    covariate_support = "matrixORNULL" # ne doit pas contenir l'intercept
  ),
  prototype = list(
    model_func = function(phi, t) numeric(0),
    q_phi = integer(0),
    index_select = c(), # NULL?
    index_fixed = c(), # NULL?
    covariate_support = matrix(numeric(0), nrow = 0, ncol = 0) # NULL?
  ),
  validity = function(object) {
    args <- names(formals(object@model_func))

    if (!all(c("phi", "t") %in% args)) {
      return("The model function must have arguments 'phi' and 't'.")
    }

    if (!(is.integer(object@q_phi) && length(object@q_phi) == 1 &&
      object@q_phi >= 1)) {
      return("q_phi must be a positive integer value.")
    }

    if (!(is.integer(object@index_select) && all(object@index_select >= 1) &&
      all(object@index_select <= object@q_phi))) {
      return("index_select must be a vector of integer values between 1
      and q_phi.")
    }

    if (!(is.integer(object@index_fixed) && all(object@index_fixed >= 1) &&
      all(object@index_fixed <= object@q_phi))) {
      return("index_fixed must be a vector of integer values between 1
      and q_phi.")
    }

    if (length(intersect(object@index_select, object@index_fixed)) > 0) {
      return("index_select and index_fixed must not have any elements in
      common.")
    }

    supp <- object@covariate_support
    if (!is.null(supp)) {
      if (!is.matrix(supp)) {
        return("covariate_support must be a matrix or NULL.")
      }
      if (!all(dim(supp) == c(0, 0))) {
        nrows_expected <- object@q_phi - length(object@index_select)
        if (ncol(supp) != nrows_expected) {
          return(sprintf("covariate_support must have %d rows
          (q_phi - length(index_select)).", nrows_expected))
        }
        if (!all(supp %in% c(0, 1))) {
          return("covariate_support must contain only 0s and 1s.")
        }
      }
    }

    TRUE
  }
)

#' @export
modelC <- function(
    g, nphi, phi_select = c(), phi_fixed = c(), support = c()) {
  new("modelC",
    model_func = g, q_phi = as.integer(nphi),
    index_select = as.integer(phi_select),
    index_fixed = as.integer(phi_fixed),
    covariate_support = support
  )
}


## -- Initilialization of unknown parameters

#' @exportClass initC
setClass(
  "initC",
  slots = list(
    beta_hdim = "matrixORNULL",
    gamma_hdim = "matrixORNULL",
    beta_ldim = "matrixORNULL",
    gamma_ldim = "matrixORNULL",
    sigma2 = "numeric",
    alpha = "numericORNULL"
  ),
  prototype = list(
    beta_hdim = NULL,
    gamma_hdim = NULL,
    beta_ldim = NULL,
    gamma_ldim = NULL,
    sigma2 = 10,
    alpha = NULL
  ),
  validity = function(object) {
    if ((!is.null(object@beta_hdim)) && (!is.null(object@gamma_hdim))) {
      if (ncol(object@gamma_hdim) != ncol(object@beta_hdim)) {
        return(
          paste0(
            "The initialization of matrices 'beta_s' and ",
            "'gamma_s' must have the same number of columns."
          )
        )
      }
    }

    if ((!is.null(object@beta_ldim)) && (!is.null(object@gamma_ldim))) {
      if (ncol(object@gamma_ldim) != ncol(object@beta_ldim)) {
        return(
          paste0(
            "The initialization of matrices 'beta_ns' and ",
            "'gamma_ns' must have the same number of columns."
          )
        )
      }
    }

    if (!is.null(object@gamma_hdim)) {
      check_covariance(object@gamma_hdim, "The initial value for 'gamma_s' ")
    }
    if (!is.null(object@gamma_ldim)) {
      check_covariance(object@gamma_ldim, "The initial value for 'gamma_ns' ")
    }

    check_positive_slot(object@sigma2, "The value for 'sigma2' ")

    if ((!is.null(object@beta_hdim)) && (!is.null(object@alpha))) {
      if (ncol(object@beta_hdim) != length(object@alpha)) {
        return(
          paste0(
            "The number of elements in 'alpha' must be equal to the",
            "number of columns in 'beta_s'."
          )
        )
      }
    }


    TRUE
  }
)

#' @export
initC <- function(
    beta_s = NULL, beta_ns = NULL, gamma_s = NULL,
    gamma_ns = NULL, sigma2 = 10, alpha = NULL) {
  new("initC",
    beta_hdim = beta_s, beta_ldim = beta_ns,
    gamma_hdim = gamma_s, gamma_ldim = gamma_ns, sigma2 = sigma2,
    alpha = alpha
  )
}


## -- Hyperparameters

#' @exportClass hyperC
setClass(
  "hyperC",
  slots = list(
    # nu0 = "numericORNULL",
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
    param = "list",
    support = "list",
    unique_support = "list",
    map_to_unique_support = "numeric",
    nu0_grid = "numeric"
  ),
  prototype = list(
    pen = character(0),
    crit_values = numeric(0),
    thresholds = list(),
    param = list(),
    support = list(),
    unique_support = list(),
    map_to_unique_support = numeric(0),
    nu0_grid = numeric(0)
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
