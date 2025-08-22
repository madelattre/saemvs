setClassUnion("matrixORNULL", c("matrix", "NULL"))
setClassUnion("numericORNULL", c("numeric", "NULL"))
setClassUnion("listORNULL", c("list", "NULL"))

## -- Data

#' @exportClass dataC
setClass(
  "dataC",
  slots = list(
    y_list = "list",
    t_list = "list",
    x_sel = "matrixORNULL",
    x_forced = "matrixORNULL" # ,
    # tx_x = "matrixORNULL",
    # kron_tx_x = "matrixORNULL",
    # x_forced_sel = "matrixORNULL"
  ),
  prototype = list(
    y_list = list(),
    t_list = list(),
    x_sel = NULL,
    x_forced = NULL # ,
    # tx_x = NULL,
    # kron_tx_x = NULL,
    # x_forced_sel = NULL
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

    if (!is.null(object@x_sel) && nrow(object@x_sel) != n) {
      return(
        paste0(
          "The matrix 'x_sel' must have exactly", n,
          "rows, to match the", n, "individual sequences in 'y'.
          Currently, 'x_sel' has,", nrow(object@x_sel), "rows."
        )
      )
    }

    if (!is.null(object@x_forced) && nrow(object@x_forced) != n) {
      return(
        paste0(
          "The matrix 'x_forced' must have exactly", n,
          "rows, to match the", n, "individual sequences in 'y'.
          Currently, 'x_forced' has,", nrow(object@x_forced), "rows."
        )
      )
    }

    TRUE
  }
)

#' @export
dataC <- function(
    y, t, x_sel = NULL, x_forced = NULL) {
  new("dataC",
    y_list = y, t_list = t, x_sel = x_sel, x_forced = x_forced
  )
}

#' @exportClass dataAlgo
setClass(
  "dataAlgo",
  slots = list( # On déclare ici les slots supplémentaires par rapport à dataC
    # Ce sont les slots qui vont être utilisés dans les algos
    # (pas les autres déjà présents dans dataC)
    x_phi_sel = "matrixORNULL",
    x_phi_insel = "matrixORNULL",
    tx_x_phi_sel = "matrixORNULL",
    kron_tx_x_phi_sel = "matrixORNULL",
    x_phi_insel_list = "listORNULL"
  ),
  prototype = list(
    x_phi_sel = NULL,
    x_phi_insel = NULL,
    tx_x_phi_sel = NULL,
    kron_tx_x_phi_sel = NULL,
    x_phi_insel_list = NULL
  ),
  contains = "dataC"
)



## -- Model

#' @exportClass modelC
setClass(
  "modelC",
  slots = list(
    model_func = "function",
    phi_dim = "numeric",
    phi_sel_idx = "numericORNULL",
    phi_fixed_idx = "numericORNULL",
    x_forced_support = "matrixORNULL"
  ),
  prototype = list(
    model_func = function(phi, t) numeric(0),
    phi_dim = integer(0),
    phi_sel_idx = c(),
    phi_fixed_idx = c(),
    x_forced_support = matrix(numeric(0), nrow = 0, ncol = 0)
  ),
  validity = function(object) {
    args <- names(formals(object@model_func))

    if (!all(c("phi", "t") %in% args)) {
      return("The model function must have arguments 'phi' and 't'.")
    }

    if (!(is.integer(object@phi_dim) && length(object@phi_dim) == 1 &&
      object@phi_dim >= 1)) {
      return("'phi_dim' must be a positive integer value.")
    }

    if (!(is.integer(object@phi_sel_idx) && all(object@phi_sel_idx >= 1) &&
      all(object@phi_sel_idx <= object@phi_dim))) {
      return("'phi_sel_idx' must be a vector of integer values between 1
      and 'phi_dim'.")
    }

    if (!(is.integer(object@phi_fixed_idx) && all(object@phi_fixed_idx >= 1) &&
      all(object@phi_fixed_idx <= object@phi_dim))) {
      return("'phi_fixed_idx' must be a vector of integer values between 1
      and 'phi_dim'.")
    }

    if (length(intersect(object@phi_sel_idx, object@phi_fixed_idx)) > 0) {
      return("'phi_sel_idx' and 'phi_fixed_idx' must not have any elements in
      common.")
    }

    supp <- object@x_forced_support

    if (!is.null(supp)) {
      if (!is.matrix(supp)) {
        return("'x_forced_support' must be a matrix or NULL.")
      }
      if (!all(dim(supp) == c(0, 0))) {
        if (ncol(supp) != object@phi_dim) {
          return(
            paste0(
              "'x_forced_support' must have ", object@phi_dim,
              " columns: one for each component in 'phi'."
            )
          )
        }
        # On vérifie le nombre de lignes dans 'x_forced_support' en utilisant les information de l'objet data
      }
    }
    TRUE
  }
)

#' @export
modelC <- function(
    g, phi_dim, phi_sel_idx = c(), phi_fixed_idx = c(),
    x_forced_support = matrix(numeric(0), nrow = 0, ncol = 0)) {
  new("modelC",
    model_func = g,
    phi_dim = as.integer(phi_dim),
    phi_sel_idx = as.integer(phi_sel_idx),
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
