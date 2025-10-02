# ---------------------------------------------
# -- check_data: verify that the user data fits the model
# ---------------------------------------------

#' Check Data Consistency with Model
#'
#' Internal function to verify that the data provided by the user
#' is consistent with the SAEMVS model specification.
#'
#' @param data An object of class \code{saemvsData}.
#' @param model An object of class \code{saemvsModel}.
#'
#' @details
#' Performs the following checks:
#' \itemize{
#'   \item If \code{phi_to_select_idx} is non-empty, \code{x_candidates} must be provided.
#'   \item \code{x_candidates}, if present, must have one row per sequence in \code{y_series}.
#'   \item \code{x_forced} and \code{x_forced_support}, if present, must be numeric matrices
#'         with compatible dimensions.
#' }
#'
#' @keywords internal
setGeneric(
  "check_data",
  function(data, model) {
    standardGeneric("check_data")
  }
)

setMethod(
  "check_data",
  signature(data = "saemvsData", model = "saemvsModel"),
  function(data, model) {
    n_ind <- length(data@y_series)

    if (length(model@phi_to_select_idx) > 0 && is.null(data@x_candidates)) {
      stop(sprintf(
        "Parameter selection requested (%d phi_to_select_idx), but 'x_candidates' is missing.",
        length(model@phi_to_select_idx)
      ))
    }

    if (!is.null(data@x_candidates) && nrow(data@x_candidates) != n_ind) {
      stop(sprintf(
        "'x_candidates' must have %d rows (matching y_series); got %d.",
        n_ind, nrow(data@x_candidates)
      ))
    }

    if (!is.null(data@x_forced) && !is.null(model@x_forced_support)) {
      if (!is.matrix(data@x_forced)) stop("'x_forced' must be a matrix.")
      if (!is.numeric(data@x_forced)) stop("'x_forced' must be numeric.")
      if (ncol(data@x_forced) != nrow(model@x_forced_support)) {
        stop(sprintf(
          "Number of columns in 'x_forced' (%d) must equal number of rows in 'x_forced_support' (%d).",
          ncol(data@x_forced), nrow(model@x_forced_support)
        ))
      }
    }
  }
)

# ---------------------------------------------
# -- check_init: verify that initial parameters fit the model
# ---------------------------------------------

#' Check Initial Values Consistency
#'
#' Internal function to verify that the initial values provided
#' by the user match the SAEMVS model and data dimensions.
#'
#' @param init An object of class \code{saemvsInit}.
#' @param data An object of class \code{saemvsData}.
#' @param model An object of class \code{saemvsModel}.
#'
#' @details
#' Performs the following checks:
#' \itemize{
#'   \item Length of \code{intercept} must match \code{phi_dim} in the model.
#'   \item \code{beta_candidates} and \code{beta_forced}, if present, must have
#'         the same number of columns as \code{intercept}.
#'   \item \code{cov_re} must be a numeric square matrix of size \code{phi_dim}.
#'   \item \code{sigma2} must be a single strictly positive number.
#' }
#'
#' @keywords internal
setGeneric(
  "check_init",
  function(init, data, model) {
    standardGeneric("check_init")
  }
)

setMethod(
  "check_init",
  signature(init = "saemvsInit", data = "saemvsData", model = "saemvsModel"),
  function(init, data, model) {
    n_phi <- length(init@intercept)
    if (n_phi != model@phi_dim) {
      stop(sprintf(
        "Length of 'intercept' (%d) must match 'phi_dim' (%d).",
        n_phi, model@phi_dim
      ))
    }

    if (!is.null(init@beta_candidates)) {
      if (!is.matrix(init@beta_candidates) || ncol(init@beta_candidates) != n_phi) {
        stop(sprintf(
          "'beta_candidates' must be a matrix with %d columns (matching intercept); got %d.",
          n_phi, ifelse(is.null(ncol(init@beta_candidates)), NA, ncol(init@beta_candidates))
        ))
      }
    }

    if (!is.null(init@beta_forced)) {
      if (!is.matrix(init@beta_forced) || ncol(init@beta_forced) != n_phi) {
        stop(sprintf(
          "'beta_forced' must be a matrix with %d columns (matching intercept); got %d.",
          n_phi, ifelse(is.null(ncol(init@beta_forced)), NA, ncol(init@beta_forced))
        ))
      }
    }

    if (!is.matrix(init@cov_re) || nrow(init@cov_re) != n_phi || ncol(init@cov_re) != n_phi) {
      stop(sprintf(
        "'cov_re' must be a %dx%d numeric matrix; got %dx%d.",
        n_phi, n_phi, nrow(init@cov_re), ncol(init@cov_re)
      ))
    }

    if (!is_empty_support(model@x_forced_support)) {
      xfs <- model@x_forced_support
      bfs <- init@beta_forced

      if (is_empty_matrix(bfs) && isFALSE(init@default)) {
        stop(
          "Non-empty 'x_forced_support' requires a non-empty initialisation of 'beta_forced'. ",
          "Alternatively, set 'init@default = TRUE'."
        )
      }

      if (!all(dim(xfs) == dim(bfs))) {
        stop(
          "'x_forced_support' and 'beta_forced' must have the same dimensions: ",
          "got ", paste0(nrow(xfs), "x", ncol(xfs)), " vs ",
          paste0(nrow(bfs), "x", ncol(bfs)), "."
        )
      }

      if (isFALSE(init@default)) {
        mask_violation <- any(bfs[xfs == 0] != 0)
        if (mask_violation) {
          stop(
            "Inconsistent forced values: where 'x_forced_support' has 0, ",
            "'beta_forced' must also be 0 (since init@default = FALSE)."
          )
        }
      }
    }


    # phi_idx <- model@phi_fixed_idx
    # if (length(phi_idx) > 0) {
    #   diag_vals <- diag(init@cov_re[phi_idx, phi_idx, drop = FALSE])

    #   if (any(diag_vals != 0)) {
    #     warning(sprintf(
    #       "Diagonal components of 'cov_re' at phi_fixed_idx (%s) are not zero (%s).
    #    According to the model specification they should be zero, and will be ignored in the procedure.",
    #       paste(phi_idx[diag_vals != 0], collapse = ", "),
    #       paste(diag_vals[diag_vals != 0], collapse = ", ")
    #     ))
    #   }
    # }

    phi_idx <- model@phi_fixed_idx
    diag_vals <- diag(init@cov_re)

    zero_diag_idx <- which(diag_vals == 0)

    if (!setequal(zero_diag_idx, phi_idx)) {
      stop(sprintf(
        "Mismatch between zero diagonals in 'cov_re' (%s) and 'phi_fixed_idx' (%s).",
        paste(zero_diag_idx, collapse = ", "),
        paste(phi_idx, collapse = ", ")
      ))
    }

    if (!is.numeric(init@sigma2) || length(init@sigma2) != 1 || init@sigma2 <= 0) {
      stop("'sigma2' must be a single strictly positive number.")
    }
  }
)

# ---------------------------------------------
# -- check_hyper: verify that hyperparameters fit the model and tuning
# ---------------------------------------------

#' Check Hyperparameters Consistency
#'
#' Internal function to verify that the hyperparameters for the spike-and-slab
#' prior are consistent with the SAEMVS model and tuning parameters.
#'
#' @param hyper An object of class \code{saemvsHyperSpikeAndSlab}.
#' @param model An object of class \code{saemvsModel}.
#' @param tuning An object of class \code{saemvsTuning}.
#'
#' @details
#' Performs the following checks:
#' \itemize{
#'   \item \code{inclusion_prob_prior_a} and \code{inclusion_prob_prior_b} lengths
#'         must match the number of phi parameters subject to selection.
#'   \item \code{cov_re_prior_scale} must be a numeric square matrix with
#'         dimensions equal to the number of phi parameters selected.
#'   \item All values in \code{spike_values_grid} must be smaller than \code{slab_parameter}.
#' }
#'
#' @keywords internal
setGeneric(
  "check_hyper",
  function(hyper, model, tuning) {
    standardGeneric("check_hyper")
  }
)

setMethod(
  "check_hyper",
  signature(
    hyper = "saemvsHyperSpikeAndSlab", model = "saemvsModel",
    tuning = "saemvsTuning"
  ),
  function(hyper, model, tuning) {
    nb_selected <- length(model@phi_to_select_idx)

    if (nb_selected > 0) {
      if (!is.null(hyper@inclusion_prob_prior_a) && length(hyper@inclusion_prob_prior_a) != nb_selected) {
        stop(sprintf(
          "Hyperparameter 'inclusion_prob_prior_a' must have %d components; got %d.",
          nb_selected, length(hyper@inclusion_prob_prior_a)
        ))
      }
      if (!is.null(hyper@inclusion_prob_prior_b) && length(hyper@inclusion_prob_prior_b) != nb_selected) {
        stop(sprintf(
          "Hyperparameter 'inclusion_prob_prior_b' must have %d components; got %d.",
          nb_selected, length(hyper@inclusion_prob_prior_b)
        ))
      }
      if (!is.null(hyper@cov_re_prior_scale)) {
        if (!is.matrix(hyper@cov_re_prior_scale) || !is.numeric(hyper@cov_re_prior_scale)) {
          stop("'cov_re_prior_scale' must be a numeric matrix.")
        }
        if (any(dim(hyper@cov_re_prior_scale) != nb_selected)) {
          stop(sprintf(
            "'cov_re_prior_scale' must be a %dx%d matrix; got %dx%d.",
            nb_selected, nb_selected,
            nrow(hyper@cov_re_prior_scale), ncol(hyper@cov_re_prior_scale)
          ))
        }
      }
    }

    if (!is.null(tuning@spike_values_grid) && length(tuning@spike_values_grid) > 0) {
      if (!all(tuning@spike_values_grid < hyper@slab_parameter)) {
        stop("'spike_values_grid' values must all be smaller than 'slab_parameter'.")
      }
    }
  }
)
