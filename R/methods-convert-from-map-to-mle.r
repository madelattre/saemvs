## ===============================================================
## MAP → MLE conversion methods (internal)
## ===============================================================

#' Transform data object from MAP to MLE representation
#'
#' Internal method that restructures a \code{saemvsData} object when moving
#' from a MAP (maximum a posteriori) representation to an MLE (maximum likelihood
#' estimation) representation.  
#'
#' The function selects only the candidate covariates that are active in the
#' support matrix, and appends them to the forced covariates. The resulting
#' dataset is compatible with an MLE-based model fit.
#'
#' @param data An object of class \code{saemvsData}, containing the response
#'   vector, time series, forced covariates, and candidate covariates.
#' @param cand_support A binary matrix of dimension \code{(p + 1) x q}, where:
#'   \itemize{
#'     \item The first row corresponds to the intercept and must contain only ones.
#'     \item The following \code{p} rows correspond to the candidate covariates
#'           in \code{data@x_candidates}.
#'     \item Each column corresponds to one parameter.
#'     \item An entry equal to 1 indicates inclusion of the covariate in the support.
#'   }
#'
#' @return A new object of class \code{saemvsData}, with updated forced covariates
#'   that include both the original forced covariates and the subset of active
#'   candidate covariates.
#'
#' @details
#' Robustness checks are performed to ensure:
#' \itemize{
#'   \item The number of rows in \code{cand_support} matches
#'         \code{ncol(data@x_candidates) + 1}.
#'   \item The first row of \code{cand_support} (intercept) contains only ones.
#' }
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'   new_data <- from_map_to_mle_data(saemvsData_obj, cand_support)
#' }
setGeneric(
  "map_to_mle_data",
  function(data, cand_support) {
    standardGeneric("map_to_mle_data")
  }
)

setMethod(
  "map_to_mle_data",
  signature(data = "saemvsData", cand_support = "matrix"),
  function(data, cand_support) {
    if (!is.matrix(cand_support)) {
      stop("'cand_support' must be a matrix.")
    }

    expected_nrow <- ncol(data@x_candidates) + 1
    if (nrow(cand_support) != expected_nrow) {
      stop(sprintf(
        "Dimension mismatch: 'cand_support' has %d rows, but should have %d (1 intercept + %d candidates).",
        nrow(cand_support), expected_nrow, ncol(data@x_candidates)
      ))
    }

    # Sanity check: first row must be intercept (all ones)
    if (!all(cand_support[1, ] == 1)) {
      stop("Internal error: first row of 'cand_support' must correspond to the intercept (all ones).")
    }

    # Identify candidate covariates with at least one inclusion (excluding intercept row)
    active_candidate_idx <- which(rowSums(cand_support[-1, , drop = FALSE]) > 0)

    # Restricted candidate design matrix
    restricted_candidates <- data@x_candidates[, active_candidate_idx, drop = FALSE]


    # New forced design matrix = old forced + restricted candidates
    new_x_forced <- cbind(data@x_forced, restricted_candidates)

    saemvsData(
      y = data@y_series,
      t = data@t_series,
      x_forced = new_x_forced
    )
  }
)


#' Convert model from MAP to MLE representation (internal)
#'
#' Internal method that transforms a \code{saemvsModel} object
#'   when moving from the MAP setting to the MLE setting, updating the support
#'   of forced and selected covariates accordingly.
#'
#' @param model An object of class \code{saemvsModel}.
#' @param cand_support A binary matrix indicating candidate covariate
#'   support (rows = covariates, columns = parameters).
#'
#' @return A \code{saemvsModel} object with updated support matrix.
#'
#' @keywords internal
setGeneric(
  "map_to_mle_model",
  function(model, cand_support) {
    standardGeneric("map_to_mle_model")
  }
)

setMethod(
  "map_to_mle_model",
  signature(model = "saemvsModel", cand_support = "matrix"),
  function(model, cand_support) {
    if (!is.matrix(cand_support)) {
      stop("'cand_support' must be a matrix.")
    }
    nb_phi_select <- length(model@phi_to_select_idx)

    if (ncol(cand_support) != nb_phi_select) {
      stop("Number of columns in 'cand_support' must match number of phi-to-select parameters.")
    }

    # Extract covariates contributing to phi_to_select
    active_rows <- which(rowSums(cand_support) > 0)
    covariate_restricted_support <- cand_support[active_rows, , drop = FALSE]

    # If more than intercept, build the restricted support matrix
    if (nrow(covariate_restricted_support) > 1) {
      selected_support_matrix <- matrix(
        covariate_restricted_support[-1, ],
        ncol = nb_phi_select
      )
    } else {
      selected_support_matrix <- NULL
    }

    # Build new full support by merging forced support with restricted selection
    all_phi_idx <- seq_len(model@phi_dim)
    phi_unselect_idx <- setdiff(all_phi_idx, model@phi_to_select_idx)
    perm <- c(model@phi_to_select_idx, phi_unselect_idx)
    inv_perm <- match(seq_along(perm), perm)

    new_support <- merge_support(
      model@x_forced_support,
      selected_support_matrix,
      nb_phi_select,
      model@phi_dim - nb_phi_select,
      inv_perm
    )

    saemvsModel(
      g = model@model_func,
      phi_dim = model@phi_dim,
      phi_fixed_idx = model@phi_fixed_idx,
      x_forced_support = new_support
    )
  }
)


#' Convert initialization from MAP to MLE representation (internal)
#'
#' Internal method that transforms a \code{saemvsInit} object
#'   when moving from the MAP setting to the MLE setting, restricting
#'   beta initialization to forced and selected candidate covariates.
#'
#' @param init An object of class \code{saemvsInit}.
#' @param model An object of class \code{saemvsModel}.
#' @param cand_support A binary matrix indicating candidate covariate
#'   support (rows = covariates, columns = parameters).
#'
#' @return A \code{saemvsInit} object with restricted initialization.
#'
#' @keywords internal
setGeneric(
  "map_to_mle_init",
  function(init, model, cand_support) {
    standardGeneric("map_to_mle_init")
  }
)

setMethod(
  "map_to_mle_init",
  signature(init = "saemvsInit", model = "saemvsModel", cand_support = "matrix"),
  function(init, model, cand_support) {
    if (!is.matrix(cand_support)) {
      stop("'cand_support' must be a matrix.")
    }

    if (ncol(cand_support) != length(model@phi_to_select_idx)) {
      stop("Number of columns in 'cand_support' must match number of phi-to-select parameters.")
    }

    # Identify rows with at least one active covariate (excluding intercept)
    active_candidate_idx <- which(rowSums(cand_support[-1, , drop = FALSE]) > 0)

    # Build restricted beta initialization
    new_beta_init <- rbind(
      init@beta_forced,
      init@beta_candidates[active_candidate_idx, , drop = FALSE]
    )

    saemvsInit(
      intercept = init@intercept,
      beta_forced = new_beta_init,
      cov_re = init@cov_re,
      sigma2 = init@sigma2
    )
  }
)



# ## -- Prepare data, model and initialization from map to mle

# setGeneric(
#   "from_map_to_mle_data",
#   function(data, cand_support) {
#     standardGeneric("from_map_to_mle_data")
#   }
# )

# setMethod(
#   "from_map_to_mle_data",
#   signature(
#     data = "saemvsData", cand_support = "matrix"
#   ),
#   function(data, cand_support) {
#     n <- length(data@y_series)

#     lines_with_ones <- which(rowSums(cand_support) > 0)

#     v_restricted <- matrix(
#       cbind(1, data@x_candidates)[, lines_with_ones],
#       nrow = n
#     )[, -1] # -1 pour ne pas contenir l'intercept

#     # Si on ne force pas l'inclusion de variables dans phi_sel
#     new_x <- matrix(
#       cbind(data@x_forced, v_restricted),
#       nrow = n
#     )

#     new_data <- saemvsData(
#       y = data@y_series,
#       t = data@t_series,
#       x_forced = new_x
#     )

#     return(new_data)
#   }
# )

# setGeneric(
#   "from_map_to_mle_model",
#   function(model, cand_support) {
#     standardGeneric("from_map_to_mle_model")
#   }
# )

# setMethod(
#   "from_map_to_mle_model",
#   signature(
#     model = "saemvsModel", cand_support = "matrix"
#   ),
#   function(model, cand_support) {
#     # Ici il faut tenir compte des indices des phi_s et des phi_ns

#     all_phi_index <- seq(1, model@phi_dim)
#     index_unselect <- setdiff(all_phi_index, model@phi_to_select_idx)
#     perm <- c(model@phi_to_select_idx, index_unselect)
#     inv_perm <- match(seq_along(perm), perm)

#     nb_phi_s <- length(model@phi_to_select_idx)

#     lines_with_ones <- which(rowSums(cand_support) > 0)

#     covariate_restricted_support <- cand_support[lines_with_ones, ]

#     if (is.matrix(covariate_restricted_support)) {
#       selected_support_matrix <- matrix(
#         covariate_restricted_support[-1, ],
#         ncol = nb_phi_s,
#         nrow = length(covariate_restricted_support[-1, ]) / nb_phi_s
#       )
#     } else {
#       selected_support_matrix <- NULL
#     }

#     new_covariate_support <- merge_support(
#       model@x_forced_support,
#       selected_support_matrix,
#       # data@w,
#       nb_phi_s,
#       model@phi_dim - nb_phi_s,
#       inv_perm
#     )

#     new_model <- saemvsModel(
#       g = model@model_func,
#       phi_dim = model@phi_dim,
#       phi_fixed_idx = model@phi_fixed_idx,
#       x_forced_support = new_covariate_support
#     )

#     return(new_model)
#   }
# )


# setGeneric(
#   "from_map_to_mle_init",
#   function(init, model, cand_support) {
#     standardGeneric("from_map_to_mle_init")
#   }
# )

# setMethod(
#   "from_map_to_mle_init",
#   signature(
#     init = "saemvsInit", model = "saemvsModel", cand_support = "matrix"
#   ),
#   # Ici il faut tenir compte des indices des phi_s et des phi_ns
#   function(init, model, cand_support) {
#     lines_with_ones <- which(rowSums(cand_support[-1, ]) > 0)

#     new_beta_init <- rbind(
#       init@beta_forced,
#       init@beta_candidates[lines_with_ones, ]
#     )

#     # if (is.null(init@gamma_ldim)) {
#     #   new_gamma_init <- init@gamma_hdim
#     # } else {
#     #   new_gamma_init <- as.matrix(Matrix::bdiag(
#     #     init@gamma_hdim,
#     #     init@gamma_ldim
#     #   ))
#     # }

#     new_init <- saemvsInit(
#       intercept = init@intercept,
#       beta_forced = new_beta_init,
#       cov_re = init@cov_re,
#       sigma2 = init@sigma2
#     )

#     return(new_init)
#   }
# )
