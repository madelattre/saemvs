#' Transform data object from MAP to MLE representation
#'
#' Internal method that restructures a \code{saemvsData} object when moving
#' from a MAP (maximum a posteriori) representation to an MLE (maximum
#' likelihood estimation) representation.
#'
#' The function selects only the candidate covariates that are active in the
#' support matrix, and appends them to the forced covariates. The resulting
#' dataset is compatible with an MLE-based model fit.
#'
#' @param data An object of class \code{saemvsData}, containing the response
#'   vector, time series, forced covariates, and candidate covariates.
#' @param cand_support A binary matrix of dimension \code{(p + 1) x q}, where:
#'   \itemize{
#'     \item The first row corresponds to the intercept and must contain only
#'     ones.
#'     \item The following \code{p} rows correspond to the candidate covariates
#'     in \code{data@x_candidates}.
#'     \item Each column corresponds to one parameter.
#'     \item An entry equal to 1 indicates inclusion of the covariate in the
#'     support.
#'   }
#'
#' @return A new object of class \code{saemvsData}, with updated forced
#' covariates that include both the original forced covariates and the subset of
#' active candidate covariates.
#'
#' @details
#' Robustness checks are performed to ensure:
#' \itemize{
#'   \item The number of rows in \code{cand_support} matches
#'   \code{ncol(data@x_candidates) + 1}.
#'   \item The first row of \code{cand_support} (intercept) contains only ones.
#' }
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' new_data <- from_map_to_mle_data(saemvsData_obj, cand_support)
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
      stop(
        sprintf(
          paste0(
            "Dimension mismatch: 'cand_support' has %d rows, but should have ",
            "%d (1 intercept + %d candidates)."
          ),
          nrow(cand_support),
          expected_nrow,
          ncol(data@x_candidates)
        )
      )
    }

    if (!all(cand_support[1, ] == 1)) {
      stop(
        paste0(
          "Internal error: first row of 'cand_support' must correspond ",
          "to the intercept (all ones)."
        )
      )
    }

    active_candidate_idx <- which(rowSums(cand_support[-1, , drop = FALSE]) > 0)

    restricted_candidates <-
      data@x_candidates[, active_candidate_idx, drop = FALSE]

    if (dim(restricted_candidates)[2] == 0) {
      restricted_candidates <- NULL
    }


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
#' Internal method that transforms a \code{saemvsModel} object when moving from
#' the MAP setting to the MLE setting, updating the support of forced and
#' selected covariates accordingly.
#'
#' @param model An object of class \code{saemvsProcessedModel}.
#' @param cand_support A binary matrix indicating candidate covariate support
#' (rows = covariates, columns = parameters).
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
  signature(model = "saemvsProcessedModel", cand_support = "matrix"),
  function(model, cand_support) {
    if (!is.matrix(cand_support)) {
      stop("'cand_support' must be a matrix.")
    }
    nb_phi_select <- length(model@phi_to_select_idx)

    if (ncol(cand_support) != nb_phi_select) {
      stop(
        paste0(
          "Number of columns in 'cand_support' must match ",
          "number of phi-to-select parameters."
        )
      )
    }

    active_rows <- which(rowSums(cand_support) > 0)
    covariate_restricted_support <- cand_support[active_rows, , drop = FALSE]

    if (nrow(covariate_restricted_support) > 1) {
      selected_support_matrix <- matrix(
        covariate_restricted_support[-1, ],
        ncol = nb_phi_select
      )
    } else {
      selected_support_matrix <- NULL
    }

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

    methods::new("saemvsProcessedModel",
      model_func = model@model_func,
      phi_dim = model@phi_dim,
      phi_fixed_idx = model@phi_fixed_idx,
      x_forced_support = new_support
    )
  }
)


#' Convert initialization from MAP to MLE representation (internal)
#'
#' Internal method that transforms a \code{saemvsInit} object when moving from
#' the MAP setting to the MLE setting, restricting beta initialization to forced
#' and selected candidate covariates.
#'
#' @param init An object of class \code{saemvsInit}.
#' @param model An object of class \code{saemvsModel}.
#' @param cand_support A binary matrix indicating candidate covariate support
#' (rows = covariates, columns = parameters).
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
  signature(
    init = "saemvsProcessedInit", model = "saemvsProcessedModel",
    cand_support = "matrix"
  ),
  function(init, model, cand_support) {
    if (!is.matrix(cand_support)) {
      stop("'cand_support' must be a matrix.")
    }

    if (ncol(cand_support) != length(model@phi_to_select_idx)) {
      stop(
        paste0(
          "Number of columns in 'cand_support' must match ",
          "number of phi-to-select parameters."
        )
      )
    }

    active_candidate_idx <- which(rowSums(cand_support[-1, , drop = FALSE]) > 0)

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
