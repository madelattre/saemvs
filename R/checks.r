check_covariance <- function(mat, name_mat) {
  if (ncol(mat) != nrow(mat)) {
    return(
      paste0(name_mat, "must be a squared matrix.")
    )
  }

  if (!all.equal(mat, t(mat), tolerance = 1e-8)) {
    return(
      paste0(name_mat, "must be a symmetric matrix.")
    )
  }

  eigenvalues <- eigen(mat, symmetric = TRUE)$values
  if (min(eigenvalues) < 1e-8) {
    return(
      paste0(name_mat, "must be a positive definite matrix.")
    )
  }
}

check_beta_gamma <- function(beta, gamma, q, case = "ns") {
  if (case == "ns") {
    start_case <- paste0(
      "As ", q, " parameters are not subject to selection",
      "(", q, " missing parameter indices in 'index_select'),"
    )
  } else {
    start_case <- paste0(
      "As ", q, " parameters are subject to selection",
      "(", q, " indices in 'index_select'),"
    )
  }

  if (q != 0) {
    if (is.null(beta) || is.null(gamma)) {
      stop(paste0(
        start_case, " 'beta_", case, "' and 'gamma_", case,
        "' should both be non empty matrices."
      ))
    }
  } else {
    if (ncol(beta) != ncol(gamma)) {
      stop(paste0(
        start_case, " 'beta_", case, "' and 'gamma_", case,
        "' must both have", q, "columns."
      ))
    }
  }
}

check_beta_support <- function(beta, support) {
  # A appeler uniquement si le support n'est pas NULL
  if ((is.null(support)) && (nrow(beta) > 1)) {
    stop(
      paste0(
        "'beta_ns' must have only one row (intercepts) ",
        "as no 'covariate_support' is provided in the model."
      )
    )
  }
  if (!is.null(support) && (nrow(support) != (nrow(beta) - 1))) {
    stop(
      paste0(
        "'beta_ns' must have ", nrow(support) + 1, " rows ",
        "(intercept + ", nrow(support), " covariates) as 'support '",
        "is designed for ", nrow(support), " covariates (",
        nrow(support), " rows)."
      )
    )
  }
}

check_beta_hdim <- function(beta, pv) {
  if (pv != (nrow(beta))) {
    stop(
      paste0(
        "'beta_s' must have ", pv, " rows ",
        "(intercept + ", pv - 1, " covariates) as matrix 'v' ",
        "contains ", pv - 1, " covariates (",
        pv - 1, " columns)."
      )
    )
  }
}

check_positive_slot <- function(slot, name_slot) {
  if (slot <= 0) {
    return(
      paste0(name_slot, "must be strictly positive.")
    )
  }
}

check_positive_or_null_slot <- function(slot, name_slot) {
  if (!is.null(slot)) {
    if (slot <= 0) {
      return(
        paste0(name_slot, "must be NULL or strictly positive.")
      )
    }
  }
}

check_positive_integer_slot <- function(slot, name_slot) {
  if ((slot <= 0) || !is.integer(slot)) {
    return(
      paste0(name_slot, "must be a strictly positive integer.")
    )
  }
}
