# Internal helper: generate convergence diagnostics for SAEM
#'
#' @param res_saem A \code{saemResults} object.
#' @param component Character string specifying which element to plot. One of:
#'   \itemize{
#'     \item `"sigma2"`: residual variance
#'     \item `"coef_phi_sel"`: regression coefficients for parameters subject
#'  to selection
#'     \item `"coef_phi_non_sel"`: regression coefficients for parameters not
#'  subject to selection
#'     \item `"variance_phi_sel"`: covariance matrix of parameters subject to
#'  selection
#'     \item `"variance_phi_non_sel"`: covariance matrix of parameters not
#' subject to selection
#'   }
#' @param sel_components Optional character vector specifying which components
#'  to display.
#'   Can be `"random"` or `"top:n"` to select up to 16 components. Defaults to
#'  the first 16.
#' @param phi Optional character or integer specifying a parameter to focus on.
#' @param ... Additional arguments (not used).
#'
#' @return A \code{ggplot2} object showing the evolution of the selected
#'  component across iterations.
#'   Returns \code{invisible(NULL)} if no data is available.
#'
#' @keywords internal
#' @noRd
.convergence_plot_backend <- function(
    res_saem, component, sel_components, phi = NULL) {
  max_ticks <- 4

  phi_to_select_idx <- res_saem@phi_to_select_idx
  phi_not_to_select_idx <- res_saem@phi_not_to_select_idx

  phi_names <- res_saem@phi_names
  x_candidates_names <- res_saem@x_candidates_names
  x_forced_names <- res_saem@x_forced_names

  if (!is.null(phi)) {
    if (!(phi %in% phi_names)) {
      stop(
        sprintf("Parameter '%s' is not present in phi_names of the SAEM object.", phi) # nolint: line_length_linter
      )
    }

    phi_idx <- which(phi_names == phi)

    if (component %in% c("coef_phi_sel", "variance_phi_sel")) {
      if (is.null(phi_to_select_idx) || !(phi_idx %in% phi_to_select_idx)) {
        stop(
          sprintf(
            paste0(
              "Parameter %d is not subject to selection, thus not compatible",
              " with component '%s'."
            ),
            phi,
            component
          )
        )
      }
    } else if (component %in% c("coef_phi_non_sel", "variance_phi_non_sel")) {
      if (is.null(phi_not_to_select_idx) ||
        !(phi_idx %in% phi_not_to_select_idx)) {
        stop(
          sprintf(
            paste0(
              "Parameter %d is subject to selection, thus not compatible ",
              "with component '%s'."
            ),
            phi,
            component
          )
        )
      }
    }
  }
  # --- Handle case: residual variance (sigma2) ---
  if (component == "sigma2") {
    df <- data.frame(
      iteration = seq_along(res_saem@sigma2),
      value = res_saem@sigma2
    )

    n_iter <- length(unique(df$iteration))
    step <- max(1, floor(n_iter / max_ticks))
    breaks <- seq(1, n_iter, by = step)

    g <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = value)) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = "Iteration", y = "Estimated value",
        title = "Evolution of residual variance (sigma2)"
      ) +
      ggplot2::scale_x_continuous(breaks = breaks) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1
      ))

    # --- Handle case: beta/gamma matrices ---
  } else if (component %in% c(
    "coef_phi_sel", "coef_phi_non_sel",
    "variance_phi_sel", "variance_phi_non_sel"
  )) {
    list_est <- switch(component,
      coef_phi_sel = res_saem@beta_to_select,
      coef_phi_non_sel = res_saem@beta_not_to_select,
      variance_phi_sel = res_saem@gamma_to_select,
      variance_phi_non_sel = res_saem@gamma_not_to_select
    )

    if (length(list_est) == 0 || all(sapply(list_est, is.null))) {
      message(
        sprintf(
          paste0(
            "No data available for component '%s'. ",
            "Plot cannot be generated."
          ),
          component
        )
      )
      return(invisible(NULL))
    }

    df <- do.call(rbind, lapply(seq_along(list_est), function(iter) {
      mat <- list_est[[iter]]
      data.frame(
        iteration = iter,
        i = rep(seq_len(nrow(mat)), ncol(mat)),
        j = rep(seq_len(ncol(mat)), each = nrow(mat)),
        value = as.vector(mat)
      )
    }))

    # ---- Filter according to phi if provided ----
    if (!is.null(phi)) {
      if (component %in% c("coef_phi_sel", "coef_phi_non_sel")) {
        # coefficients -> phi correspond aux colonnes
        df <- base::subset(df, j == phi_idx)
      } else if (component %in% c("variance_phi_sel", "variance_phi_non_sel")) {
        # covariance -> garder lignes ou colonnes liées à phi
        df <- base::subset(df, i == phi_idx | j == phi_idx)
      }
    }


    n_iter <- length(unique(df$iteration))
    step <- max(1, floor(n_iter / max_ticks))
    breaks <- seq(1, n_iter, by = step)

    # --- Create facet labels according to component ---
    make_facet_label <- function(comp, i, j) {
      # intercept
      var_label <- if (i == 1) {
        "intercept"
      } else if (comp == "coef_phi_sel") {
        x_candidates_names[i - 1]
      } else if (comp == "coef_phi_non_sel") {
        x_forced_names[i - 1]
      } else {
        NA
      }

      if (comp == "coef_phi_sel") {
        paste0(
          phi_names[phi_to_select_idx[j]],
          ", ",
          var_label
        )
      } else if (comp == "coef_phi_non_sel") {
        paste0(
          phi_names[phi_not_to_select_idx[j]],
          ", ",
          var_label
        )
      } else if (comp == "variance_phi_sel") {
        paste0(
          "cov(",
          phi_names[phi_to_select_idx[i]],
          ", ",
          phi_names[phi_to_select_idx[j]],
          ")"
        )
      } else if (comp == "variance_phi_non_sel") {
        paste0(
          "cov(",
          phi_names[phi_not_to_select_idx[i]],
          ", ",
          phi_names[phi_not_to_select_idx[j]],
          ")"
        )
      } else {
        NA
      }
    }



    df$facet_label <- mapply(make_facet_label, component, df$i, df$j)
    possible_components <- unique(df$facet_label)

    # --- Default / random / top:n selection ---
    if (missing(sel_components)) {
      sel_components <- utils::head(possible_components, 16)
    } else if (length(sel_components) == 1 && sel_components == "random") {
      sel_components <- sample(
        possible_components,
        min(16, length(possible_components))
      )
    } else if (length(sel_components) == 1 && grepl(
      "^top:[0-9]+$",
      sel_components
    )) {
      n <- as.numeric(sub("top:", "", sel_components))

      if (!is.null(phi)) {
        df_col <- base::subset(df, j == phi_idx)
      } else {
        df_col <- df
      }

      var_df <- stats::aggregate(value ~ facet_label, df_col, var)

      top_sel <- utils::head(var_df[order(-var_df$value), "facet_label"], n)
      sel_components <- top_sel
    }

    invalid <- setdiff(sel_components, possible_components)
    if (length(invalid) > 0) {
      stop(
        "The following components do not exist in the matrix: ",
        paste(invalid, collapse = ", ")
      )
    }

    if (length(sel_components) > 16) {
      sel_components <- sel_components[1:16]
      warning(
        paste0(
          "Only the first 16 components of 'sel_components' are displayed."
        )
      )
    }

    df <- base::subset(df, facet_label %in% sel_components)

    # --- Build plot ---
    g <- ggplot2::ggplot(df, ggplot2::aes(
      x = iteration, y = value,
      color = facet_label
    )) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~facet_label, scales = "free_y", ncol = 4) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = "Iteration", y = "Estimated value" # ,
      ) +
      ggplot2::scale_x_continuous(breaks = breaks) +
      ggplot2::scale_color_discrete(name = "component") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1
      ))
  } else {
    stop(
      sprintf(
        "'component' must be one of: %s.",
        paste(
          shQuote(c(
            "sigma2",
            "coef_phi_sel",
            "coef_phi_non_sel",
            "variance_phi_sel",
            "variance_phi_non_sel"
          )),
          collapse = ", "
        )
      )
    )
  }

  print(g)
}

#' Internal helper to generate SAEMVS diagnostic plots
#'
#' Generates either the selection criterion plot or the coefficient path
#' plot for a given parameter from a \code{saemvsResults} object.
#'
#' @param res_saemvs A \code{saemvsResults} object.
#' @param type Character string specifying which plot to generate.
#'   Must be one of \code{"criterion"} or \code{"coefficients"}.
#' @param param Parameter for which to display the coefficient path when
#'   \code{type = "coefficients"}. Can be either:
#'   \itemize{
#'     \item an integer index,
#'     \item a character string matching a name in
#'       \code{res_saemvs@phi_names}.
#'   }
#'   Must not be \code{NULL} when \code{type = "coefficients"}.
#'
#' @return A single \pkg{ggplot2} object:
#' \itemize{
#'   \item the selection criterion plot if \code{type = "criterion"},
#'   \item the coefficient path plot for the selected parameter if
#'     \code{type = "coefficients"}.
#' }
#'
#' @details
#' This function is an internal backend used by the \code{plot} method for
#' \code{saemvsResults}. It does not handle multiple parameters; looping and
#' layout (e.g., for \code{param = "all"}) are managed at a higher level.
#'
#' @keywords internal
#' @noRd
.prepare_grid_plot_backend <- function(res_saemvs, type, param = NULL) {
  type <- match.arg(type, c("criterion", "coefficients"))

  # --- Shared extraction ---
  ebic <- res_saemvs@criterion_values
  map_to_unique_support <- res_saemvs@support_mapping
  nu0_grid <- res_saemvs@spike_values_grid
  pen <- res_saemvs@criterion

  data_crit <- data.frame(
    nu0 = nu0_grid,
    crit = ebic[map_to_unique_support]
  )

  x_min <- log(data_crit$nu0[which.min(data_crit$crit)])

  # ============================================================
  # 1. CRITERION
  # ============================================================
  if (type == "criterion") {
    return(
      ggplot2::ggplot(data_crit, ggplot2::aes(x = log(nu0), y = crit)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::xlab(expression(paste("log(", nu[0], ")"))) +
        ggplot2::ylab(pen) +
        ggplot2::ggtitle(pen) +
        ggplot2::geom_vline(
          xintercept = x_min,
          color = "red",
          linetype = "dashed"
        )
    )
  }

  # ============================================================
  # 2. COEFFICIENTS
  # ============================================================

  if (is.null(param)) {
    stop("Argument 'param' must be provided when type = 'coefficients'")
  }

  phi_names <- res_saemvs@phi_names[res_saemvs@phi_to_select_idx]

  m <- if (is.character(param)) match(param, phi_names) else param

  if (is.na(m) || m < 1 || m > length(phi_names)) {
    stop(sprintf("Parameter '%s' not found.", param))
  }

  support <- res_saemvs@support
  p <- dim(support[[1]])[1] - 1
  q <- dim(support[[1]])[2]

  threshold <- matrix(simplify2array(res_saemvs@thresholds), nrow = q)
  beta <- simplify2array(res_saemvs@beta_map)

  # nombre de covariables forcées
  n_forced <- length(res_saemvs@x_forced_names)

  # ============================================================
  # Construction tidy des données
  # ============================================================

  # --- seuil bas
  df_low <- data.frame(
    nu0 = nu0_grid,
    value = -threshold[m, ],
    var = "threshold_low",
    type = "threshold"
  )

  # --- coefficients beta
  df_beta <- as.data.frame(t(beta[, m, ])) # nb_nu0 x p
  colnames(df_beta) <- paste0("beta_", seq_len(p))
  df_beta$nu0 <- nu0_grid

  df_beta <- tidyr::pivot_longer(
    df_beta,
    cols = starts_with("beta_"),
    names_to = "var",
    values_to = "value"
  )

  # index numérique des beta
  df_beta$idx <- as.integer(sub("beta_", "", df_beta$var))

  # type des coefficients
  df_beta$type <- ifelse(
    df_beta$idx <= n_forced,
    "forced",
    "selected"
  )

  # --- seuil haut
  df_high <- data.frame(
    nu0 = nu0_grid,
    value = threshold[m, ],
    var = "threshold_high",
    type = "threshold"
  )

  # --- concat
  df <- dplyr::bind_rows(df_low, df_beta, df_high)

  # ordre des facteurs pour une légende propre
  df$type <- factor(df$type, levels = c("threshold", "forced", "selected"))

  # ============================================================
  # Plot
  # ============================================================

  return(
    ggplot2::ggplot(df, ggplot2::aes(
      x = log(nu0),
      y = value,
      group = var,
      color = type
    )) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::scale_color_manual(
        values = c(
          threshold = "red",
          forced = "grey60",
          selected = "black"
        ),
        labels = c(
          threshold = "threshold",
          forced = "forced covariates",
          selected = "candidate covariates"
        )
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(color = NULL) +
      ggplot2::xlab(expression(paste("log(", nu[0], ")"))) +
      ggplot2::ylab(expression(hat(beta))) +
      ggplot2::ggtitle(paste("Parameter", phi_names[m])) +
      ggplot2::geom_vline(
        xintercept = x_min,
        color = "red",
        linetype = "dashed"
      )
  )
}

#' Plot SAEM convergence diagnostics
#'
#' Produces convergence plots for a SAEM fit.
#'
#' This function is intended for objects returned by \code{test_saemvs}.
#' It visualizes the evolution of selected components (residual variance,
#' regression coefficients, or covariance matrices) across iterations of
#' the SAEM algorithm. It is useful for diagnosing convergence issues
#' in the SAEM algorithm.
#'
#' @param x A \code{saemResults} object, typically obtained from
#'   \code{\link{test_saemvs}}.
#' @param type Character string specifying which element to plot. One of:
#'   \itemize{
#'     \item `"sigma2"`: residual variance
#'     \item `"coef_phi_sel"`: regression coefficients for parameters subject
#'           to selection
#'     \item `"coef_phi_non_sel"`: regression coefficients for parameters not
#'           subject to selection
#'     \item `"variance_phi_sel"`: covariance matrix of parameters subject
#'           to selection
#'     \item `"variance_phi_non_sel"`: covariance matrix of parameters not
#'           subject to selection
#'   }
#' @param ... Further arguments passed to the internal plotting function:
#'   \itemize{
#'     \item `sel_components` Optional character vector specifying which
#'           components to display. Can be `"random"` or `"top:n"` to select
#'           up to 16 components. Defaults to the first 16.
#'     \item `phi` Optional character or integer specifying a parameter to
#'           focus on.
#'   }
#'
#' @return A \pkg{ggplot2} object showing the evolution of the selected
#'   component across iterations.
#'
#' @name plot-saemResults
#' @rdname plot-saemResults
#' @aliases plot,saemResults,missing-method
#' @exportMethod plot
setMethod(
  "plot",
  signature(x = "saemResults", y = "missing"),
  function(x, component = "coef_phi_sel", ...) {
    .convergence_plot_backend(
      res_saem = x,
      component = component,
      ...
    )
  }
)

#' Plot SAEMVS selection path
#'
#' Generates diagnostic plots for a SAEMVS fit, showing either the selection
#' criterion or the evolution of regression coefficients across the spike
#' variance grid.
#'
#' @param x A \code{saemvsResults} object, typically returned by
#'   \code{\link{saemvs}}.
#'
#' @param type Character string specifying which plot to generate.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"criterion"}: plots the selection criterion (BIC or e-BIC)
#'       along the spike variance grid. The plot highlights the grid value
#'       corresponding to the best model (minimum criterion value).
#'     \item \code{"coefficients"}: plots the regression coefficients for a
#'       given parameter along the spike variance grid, allowing visualization
#'       of how coefficients change with the spike parameter.
#'   }
#'
#' @param param Parameter for which to display the coefficient path when
#'   \code{type = "coefficients"}. Can be either:
#'   \itemize{
#'     \item an integer index,
#'     \item a character string matching a name in \code{x@phi_names},
#'     \item \code{"all"} to display all parameters in a grid of plots.
#'   }
#'   Ignored when \code{type = "criterion"}.
#'
#' @return
#' \itemize{
#'   \item If \code{type = "criterion"}: a single \pkg{ggplot2} object.
#'   \item If \code{type = "coefficients"} and \code{param} is a single
#'     parameter: a single \pkg{ggplot2} object.
#'   \item If \code{type = "coefficients"} and \code{param = "all"}:
#'     the plots are displayed in a grid layout, and an (invisible) list
#'     of \pkg{ggplot2} objects is returned.
#' }
#'
#' @examples
#' \dontrun{
#' # Fit a SAEMVS model
#' res <- saemvs(...)
#'
#' # Plot the selection criterion
#' plot(res, type = "criterion")
#'
#' # Plot coefficients for a specific parameter (by index)
#' plot(res, type = "coefficients", param = 1)
#'
#' # Plot coefficients for a specific parameter (by name)
#' plot(res, type = "coefficients", param = "CL")
#'
#' # Plot coefficients for all parameters
#' plot(res, type = "coefficients", param = "all")
#' }
#'
#' @name plot-saemvsResults
#' @rdname plot-saemvsResults
#' @aliases plot,saemvsResults,missing-method
#' @exportMethod plot
setMethod(
  "plot",
  signature(x = "saemvsResults", y = "missing"),
  function(x, type = c("criterion", "coefficients"), param = NULL) {
    type <- match.arg(type)
    phi_names <- x@phi_names[x@phi_to_select_idx]

    if (type == "coefficients" && identical(param, "all")) {
      plots <- lapply(seq_along(phi_names), function(m) {
        .prepare_grid_plot_backend(x, type = "coefficients", param = m)
      })

      n <- length(plots)

      ncol <- ceiling(sqrt(n))
      nrow <- ceiling(n / ncol)

      gridExtra::grid.arrange(grobs = plots, ncol = ncol, nrow = nrow)

      return(invisible(plots))
    }

    .prepare_grid_plot_backend(x, type = type, param = param)
  }
)
