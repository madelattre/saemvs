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
  res_saem, component, sel_components, phi = NULL
) {
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
      ggplot2::geom_point() +
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

# Internal helper: prepare grid plots for SAEMVS
#'
#' @param res_saemvs A \code{saemvsResults} object.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{reg_plot}: list of \pkg{ggplot2} objects, one per parameter
#'  block
#'   \item \code{ebic_plot}: \pkg{ggplot2} object showing the selection
#'  criterion across the grid
#' }
#'
#' @keywords internal
#' @noRd
.prepare_grid_plot_backend <- function(res_saemvs) {
  # --- Extract key results from saemvsResults ---
  ebic <- res_saemvs@criterion_values
  support <- res_saemvs@support
  map_to_unique_support <- res_saemvs@support_mapping
  nu0_grid <- res_saemvs@spike_values_grid
  nb_nu0 <- length(nu0_grid)
  pen <- res_saemvs@criterion
  phi_names <- res_saemvs@phi_names[res_saemvs@phi_to_select_idx]

  p <- dim(support[[1]])[1] - 1
  q <- dim(support[[1]])[2]

  threshold <- matrix(simplify2array(res_saemvs@thresholds), nrow = q)
  beta <- simplify2array(res_saemvs@beta_map)

  # --- Plot criterion values vs log(nu0) ---
  data2 <- data.frame(nu0_grid = nu0_grid, crit = ebic[map_to_unique_support])
  x_min <- log(data2$nu0_grid[which.min(data2$crit)])

  g2 <- ggplot2::ggplot(data2, ggplot2::aes(x = log(nu0_grid), y = crit)) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::xlab(expression(paste("log(", nu[0], " ) "))) +
    ggplot2::ylab(paste(pen)) +
    ggplot2::ggtitle(paste(pen)) +
    ggplot2::geom_vline(
      xintercept = x_min, color = "red",
      linetype = "dashed"
    )

  # --- Plot regression coefficients and thresholds ---
  id_var <- rep(c(1:(p + 2)), nb_nu0)
  g <- list()

  for (m in 1:q) {
    y <- rbind(-threshold[m, ], beta[, m, ], threshold[m, ])
    y <- c(y)

    x <- rep(nu0_grid, each = p + 2)

    data <- data.frame(id_var, y, x)

    g[[m]] <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = log(x), y = y, group = id_var,
        color = as.factor(id_var)
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::scale_color_manual(values = c("red", rep("black", p), "red")) +
      ggplot2::theme_bw() +
      ggplot2::xlab(expression(paste("log(", nu[0], " ) "))) +
      ggplot2::ylab(expression(hat(beta))) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(paste("Parameter ", phi_names[m])) +
      ggplot2::geom_vline(
        xintercept = x_min, color = "red",
        linetype = "dashed"
      )
  }

  return(list(reg_plot = g, ebic_plot = g2)) # nolint: return_linter
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
  function(x, type = "coef_phi_sel", ...) {
    .convergence_plot_backend(
      res_saem = x,
      component = type,
      ...
    )
  }
)

#' Plot SAEMVS selection path
#'
#' Generates diagnostic plots for a SAEMVS fit, showing either the selection
#' criterion or the evolution of regression coefficients across the spike
#'  variance grid.
#'
#' @param x A \code{saemvsResults} object, typically returned by
#'  \code{\link{saemvs}}.
#' @param type Character string specifying which plot to generate.
#' Must be one of:
#'   \itemize{
#'     \item \code{"criterion"}: plots the selection criterion (BIC or e-BIC)
#'           along the spike variance grid. The plot highlights the grid value
#'           corresponding to the best model (minimum criterion value).
#'     \item \code{"coefficients"}: plots the regression coefficients for each
#'           parameter along the spike variance grid, allowing visualization
#'           of how coefficients change with the spike parameter.
#'   }
#'
#' @return Either a \pkg{ggplot2} object (\code{"criterion"}) or a list of
#'   \pkg{ggplot2} objects (\code{"coefficients"}), depending on the
#'  \code{type}.
#'
#' @examples
#' \dontrun{
#' # Fit a SAEMVS model
#' res <- saemvs(...)
#' # Plot the selection criterion along the spike variance grid
#' plot(res, type = "criterion")
#' # Plot regression coefficients along the grid for first model parameter
#' plot(res, type = "coefficients")[[1]]
#' # Plot regression coefficients along the grid for second model parameter
#' plot(res, type = "coefficients")[[2]]
#' ...
#' }
#'
#' @name plot-saemvsResults
#' @rdname plot-saemvsResults
#' @aliases plot,saemvsResults,missing-method
#' @exportMethod plot
setMethod(
  "plot",
  signature(x = "saemvsResults", y = "missing"),
  function(x, type = c("criterion", "coefficients")) {
    type <- match.arg(type)
    plots <- .prepare_grid_plot_backend(x)

    if (type == "criterion") {
      plots$ebic_plot
    } else {
      plots$reg_plot
    }
  }
)
