#' Convergence Diagnostic Plot for SAEM Estimation
#'
#' The \code{convergence_plot} function provides visual diagnostics of
#' convergence for specific components estimated in a \linkS4class{saemResults}
#' object.
#'
#' Supported components are:
#' \itemize{
#'   \item \code{"sigma2"}: evolution of the residual variance.
#'   \item \code{"beta_s"}: estimated coefficients for the variables subject to
#'   selection.
#'   \item \code{"beta_ns"}: estimated coefficients for the variables not
#'   subject to selection.
#'   \item \code{"gamma_s"}: estimated random effects covariance components
#'   subject to selection.
#'   \item \code{"gamma_ns"}: estimated random effects covariance components not
#'   subject to selection.
#' }
#'
#' @param res_saem An object of class \linkS4class{saemResults}.
#' @param component A character string indicating which component to plot.
#'   Must be one of \code{"sigma2"}, \code{"beta_s"}, \code{"beta_ns"},
#'   \code{"gamma_s"}, \code{"gamma_ns"}.
#' @param sel_components A character vector specifying which elements of
#'   the selected matrix to track. Elements must be of the form \code{"(i,j)"}.
#'   At most 16 components are shown. Special values:
#'   \itemize{
#'     \item if missing, the first 16 components are chosen,
#'     \item if \code{"random"}, 16 components are chosen at random,
#'     \item if \code{"top:n"}, the \code{n} most variable components are shown.
#'   }
#'   Ignored if \code{component = "sigma2"}.
#'
#' @return A \pkg{ggplot2} object showing the evolution of the
#'   requested parameter(s) across SAEM iterations.
#'
#' @export
setGeneric(
  "convergence_plot",
  function(res_saem, component, sel_components, ...) {
    standardGeneric("convergence_plot")
  }
)

#' @rdname convergence_plot
#' @exportMethod convergence_plot
setMethod(
  "convergence_plot",
  signature(
    res_saem = "saemResults", component = "character", sel_components = "ANY"
  ),
  function(res_saem, component, sel_components, phi = NULL) {
    max_ticks <- 4

    phi_to_select_idx <- res_saem@phi_to_select_idx
    phi_not_to_select_idx <- res_saem@phi_not_to_select_idx

    if (!is.null(phi)) {
      if (component %in% c("coef_phi_sel", "variance_phi_sel")) {
        if (is.null(phi_to_select_idx) || !(phi %in% phi_to_select_idx)) {
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
          !(phi %in% phi_not_to_select_idx)) {
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
        if (comp == "coef_phi_sel") {
          paste0("\u03C6", phi_to_select_idx[j], ", var. ", i - 1)
        } else if (comp == "coef_phi_non_sel") {
          paste0("\u03C6", phi_not_to_select_idx[j], ", var. ", i - 1)
        } else if (comp == "variance_phi_sel") {
          paste0(
            "cov(", "\u03C6", phi_to_select_idx[i],
            ", \u03C6", phi_to_select_idx[j], ")"
          )
        } else if (comp == "variance_phi_non_sel") {
          paste0(
            "cov(", "\u03C6", phi_not_to_select_idx[i],
            ", \u03C6", phi_not_to_select_idx[j], ")"
          )
        } else {
          NA
        }
      }

      df$facet_label <- mapply(make_facet_label, component, df$i, df$j)
      possible_components <- unique(df$facet_label)

      # --- Default / random / top:n selection ---
      if (missing(sel_components)) {
        sel_components <- head(possible_components, 16)
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
          df_col <- subset(df, j == phi)
        } else {
          df_col <- df
        }

        var_df <- aggregate(value ~ facet_label, df_col, var)

        top_sel <- head(var_df[order(-var_df$value), "facet_label"], n)
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

      df <- subset(df, facet_label %in% sel_components)

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
        paste0(
          "'component' must be one of: 'sigma2', 'beta_s', 'beta_ns', ",
          "'gamma_s', 'gamma_ns'."
        )
      )
    }

    print(g)
  }
)


#' Prepare Grid Plots for SAEMVS Results
#'
#' The \code{prepare_grid_plot} function generates diagnostic plots to visualize
#' the SAEMVS selection process across the grid of spike variances
#' (\eqn{\nu_0}).
#'
#' Two types of plots are returned:
#' \itemize{
#'   \item A criterion plot (\code{ebic_plot}): values of the selection
#'   criterion (BIC or e-BIC) across the grid of \eqn{\nu_0}, with the optimal
#'   point highlighted.
#'   \item A list of regression coefficient plots (\code{reg_plot}):
#'   for each parameter block, the estimated regression coefficients
#'   (\eqn{\hat{\beta}}) are plotted across the grid of \eqn{\nu_0}, together
#'   with their threshold bounds.
#' }
#'
#' @param res_saemvs An object of class \linkS4class{saemvsResults},
#'   typically obtained from a call to \code{\link{saemvs}}.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{reg_plot}: A list of \pkg{ggplot2} objects, one for each
#'   parameter block, showing regression estimates versus spike variance.
#'   \item \code{ebic_plot}: A single \pkg{ggplot2} object showing the
#'   selection criterion (BIC/e-BIC) across the grid of \eqn{\nu_0}.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming res is an object of class saemvsResults
#' plots <- prepare_grid_plot(res)
#' plots$ebic_plot # Plot criterion values
#' plots$reg_plot[[1]] # Plot first parameter block
#' }
#'
#' @seealso \code{\link{saemvs}}, \linkS4class{saemvsResults}
#'
#' @export
setGeneric(
  "prepare_grid_plot",
  function(res_saemvs) standardGeneric("prepare_grid_plot")
)

#' @rdname prepare_grid_plot
#' @exportMethod prepare_grid_plot
setMethod(
  "prepare_grid_plot",
  signature(
    res_saemvs = "saemvsResults"
  ),
  function(res_saemvs) {
    # --- Extract key results from saemvsResults ---
    ebic <- res_saemvs@criterion_values
    support <- res_saemvs@support
    map_to_unique_support <- res_saemvs@support_mapping
    nu0_grid <- res_saemvs@spike_values_grid
    nb_nu0 <- length(nu0_grid)
    pen <- res_saemvs@criterion

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
        ggplot2::ggtitle(bquote("Parameter " ~ varphi[.(m)])) +
        ggplot2::geom_vline(
          xintercept = x_min, color = "red",
          linetype = "dashed"
        )
    }

    return(list(reg_plot = g, ebic_plot = g2))
  }
)
