#' Convergence Diagnostic Plot for SAEM Estimation
#'
#' The \code{convergence_plot} function provides visual diagnostics
#' of convergence for specific components estimated in a
#' \linkS4class{saemResults} object. It can be used to monitor the
#' evolution of parameter estimates across SAEM iterations.
#'
#' Supported components are:
#' \itemize{
#'   \item \code{"sigma2"}: evolution of the residual variance.
#'   \item \code{"beta_s"}: estimated coefficients for the variables
#'         subject to selection.
#'   \item \code{"beta_ns"}: estimated coefficients for the variables
#'         not subject to selection.
#'   \item \code{"gamma_s"}: estimated random effects covariance
#'         components subject to selection.
#'   \item \code{"gamma_ns"}: estimated random effects covariance
#'         components not subject to selection.
#' }
#'
#' @param res_saem An object of class \linkS4class{saemResults}.
#' @param component A character string indicating which component to plot.
#'   Must be one of \code{"sigma2"}, \code{"beta_s"}, \code{"beta_ns"},
#'   \code{"gamma_s"}, \code{"gamma_ns"}.
#' @param sel_components A character vector specifying which elements of
#'   the selected matrix (\code{"beta_*"} or \code{"gamma_*"}) to track.
#'   Elements must be of the form \code{"(i,j)"}, where \code{i} and
#'   \code{j} are row and column indices. At most 16 components are shown.
#'   Ignored if \code{component = "sigma2"}.
#'
#' @return A \pkg{ggplot2} object (printed by default), showing the
#'   evolution of the requested parameter(s) across SAEM iterations.
#'
#' @examples
#' \dontrun{
#' # Example: Plot sigma^2 evolution
#' convergence_plot(res, component = "sigma2")
#'
#' # Example: Plot selected beta coefficients
#' convergence_plot(res, component = "beta_s", sel_components = c("(1,1)", "(2,1)"))
#' }
#'
#' @seealso \linkS4class{saemResults}
#'
#' @export
setGeneric(
  "convergence_plot",
  function(res_saem, component, sel_components) {
    standardGeneric("convergence_plot")
  }
)

#' @rdname convergence_plot
#' @exportMethod convergence_plot
setMethod(
  "convergence_plot",
  signature(
    res_saem = "saemResults", component = "character",
    sel_components = "character"
  ),
  function(res_saem, component, sel_components) {
    # --- Handle case: residual variance (sigma2) ---
    if (component == "sigma2") {
      df <- data.frame(
        iteration = seq_along(res_saem@sigma2),
        value = res_saem@sigma2
      )

      g <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = value)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::labs(
          x = "Iteration", y = "Estimated value",
          title = "Evolution of the selected matrix component"
        )

      # --- Handle case: beta/gamma matrices ---
    } else if (component %in% c("beta_s", "beta_ns", "gamma_s", "gamma_ns")) {
      # Select appropriate list of matrices based on requested component
      switch(component,
        beta_s = {
          list_est <- res_saem@beta_to_select
        },
        beta_ns = {
          list_est <- res_saem@beta_not_to_select
        },
        gamma_s = {
          list_est <- res_saem@gamma_to_select
        },
        gamma_ns = {
          list_est <- res_saem@gamma_not_to_select
        }
      )

      # Build long-format dataframe of all iterations
      df <- do.call(rbind, lapply(seq_along(list_est), function(iter) {
        mat <- list_est[[iter]]
        data.frame(
          iteration = iter,
          i = rep(seq_len(nrow(mat)), ncol(mat)),
          j = rep(seq_len(ncol(mat)), each = nrow(mat)),
          value = as.vector(mat)
        )
      }))
      df$component <- paste0("(", df$i, ",", df$j, ")")

      # --- Validate sel_components ---
      valid_format <- grepl("^\\([0-9]+,[0-9]+\\)$", sel_components)
      if (any(!valid_format)) {
        stop("Some elements in 'sel_components' do not match the required format '(i,j)'.")
      }

      possible_components <- unique(df$component)
      invalid <- setdiff(sel_components, possible_components)
      if (length(invalid) > 0) {
        stop(
          paste(
            "The following components do not exist in the matrix:",
            paste(invalid, collapse = ", ")
          )
        )
      }

      # Restrict to 16 components max for readability
      if (length(sel_components) > 16) {
        sel_components <- sel_components[1:16]
        warning("Only the first 16 components of 'sel_components' have been considered.")
      }

      # Keep only selected components
      df <- subset(df, component %in% sel_components)

      # --- Build plot ---
      g <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = iteration, y = value, color = component)
      ) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~component, scales = "free_y") +
        ggplot2::theme_bw() +
        ggplot2::labs(
          x = "Iteration", y = "Estimated value",
          title = paste("Evolution of each matrix component in", component)
        )

      # --- Handle invalid input ---
    } else {
      stop(
        paste0(
          "'component' must be 'beta_s', 'beta_ns', 'gamma_s', ",
          "'gamma_ns' or 'sigma2'."
        )
      )
    }

    # Print plot by default
    print(g)
  }
)


#' Prepare Grid Plots for SAEMVS Results
#'
#' The \code{prepare_grid_plot} function generates diagnostic plots
#' to visualize the SAEMVS selection process across the grid of spike
#' variances (\eqn{\nu_0}).
#'
#' Two types of plots are returned:
#' \itemize{
#'   \item A criterion plot (\code{ebic_plot}): values of the selection
#'   criterion (BIC or e-BIC) across the grid of \eqn{\nu_0}, with the
#'   optimal point highlighted.
#'   \item A list of regression coefficient plots (\code{reg_plot}):
#'   for each parameter block, the estimated regression coefficients
#'   (\eqn{\hat{\beta}}) are plotted across the grid of \eqn{\nu_0},
#'   together with their threshold bounds.
#' }
#'
#' @param res_saemvs An object of class \linkS4class{saemvsResults},
#'   typically obtained from a call to \code{\link{saemvs}}.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{reg_plot}: A list of \pkg{ggplot2} objects, one for
#'   each parameter block, showing regression estimates versus
#'   spike variance.
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
    ebic <- res_saemvs@criterion_values # Criterion values (BIC/e-BIC)
    threshold <- as.matrix(simplify2array(res_saemvs@thresholds)) # Thresholds for selection
    beta <- simplify2array(res_saemvs@beta_map) # MAP regression estimates
    support <- res_saemvs@support # Support sets
    map_to_unique_support <- res_saemvs@support_mapping
    nu0_grid <- res_saemvs@spike_values_grid # Spike variances
    nb_nu0 <- length(nu0_grid)
    pen <- res_saemvs@criterion # Selection criterion used

    # Dimensions of the support matrix
    p <- dim(support[[1]])[1] - 1 # Number of covariates (excluding intercept)
    q <- dim(support[[1]])[2] # Number of parameter blocks

    # --- Plot criterion values vs log(nu0) ---
    data2 <- data.frame(nu0_grid = nu0_grid, crit = ebic[map_to_unique_support])
    x_min <- log(data2$nu0_grid[which.min(data2$crit)]) # Optimal nu0

    g2 <- ggplot2::ggplot(data2, ggplot2::aes(x = log(nu0_grid), y = crit)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::xlab(expression(paste("log(", nu[0], " ) "))) +
      ggplot2::ylab(paste(pen)) +
      ggplot2::ggtitle(paste(pen)) +
      ggplot2::geom_vline(xintercept = x_min, color = "red", linetype = "dashed")

    # --- Plot regression coefficients and thresholds ---
    id_var <- rep(c(1:(p + 2)), nb_nu0) # Identifiers for intercept, covariates, and thresholds
    g <- list()

    for (m in 1:q) {
      # Stack lower threshold, beta estimates, and upper threshold
      y <- rbind(-threshold[m, ], beta[, m, ], threshold[m, ])
      y <- c(y)

      # Replicate nu0 values for each variable
      x <- rep(nu0_grid, each = p + 2)

      # Build dataframe for ggplot
      data <- data.frame(id_var, y, x)

      # Plot for parameter block m
      g[[m]] <- ggplot2::ggplot(
        data,
        ggplot2::aes(x = log(x), y = y, group = id_var, color = as.factor(id_var))
      ) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_color_manual(values = c("red", rep("black", p), "red")) +
        ggplot2::theme_bw() +
        ggplot2::xlab(expression(paste("log(", nu[0], " ) "))) +
        ggplot2::ylab(expression(hat(beta))) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(bquote("Parameter " ~ varphi[.(m)])) +
        ggplot2::geom_vline(xintercept = x_min, color = "red", linetype = "dashed")
    }

    # --- Return both plots ---
    return(list(reg_plot = g, ebic_plot = g2))
  }
)
