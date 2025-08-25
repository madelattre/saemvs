#' @export
setGeneric(
  "convergence_plot",
  function(res_saem, component, sel_components) {
    standardGeneric("convergence_plot")
  }
)

#' @exportMethod convergence_plot
setMethod(
  "convergence_plot",
  signature(
    res_saem = "saemResults", component = "character",
    sel_components = "character"
  ),
  function(res_saem, component, sel_components) {
    # component = c("beta_s", "beta_ns", "gamma_s", "gamma_ns", "sigma2")




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
    } else if (component %in% c("beta_s", "beta_ns", "gamma_s", "gamma_ns")) {
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
        },
      )

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

      ## Check format of sel_components
      valid_format <- grepl("^\\([0-9]+,[0-9]+\\)$", sel_components)
      if (any(!valid_format)) {
        stop("Some elements in 'sel_components' do not match the required format '(i,j)'.")
      }

      ## Check if components exist in the matrix
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

      ## Limit to 16 components
      if (length(sel_components) > 16) {
        sel_components <- sel_components[1:16]
        warning("Only the first 16 components of 'sel_components' have been considered.")
      }

      df <- subset(df, component %in% sel_components)



      g <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = iteration, y = value, color = component)
      ) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~component, scales = "free_y") +
        # ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::labs(
          x = "Iteration", y = "Estimated value",
          title = paste("Evolution of each matrix component in ", component)
        )
    } else {
      stop(
        paste0(
          "'component' must be 'beta_s', 'beta_ns', 'gamma_s', ",
          "'gamma_ns' or 'sigma2'."
        )
      )
    }

    print(g)
  }
)

#' @export
setGeneric(
  "prepare_grid_plot",
  function(res_saemvs) standardGeneric("prepare_grid_plot")
)

#' @exportMethod prepare_grid_plot
setMethod(
  "prepare_grid_plot",
  signature(
    res_saemvs = "saemvsResults"
  ),
  function(res_saemvs) {
    ebic <- res_saemvs@criterion_values
    threshold <- simplify2array(res_saemvs@thresholds)
    beta <- simplify2array(res_saemvs@beta_map)
    support <- res_saemvs@support
    map_to_unique_support <- res_saemvs@support_mapping
    nu0_grid <- res_saemvs@spike_values_grid
    nb_nu0 <- length(nu0_grid)
    pen <- res_saemvs@criterion

    p <- dim(support[[1]])[1] - 1
    q <- dim(support[[1]])[2]


    data2 <- data.frame(nu0_grid = nu0_grid, crit = ebic[map_to_unique_support])
    x_min <- log(data2$nu0_grid[which.min(data2$crit)])

    g2 <- ggplot2::ggplot(data2, ggplot2::aes(x = log(nu0_grid), y = crit)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::xlab(expression(paste("log(", nu[0], " ) "))) +
      ggplot2::ylab(paste(pen)) +
      ggplot2::ggtitle(paste(pen)) +
      ggplot2::geom_vline(xintercept = x_min, color = "red", linetype = "dashed")

    id_var <- rep(c(1:(p + 2)), nb_nu0)

    g <- list()

    for (m in 1:q) {
      y <- rbind(-threshold[m, ], beta[, m, ], threshold[m, ])
      y <- c(y)
      x <- rep(nu0_grid, each = p + 2)
      data <- data.frame(id_var, y, x)
      g[[m]] <- ggplot2::ggplot(
        data,
        ggplot2::aes(x = log(x), y = y, group = id_var, color = as.factor(id_var))
      ) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::scale_color_manual(
          values = c("red", rep("black", p), "red")
        ) +
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
