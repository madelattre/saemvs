#' @export
setGeneric(
  "prepare_grid_plot",
  function(res_saemvs) standardGeneric("prepare_grid_plot")
)

#' @exportMethod prepare_grid_plot
setMethod(
  "prepare_grid_plot",
  signature(
    res_saemvs = "resSAEMVS"
  ),
  function(res_saemvs) {
    ebic <- res_saemvs@crit_values
    threshold <- simplify2array(res_saemvs@thresholds)
    beta <- simplify2array(res_saemvs@beta)
    support <- res_saemvs@support
    map_to_unique_support <- res_saemvs@map_to_unique_support
    nu0_grid <- res_saemvs@nu0_grid
    nb_nu0 <- length(nu0_grid)
    pen <- res_saemvs@pen

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
