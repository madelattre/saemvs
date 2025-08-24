setGeneric(
  "loglik",
  function(data, model, tuning_algo, param, pen, p) {
    standardGeneric("loglik")
  }
)


setMethod(
  "loglik",
  signature(
    data = "saemvsData", model = "saemvsModel",
    tuning_algo = "tuningC", param = "list", pen = "character",
    p = "numeric"
  ),
  function(data, model, tuning_algo, param, pen, p) {
    # check param en fonction du modèle?

    data <- prepare_data(data, model)

    nsim <- tuning_algo@nb_is

    support <- model@x_forced_support
    supp_index <- which(c(t(support)) == 1)

    yi <- data@y_series
    ti <- data@t_series

    n <- length(yi)
    ni <- lengths(yi)

    beta_x <- data@x_phi_not_to_select %*% param$beta
    beta_x_list <- split(beta_x, row(beta_x))
    gamma <- param$gamma
    sigma2 <- param$sigma2

    phi_samples <- lapply(beta_x_list, function(m) {
      mvnfast::rmvn(nsim, mu = m, sigma = gamma)
    })

    log_lik_i <- function(i) {
      errs_i <- sum(apply(phi_samples[[i]], 1, function(x) {
        exp(-sum((yi[[i]] - g_vector_cpp(x, ti[[i]]))^2) / (2 * sigma2))
      }))
      mco_i <- sum(errs_i)
      value <- log((2 * pi * sigma2)^(-ni[i] / 2) * 1 / (nsim) * mco_i)
      return(value)
    }

    loglike <- sum(sapply(seq_along(yi), log_lik_i))

    nb_beta <- length(supp_index)

    if (pen == "e-BIC") {
      biclike <- -2 * loglike + (nb_beta - model@phi_dim) * log(n) +
        2 * log(choose(p * model@phi_dim, nb_beta - model@phi_dim))
    } else if (pen == "BIC") {
      biclike <- -2 * loglike + (nb_beta) * log(n)
    } else {
      biclike <- -2 * loglike
    }

    return(biclike)
  }
)
