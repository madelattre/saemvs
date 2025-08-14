setGeneric(
  "loglik",
  function(data, model, tuning_algo, param, ebic_dim) {
    standardGeneric("loglik")
  }
)


setMethod(
  "loglik",
  signature(
    data = "dataC", model = "modelC",
    tuning_algo = "tuningC", param = "list", ebic_dim = "numeric"
  ),
  function(data, model, tuning_algo, param, ebic_dim) {
    # check param en fonction du modèle?

    data <- prepare_data(data, model)

    nsim <- tuning_algo@nb_is

    g <- model@model_func
    support <- model@covariate_support
    supp_index <- which(c(t(support)) == 1)

    yi <- data@y_list
    ti <- data@t_list

    n <- length(yi)
    ni <- lengths(yi)

    beta_x <- data@w %*% param$beta
    beta_x_list <- split(beta_x, row(beta_x))
    gamma <- param$gamma
    sigma2 <- param$sigma2

    phi_samples <- lapply(beta_x_list, function(m) {
      mvnfast::rmvn(nsim, mu = m, sigma = gamma)
    })

    log_lik_i <- function(i) {
      errs_i <- sum(apply(phi_samples[[i]], 1, function(x) {
        exp(-sum((yi[[i]] - g(x, ti[[i]]))^2) / (2 * sigma2))
      }))
      mco_i <- sum(errs_i)
      value <- log((2 * pi * sigma2)^(-ni[i] / 2) * 1 / (nsim) * mco_i)
      return(value)
    }

    loglike <- sum(sapply(seq_along(yi), log_lik_i))

    nb_beta <- length(supp_index)
    ebic <- loglike + (nb_beta - model@q_phi) * log(n) +
      2 * log(choose(ebic_dim * model@q_phi, nb_beta - model@q_phi))

    return(list(loglik = loglike, ebic = ebic))
  }
)
