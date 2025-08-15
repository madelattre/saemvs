data_map_cov <- dataC(y = y_list, t = times, v = v, w = v[, 1:2])
data_map <- dataC(y = y_list, t = times, v = v)
data_mle <- dataC(y = y_list, t = times, w = v)

# -- Prepare models

g <- function(phi, t) {
  phi[3] + phi[1] / (1 + exp(-(t - phi[2])))
}

q_phi <- 3

mod_map_full <- modelC(g = g, nphi = q_phi, phi_select = c(1, 2, 3))

mod_map_part_notfixed <- modelC(
  g = g, nphi = q_phi, phi_select = c(1, 2)
)

mod_map_part_fixed <- modelC(
  g = g, nphi = q_phi, phi_select = c(1, 2), phi_fixed = c(3)
)

mod_map_part_fixed_cov <- modelC(
  g = g, nphi = q_phi, phi_select = c(1, 2), phi_fixed = c(3),
  support = matrix(c(1, 1), ncol = 1)
)

true_support <- t(rbind(
  c(1, 1, 1, 0),
  c(1, 0, 1, 0),
  c(0, 0, 0, 0)
))


mod_mle_full <- modelC(
  g = g, nphi = q_phi, support = true_support
)

mod_mle_part_fixed <- modelC(
  g = g, nphi = q_phi, support = true_support, phi_fixed = c(3)
)

## -- Prepare initial values for each model


init_map_full_select <- initC(
  beta_s = matrix(
    c(
      1300, 80, 20, 20, 0,
      350, 15, 0, 0, 0,
      250, 5, -5, 5, -5
    ),
    ncol = 3
  ),
  gamma_s = diag(c(400, 200, 200)),
  sigma2 = 100,
  alpha = rep(0.5, 3)
)

init_map_part_select_nofixed <- initC(
  beta_s = matrix(
    c(
      1300, 80, 20, 20, 0,
      350, 15, 0, 0, 0
    ),
    ncol = 2
  ),
  gamma_s = diag(c(400, 200)),
  sigma2 = 100,
  beta_ns = matrix(250, ncol = 1, nrow = 1),
  gamma_ns = matrix(200, ncol = 1, nrow = 1),
  alpha = c(0.5, 0.5)
)

init_map_part_select_fixed <- initC(
  beta_s = matrix(
    c(
      1300, 80, 20, 20, 0,
      350, 15, 0, 0, 0
    ),
    ncol = 2
  ),
  gamma_s = diag(c(400, 200)),
  sigma2 = 100,
  beta_ns = matrix(250, 1, 1),
  gamma_ns = matrix(200, ncol = 1, nrow = 1),
  alpha = c(0.5, 0.5)
)

init_map_part_select_fixed_cov <- initC(
  beta_s = matrix(
    c(
      1300, 80, 20, 20, 0,
      350, 15, 0, 0, 0
    ),
    ncol = 2
  ),
  gamma_s = diag(c(400, 200)),
  sigma2 = 100,
  beta_ns = matrix(c(250, 5, -5), ncol = 1),
  gamma_ns = matrix(200, ncol = 1, nrow = 1),
  alpha = c(0.5, 0.5)
)


init_mle_nofixed <- initC(
  beta_ns = matrix(
    c(
      1300, 80, 20, 20, 0,
      350, 15, 0, 0, 0,
      250, 5, -5, 5, -5
    ),
    ncol = 3
  ),
  gamma_ns = diag(c(400, 200, 200)),
  sigma2 = 100
)


init_mle_fixed <- initC(
  beta_ns = matrix(
    c(
      1300, 80, 20, 20, 0,
      350, 15, 0, 0, 0,
      250, 5, -5, 5, -5
    ),
    ncol = 3
  ),
  gamma_ns = diag(c(400, 200, 200)),
  sigma2 = 100
)


## Prepare hyperparameters


hyperparam_map_2s <- hyperC(
  d = 4, sgam = diag(rep(0.2, 2))
)

full_hyperparam_map_2s <- fullHyperC(0.002, hyperparam_map_2s)

hyperparam_map_3s <- hyperC(
  d = 4, sgam = diag(rep(0.2, 3))
)

full_hyperparam_map_3s <- fullHyperC(0.002, hyperparam_map_3s)


switch(case_test,
  map_full_select = {
    model <- mod_map_full
    init <- init_map_full_select
    dat <- data_map
    hyperparam <- hyperparam_map_3s
    full_hyperparam <- full_hyperparam_map_3s
  },
  map_part_select_nofixed = {
    model <- mod_map_part_notfixed
    init <- init_map_part_select_nofixed
    dat <- data_map
    hyperparam <- hyperparam_map_2s
    full_hyperparam <- full_hyperparam_map_2s
  },
  map_part_select_fixed = {
    model <- mod_map_part_fixed
    init <- init_map_part_select_fixed
    dat <- data_map
    hyperparam <- hyperparam_map_2s
    full_hyperparam <- full_hyperparam_map_2s
  },
  map_part_select_fixed_cov = {
    model <- mod_map_part_fixed_cov
    init <- init_map_part_select_fixed_cov
    dat <- data_map_cov
    hyperparam <- hyperparam_map_2s
    full_hyperparam <- full_hyperparam_map_2s
  },
  mle_nofixed = {
    model <- mod_mle_full
    init <- init_mle_nofixed
    dat <- data_mle
    hyperparam <- hyperC(NULL, NULL, NULL) # Modifier avec les nouvelles classes
    full_hyperparam <- fullHyperC(NULL, hyperparam)
  },
  mle_part_fixed = {
    model <- mod_mle_part_fixed
    init <- init_mle_fixed
    dat <- data_mle
    hyperparam <- hyperC(NULL, NULL, NULL) # Modifier avec les nouvelles classes
    full_hyperparam <- fullHyperC(NULL, hyperparam)
  }
)


##

tuning_algo <- tuningC(nu0_grid = c(0.5, 0.002, 0.01), seed = -5)
