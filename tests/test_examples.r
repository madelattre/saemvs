library(saemvs)

data_file <- system.file(
  "examples/simulated_data_small_dim.RData",
  package = "saemvs"
)
load(data_file)



# Example 1 : variable selection on the three model parameters

case_test <- "map_full_select"

source("tests/prepare_test_examples.r")


res <- saemvs(dat, model, init, tuning_algo, hyperparam, "e-BIC")
summary(res, 3)
# plots <- prepare_grid_plot(res)
# # quartz()
# plots$reg_plot[[1]]
# plots$reg_plot[[2]]
# plots$reg_plot[[3]]

# plots$ebic_plot

# Example 2 : variable selection on the first two model parameters, third parameter fixed

case_test <- "map_part_select_fixed"

source("tests/prepare_test_examples.r")


res <- saemvs(dat, model, init, tuning_algo, hyperparam, "e-BIC")
summary(res, 3)
# plots <- prepare_grid_plot(res)
# # quartz()
# plots$reg_plot[[1]]
# plots$reg_plot[[2]]

# plots$ebic_plot


# Example 3 : variable selection on the first two model parameters, third parameter random

case_test <- "map_part_select_nofixed"

source("tests/prepare_test_examples.r")

res <- saemvs(dat, model, init, tuning_algo, hyperparam, "e-BIC")
summary(res, 3)
# plots <- prepare_grid_plot(res)
# # quartz()
# plots$reg_plot[[1]]
# plots$reg_plot[[2]]

# plots$ebic_plot

# Example 4 : variable selection on the first two model parameters, third parameter fixed and depending on fixed covariates

case_test <- "map_part_select_fixed_cov"

source("tests/prepare_test_examples.r")

res <- saemvs(dat, model, init, tuning_algo, hyperparam, "e-BIC")
summary(res, 3)


# plots <- prepare_grid_plot(res)
# # quartz()
# plots$reg_plot[[1]]
# plots$reg_plot[[2]]

# plots$ebic_plot

# pose problème si trop d'itérations à cause de la variance nulle?

# case_test <- "map_part_select_fixed_cov"




# # saem

res_test <- test_saemvs(dat, model, init, tuning_algo, full_hyperparam)
quartz()
convergence_plot(res_test, "beta_s",
  sel_components = c("(1,1)", "(1,2)", "(4,2)")
)


convergence_plot(res_test, "gamma_s",
  sel_components = c("(1,1)", "(1,2)", "(2,2)")
)



# state$beta_hdim[[tuning_algo@niter + 1]]
# state$gamma_hdim[[tuning_algo@niter + 1]]
# state$beta_ldim[[tuning_algo@niter + 1]]
# state$gamma_ldim[[tuning_algo@niter + 1]]
# state$sigma2[tuning_algo@niter + 1]

# # grille qui permet d'obtenir des supports différents sur l'exemple en petite dimension
