library(saemvs)

data_file <- system.file(
  "examples/simulated_data_small_dim.RData",
  package = "saemvs"
)
load(data_file)
case_test <- "map_part_select_fixed"

source("tests/prepare_test_examples.r")


# saem

# g_cpp <- compile_model(model@model_func)
# state <- run_saem(dat, model, init, tuning_algo, full_hyperparam)

# state$beta_hdim[[tuning_algo@niter + 1]]
# state$gamma_hdim[[tuning_algo@niter + 1]]
# state$beta_ldim[[tuning_algo@niter + 1]]
# state$gamma_ldim[[tuning_algo@niter + 1]]
# state$sigma2[tuning_algo@niter + 1]

# # première étape de saemvs

# res1 <- saemvs_one_map_run(dat, model, init, tuning_algo, full_hyperparam)

# print(res1)


# grille qui permet d'obtenir des supports différents sur l'exemple en petite dimension


res <- saemvs(dat, model, init, tuning_algo, hyperparam, "BIC")

plots <- prepare_grid_plot(res)
plots$reg_plot[[1]]
plots$reg_plot[[2]]
plots$ebic_plot

which.min(saemvs_res$ebic)

# quartz()

state <- run_saem(dat, model, init, tuning_algo, full_hyperparam)


# pose problème si trop d'itérations à cause de la variance nulle?

# case_test <- "map_part_select_fixed_cov"
# case_test <- "map_part_select_fixed"
# case_test <- "map_part_select_nofixed"
# case_test <- "map_part_select_fixed"
# case_test <- "mle_part_fixed"
# case_test <- "mle_nofixed"
