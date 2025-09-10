# ------------------------------------------------------------
# Helper for saemvsModel tests
# ------------------------------------------------------------

# Fonction simple valide
valid_model_func <- function(phi, t) phi * t

# Arguments valides de base
valid_args_model <- list(
  g = valid_model_func,
  phi_dim = 3,
  phi_to_select_idx = c(1, 2),
  phi_fixed_idx = c(3),
  x_forced_support = matrix(rnorm(6), nrow = 2, ncol = 3) # phi_dim = 3
)

# Cas invalides
invalid_cases_model <- list(
  "phi_dim not positive" = list(
    args = modifyList(valid_args_model, list(phi_dim = 0)),
    msg = "'phi_dim' must be a positive integer"
  ),
  "phi_to_select_idx out of range" = list(
    args = modifyList(valid_args_model, list(phi_to_select_idx = c(0, 2))),
    msg = "'phi_to_select_idx' must contain integers between 1 and 'phi_dim'"
  ),
  "phi_fixed_idx out of range" = list(
    args = modifyList(valid_args_model, list(phi_fixed_idx = c(3, 4))),
    msg =
      "'phi_not_to_select_idx' must contain integers between 1 and 'phi_dim'"
  ),
  "overlapping indices" = list(
    args = modifyList(
      valid_args_model,
      list(phi_to_select_idx = c(1, 2), phi_fixed_idx = c(2, 3))
    ),
    msg = "'phi_to_select_idx' and 'phi_fixed_idx' must not overlap"
  ),
  "x_forced_support wrong ncol" = list(
    args = modifyList(valid_args_model, list(
      x_forced_support = matrix(rnorm(6), nrow = 3, ncol = 2)
    )),
    msg = "'x_forced_support' must have 3 columns"
  ),
  "invalid model function" = list(
    args = modifyList(valid_args_model, list(g = function(a) a)),
    msg = "The model function must have arguments 'phi' and 't'"
  )
)
