# ------------------------------------------------------------
# Helper for saemvsProcessedData tests
# ------------------------------------------------------------

# Valid arguments for saemvsProcessedData
valid_args_processed <- list(
  y_series = list(rnorm(5), rnorm(3)), # 2 individuals
  t_series = list(1:5, 1:3),
  x_phi_to_select = matrix(rnorm(4), nrow = 2), # 2 rows = n_ind
  x_phi_not_to_select = matrix(rnorm(2), nrow = 2),
  x_phi_not_to_select_list = list(rnorm(5), rnorm(3))
)

# Invalid cases
invalid_cases_processed <- list(
  "non-numeric x_phi_to_select" = list(
    args = modifyList(valid_args_processed, list(
      x_phi_to_select = matrix(c("a", "b", "c", "d"), nrow = 2)
    )),
    msg = "All entries in 'x_phi_to_select' must be numeric"
  ),
  "non-numeric x_phi_not_to_select" = list(
    args = modifyList(
      valid_args_processed,
      list(x_phi_not_to_select = matrix(c("a", "b"), nrow = 2))
    ),
    msg = "All entries in 'x_phi_not_to_select' must be numeric"
  ),
  "wrong rows in x_phi_to_select" = list(
    args = modifyList(
      valid_args_processed,
      list(x_phi_to_select = matrix(rnorm(6), nrow = 3)) # should have 2 rows
    ),
    msg = "'x_phi_to_select' must have 2 rows"
  ),
  "wrong length in x_phi_not_to_select_list" = list(
    args = list(
      y_series = valid_args_processed$y_series,
      t_series = valid_args_processed$t_series,
      x_phi_to_select = valid_args_processed$x_phi_to_select,
      x_phi_not_to_select = valid_args_processed$x_phi_not_to_select,
      x_phi_not_to_select_list = list(rnorm(5)) # nouvelle valeur
    ),
    msg = "'x_phi_not_to_select_list' must have 2 elements"
  )
)
