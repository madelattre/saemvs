threshold <- function(nu1, nu0, alpha) {
  value <- 2 * nu0 * nu1 * log(sqrt(nu1 / nu0) *
    (1 - alpha) / alpha) / (nu1 - nu0)
  return(sqrt(value))
}

p_star_fast <- function(beta, alpha, nu0, nu1) {
  norm1 <- dnorm(beta, mean = 0, sd = sqrt(nu1))
  norm0 <- dnorm(beta, mean = 0, sd = sqrt(nu0))
  num <- norm1 * rep(alpha, each = nrow(beta))
  denom <- num + norm0 * rep(1 - alpha, each = nrow(beta))

  p_star <- num / denom
  return(p_star)
}


get_case <- function(config) {
  if (config$method_type == "map") {
    if (config$num_parameters_not_to_select == 0) {
      "map_full_select"
    } else if (length(config$fixed_parameters_indices) == 0) {
      "map_part_select_nofixed"
    } else {
      "map_part_select_fixed"
    }
  } else {
    if (length(config$fixed_parameters_indices) > 0) "mle_fixed" else "mle_nofixed"
  }
}
