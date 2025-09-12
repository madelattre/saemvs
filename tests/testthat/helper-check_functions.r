# helper-check_functions.R

# Minimal valid saemvsData object
create_test_data <- function(n = 10, p = 2) {
  y <- list(rnorm(n))
  t <- list(seq_len(n))
  x_candidates <- matrix(rnorm(p), ncol = p)
  saemvsData(y = y, t = t, x_candidates = x_candidates)
}

# Minimal valid saemvsInit object
# create_test_init <- function(phi_dim = 1) {
#   saemvsInit(
#     intercept = rep(0, phi_dim),
#     beta_forced = NULL,
#     beta_candidates = NULL,
#     cov_re = diag(phi_dim),
#     sigma2 = 1
#   )
# }

create_test_init <- function(model, sigma2 = 1) {
  intercept <- rep(0, model@phi_dim)  # matches phi_dim
  saemvsInit(
    intercept = intercept,
    beta_forced = NULL,
    beta_candidates = NULL,
    cov_re = diag(1, nrow = model@phi_dim),
    sigma2 = sigma2
  )
}

# Minimal valid saemvsHyper object
create_test_hyper <- function(n_selected = 2, n_candidates = 2) {
  hyper_slab <- saemvsHyperSlab(
    slab_parameter = 1,
    cov_re_prior_scale = diag(1, n_selected),
    cov_re_prior_df = n_selected + 1 #,
    # residual_variance_prior_shape = 1,
    # residual_variance_prior_rate = 1,
    # phi_intercept_prior_variance = rep(1, n_selected),
    # inclusion_prob_prior_a = rep(1, n_selected),
    # inclusion_prob_prior_b = rep(1, n_selected)
  )
  
  spike_parameter <- 0.1 #seq(0.1, 1, length.out = n_candidates)
  
  saemvsHyperSpikeAndSlab(spike_parameter, hyper_slab)
}

create_test_hyper_slab <- function(n_selected = 2, n_candidates = 2) {
  saemvsHyperSlab(
    slab_parameter = 1,
    cov_re_prior_scale = diag(1, n_selected),
    cov_re_prior_df = n_selected + 1 #,
    # residual_variance_prior_shape = 1,
    # residual_variance_prior_rate = 1,
    # phi_intercept_prior_variance = rep(1, n_selected),
    # inclusion_prob_prior_a = rep(1, n_selected),
    # inclusion_prob_prior_b = rep(1, n_selected)
  )
}

# Create a minimal valid saemvsModel object
create_test_model <- function(phi_dim = 2,
                              phi_to_select_idx = 1,
                              phi_fixed_idx = integer(0),
                              n_forced = 0) {
  # Dummy model function: linear in phi
  g <- function(phi, t) {
    # Just return a linear combination to keep it valid
    as.numeric(t %*% phi[seq_len(min(length(phi), ncol(t)))])
  }

  # Forced support: matrix with n_forced rows and phi_dim columns
  if (n_forced > 0) {
    x_forced_support <- matrix(
      1, nrow = n_forced, ncol = phi_dim
    )
  } else {
    x_forced_support <- matrix(numeric(0), nrow = 0, ncol = phi_dim)
  }

  saemvsModel(
    g = g,
    phi_dim = phi_dim,
    phi_to_select_idx = phi_to_select_idx,
    phi_fixed_idx = phi_fixed_idx,
    x_forced_support = x_forced_support
  )
}
