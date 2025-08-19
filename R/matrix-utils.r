# Creation of the design matrix x

# Devrait plutôt s'appeler from_x_to_list
from_v_to_x <- function(v, support, q) {
  # supp_index <- apply(support, 1, function(x) which(x == 1))
  # ne donne pas une liste si le support est complet
  split_cols <- split(support, col(support))
  supp_index <- lapply(split_cols, function(x) which(x == 1))
  v_list <- lapply(split(v, row(v)), matrix, nrow = 1)
  v_q_list <- lapply(v_list, function(x) replicate(q, x, simplify = FALSE))
  x <- lapply(v_q_list, x_per_indiv_rcpp, supp_index)
  return(x)
}

# shrink_covariance_matrix <- function(mat_old, mat_temp, indices, alphas) {
#   d <- nrow(mat_old)

#   mat_new <- mat_temp

#   shrinked <- rep(FALSE, d)
#   shrinked[indices] <- TRUE

#   diag_mask <- logical(d)
#   diag_mask[indices] <- TRUE
#   mat_new[cbind(indices, indices)] <- diag(mat_old)[indices] * alphas

#   shrinked_matrix <- outer(shrinked, shrinked, `|`) & !diag(diag(1, d))
#   mat_new[shrinked_matrix] <- 0

#   return(mat_new)
# }

shrink_covariance_matrix <- function(Sigma_old, Sigma_temp, indices, alphas) {
  # Essayer de vectoriser la fonction

  d <- nrow(Sigma_old)

  Sigma_new <- matrix(0, nrow = d, ncol = d)
  shrinked <- rep(FALSE, d)
  shrinked[indices] <- TRUE
  sqrt_alphas <- rep(1, d)
  sqrt_alphas[indices] <- sqrt(alphas)

  for (i in 1:d) {
    for (j in 1:d) {
      if (i == j && shrinked[i]) {
        Sigma_new[i, j] <- Sigma_old[i, j] * alphas[which(indices == i)]
      } else if (shrinked[i] && shrinked[j]) {
        Sigma_new[i, j] <- 0 # Sigma_old[i, j] * sqrt_alphas[i] * sqrt_alphas[j]
      } else if (shrinked[i]) {
        Sigma_new[i, j] <- 0 # Sigma_old[i, j] * sqrt_alphas[i]
      } else if (shrinked[j]) {
        Sigma_new[i, j] <- 0 # Sigma_old[i, j] * sqrt_alphas[j]
      } else {
        Sigma_new[i, j] <- Sigma_temp[i, j]
      }
    }
  }

  return(Sigma_new)
}

zero_out_shrinked <- function(mat, indices) {
  mat_zeroed <- mat
  mat_zeroed[indices, ] <- 0
  mat_zeroed[, indices] <- 0
  return(mat_zeroed)
}

merge_support <- function(
    fixed_support, selected_support, w,
    nb_phi_s, nb_phi_ns) {
  # aucun des deux supports mergés ne doit contenir l'intercept
  # la fonction génère un support sans intercept
  # -> plus pratique pour appliquer saem dessus ensuite

  if (!is.null(selected_support)) {
    nb_selec <- dim(selected_support)[1]
    v_supp <- cbind(
      selected_support, matrix(0, nrow = nb_selec, ncol = nb_phi_ns)
    )
  } else {
    v_supp <- NULL
  }

  dim_fixed <- dim(fixed_support)[1]

  if (is.null(w)) {
    w_supp <- NULL
  } else {
    phis_w <- matrix(0, nrow = dim_fixed, ncol = nb_phi_s)
    w_supp <- cbind(phis_w, fixed_support)
  }


  new_covariate_support <- rbind(
    w_supp,
    v_supp
  )

  return(new_covariate_support)
}


merge_init_beta <- function(beta_hdim, beta_ldim, lines_with_ones) {
  # a-t-on vraiment besoin de w
  nb_phi_s <- dim(beta_hdim)[2]

  if (is.null(beta_ldim)) {
    nb_phi_ns <- 0
  } else {
    nb_phi_ns <- dim(beta_ldim)[2]
  }


  intercepts <- matrix(
    c(beta_hdim[1, ], beta_ldim[1, ]),
    nrow = 1
  )

  if (!is.null(beta_ldim) && (nrow(beta_ldim) > 1)) { ## replace w by beta_ldim
    top_s <- matrix(0, ncol = nb_phi_s, nrow = dim(beta_ldim)[1] - 1)
    top <- cbind(top_s, beta_ldim[-1, ])
  } else {
    top <- NULL
  }

  if (length(lines_with_ones) > 1) {
    bottom_ns <- matrix(0,
      ncol = nb_phi_ns,
      nrow = length(lines_with_ones) - 1
    )
    bottom_s <- matrix(
      beta_hdim[lines_with_ones[-1], ],
      ncol = nb_phi_s
    )
    bottom <- cbind(bottom_s, bottom_ns)
  } else {
    bottom_ns <- NULL
    bottom_s <- NULL
    bottom <- NULL
  }


  final_beta_init <- rbind(
    intercepts,
    top,
    bottom
  )

  return(final_beta_init)
}
