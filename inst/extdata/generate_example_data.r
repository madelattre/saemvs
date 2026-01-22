# ============================================================
#  Simulation de données longitudinales avec effets aléatoires
# ============================================================

library(mvtnorm)

set.seed(1123)

# ------------------------------------------------------------
# Model function
# ------------------------------------------------------------
g <- function(phi, t) {
  # Fonction logistique décalée
  phi[2] + phi[1] / (1 + exp(-(t - phi[3])))
}

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
n     <- 100
jmin  <- 25
jmax  <- 30
tmin  <- 100
tmax  <- 3000

sigma2 <- 30

# ------------------------------------------------------------
# Random effects
# ------------------------------------------------------------
gamma2 <- diag(c(200, 10, 50))  # Covariance matrix

# ------------------------------------------------------------
# Parameters and covariates
# ------------------------------------------------------------
mu <- matrix(c(1200, 200, 500), ncol = 1)

beta <- rbind(
  c(100, 50, 20, rep(0, 47)),
  rep(0, 50),
  c(20,  0, 10, rep(0, 47))
)

beta_tilde <- cbind(mu, beta)
p <- ncol(beta)

covariates <- cbind(
  1,
  scale(matrix(rnorm(n * p), nrow = n, ncol = p))
)

x <- covariates[, -1]

# ------------------------------------------------------------
# Random effects
# ------------------------------------------------------------
phi_mean <- beta_tilde %*% t(covariates)

phi <- t(
  phi_mean +
    t(rmvnorm(n, mean = c(0, 0, 0), sigma = gamma2))
)

# ------------------------------------------------------------
# Longitudinal data
# ------------------------------------------------------------
j_vec <- sample(jmin:jmax, n, replace = TRUE)

time_all <- vector("list", n)
y_all    <- vector("list", n)
id_all   <- vector("list", n)

for (i in seq_len(n)) {
  t_i <- seq(tmin, tmax, length.out = j_vec[i])

  y_i <- g(phi[i, ], t_i) +
    rnorm(j_vec[i], sd = sqrt(sigma2))

  time_all[[i]] <- t_i
  y_all[[i]]    <- y_i
  id_all[[i]]   <- rep(i, j_vec[i])
}

# ------------------------------------------------------------
# Dataframes
# ------------------------------------------------------------
df_long <- data.frame(
  id   = unlist(id_all),
  time = unlist(time_all),
  y    = unlist(y_all)
)

df_cov <- as.data.frame(x)
colnames(df_cov) <- paste0("x", seq_len(p))
df_cov$id <- seq_len(n)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
save(df_long, file = "data/df_long.rda")
save(df_cov,  file = "data/df_cov.rda", compress = "xz")
