library(mvtnorm)

set.seed(1123)

# ----- Define the function -----
g <- function(phi, psi, t) {
  # Logistic-like function: psi[1] + phi[1] / (1 + exp(-(t - phi[2])))
  psi[1] + phi[1] / (1 + exp(-(t - phi[2])))
}

# ----- Parameters -----
n <- 100 # Number of individuals
jmin <- 25 # Minimum number of observations per individual
jmax <- 30 # Maximum number of observations per individual
tmin <- 100 # Minimum time
tmax <- 3000 # Maximum time
sigma2 <- 30 # Observation variance
gamma2 <- diag(c(200, 10)) # Covariance for phi random effects
mu <- matrix(c(1200, 500), ncol = 1) # Intercepts for phi
beta <- rbind(
  c(100, 50, 20, rep(0, 47)),
  c(20, 0, 10, rep(0, 47))
)
beta_tilde <- cbind(mu, beta) # Combine mu and beta
psi <- c(200) # Constant for g function


p <- dim(beta)[2] # Number of covariates

# ----- Vectorized generation of v_tilde and v -----
v_tilde <- cbind(
  1, # Column of ones for intercept
  scale(matrix(rnorm(n * p), nrow = n, ncol = p)) # Standardized random covariates
)
x <- v_tilde[, -1] # Remove intercept column

# ----- Vectorized calculation of phi -----∏
phi_mean <- beta_tilde %*% t(v_tilde) # (2 x n) matrix
# Add multivariate normal noise for each individual
phi <- t(phi_mean + t(rmvnorm(n, mean = c(0, 0), sigma = gamma2))) # Result: n x 2

# ----- Generate times and observations -----
t_list <- vector("list", n)
y_list <- vector("list", n)

# Sample number of observations per individual
j_vec <- sample(jmin:jmax, n, replace = TRUE)

# Pré-allocation
y_all <- vector("list", n)
time_all <- vector("list", n)
id_all <- vector("list", n)

for (i in seq_len(n)) {
  t_i <- seq(tmin, tmax, length.out = j_vec[i])
  y_i <- g(phi[i, ], psi, t_i) + rnorm(j_vec[i], sd = sqrt(sigma2))

  time_all[[i]] <- t_i
  y_all[[i]] <- y_i
  id_all[[i]] <- rep(i, j_vec[i])
}

# Combine into a single long dataframe
df_long <- data.frame(
  id = unlist(id_all),
  time = unlist(time_all),
  y = unlist(y_all)
)

# Covariates dataframe
df_cov <- as.data.frame(x)
colnames(df_cov) <- paste0("x", 1:p)
df_cov$id <- 1:n

example_data_list <- list(
  y_list = y_all,
  t_list = time_all,
  x = x
)

# ----- Save results -----

save(df_long, file = "data/df_long.rda")
save(df_cov, file = "data/df_cov.rda",
     compress = "xz")