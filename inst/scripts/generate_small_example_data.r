library(mvtnorm)

set.seed(1123)

# ----- Define the function -----
g <- function(phi, psi, t) {
  # Logistic-like function: psi[1] + phi[1] / (1 + exp(-(t - phi[2])))
  psi[1] + phi[1] / (1 + exp(-(t - phi[2])))
}

# ----- Parameters -----
n <- 200 # Number of individuals
jmin <- 25 # Minimum number of observations per individual
jmax <- 30 # Maximum number of observations per individual
tmin <- 100 # Minimum time
tmax <- 3000 # Maximum time
sigma2 <- 30 # Observation variance
gamma2 <- diag(c(200, 10)) # Covariance for phi random effects
mu <- matrix(c(1200, 500), ncol = 1) # Intercepts for phi
beta <- rbind(
  c(100, 50, 20, 0),
  c(20, 0, 10, 0)
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

for (i in 1:n) {
  t_i <- sort(runif(j_vec[i], tmin, tmax)) # Random times
  y_i <- g(phi[i, ], psi, t_i) + rnorm(j_vec[i], sd = sqrt(sigma2)) # Observations with noise
  t_list[[i]] <- t_i
  y_list[[i]] <- y_i
}

small_example_data <- list(
  y_list = y_list,
  t_list = t_list,
  x = x
)

# ----- Save results -----
save(small_example_data, file = "data/small_example_data.rda")
