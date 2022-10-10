library(CVtreeMLE)
library(tidyr)

# 1.  All mixtures lower than specified thresholds
# 2.  M1 is higher but M2 and M3 are lower
# 3.  M2 is higher but M1 and M3 are lower
# 4.  M1 and M2 are higher and M3 is lower
# 5.  M3 is higher and M1 and M2 are lower
# 6.  M1 and M3 are higher and M2 is lower
# 7.  M2 and M3 are higher and M1 is lower
# 8.  All mixtures are higher than thresholds


ate <- 6

n_obs <- 10000
splits <- c(0.99, 2.0, 2.5)
mins <- c(0, 0, 0)
maxs <- c(3, 4, 5)
mu <- c(0, 0, 0)
sigma <- matrix(c(1, 0.5, 0.8, 0.5, 1, 0.7, 0.8, 0.7, 1), nrow = 3, ncol = 3)
w1_betas <- c(0.0, 0.01, 0.03, 0.06, 0.1, 0.05, 0.2, 0.04)
w2_betas <- c(0.0, 0.04, 0.01, 0.07, 0.15, 0.1, 0.1, 0.04)
mix_subspace_betas <- c(0.00, 0.08, 0.05, 0.01, 0.05, 0.033, 0.07, 0.09)
subspace_assoc_strength_betas <- c(0, 0, 0, 0, 0, 0, ate, 0)
marginal_impact_betas <- c(0, 0, 0)
eps_sd <- 0.01
binary <- FALSE


sim_data_exp <- simulate_mixture_cube(
  n_obs = n_obs,
  splits = splits,
  mins = mins,
  maxs = maxs,
  mu = mu,
  sigma = sigma,
  w1_betas = w1_betas,
  w2_betas = w2_betas,
  mix_subspace_betas = mix_subspace_betas,
  subspace_assoc_strength_betas = subspace_assoc_strength_betas,
  marginal_impact_betas = marginal_impact_betas,
  eps_sd = eps_sd,
  binary = binary
)

expected_y_exp <- subset(x = sim_data_exp, M1 < 0.99 & M2 > 2.0 & M3 > 2.5, select = y)

ate <- 0

n_obs <- 10000
splits <- c(0.4, 2.2, 4)
mins <- c(0, 0, 0)
maxs <- c(3, 4, 5)
mu <- c(0, 0, 0)
sigma <- matrix(c(1, 0.5, 0.8, 0.5, 1, 0.7, 0.8, 0.7, 1), nrow = 3, ncol = 3)
w1_betas <- c(0.0, 0.01, 0.03, 0.06, 0.1, 0.05, 0.2, 0.04)
w2_betas <- c(0.0, 0.04, 0.01, 0.07, 0.15, 0.1, 0.1, 0.04)
mix_subspace_betas <- c(0.00, 0.08, 0.05, 0.01, 0.05, 0.033, 0.07, 0.09)
subspace_assoc_strength_betas <- c(0, ate, 0, 0, 0, 0, 0, 0)
marginal_impact_betas <- c(0, 0, 0)
eps_sd <- 0.01
binary <- FALSE


sim_data_unexp <- simulate_mixture_cube(
  n_obs = n_obs,
  splits = splits,
  mins = mins,
  maxs = maxs,
  mu = mu,
  sigma = sigma,
  w1_betas = w1_betas,
  w2_betas = w2_betas,
  mix_subspace_betas = mix_subspace_betas,
  subspace_assoc_strength_betas = subspace_assoc_strength_betas,
  marginal_impact_betas = marginal_impact_betas,
  eps_sd = eps_sd,
  binary = binary
)

expected_y_unexp <- subset(x = sim_data_unexp, M1 < 0.99 & M2 > 2.0 & M3 > 2.5, select = y)


expect_equal(mean(expected_y_exp$y) - mean(expected_y_unexp$y) , 6, tolerance = 0.07)
