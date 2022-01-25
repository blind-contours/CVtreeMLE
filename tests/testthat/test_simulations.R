library(CVtreeMLE)

ATE <- 6

n_obs <- 500 # number of observations we want to simulate
splits <- c(0.99, 2.0, 2.5) # split points for each mixture
mins <- c(0, 0, 0) # minimum values for each mixture
maxs <- c(3, 4, 5) # maximum value for each mixture
mu <- c(0, 0, 0) # mu for each mixture
sigma <- matrix(c(1, 0.5, 0.8, 0.5, 1, 0.7, 0.8, 0.7, 1), nrow = 3, ncol = 3) # variance/covariance of mixture variables
w1_betas <- c(0.0, 0.01, 0.03, 0.06, 0.1, 0.05, 0.2, 0.04) # subspace probability relationship with covariate W1
w2_betas <- c(0.0, 0.04, 0.01, 0.07, 0.15, 0.1, 0.1, 0.04) # subspace probability relationship with covariate W2
mix_subspace_betas <- c(0.00, 0.08, 0.05, 0.01, 0.05, 0.033, 0.07, 0.09) # probability of mixture subspace (for multinomial outcome generation)
subspace_assoc_strength_betas <- c(0, 0, 0, 0, 0, 0, ATE, 0) # mixture subspace impact on outcome Y, here the subspace where M1 is lower and M2 and M3 are higher based on values in splits
marginal_impact_betas <- c(0, 0, 0) # marginal impact of mixture component on Y
eps_sd <- 0.01 # random error
binary <- FALSE # if outcome is binary


sim_data <- simulate_mixture_cube(
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

expected_y <- subset(x = sim_data, M1 < 0.99 & M2 > 2.0 & M3 > 2.5, select = y)


# not going to be exactly equal because of confounding and random error so tolerance is set high
expect_equal(mean(expected_y$y), ATE, tolerance = 0.07)


ATE <- -2

n_obs <- 1000 # number of observations we want to simulate
splits <- c(0.4, 2.2, 4) # split points for each mixture
mins <- c(0, 0, 0) # minimum values for each mixture
maxs <- c(3, 4, 5) # maximum value for each mixture
mu <- c(0, 0, 0) # mu for each mixture
sigma <- matrix(c(1, 0.5, 0.8, 0.5, 1, 0.7, 0.8, 0.7, 1), nrow = 3, ncol = 3) # variance/covariance of mixture variables
w1_betas <- c(0.0, 0.01, 0.03, 0.06, 0.1, 0.05, 0.2, 0.04) # subspace probability relationship with covariate W1
w2_betas <- c(0.0, 0.04, 0.01, 0.07, 0.15, 0.1, 0.1, 0.04) # subspace probability relationship with covariate W2
mix_subspace_betas <- c(0.00, 0.08, 0.05, 0.01, 0.05, 0.033, 0.07, 0.09) # probability of mixture subspace (for multinomial outcome generation)
subspace_assoc_strength_betas <- c(0, ATE, 0, 0, 0, 0, 0, 0) # mixture subspace impact on outcome Y, here the subspace where M1 is lower and M2 and M3 are higher based on values in splits
marginal_impact_betas <- c(0, 0, 0) # marginal impact of mixture component on Y
eps_sd <- 0.01 # random error
binary <- FALSE # if outcome is binary


sim_data <- simulate_mixture_cube(
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

# 1.  All mixtures lower than specified thresholds
# 2.  M1 is higher but M2 and M3 are lower
# 3.  M2 is higher but M1 and M3 are lower
# 4.  M1 and M2 are higher and M3 is lower
# 5.  M3 is higher and M1 and M2 are lower
# 6.  M1 and M3 are higher and M2 is lower
# 7.  M2 and M3 are higher and M1 is lower
# 8.  All mixtures are higher than thresholds

expected_y <- subset(x = sim_data, M1 >= 0.99 & M2 <= 2.0 & M3 <= 2.5, select = y)


# not going to be exactly equal because of confounding and random error so tolerance is set high
expect_equal(mean(expected_y$y), ATE, tolerance = 0.07)
