
generate_P_0 <- function(truth, n) {
  # number of observations we want to simulate
  n_obs <- n
  # split points for each mixture
  splits <- c(0.8, 2.5, 3.6)
  # minimum values for each mixture
  mins <- c(0, 0, 0)
  # maximum value for each mixture
  maxs <- c(3, 4, 5)
  # mu for each mixture
  mu <- c(0, 0, 0)
  # variance/covariance of mixture variables
  sigma <- matrix(c(1, 0.5, 0.8, 0.5, 1, 0.7, 0.8, 0.7, 1), nrow = 3, ncol = 3)
  # subspace probability relationship with covariate W1
  w1_betas <- c(0.0, 0.01, 0.03, 0.06, 0.1, 0.05, 0.2, 0.04)
  # subspace probability relationship with covariate W2
  w2_betas <- c(0.0, 0.04, 0.01, 0.07, 0.15, 0.1, 0.1, 0.04)
  # probability of mixture subspace (for multinomial outcome generation)
  mix_subspace_betas <- c(0.00, 0.08, 0.05, 0.01, 0.05, 0.033, 0.07, 0.09)
  # mixture subspace impact on outcome Y, here the subspace where M1 is lower and M2 and M3 are higher based on values in splits
  subspace_assoc_strength_betas <- truth
  # marginal impact of mixture component on Y
  marginal_impact_betas <- c(0, 0, 0)
  # random error
  eps_sd <- 0.01
  # if outcome is binary
  binary <- FALSE

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

  return(sim_data)
}
