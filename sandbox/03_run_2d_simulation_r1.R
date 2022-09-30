# packages
library(here)
library(data.table)
library(tidyverse)
library(origami)
library(sl3)
library(dplyr)
library(foreach)

library(CVtreeMLE)
source(here("sandbox", "02_fit_estimators.R"))
source(here("sandbox", "simulate_2d_data.R"))

# Sys.unsetenv("GITHUB_PAT")
devtools::install(upgrade = "never")

# simulation parameters
seed <- set.seed(7259)
n_sim <- 5 # number of simulations
n_obs <- (cumsum(rep(sqrt(40), 7))^2)[-1] # sample sizes at root-n scale
# n_obs <- n_obs[1:3]
true_rule <- "m1 > 4.15 & m2 > 5.17"


# Establish globals ---------------------------

# number of observations to create
n <- 500000

# means of the exposures
mu <- c(1,2)

# covariance of the exposures
sigma <- matrix(c(1, 0.1,
                  0.1, 1),
                nrow = 2,
                ncol = 2)

# number of cuts per exposure
n_cuts <- 5

# exposure grid
exposure_grid <- expand.grid(seq(5), seq(5))
labels <- apply(exposure_grid, 1, paste, collapse = " ")
exposure_grid <- cbind.data.frame(exposure_grid, labels)

# beta matrix
c_matrix <- matrix(c(0.2,0.3,0.2,0.5),
                   ncol  = 2,
                   nrow = 2)

# Generate exposures -----------------

exposure_results <- gen_exposures(mu = mu,
                                  sigma = sigma,
                                  n = 5000000,
                                  n_cuts)

exposures <- exposure_results$exposures_labeled

# Generate simulated data -----------------
P_0_sim <- gen_covariates(n) %>% # gen covariates
  gen_multinom(exposure_grid = exposure_grid) %>% # assign obs to exposure cube
  assign_exposures(exposure_grid = exposure_grid,
                   exposures = exposures) %>% # put continuous exposure in cube
  assign_outcomes(exposure_grid = exposure_grid,
                  c_matrix = c_matrix)  # assign outcome based on cube

P_0_data <- P_0_sim$data

# perform simulation across sample sizes
sim_results <- lapply(n_obs, function(sample_size) {
  # get results in parallel
  results <- list()

  for(this_iter in seq_len(n_sim)) {

    data_sim <-  P_0_data %>%
      slice_sample(n = sample_size)

    # print(mean(data_sim$outcome_obs))

    # data_sim <- P_0_data %>%
    #   group_by(region_label) %>%
    #   sample_n(sample_size / 25)

    est_out <- fit_estimators(data = as.data.frame(data_sim),
                              true_rule,
                              covars = c("age", "sex", "bmi"),
                              exposures = c("m1", "m2"),
                              outcome = "outcome_obs",
                              P_0_data)

    results[[this_iter]] <- est_out

  }
  # concatenate iterations
  results_out <- bind_rows(results, .id = "sim_iter")
  return(results_out)
})


# save results to file
names(sim_results) <- paste("n", n_obs, sep = "_")
timestamp <- str_replace_all(Sys.time(), " ", "_")
saveRDS(object = sim_results,
        file = here("sandbox/data", paste0("CVtreeMLE_", "run_1", ".rds")))
