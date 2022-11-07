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

Sys.unsetenv("GITHUB_PAT")
devtools::install(upgrade = "never")

# simulation parameters
n_sim <- 10 # number of simulations
n_obs <- c(200, 350, 500, 750, 1000, 1500, 2000, 3000, 5000)

# Establish globals ---------------------------

# number of observations to create
n <- 500000


P_0_data <- simulate_mixture_cube(n_obs = n,
                                  subspace_assoc_strength_betas = c(0.1, 0.2, 0.2, 5,
                                                                    0.5, 0.4, 0.3, 0.1))
true_rule <- "M3 <= 2.5 & M1 >= 0.99 & M2 >= 2.0"

P_0_data_ate_eval <- P_0_data %>%
  dplyr::mutate(A = ifelse(eval(parse(text = true_rule)), 1, 0))

region_ave_outcome <- mean(subset(P_0_data_ate_eval, A == 1)$y)
nonregion_ave_outcome <- mean(subset(P_0_data_ate_eval, A == 0)$y)
true_ate <- region_ave_outcome - nonregion_ave_outcome

exposure_dim <- 3
sim_results_df <- data.frame()

for (sample_size in n_obs) {
  # get results in parallel
  results <- list()
  print(sample_size)

  for(this_iter in seq_len(n_sim)) {
    seed <- sample(1:10000,1)
    set.seed(seed)
    print(this_iter)

    data_sim <-  P_0_data %>%
      slice_sample(n = sample_size)

    est_out <- fit_estimators(data = as.data.frame(data_sim),
                              covars = c("age", "sex", "bmi"),
                              exposures = c("M1", "M2", "M3"),
                              outcome = "y",
                              seed = seed,
                              P_0_data = P_0_data,
                              true_rule = true_rule,
                              true_ate = true_ate,
                              target_var_set = "M1-M2-M3",
                              exposure_dim = exposure_dim)

    est_out$n_obs <- sample_size

    results[[this_iter]] <- est_out

  }
  # concatenate iterations
  results_out <- bind_rows(results, .id = "sim_iter")
  sim_results_df <- rbind(sim_results_df, results_out)
}

# save results to file
timestamp <- str_replace_all(Sys.time(), " ", "_")
saveRDS(object = sim_results_df,
        file = here("sandbox/data", paste0("3D_CVtreeMLE_", "run_3", ".rds")))
