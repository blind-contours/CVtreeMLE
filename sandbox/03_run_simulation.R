# # be sure to set this env variable by `export R_LIBDIR=/path/to/your/R/libs`
# r_libdir <- Sys.getenv("R_LIBDIR")
#
# # set user-specific package library
# if (grepl("savio2", Sys.info()["nodename"])) {
#   .libPaths(r_libdir)
#   Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
# }

# packages
library(here)
library(data.table)
library(tidyverse)
library(origami)
library(sl3)
library(dplyr)
library(foreach)

Sys.unsetenv("GITHUB_PAT")
devtools::install(upgrade = "never")
library(CVtreeMLE)
# simulation parameters
seed <- set.seed(7259)
n_sim <- 25 # number of simulations
n_obs <- (cumsum(rep(sqrt(40), 6))^2)[-1] # sample sizes at root-n scale
# n_obs <- n_obs[1:3]
truth <- c(1, 1, 1, 1, 1, 1, 1, 7)
true_rule <- "M1 > 1.0 & M2 > 2.5 & M3 > 3.5"

source(here("sandbox", "02_fit_estimators.R"))

P_0_data <- simulate_mixture_cube(n_obs = 100000,
                                  subspace_assoc_strength_betas = truth)

# perform simulation across sample sizes
sim_results <- lapply(n_obs, function(sample_size) {
  # get results in parallel
  results <- foreach(this_iter = seq_len(n_sim)) %do% {
    gc()
    data_sim <-  P_0_data[sample(nrow(P_0_data), sample_size), ]
    est_out <- fit_estimators(data = data_sim, true_rule)
    return(est_out)
  }
  # concatenate iterations
  results_out <- bind_rows(results, .id = "sim_iter")
  return(results_out)
})


# save results to file
names(sim_results) <- paste("n", n_obs, sep = "_")
timestamp <- str_replace_all(Sys.time(), " ", "_")
saveRDS(object = sim_results,
        file = here("sandbox/data", paste0("CVtreeMLE_", timestamp, ".rds")))
