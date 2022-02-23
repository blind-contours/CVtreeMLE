# be sure to set this env variable by `export R_LIBDIR=/path/to/your/R/libs`
r_libdir <- Sys.getenv("R_LIBDIR")

# set user-specific package library
if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths(r_libdir)
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
}

# packages
library(here)
library(foreach)
library(future)
library(doFuture)
library(doRNG)
library(data.table)
library(tidyverse)
library(origami)
library(sl3)
library(dplyr)

devtools::load_all(here())

registerDoFuture()
plan(multiprocess)

# simulation parameters
set.seed(7259)
n_sim <- 2 # number of simulations
n_obs <- (cumsum(rep(sqrt(40), 6))^2)[-1] # sample sizes at root-n scale
# n_obs <- n_obs[1:3]
truth <- c(1, 1, 1, 1, 1, 1, 6, 7)
true_rule <- "M2 > 2.5 & M3 > 3.6"

source(here("sandbox", "01_setup_data.R"))
source(here("sandbox", "02_fit_estimators.R"))

P_0_data <- generate_P_0(truth = truth, n = 100000)

# perform simulation across sample sizes
sim_results <- lapply(n_obs, function(sample_size) {
  # get results in parallel
  results <- foreach(this_iter = seq_len(n_sim),
                     .options.multicore = list(preschedule = FALSE),
                     .errorhandling = "remove") %do% {
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
