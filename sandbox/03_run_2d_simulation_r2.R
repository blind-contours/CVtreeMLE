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
n_sim <- 25 # number of simulations
n_obs <- c(200, 350, 500, 750, 1000, 1500, 2000, 3000)

# Establish globals ---------------------------

# number of observations to create
n <- 500000

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

# Generate simulated data -----------------
P_0_sim <- gen_covariates(n) %>% # gen covariates
  gen_multinom(exposure_grid = exposure_grid) %>% # assign obs to exposure cube
  assign_outcomes(exposure_grid = exposure_grid,
                  c_matrix = c_matrix)  # assign outcome based on cube

P_0_data <- P_0_sim$data
P_0_data_filt <- P_0_data[complete.cases(P_0_data) ,]

# Calculate Empirical Truth -------------------------

data <- P_0_data_filt
m1_seq <- seq(1:5)
m2_seq <- seq(1:5)
m1_dir <- c(">=", "<=")
m2_dir <- c(">=", "<=")

result_list <- list()
iter <- 1
for (i in m1_seq) {
  for (j in m2_seq) {
    for (k in m1_dir) {
      for (l in m2_dir) {
        rule <- (paste("region1 ",k, i, "&", "region2", l, j))
        ate_res <- calc_empir_truth(data, rule)
        result_list[[iter]] <- c(rule, ate_res)
        iter <- iter + 1

      }
    }
  }
}

result_df <- as.data.frame(do.call(rbind,result_list))
colnames(result_df) <- c("Rule", "EY_region", "EY_compl", "ATE")
result_df$ATE <- as.numeric(result_df$ATE)
p0_max_ate <- result_df[which.max(result_df$ATE),]
true_rule <- p0_max_ate$Rule
true_ate <- p0_max_ate$ATE
# perform simulation across sample sizes
sim_results_df <- data.frame()

for (sample_size in n_obs) {
  # get results in parallel
  results <- list()

  for(this_iter in seq_len(n_sim)) {

    seed <- sample(1:10000,1)
    set.seed(seed)

    data_sim <-  P_0_data_filt %>%
      slice_sample(n = sample_size)

    est_out <- fit_estimators(data = as.data.frame(data_sim),
                              covars = c("age", "sex", "bmi"),
                              exposures = c("region1", "region2"),
                              outcome = "outcome_obs",
                              seed = seed,
                              P_0_data = P_0_data_filt,
                              true_rule = true_rule,
                              true_ate = true_ate)

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
        file = here("sandbox/data", paste0("CVtreeMLE_", "run_2", ".rds")))
