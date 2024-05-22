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
n_obs <- c(500, 750, 1000, 1500, 2000, 5000)
exposure_dim <- 2

# Establish globals ---------------------------

# number of observations to create
n <- 1000000

# number of cuts per exposure
n_cuts <- 5

# exposure grid
exposure_grid <- expand.grid(seq(5), seq(5))
labels <- apply(exposure_grid, 1, paste, collapse = " ")
exposure_grid <- cbind.data.frame(exposure_grid, labels)

# beta matrix
c_matrix <- matrix(c(0.1,0.7,0.8,0.9),
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
        region_res <- calc_empir_truth(data, rule, exposure_dim = 2)
        result_list[[iter]] <- c(rule, region_res)
        iter <- iter + 1

      }
    }
  }
}

result_df <- as.data.frame(do.call(rbind,result_list))
colnames(result_df) <- c("Rule", "EY_PIE", "EY_region")
result_df$EY_PIE <- as.numeric(result_df$EY_PIE)
result_df$EY_region <- as.numeric(result_df$EY_region)

p0_min_ate <- result_df[which.max(result_df$EY_PIE),]
true_rule <- "region1  >= 5 & region2 >= 5"
true_ate <- 38.00052
true_min_ave <- 70.10784

# perform simulation across sample sizes
sim_results_df <- data.frame()
cross_validations <- c(2,4,5,10,10,10, 10)

data$A <- ifelse(data$region1 >= 3 & data$region2 >= 4, 1, 0)
data_1 <- data
data_1$A <- 1

glm_fit <- glm(outcome_true ~ A*sex*age*bmi, data = data)
preds <- predict(glm_fit, newdata = data_1)
true_ate <- mean(preds - data$outcome_true )


for (i in seq_along(n_obs)) {
  sample_size <- n_obs[i]
  cv <- cross_validations[i]
  # get results in parallel
  results <- list()
  print(sample_size)

  for(this_iter in seq_len(n_sim)) {
    print(this_iter)
    stop <- FALSE
    while (stop == FALSE) {

      seed <- sample(1:10000,1)
      set.seed(seed)

      data_sim <-  P_0_data_filt %>%
        slice_sample(n = sample_size)

      n_true <- dim(data_sim[data_sim$Label == "1 1",])
      if (n_true[1] != 0) {
        stop <- TRUE
      }

    }

    est_out <- fit_estimators(data = as.data.frame(data_sim),
                              covars = c("age", "sex", "bmi"),
                              exposures = c("region1", "region2"),
                              outcome = "outcome_true",
                              seed = seed,
                              P_0_data = P_0_data_filt,
                              true_rule = true_rule,
                              true_ate = true_ate,
                              true_region = true_min_ave,
                              exposure_dim = exposure_dim,
                              cv = cv,
                              region = true_rule,
                              min_max = "max")

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
        file = here("sandbox/data", paste0("CVtreeMLE_", "run_13", ".rds")))
