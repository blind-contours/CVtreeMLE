library(CVtreeMLE)
library(testthat)
library(sl3)
library(partykit)
library(pre)
library(tidyr)


seed <- 6442
set.seed(seed)
data <- simulate_mixture_cube(subspace_assoc_strength_betas = c(0.1, 0.1, 0.2,
                                                                0.3, 0.4, 0.7,
                                                                0.8, 2.5))


data$M4 <- rnorm(nrow(data))

sls <- create_sls()
w_stack <- sls$W_stack
tree_stack <- sls$A_stack
mix_comps <- c("M1", "M2", "M3", "M4")
w <- c("sex", "age", "bmi")
data$y_scaled <- data$y

example_output <- fit_marg_rule_backfitting(mix_comps = mix_comps,
                                                      at = data,
                                                      w = w,
                                                      w_stack = w_stack,
                                                      tree_stack = tree_stack,
                                                      fold = 1,
                                                      max_iter = 1,
                                                      verbose = FALSE,
                                                      parallel_cv = FALSE,
                                                      seed = 6442)

marginal_df <- example_output$marginal_df

# Mixture variable 4 has no impact and so we expect no rules are found for it:
expect_true(class(marginal_df) == "data.frame")
