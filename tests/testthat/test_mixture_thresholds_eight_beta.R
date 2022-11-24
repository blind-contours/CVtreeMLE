library(CVtreeMLE)
library(testthat)
library(sl3)
library(pre)
library(partykit)


seed <- 6442
set.seed(seed)
data <- simulate_mixture_cube(
  n_obs = 1000,
  subspace_assoc_strength_betas = c(
    0.1, 0.1, 0.2,
    0.3, 0.4, 0.7,
    0.8, 10
  )
)

sls <- create_sls()
w_stack <- sls$W_stack
tree_stack <- sls$A_stack
mix_comps <- c("M1", "M2", "M3")
w <- c("age", "sex", "bmi")
data$y_scaled <- data$y

example_output <- fit_mix_rule_backfitting(
  at = data,
  a = mix_comps,
  w = w,
  y = "y_scaled",
  direction = "positive",
  w_stack = w_stack,
  fold = 1,
  max_iter = 4,
  verbose = FALSE,
  parallel = FALSE,
  seed = seed
)

mixture_rules <- example_output$rules

expect_true(mixture_rules[
  which.max(mixture_rules$coefficient), "test"
] ==
  "M1-M2-M3")
