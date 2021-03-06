library(CVtreeMLE)
library(testthat)
library(sl3)



data <- simulate_mixture_cube()

sls <- create_sls()
w_stack <- sls$W_stack
tree_stack <- sls$A_stack
mix_comps <- c("M1", "M2", "M3")
w <- c("w", "w2")
data$y_scaled <- data$y

example_output <- fit_mix_rule_backfitting(at = data,
                                             a = mix_comps,
                                             w = w,
                                             y = "y_scaled",
                                             w_stack = w_stack,
                                             fold = 1,
                                             max_iter = 1,
                                             verbose = FALSE,
                                             parallel = FALSE,
                                             seed = 6442)

expect_true(example_output[
  which.max(example_output$coefficient), "test"] == "M1M2M3")
