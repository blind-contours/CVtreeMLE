library(CVtreeMLE)
library(testthat)
library(sl3)


seed <- 6442
set.seed(seed)
data <- simulate_mixture_cube(subspace_assoc_strength_betas = c(0.1, 0.1, 0.2,
                                                                0.3, 0.4, 0.7,
                                                                0.8, 2.5))

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
                                             max_iter = 2,
                                             verbose = FALSE,
                                             parallel = FALSE,
                                             seed = seed)

expect_true(example_output[
  which.max(example_output$coefficient), "description"] ==
  "M3 > 2.50969709970333 & M1 > 0.947419263385584 & M2 > 1.99541667300285")


