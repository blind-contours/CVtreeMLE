library(CVtreeMLE)
library(testthat)
library(sl3)
library(pre)
library(partykit)

sim_data <- simulate_mixture_cube(n_obs = 200)

ptm <- proc.time()

sim_results <- CVtreeMLE(
  data = sim_data,
  w = c("age", "sex", "bmi"),
  a = c(paste("M", seq(3), sep = "")),
  y = "y",
  n_folds = 2,
  parallel_cv = TRUE,
  seed = 2333,
  parallel_type = "multi_session",
  family = "continuous",
  num_cores = 2
)
proc.time() - ptm

## the fit marginal param is set to false so this should be NULL
expect_true(is.null(sim_results$`Pooled TMLE Marginal Results`))

## the true data has a triple interaction, make sure we find it in small samples
expect_true("M1-M2-M3" %in% sim_results$`Pooled TMLE Mixture Results`$Vars)

## the true ARE is 6 - ensure our max region estimate is close to truth for
## the simulation
expect_equal(max(sim_results$`Pooled TMLE Mixture Results`$`Mixture ATE`),
  6,
  tolerance = 0.2
)
