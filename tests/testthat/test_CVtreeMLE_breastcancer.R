library(testthat)
library(CVtreeMLE)
library(xml2)
library(sl3)

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Use multicore parallelization to speed up processing.
  future::plan("multiprocess", workers = 2)
}
#################################
# mlbench BreastCancer dataset.


data(BreastCancer, package = "mlbench")
data <- BreastCancer
names(data) <- tolower(names(data))

set.seed(3, "L'Ecuyer-CMRG")

# Reduce to a dataset of 200 observations to speed up testing.

# Create a numeric outcome variable.
data$y <- as.numeric(data$class == "malignant")
table(data$y)

#################################
# set up learners
lrnr_glm <- Lrnr_glm$new()
lrnr_gam <- Lrnr_gam$new()
lrnr_ranger <- Lrnr_ranger$new()

# put all the learners together (this is just one way to do it)
learners <- c(
  lrnr_glm,
  lrnr_gam,
  lrnr_ranger
)

Q1_stack <- make_learner(Stack, learners)

lrnr_glmtree_001 <- Lrnr_glmtree$new(alpha = 0.5, maxdepth = 3)
lrnr_glmtree_002 <- Lrnr_glmtree$new(alpha = 0.6, maxdepth = 4)
lrnr_glmtree_003 <- Lrnr_glmtree$new(alpha = 0.7, maxdepth = 2)
lrnr_glmtree_004 <- Lrnr_glmtree$new(alpha = 0.8, maxdepth = 1)

learners <- c(lrnr_glmtree_001, lrnr_glmtree_002, lrnr_glmtree_003, lrnr_glmtree_004)
discrete_sl_metalrn <- Lrnr_cv_selector$new()

tree_stack <- make_learner(Stack, learners)

discrete_tree_sl <- Lrnr_sl$new(
  learners = tree_stack,
  metalearner = discrete_sl_metalrn
)

# we will arbitrarily investigate cell qualities as a mixture and the other variables as covariates

W <- c("marg.adhesion", "epith.c.size", "bare.nuclei", "bl.cromatin", "normal.nucleoli", "mitoses")
A <- c("cl.thickness", "cell.size", "cell.shape")

ptm <- proc.time()

breast_cancr_results <- CVtreeMLE(
  data = data,
  W = W,
  Y = "y",
  A = A,
  back_iter_SL = Q1_stack,
  tree_SL = discrete_tree_sl,
  n_folds = 2,
  family = "binomial",
  H.AW_trunc_lvl = 10,
  parallel = TRUE,
  num_cores = 2,
  max_iter = 3,
  verbose = FALSE
)

proc.time() - ptm
