library(testthat)
library(CVtreeMLE)
library(xml2)
library(sl3)

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Use multicore parallelization to speed up processing.
  future::plan("multiprocess", workers = 2)
}

data(BreastCancer, package = "mlbench")
data <- BreastCancer
names(data) <- tolower(names(data))

set.seed(3, "L'Ecuyer-CMRG")

# Create a numeric outcome variable.
data$y <- as.numeric(data$class == "malignant")
table(data$y)

y <- "y"

w <- c("marg.adhesion", "epith.c.size", "bare.nuclei",
       "bl.cromatin", "normal.nucleoli", "mitoses")

a <- c("cl.thickness", "cell.size", "cell.shape")

ptm <- proc.time()

breast_cancr_results <- CVtreeMLE(
  data = data,
  w = w,
  a = a,
  y = y,
  n_folds = 4,
  family = "binomial",
  parallel = TRUE,
  num_cores = 2,
  max_iter = 3
)

proc.time() - ptm
