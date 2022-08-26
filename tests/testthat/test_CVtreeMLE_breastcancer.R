library(testthat)
library(CVtreeMLE)
library(xml2)
library(sl3)
library(dplyr)

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Use multicore parallelization to speed up processing.
  future::plan("multiprocess", workers = 2)
}

data(BreastCancer, package = "mlbench")
data <- BreastCancer
names(data) <- tolower(names(data))

seed <- 3349992

# Create a numeric outcome variable.
data$y <- as.numeric(data$class == "malignant")
table(data$y)

y <- "y"

w <- c("marg.adhesion", "epith.c.size", "bare.nuclei",
       "bl.cromatin", "normal.nucleoli", "mitoses")

a <- c("cl.thickness", "cell.size", "cell.shape")

data <- data.frame(lapply(data, function(x) as.numeric(as.character(x))))

data <- data %>%
  dplyr::mutate_if(is.numeric, ~tidyr::replace_na(., mean(., na.rm = TRUE)))

data <- data[colnames(data) != "class"]

ptm <- proc.time()

breast_cancr_results <- CVtreeMLE(
  data = data,
  w = w,
  a = a,
  y = y,
  n_folds = 4,
  seed = seed,
  family = "binomial",
  parallel = TRUE,
  num_cores = 2,
  max_iter = 3
)

## test mixture result outputs given interaction exist in this data:
expect_true(
  class(breast_cancr_results$`V-Specific Mix Results`) == "list")

expect_true(
  class(breast_cancr_results$`Pooled TMLE Mixture Results`) == "data.frame")

expect_true(
  dim(breast_cancr_results$`Pooled TMLE Mixture Results`)[1] ==
    length(breast_cancr_results$`Mixture Data List`)
  )


proc.time() - ptm
