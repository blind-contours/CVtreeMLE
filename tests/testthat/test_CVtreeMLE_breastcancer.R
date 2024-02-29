library(testthat)
library(CVtreeMLE)
library(sl3)
library(partykit)
library(pre)
library(xml2)
library(sl3)
library(dplyr)

data(BreastCancer, package = "mlbench")
data <- BreastCancer
names(data) <- tolower(names(data))

seed <- 3349992

# Create a numeric outcome variable.
data$y <- as.numeric(data$class == "malignant")
table(data$y)

y <- "y"

w <- c(
  "marg.adhesion", "epith.c.size",
  "bl.cromatin", "normal.nucleoli", "mitoses"
)

a <- c("cl.thickness", "cell.size", "cell.shape")

data <- suppressWarnings(data.frame(lapply(
  data,
  function(x) {
    as.numeric(
      as.character(x)
    )
  }
)))

data <- data[colnames(data) != "class"]

ptm <- proc.time()
breast_cancr_results <- CVtreeMLE(
  data = data,
  w = w,
  a = a,
  y = y,
  n_folds = 5,
  seed = seed,
  family = "binomial",
  parallel = TRUE,
  num_cores = 2,
)
  proc.time() - ptm

## test to make sure the RMSE table for ensemble specific trees has the same
## trees as the v-fold results, basically that v-fold results are estimated
## for each data-adaptively identified tree

expect_true(
  all(names(
    breast_cancr_results$`V-Specific Mix Results`
  )
  %in% breast_cancr_results$`Model RMSEs`$`Var(s)`)
)

## test mixture result outputs given interaction exist in this data:
expect_true(
  nrow(breast_cancr_results$`Oracle Region Results`) == 1
)

expect_true(
  class(breast_cancr_results$`Pooled TMLE Mixture Results`) == "data.frame"
)

