library(testthat)
library(CVtreeMLE)
library(sl3)
library(partykit)
library(pre)
library(xml2)
library(sl3)
library(dplyr)
library(tidyr)

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
  dplyr::mutate_if(is.numeric, ~ replace_na(., mean(., na.rm = TRUE)))

data <- data[colnames(data) != "class"]

ptm <- proc.time()

breast_cancr_results <- CVtreeMLE(
  data = data,
  w = w,
  a = a,
  y = y,
  n_folds = 2,
  seed = seed,
  family = "binomial",
  parallel = TRUE,
  num_cores = 2,
  max_iter = 3
)

proc.time() - ptm

## test mixture result outputs given interaction exist in this data:
expect_true(
  class(breast_cancr_results$`V-Specific Mix Results`) == "list")

expect_true(
  class(breast_cancr_results$`Pooled TMLE Mixture Results`) == "data.frame")


