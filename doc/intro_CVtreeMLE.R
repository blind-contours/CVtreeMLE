## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/vignette_"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(data.table)
library(CVtreeMLE)
library(sl3)
library(kableExtra)
library(dplyr)
library(ggplot2)
library(purrr)

seed <- 5454433
set.seed(seed)

## ----load_NIEHS_data, warning=FALSE-------------------------------------------
niehs_data <- NIEHS_data_1

head(niehs_data) %>%
  kableExtra::kbl(caption = "NIEHS Data") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----add_covariates, warning=FALSE--------------------------------------------
niehs_data$Z2 <- rbinom(nrow(niehs_data),
  size = 1,
  prob = 0.3
)

niehs_data$Z3 <- rbinom(nrow(niehs_data),
  size = 1,
  prob = 0.1
)

## ----run_simulation, warnings = FALSE-----------------------------------------
ptm <- proc.time()

niehs_results <- CVtreeMLE(
  data = as.data.frame(niehs_data),
  w = c("Z", "Z2", "Z3"),
  a = c(paste("X", seq(7), sep = "")),
  y = "Y",
  n_folds = 5,
  seed = seed,
  parallel_cv = TRUE,
  parallel = TRUE,
  family = "continuous",
  num_cores = 8,
  min_max = "min",
  min_obs = 25
)
proc.time() - ptm

## ----pooled_mixture_results---------------------------------------------------
pooled_mixture_results <- niehs_results$`Oracle Region Results`

pooled_mixture_results %>%
  kableExtra::kbl(caption = "Oracle Mixture Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----region_specific_pooling--------------------------------------------------
region_specific_pooling <- niehs_results$`Pooled TMLE Mixture Results`

region_specific_pooling %>%
  kableExtra::kbl(caption = "Region Specific Mixture Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----k_fold_results-----------------------------------------------------------
k_fold_results <- niehs_results$`V-Specific Mix Results`

k_fold_results %>%
  kableExtra::kbl(caption = "K-fold Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

