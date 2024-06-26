## ----include = FALSE----------------------------------------------------------
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

## ----NHANES data and variables------------------------------------------------
data("NHANES_eurocim")

exposures <- c("LBX074LA", # PCB74 Lipid Adj (ng/g)
               "LBX099LA", # PCB99 Lipid Adj (ng/g)
               "LBX118LA", # PCB118 Lipid Adj (ng/g)
               "LBX138LA", # PCB138 Lipid Adj (ng/g)
               "LBX153LA", # PCB153 Lipid Adj (ng/g)
               "LBX170LA", # PCB170 Lipid Adj (ng/g)
               "LBX180LA", # PCB180 Lipid Adj (ng/g)
               "LBX187LA", # PCB187 Lipid Adj (ng/g)
               "LBX194LA", # PCB194 Lipid Adj (ng/g)
               "LBXD03LA", # 1,2,3,6,7,8-hxcdd Lipid Adj (pg/g)
               "LBXD05LA", # 1,2,3,4,6,7,8-hpcdd Lipid Adj (pg/g)
               "LBXD07LA", # 1,2,3,4,6,7,8,9-ocdd Lipid Adj (pg/g)
               "LBXF03LA", # 2,3,4,7,8-pncdf Lipid Adj (pg/g)
               "LBXF04LA", # 1,2,3,4,7,8-hxcdf Lipid Adj (pg/g)
               "LBXF05LA", # 1,2,3,6,7,8-hxcdf Lipid Adj (pg/g)
               "LBXF08LA", # 1,2,3,4,6,7,8-hxcdf Lipid Adj (pg/g)
               "LBXHXCLA", # 3,3',4,4',5,5'-hxcb Lipid Adj (pg/g)
               "LBXPCBLA") # 3,3',4,4',5-pcnb Lipid Adj (pg/g)

NHANES_eurocim <- NHANES_eurocim[complete.cases(NHANES_eurocim[, exposures]), ]

outcome <- "TELOMEAN"

covariates <- c("LBXWBCSI", # White blood cell count (SI)
                "LBXLYPCT", # Lymphocyte percent (%)
                "LBXMOPCT", # Monocyte percent (%)
                "LBXEOPCT", # Eosinophils percent (%)
                "LBXBAPCT", # Basophils percent (%)
                "LBXNEPCT", # Segmented neutrophils percent (%)
                "male", # Sex
                "age_cent", # Age at Screening, centered
                "race_cat", # race
                "bmi_cat3", # Body Mass Index (kg/m**2)
                "ln_lbxcot", # Cotinine (ng/mL), log-transformed
                "edu_cat") # Education Level - Adults 20+



## ----remove correlated exposures----------------------------------------------
# Calculate the correlation matrix for the exposures
cor_matrix <- cor(NHANES_eurocim[, exposures], use = "complete.obs")

# Set a threshold for high correlation
threshold <- 0.8

# Find pairs of highly correlated exposures
highly_correlated_pairs <- which(abs(cor_matrix) > threshold & lower.tri(cor_matrix), arr.ind = TRUE)

# Initiate a vector to keep track of exposures to remove
exposures_to_remove <- c()

# Loop through the highly correlated pairs and decide which exposure to remove
for (pair in seq_len(nrow(highly_correlated_pairs))) {
  row <- highly_correlated_pairs[pair, "row"]
  col <- highly_correlated_pairs[pair, "col"]

  if (!(colnames(cor_matrix)[row] %in% exposures_to_remove)) {
    exposures_to_remove <- c(exposures_to_remove, colnames(cor_matrix)[row])
  }
}

# Keep only uncorrelated exposures
exposures_to_keep <- setdiff(exposures, exposures_to_remove)


## ----run CVtreeMLE for NHANES-------------------------------------------------

# Convert continuous X variables to their corresponding deciles
NHANES_eurocim[, exposures_to_keep] <- apply(NHANES_eurocim[, exposures_to_keep], 2, round, 1)

nhanes_results <- CVtreeMLE(
  data = NHANES_eurocim, ## dataframe
  w = covariates,
  a = exposures_to_keep,
  y = outcome,
  n_folds = 5,
  seed = 34421,
  parallel_cv = TRUE,
  parallel = TRUE,
  family = "continuous",
  num_cores = 8,
  min_max = "max",
  min_obs = 20,
  p_val_thresh = 1
)

## ----k_fold_results-----------------------------------------------------------
k_fold_results <- nhanes_results$`V-Specific Mix Results`

k_fold_results %>%
  kableExtra::kbl(caption = "K-fold Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----pooled_mixture_results---------------------------------------------------
pooled_mixture_results <- nhanes_results$`Oracle Region Results`

pooled_mixture_results %>%
  kableExtra::kbl(caption = "Oracle Mixture Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----region_specific_pooling--------------------------------------------------
region_specific_pooling <- nhanes_results$`Pooled TMLE Mixture Results`

region_specific_pooling %>%
  kableExtra::kbl(caption = "Region Specific Mixture Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

