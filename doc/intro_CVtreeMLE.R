## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/vignette-"
)

## ----setup--------------------------------------------------------------------
library(data.table)
library(CVtreeMLE)
library(sl3)
library(kableExtra)
library(dplyr)

seed <- 5454432
set.seed(seed)

## ----simulate data, warning=FALSE---------------------------------------------
NIEHS_data <- NIEHS_data_1

head(NIEHS_data) %>%
  kableExtra::kbl(caption = "NIEHS Data") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----add covariates, warning=FALSE--------------------------------------------
NIEHS_data$Z2 <- rbinom(nrow(NIEHS_data),
                                        size = 1,
                                        prob = 0.3)

NIEHS_data$Z3 <- rbinom(nrow(NIEHS_data),
                                        size = 1,
                                        prob = 0.1)

## ----run simulation, warnings = FALSE-----------------------------------------
ptm <- proc.time()

NIEHS_results <- CVtreeMLE(data = as.data.frame(NIEHS_data),
                         w = c("Z", "Z2", "Z3"),
                         a = c(paste("X", seq(7), sep = "")),
                         y = "Y",
                         n_folds = 5,
                         seed = seed,
                         parallel_cv = TRUE,
                         family = "gaussian",
                         num_cores = 6,
                         max_iter = 1)
proc.time() - ptm

## ----pooled mixture results---------------------------------------------------
pooled_mixture_results <- NIEHS_results$`Pooled TMLE Mixture Results`

pooled_mixture_results %>%
  dplyr::filter(Proportion_Folds == 1.0)  %>%
  dplyr::arrange(desc(`Mixture ATE`)) %>%
  kableExtra::kbl(caption = "Pooled TMLE Mixture Results") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----fold specific mixture results--------------------------------------------
vfold_mixture_results <- NIEHS_results$`V-Specific Mix Results`

vfold_mixture_results$X5X7 %>%
  kableExtra::kbl(caption = "X5 and X7 Interaction") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----plot sim_mixture_results, fig.height = 3, fig.width = 8------------------
mixture_plots <- plot_mixture_results(
  v_intxn_results = NIEHS_results$`V-Specific Mix Results`,hjust = 0.8)
mixture_plots$X5X7

## ----marginal results---------------------------------------------------------
pooled_marginal_results <- NIEHS_results$`Pooled TMLE Marginal Results`

pooled_marginal_results %>%
  dplyr::filter(`Proportion in Fold` == 1)  %>%
  kableExtra::kbl(caption = "Pooled Marginal Results") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----marginal refs------------------------------------------------------------
pooled_marginal_refs<- NIEHS_results$`Pooled Marginal Refs`

pooled_marginal_refs %>%
  kableExtra::kbl(caption = "Pooled Marginal Reference Regions") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----marginal_results_example-------------------------------------------------
pooled_marginal_results %>%
  dplyr::filter(`Proportion in Fold` == 1)  %>%
  dplyr::filter(var == "X1")  %>%
  kableExtra::kbl(caption = "Pooled Marginal Results") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----marginal refs example----------------------------------------------------
pooled_marginal_refs %>%
  dplyr::filter(`Variable Quantile` == "X1_1") %>%
  kableExtra::kbl(caption = "Pooled Marginal Reference Regions X1") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----fold marginal results----------------------------------------------------
fold_marginal_results<- NIEHS_results$`V-Specific Marg Results`

fold_marginal_results$`X1_2-X1_1` %>%
  kableExtra::kbl(caption = "V-fold specific results for X1 regions 1 and 2") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

