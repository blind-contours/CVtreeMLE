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

set.seed(11249)

## ----simulate data, warning=FALSE---------------------------------------------
NIEHS_data <- NIEHS_data_1

head(NIEHS_data) %>%
  kbl(caption = "NIEHS Data") %>%
  kable_classic(full_width = F, html_font = "Cambria")

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
                         family = "gaussian",
                         num_cores = 6,
                         max_iter = 1)
proc.time() - ptm

## ----pooled mixture results---------------------------------------------------
pooled_mixture_results <- NIEHS_results$`Pooled TMLE Mixture Results`

pooled_mixture_results %>%
  dplyr::filter(Proportion_Folds == 1.0)  %>%
  kbl(caption = "Pooled TMLE Mixture Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")

## ----fold specific mixture results--------------------------------------------
vfold_mixture_results <- NIEHS_results$`V-Specific Mix Results`

vfold_mixture_results$X1X5 %>%
  kbl(caption = "X1 and X5 Interaction") %>%
  kable_classic(full_width = F, html_font = "Cambria")


## ----plot sim_mixture_results, fig.height = 3, fig.width = 8------------------
mixture_plots <- plot_mixture_results(v_intxn_results = NIEHS_results$`V-Specific Mix Results`,hjust = 0.8)
mixture_plots$X1X5

## ----marginal results---------------------------------------------------------
pooled_marginal_results <- NIEHS_results$`Pooled TMLE Marginal Results`

pooled_marginal_results %>%
  dplyr::filter(`Proportion in Fold` == 1)  %>%
  kbl(caption = "Pooled Marginal Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")

## ----marginal refs------------------------------------------------------------
pooled_marginal_refs<- NIEHS_results$`Pooled Marginal Refs`

pooled_marginal_refs %>%
  kbl(caption = "Pooled Marginal Reference Regions") %>%
  kable_classic(full_width = F, html_font = "Cambria")

