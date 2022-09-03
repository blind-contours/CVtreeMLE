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
library(pre)
library(partykit)
library(kableExtra)
library(dplyr)
library(ggplot2)
library(purrr)

seed <- 5454432
set.seed(seed)

## ----load_NIEHS_data, warning=FALSE-------------------------------------------
niehs_data <- NIEHS_data_1

head(niehs_data) %>%
  kableExtra::kbl(caption = "NIEHS Data") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----add_covariates, warning=FALSE--------------------------------------------
niehs_data$Z2 <- rbinom(nrow(niehs_data),
                                        size = 1,
                                        prob = 0.3)

niehs_data$Z3 <- rbinom(nrow(niehs_data),
                                        size = 1,
                                        prob = 0.1)

## ----run_simulation, warnings = FALSE-----------------------------------------
ptm <- proc.time()

niehs_results <- CVtreeMLE(data = as.data.frame(niehs_data),
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

## ----pooled_mixture_results---------------------------------------------------
pooled_mixture_results <- niehs_results$`Pooled TMLE Mixture Results`

pooled_mixture_results %>%
  dplyr::filter(Proportion_Folds == 1.0)  %>%
  dplyr::arrange(desc(`Mixture ATE`)) %>%
  kableExtra::kbl(caption = "Pooled TMLE Mixture Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----fold_specific_mixture_results--------------------------------------------
vfold_mixture_results <- niehs_results$`V-Specific Mix Results`

vfold_mixture_results$X5X7 %>%
  kableExtra::kbl(caption = "X5 and X7 Interaction") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----mixture_plots,  fig.height = 7, fig.width = 6----------------------------
fold_1_model <- niehs_results$`Mixture Models`[[1]]
plot(fold_1_model, nterms = 9, cex = .5)

## ----consistent_mixture_results-----------------------------------------------
pooled_mixture_results %>%
  dplyr::filter(Proportion_Folds >= 0.8)  %>%
  dplyr::arrange(desc(`Mixture ATE`)) %>%
  kableExtra::kbl(caption = "Consistent Pooled TMLE Mixture Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----plot_mixture_results_dot_whisker, fig.height = 3, fig.width = 8----------
mixture_plots <- plot_mixture_results(
  v_intxn_results = niehs_results$`V-Specific Mix Results`, 
  hjust = 0.8)
mixture_plots$X5X7

## ----marginal_tree_plots, fig.height = 4, fig.width = 7-----------------------
x1_models <- map(niehs_results$`Marginal Models`, "X1")
sapply(x1_models, plot)

## ----fold_marginal_results----------------------------------------------------
fold_marginal_results <- niehs_results$`V-Specific Marg Results`

fold_marginal_results %>%
  filter(Levels == "X1_2-X1_1") %>%
  kableExtra::kbl(
    caption = "V-fold specific results for X1 regions 1 and 2") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----X1_fold_1_marginal_results-----------------------------------------------
fold_marginal_results <- niehs_results$`V-Specific Marg Results`

fold_marginal_results %>%
  filter(fold == "1") %>%
  filter(var == "X1") %>%
  kableExtra::kbl(
    caption = "V-fold specific results for X1 in Fold 1") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----fold_1_X1_tree, fig.height = 4, fig.width = 7----------------------------
plot(x1_models[[1]])

## ----X1_higher_thresholds_reduce_example--------------------------------------
fold_marginal_results %>%
  filter(var == "X1") %>%
  kableExtra::kbl(
    caption = "V-fold specific results for X1") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----pooled_marginal_results--------------------------------------------------
pooled_marginal_results <- niehs_results$`Pooled TMLE Marginal Results`

pooled_marginal_results %>%
  dplyr::filter(`Proportion in Fold` == 1)  %>%
  kableExtra::kbl(caption = "Pooled Marginal Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----marginal_refs------------------------------------------------------------
pooled_marginal_refs <- niehs_results$`Pooled Marginal Refs`

pooled_marginal_refs %>%
  kableExtra::kbl(caption = "Pooled Marginal Reference Regions") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----marginal_results_example-------------------------------------------------
pooled_marginal_results %>%
  dplyr::filter(`Proportion in Fold` == 1)  %>%
  dplyr::filter(var == "X1")  %>%
  kableExtra::kbl(caption = "Pooled Marginal Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----marginal_refs_example----------------------------------------------------
pooled_marginal_refs %>%
  dplyr::filter(`Variable Quantile` == "X1_1") %>%
  kableExtra::kbl(caption = "Pooled Marginal Reference Regions X1") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")

## ----plot_sim_marginal_results_dot_whisker, fig.height = 4, fig.width = 12----
marginal_plots <- plot_marginal_results(
  v_marginal_results = niehs_results$`V-Specific Marg Results`,
  mix_comps = c("X1", "X2", "X3", "X4", "X5", "X6", "X7"),
  hjust = 0.35)

marginal_plots$X1

