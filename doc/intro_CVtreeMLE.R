## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(data.table)
library(CVtreeMLE)
library(sl3)
library(kableExtra)
library(qgcomp)
library(dplyr)

set.seed(11249)

## ----simulation inputs--------------------------------------------------------
 # number of observations we want to simulate
n_obs <- 300
# split points for each mixture
splits <- c(0.8, 2.5, 3.6) 
# minimum values for each mixture
mins <- c(0, 0, 0)
# maximum value for each mixture
maxs <- c(3, 4, 5) 
# mu for each mixture
mu <- c(0, 0, 0)
# variance/covariance of mixture variables
sigma <- matrix(c(1, 0.5, 0.8, 0.5, 1, 0.7, 0.8, 0.7, 1), nrow = 3, ncol = 3) 
# subspace probability relationship with covariate W1
w1_betas <- c(0.0, 0.01, 0.03, 0.06, 0.1, 0.05, 0.2, 0.04)
# subspace probability relationship with covariate W2
w2_betas <- c(0.0, 0.04, 0.01, 0.07, 0.15, 0.1, 0.1, 0.04)
# probability of mixture subspace (for multinomial outcome generation)
mix_subspace_betas <- c(0.00, 0.08, 0.05, 0.01, 0.05, 0.033, 0.07, 0.09)
# mixture subspace impact on outcome Y, here the subspace where M1 is lower and M2 and M3 are higher based on values in splits
subspace_assoc_strength_betas <- c(0, 3, 0, 0, 0, 0, 0, 0)
# marginal impact of mixture component on Y
marginal_impact_betas <- c(0, 0, 0) 
# random error
eps_sd <- 0.01 
# if outcome is binary
binary <- FALSE

## ----simulate data, warning=FALSE---------------------------------------------
sim_data <- simulate_mixture_cube(
  n_obs = n_obs, 
  splits = splits,
  mins = mins,
  maxs = maxs,
  mu = mu,
  sigma = sigma,
  w1_betas = w1_betas,
  w2_betas = w2_betas,
  mix_subspace_betas = mix_subspace_betas,
  subspace_assoc_strength_betas = subspace_assoc_strength_betas,
  marginal_impact_betas = marginal_impact_betas,
  eps_sd = eps_sd,
  binary = binary
)

head(sim_data) %>%
  kbl(caption = "Simulated Data") %>%
  kable_classic(full_width = F, html_font = "Cambria")

## ----setup first stack learners-----------------------------------------------
lrnr_glm <- Lrnr_glm$new()
lrnr_bayesglm <- Lrnr_bayesglm$new()
lrnr_gam <- Lrnr_gam$new()

# put all the learners together (this is just one way to do it)
learners <- c(lrnr_glm, lrnr_bayesglm, lrnr_gam)

Q1_stack <- make_learner(Stack, learners)

## ----setup tree stack learners------------------------------------------------
lrnr_glmtree_001 <- Lrnr_glmtree$new(alpha = 0.5, maxdepth = 3)
lrnr_glmtree_002 <- Lrnr_glmtree$new(alpha = 0.6,  maxdepth = 4)
lrnr_glmtree_003 <- Lrnr_glmtree$new(alpha = 0.7, maxdepth = 2)
lrnr_glmtree_004 <- Lrnr_glmtree$new(alpha = 0.8, maxdepth = 1)

learners <- c( lrnr_glmtree_001, lrnr_glmtree_002, lrnr_glmtree_003, lrnr_glmtree_004)
discrete_sl_metalrn <- Lrnr_cv_selector$new()

tree_stack <- make_learner(Stack, learners)

discrete_tree_sl <- Lrnr_sl$new(
  learners = tree_stack,
  metalearner = discrete_sl_metalrn
)


## ----run simulation-----------------------------------------------------------
ptm <- proc.time()

sim_results <- CVtreeMLE(data = sim_data,
                         W = c("W", "W2"),
                         Y = "y",
                         A = c(paste("M", seq(3), sep = "")),
                         back_iter_SL = Q1_stack,
                         tree_SL = discrete_tree_sl, 
                         n_folds = 2,
                         family = "gaussian")
proc.time() - ptm

## ----model RMSE---------------------------------------------------------------
RMSE_results <- sim_results$`Model RMSEs`
RMSE_results %>%
kbl(caption = "Model Fit Results") %>%
kable_classic(full_width = F, html_font = "Cambria")

## ----additive quantile g-comp-------------------------------------------------
 qgcomp_additive_model <- qgcomp(y~M1+M2+M3+W+W2, expnms=c("M1", "M2", "M3"), data=sim_data)
 sqrt(mean((predict(qgcomp_additive_model$fit) - sim_data$y)^2))

## ----multi quantile g-comp----------------------------------------------------
qgcomp_multi_model <- qgcomp(y~M1*M2*M3+W+W2, expnms=c("M1", "M2", "M3"), data=sim_data)
sqrt(mean((predict(qgcomp_multi_model$fit) - sim_data$y)^2))

## ----mixture results----------------------------------------------------------
pooled_mixture_results <- sim_results$`Pooled TMLE Mixture Results`
pooled_mixture_results %>%
kbl(caption = "Pooled TMLE Mixture Results") %>%
kable_classic(full_width = F, html_font = "Cambria")

## ----fold specific results----------------------------------------------------
mixture_v_results <- sim_results$`V-Specific Mix Results`
mixture_v_results$M1M2M3 %>%
kbl(caption = "V-Fold Mixture Results") %>%
kable_classic(full_width = F, html_font = "Cambria")

## ----plot sim_mixture_results, fig.height = 3, fig.width = 8------------------
mixture_plots <- plot_mixture_results(v_intxn_results = sim_results$`V-Specific Mix Results`,hjust = 0.8)
mixture_plots$M1M2M3

## ----marginal results---------------------------------------------------------
pooled_marginal_results <- sim_results$`Pooled TMLE Marginal Results`
pooled_marginal_results %>%
kbl(caption = "Pooled Marginal Results") %>%
kable_classic(full_width = F, html_font = "Cambria")

## ----plot sim marginal results, fig.height = 3, fig.width = 8-----------------
marginal_plots <- plot_marginal_results(v_marginal_results =  sim_results$`V-Specific Marg Results`, mix_comps = c(paste("M", seq(3), sep = "")),hjust = 0.8 )
marginal_plots$M2

