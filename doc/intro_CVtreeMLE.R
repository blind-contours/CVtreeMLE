## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(data.table)
library(devtools)
library(CVtreeMLE)
devtools::load_all('~/sl3')
library(kableExtra)
library(qgcomp)
library(dplyr)

set.seed(11249)

## ----simulation inputs--------------------------------------------------------
n_obs <- 500 # number of observations we want to simulate
splits <- c(0.8, 2.5, 3.6) # split points for each mixture
mins <- c(0, 0, 0) # minimum values for each mixture
maxs <- c(3, 4, 5) # maximum value for each mixture
mu <- c(0, 0, 0) # mu for each mixture
sigma <- matrix(c(1, 0.5, 0.8, 0.5, 1, 0.7, 0.8, 0.7, 1), nrow = 3, ncol = 3) # variance/covariance of mixture variables
w1_betas <- c(0.0, 0.01, 0.03, 0.06, 0.1, 0.05, 0.2, 0.04) # subspace probability relationship with covariate W1
w2_betas <- c(0.0, 0.04, 0.01, 0.07, 0.15, 0.1, 0.1, 0.04) # subspace probability relationship with covariate W2
mix_subspace_betas <- c(0.00, 0.08, 0.05, 0.01, 0.05, 0.033, 0.07, 0.09) # probability of mixture subspace (for multinomial outcome generation)
subspace_assoc_strength_betas <- c(0, 3, 0, 0, 0, 0, 0, 0) # mixture subspace impact on outcome Y, here the subspace where M1 is lower and M2 and M3 are higher based on values in splits
marginal_impact_betas <- c(0, 0, 0) # marginal impact of mixture component on Y
eps_sd <- 0.01 # random error
binary <- FALSE # if outcome is binary

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
  kableExtra::kbl(caption = "Simulated Data") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----setup first stack learners-----------------------------------------------
lrnr_glm <- Lrnr_glm$new()
lrnr_bayesglm <- Lrnr_bayesglm$new()
lrnr_gam <- Lrnr_gam$new()
lrnr_ranger <- Lrnr_ranger$new()

# put all the learners together (this is just one way to do it)
learners <- c(lrnr_glm, lrnr_bayesglm, lrnr_gam, lrnr_ranger)

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
                         n_folds = 3,
                         family = "gaussian",
                         H.AW_trunc_lvl = 10,
                         parallel = TRUE,
                         num_cores = 2,
                         max_iter = 5,
                         verbose = TRUE)
proc.time() - ptm

## ----model RMSE---------------------------------------------------------------
RMSE_results <- sim_results$`Model RMSEs`
RMSE_results %>%
kableExtra::kbl(caption = "Model Fit Results") %>%
kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----additive quantile g-comp-------------------------------------------------
 qgcomp_additive_model <- qgcomp(y~M1+M2+M3+W+W2, expnms=c("M1", "M2", "M3"), data=sim_data)
 sqrt(mean((predict(qgcomp_additive_model$fit) - sim_data$y)^2))

## ----multi quantile g-comp----------------------------------------------------
qgcomp_multi_model <- qgcomp(y~M1*M2*M3+W+W2, expnms=c("M1", "M2", "M3"), data=sim_data)
sqrt(mean((predict(qgcomp_multi_model$fit) - sim_data$y)^2))

## ----mixture results----------------------------------------------------------
pooled_mixture_results <- sim_results$`Pooled TMLE Mixture Results`
pooled_mixture_results %>%
  kableExtra::kbl(caption = "Pooled TMLE Mixture Results") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----fold specific results----------------------------------------------------
mixture_v_results <- sim_results$`V-Specific Mix Results`
mixture_v_results$M1M2M3 %>%
  kableExtra::kbl(caption = "V-Fold Mixture Results") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----plot sim_mixture_results, fig.height = 3, fig.width = 8------------------
mixture_plots <- plot_mixture_results(v_intxn_results = sim_results$`V-Specific Mix Results`)
mixture_plots$M1M2M3

## ----marginal results---------------------------------------------------------
pooled_marginal_results <- sim_results$`Pooled TMLE Marginal Results`
pooled_marginal_results %>%
  kableExtra::kbl(caption = "Pooled Marginal Results") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----plot sim marginal results, fig.height = 3, fig.width = 8-----------------
marginal_plots <- plot_marginal_results(v_marginal_results =  sim_results$`V-Specific Marg Results`, mix_comps = c(paste("M", seq(3), sep = "")))
marginal_plots$M1

## ----metals example-----------------------------------------------------------
data("NIEHS_data_1", package="CVtreeMLE")

## ----run NIEHS----------------------------------------------------------------
ptm <- proc.time()
NIEHS_data_1 <- as.data.frame(NIEHS_data_1)
NIEH_1_results <- CVtreeMLE(data = NIEHS_data_1,
                            W = "Z",
                            Y = "Y",
                            A = c("X1", "X2", "X3", "X4", "X5", "X6", "X7"),
                            back_iter_SL = Q1_stack,
                            tree_SL = discrete_tree_sl, 
                            n_folds = 2,
                            family = "gaussian",
                            H.AW_trunc_lvl = 10,
                            parallel = TRUE,
                            num_cores = 2,
                            max_iter = 5,
                            verbose = FALSE)
proc.time() - ptm

## ----NIEHS RMSE---------------------------------------------------------------
NIEH_1_RMSE <- NIEH_1_results$`Model RMSEs`
NIEH_1_RMSE %>%
  kableExtra::kbl(caption = "NIEH Data 1 Model Fit Results") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----quant_sum NIEHS----------------------------------------------------------
NIEHS_data_1_model <- qgcomp(Y~., expnms=c("X1", "X2", "X3", "X4", "X5", "X6", "X7"), data=NIEHS_data_1)
sqrt(mean((predict(NIEHS_data_1_model$fit) - NIEHS_data_1$Y)^2))

## ----NIEH mixture results-----------------------------------------------------
NIEH_1_intxn_results <- NIEH_1_results$`Pooled TMLE Mixture Results`
NIEH_1_intxn_results %>%
  kableExtra::kbl(caption = "NIEH Mixture Results") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----NIEH v-fold mixture results----------------------------------------------
NIEH_v_intxn_results <- NIEH_1_results$`V-Specific Mix Results`
NIEH_v_intxn_results[[1]] %>%
  kableExtra::kbl(caption = "v-fold NIEH Mixture Results") %>%
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

## ----plot marginal results----------------------------------------------------
marginal_plots <- plot_marginal_results(v_marginal_results = NIEH_1_results$`V-Specific Marg Results`, mix_comps =  c("X1", "X2", "X3", "X4", "X5", "X6", "X7"))
names(marginal_plots)


## ----NIEH marginal plot example, fig.height = 3, fig.width = 8----------------
marginal_plots$X1

## ----plot interaction results, fig.height = 3, fig.width = 8------------------
mixture_plots <- plot_mixture_results(v_intxn_results = NIEH_1_results$`V-Specific Mix Results`)
mixture_plots$X2X5
