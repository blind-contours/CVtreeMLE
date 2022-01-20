
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `CVtreeMLE` <img src="inst/figures/CVtreeMLE_sticker.png" height="300" align="right"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/blind-contours/CVtreeMLE/workflows/R-CMD-check/badge.svg)](https://github.com/blind-contours/CVtreeMLE/actions)
[![Coverage
Status](https://codecov.io/gh/blind-contours/CVtreeMLE/branch/main/graph/badge.svg?token=HJP5PYQSG4)](https://codecov.io/github/blind-contours/CVtreeMLE?branch=master)
[![CRAN](https://www.r-pkg.org/badges/version/CVtreeMLE)](https://www.r-pkg.org/pkg/CVtreeMLE)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/CVtreeMLE)](https://CRAN.R-project.org/package=CVtreeMLE)
[![CRAN total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/CVtreeMLE)](https://CRAN.R-project.org/package=CVtreeMLE)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070042.svg)](https://doi.org/10.5281/zenodo.4070042)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02447/status.svg)](https://doi.org/10.21105/joss.02447)
<!-- badges: end -->

> Efficient Estimation of the Causal Effects of Joint Exposure using
> Data Adaptive Decision Trees and Cross-Validated Targeted Maximum
> Likelihood Estimation **Authors:** [David
> McCoy](https://davidmccoy.org)

------------------------------------------------------------------------

## What is `CVtreeMLE`?

The `CVtreeMLE` (Cross-Validated Decision Trees with Targeted Maximum
Likelihood Estimation) R package is designed to provide statistical
software for the construction of efficient estimators from data adaptive
decision trees. The target parameter is the average treatment effect
(ATE) causal parameter defined as the counterfactual mean outcome if all
individuals were jointly exposed to a combination of exposure levels in
a mixed exposure compared to if all individuals were not exposed. Here,
the levels of a joint exposure are data-adaptively identified based on
decision trees applied to a set of exposure variables while flexibly
controlling for covariates non-parametrically. For more information on
data- adaptive parameters see (Hubbard, Kherad-Pajouh, and Van Der Laan
2016). `CVtreeMLE` uses data-adaptive parameters by implementing V-fold
cross-validation (CV), that is, in 10-fold CV, the data is split 10
times (folds), where 90% of the data is used to determine rules in a
mixture, and the *g* and *Q* estimators needed for the ATE. These rules
and estimators created in training data are applied to the validation
data in order to calculate the final ATE target parameter. In order to
optimize the optimum bias-variance trade-off for our causal parameter of
interest we use cross-validated targeted minimum loss based estimation
(CV-TMLE). `CVtreeMLE` builds off of the CV-TMLE general theorem of
cross-validated minimum loss based estimation Zheng and Laan (2010)
which allows the full utilization of loss based super learning to obtain
the initial estimators needed for our target parameter without risk of
overfitting. Thus, `CVtreeMLE` makes possible the non-parametric
estimation of the causal effects of a mixed exposure that both results
in interpretable results which are useful for public policy and is
asymptotically efficient.

`CVtreeMLE` integrates with the
[`sl3`package](https://github.com/tlverse/sl3) (Coyle et al. 2021) to
allow for ensemble machine learning to be leveraged in the estimation
procedure. `sl3` is used to create ensemble machine learning estimators
for the *Q* and *g* mechanisms for the average treatment effect (ATE)
target parameter and is also used in the iterative backfitting procedure
In the iterative backfitting procedure, for an enselbe of decision trees
are fit on the full mixture modeled together, the [`pre`
package](https://github.com/marjoleinF/pre)(Fokkema 2020) is used to fit
rule ensembles. In backfitting procedure to find quantiles in each
mixture component individually, an Super Learner of decision trees
generated from the [`partykit`
package](http://partykit.r-forge.r-project.org/partykit/)\[partykit2015\]
is created. In each case, the goal is to find the best fitting decision
tree from which we extract decision tree rules, we then calculate the
ATE for these rules.

------------------------------------------------------------------------

## Installation

For standard use, we recommend installing the package from
[CRAN](https://CRAN.R-project.org/package=CVtreeMLE) via

``` r
install.packages("CVtreeMLE")
```

*Note:* If `CVtreeMLE` is installed from
[CRAN](https://CRAN.R-project.org/package=CVtreeMLE), the `sl3`, an
enhancing dependency that allows ensemble machine learning to be used
for nuisance parameter estimation, won’t be included. We highly
recommend additionally installing `sl3` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("tlverse/sl3@master")
```

For the latest features, install the most recent *stable version* of
`CVtreeMLE` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("blind-contours/CVtreeMLE@master")
```

To contribute, install the *development version* of `CVtreeMLE` from
GitHub via [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("blind-contours/CVtreeMLE@devel")
```

------------------------------------------------------------------------

## Example

To illustrate how `CVtreeMLE` may be used to ascertain the effect of a
joint exposure, consider the following example:

First load the package and other packages needed

``` r
library(CVtreeMLE)
library(devtools)
#> Loading required package: usethis
load_all("~/sl3")
#> ℹ Loading sl3
library(kableExtra)
library(ggplot2)
library(jtools)


set.seed(429153)
```

Use the `simulate_mixture_cube` function to generate simulated data that
represents ground-truth. Here, we create three continuous mixture
variables, *A*, that are correlated and baseline covariates, *W*, that
are potential confounders. Our outcome will be generated such that
individuals with a specific set of exposures have a different outcome
compared to individuals who are not exposed to this combination of
exposure levels.

![](inst/figures/The_Cube.png) The above figure illustrates the data we
will generate using this function. Here, individuals exposed to
*M*<sub>1</sub> at values less than 1.0, *M*<sub>2</sub> at levels more
than 2.0, and *M*<sub>3</sub> at levels at or greater than 2.5 have an
outcome of 6, compared to individuals not exposed to this combination of
thresholds who have an expected outcome of 0 - thus our ATE is 6. Two
covariates *W* confound this relationship. Let’s simulate this scenario.

## Simulate Data

``` r
n_obs <- 500 # number of observations we want to simulate
splits <- c(0.99, 2.0, 2.5) # split points for each mixture
mins <- c(0, 0, 0) # minimum values for each mixture
maxs <- c(3, 4, 5) # maximum value for each mixture
mu <- c(0, 0, 0) # mu for each mixture
sigma <- matrix(c(1, 0.5, 0.8, 0.5, 1, 0.7, 0.8, 0.7, 1), nrow = 3, ncol = 3) # variance/covariance of mixture variables
w1_betas <- c(0.0, 0.01, 0.03, 0.06, 0.1, 0.05, 0.2, 0.04) # subspace probability relationship with covariate W1
w2_betas <- c(0.0, 0.04, 0.01, 0.07, 0.15, 0.1, 0.1, 0.04) # subspace probability relationship with covariate W2
mix_subspace_betas <- c(0.00, 0.08, 0.05, 0.01, 0.05, 0.033, 0.07, 0.09) # probability of mixture subspace (for multinomial outcome generation)
subspace_assoc_strength_betas <- c(0, 0, 0, 0, 0, 0, 6, 0) # mixture subspace impact on outcome Y, here the subspace where M1 is lower and M2 and M3 are higher based on values in splits
marginal_impact_betas <- c(0, 0, 0) # marginal impact of mixture component on Y
eps_sd <- 0.01 # random error
binary <- FALSE # if outcome is binary
```

Above, the `subspace_assoc_strength_betas` parameter is used to indicate
the subspace we want to use and the expected outcome in that subspace.
The indices correspond to an area in the cube:

1.  All mixtures lower than specified thresholds
2.  M1 is higher but M2 and M3 are lower
3.  M2 is higher but M1 and M3 are lower
4.  M1 and M2 are higher and M3 is lower
5.  M3 is higher and M1 and M2 are lower
6.  M1 and M3 are higher and M2 is lower
7.  M2 and M3 are higher and M1 is lower
8.  All mixtures are higher than thresholds

``` r
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
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Simulated Data
</caption>
<thead>
<tr>
<th style="text-align:right;">
W
</th>
<th style="text-align:right;">
W2
</th>
<th style="text-align:right;">
M1
</th>
<th style="text-align:right;">
M2
</th>
<th style="text-align:right;">
M3
</th>
<th style="text-align:right;">
y
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.3125971
</td>
<td style="text-align:right;">
0.0281233
</td>
<td style="text-align:right;">
0.0918698
</td>
<td style="text-align:right;">
2.1287872
</td>
<td style="text-align:right;">
0.5587179
</td>
<td style="text-align:right;">
0.3523733
</td>
</tr>
<tr>
<td style="text-align:right;">
0.4677867
</td>
<td style="text-align:right;">
0.1344452
</td>
<td style="text-align:right;">
2.2546606
</td>
<td style="text-align:right;">
0.0433708
</td>
<td style="text-align:right;">
1.2521522
</td>
<td style="text-align:right;">
0.6081815
</td>
</tr>
<tr>
<td style="text-align:right;">
-0.7798901
</td>
<td style="text-align:right;">
-0.5718299
</td>
<td style="text-align:right;">
0.1683829
</td>
<td style="text-align:right;">
2.3168107
</td>
<td style="text-align:right;">
0.1473310
</td>
<td style="text-align:right;">
-1.3717913
</td>
</tr>
<tr>
<td style="text-align:right;">
-0.0760508
</td>
<td style="text-align:right;">
-0.6356721
</td>
<td style="text-align:right;">
1.3832386
</td>
<td style="text-align:right;">
0.9429022
</td>
<td style="text-align:right;">
2.9218706
</td>
<td style="text-align:right;">
-0.7252087
</td>
</tr>
<tr>
<td style="text-align:right;">
-0.1238976
</td>
<td style="text-align:right;">
-0.3105393
</td>
<td style="text-align:right;">
0.4341428
</td>
<td style="text-align:right;">
3.0653637
</td>
<td style="text-align:right;">
1.6723064
</td>
<td style="text-align:right;">
-0.4330232
</td>
</tr>
<tr>
<td style="text-align:right;">
0.2149969
</td>
<td style="text-align:right;">
0.1632984
</td>
<td style="text-align:right;">
1.0384095
</td>
<td style="text-align:right;">
2.2392699
</td>
<td style="text-align:right;">
0.1746694
</td>
<td style="text-align:right;">
0.3813563
</td>
</tr>
</tbody>
</table>

## Set up Estimators used in Super Learners

Here, we set up our Super Learner using `SL3` for the iterative
backfitting procedure and for our *Q* and *g* mechanisms. These learners
will fit *Y*\|*W* offset by *Y*\|*A* as we fit decision trees to the
exposure variables both jointly and individially. Once rules are
established, this Super Learner will also estimate the the propensity of
being exposed to the determined rule as well as estimating the outcome
model.

``` r
lrnr_glm <- Lrnr_glm$new()
lrnr_bayesglm <- Lrnr_bayesglm$new()
lrnr_gam <- Lrnr_gam$new()
lrnr_lasso <- Lrnr_glmnet$new(alpha = 1)
lrnr_earth <- Lrnr_earth$new()
lrnr_ranger <- Lrnr_ranger$new()
# put all the learners together (this is just one way to do it)
learners <- c(lrnr_glm, lrnr_bayesglm, lrnr_gam, lrnr_ranger)

Q1_stack <- make_learner(Stack, learners)
```

Here, we set up a Super Learner of decision trees. We use a new learner
developed for this package in the `sl3` ecosystem, `Lrnr_glmtree`.

``` r
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
```

## Run `CVtreeMLE`

We will now pass the simulated data, learners, and variable names for
each node in *O* = *W*, *A*, *Y* to the `CVtreeMLE` function:

``` r
ptm <- proc.time()

sim_results <- CVtreeMLE(data = sim_data,
                         W = c("W", "W2"),
                         Y = "y",
                         A = c(paste("M", seq(3), sep = "")),
                         back_iter_SL = Q1_stack,
                         tree_SL = discrete_tree_sl, 
                         n_folds = 2,
                         family = "gaussian",
                         H.AW_trunc_lvl = 10,
                         parallel = TRUE,
                         num_cores = 8,
                         max_iter = 5,
                         verbose = TRUE)

proc.time() - ptm 
#>    user  system elapsed 
#> 414.331  17.377 238.249
```

Let’s first look at the RMSE for the iterative back-fitting models.
Because our rules are determined in these models, from which our target
parameter is derived it’s important that our models fit well. Given our
simulated data, we would expect the mixture model to have the lowest
RMSE.

``` r
RMSE_results <- sim_results$`Model RMSEs`
head(RMSE_results) %>%
  kbl(caption = "Model Fit Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Model Fit Results
</caption>
<thead>
<tr>
<th style="text-align:left;">
Var(s)
</th>
<th style="text-align:right;">
RMSE
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
0.0654085
</td>
</tr>
<tr>
<td style="text-align:left;">
M2
</td>
<td style="text-align:right;">
0.0638730
</td>
</tr>
<tr>
<td style="text-align:left;">
M3
</td>
<td style="text-align:right;">
0.0633490
</td>
</tr>
<tr>
<td style="text-align:left;">
M1M2M3
</td>
<td style="text-align:right;">
0.0409719
</td>
</tr>
</tbody>
</table>

Our mixture decision tree model has the lowest RMSE.

We can look at the pooled TMLE results for this model:

``` r
mixture_results <- sim_results$`Pooled TMLE Mixture Results`
head(mixture_results) %>%
  kbl(caption = "Mixture Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Mixture Results
</caption>
<thead>
<tr>
<th style="text-align:right;">
Mixture ATE
</th>
<th style="text-align:right;">
Standard Error
</th>
<th style="text-align:right;">
Lower CI
</th>
<th style="text-align:right;">
Upper CI
</th>
<th style="text-align:right;">
P-value
</th>
<th style="text-align:right;">
P-value Adj
</th>
<th style="text-align:left;">
Vars
</th>
<th style="text-align:right;">
RMSE
</th>
<th style="text-align:left;">
Mixture Interaction Rules
</th>
<th style="text-align:right;">
Fraction Covered
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
5.998892
</td>
<td style="text-align:right;">
0.0142107
</td>
<td style="text-align:right;">
5.971039
</td>
<td style="text-align:right;">
6.026744
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
M1M2M3
</td>
<td style="text-align:right;">
0.2681823
</td>
<td style="text-align:left;">
M1 &gt; 0.021 & M1 &lt; 0.953 & M2 &gt; 2.011 & M2 &lt; 3.994 & M3 &gt;
2.512 & M3 &lt; 4.98
</td>
<td style="text-align:right;">
0.9833333
</td>
</tr>
</tbody>
</table>

\*Note - results in explanations below may change slightly based on
runs:

Above, the mixture ATE for this rule is 2.92 (2.73 - 3.10), which covers
our true ATE used to generate the data which was 3. The mixture ATE is
interpreted as: the average counterfactual mean outcome if all
individuals were exposed to the rule shown in
`Mixture Interaction Rules` compared to if all individuals were
unexposed is 2.92. That is, those individuals who are exposed to this
rule have an outcome that is 2.92 higher compared to those that are not
exposed to this rule. The standard error, confidence intervals and
p-values are derived from the influence curve of this estimator as
described above.

``` r
mixture_v_results <- sim_results$`V-Specific Mix Results`
head(mixture_v_results) %>%
  kbl(caption = "V-Fold Mixture Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class="kable_wrapper lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
V-Fold Mixture Results
</caption>
<tbody>
<tr>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
Mixture ATE
</th>
<th style="text-align:left;">
Standard Error
</th>
<th style="text-align:left;">
Lower CI
</th>
<th style="text-align:left;">
Upper CI
</th>
<th style="text-align:left;">
P-value
</th>
<th style="text-align:left;">
P-value Adj
</th>
<th style="text-align:left;">
RMSE
</th>
<th style="text-align:left;">
Mixture Interaction Rules
</th>
<th style="text-align:left;">
Variables
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
5.76002934560639
</td>
<td style="text-align:left;">
0.0743849818229526
</td>
<td style="text-align:left;">
5.61423746024273
</td>
<td style="text-align:left;">
5.90582123097004
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.151148247583847
</td>
<td style="text-align:left;">
M2 &gt; 1.96611335121485 & M1 &lt;= 0.952865145261606 & M3 &gt;
2.48419675826922
</td>
<td style="text-align:left;">
M1M2M3
</td>
</tr>
<tr>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
5.8622877267532
</td>
<td style="text-align:left;">
0.0282004749837316
</td>
<td style="text-align:left;">
5.80701581143816
</td>
<td style="text-align:left;">
5.91755964206823
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.408007740468045
</td>
<td style="text-align:left;">
M1 &lt;= 0.934177674302425 & M2 &gt; 2.00708027766665 & M3 &gt;
2.47762214185431
</td>
<td style="text-align:left;">
M1M2M3
</td>
</tr>
<tr>
<td style="text-align:left;">
Pooled
</td>
<td style="text-align:left;">
5.8494372758852
</td>
<td style="text-align:left;">
0.0795511930149954
</td>
<td style="text-align:left;">
5.6935
</td>
<td style="text-align:left;">
6.0054
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0.375729111995058
</td>
<td style="text-align:left;">
M1 &gt; 0.021 & M1 &lt; 0.953 & M2 &gt; 2.011 & M2 &lt; 3.994 & M3 &gt;
2.512 & M3 &lt; 4.98
</td>
<td style="text-align:left;">
M1M2M3
</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>

We can plot our v-fold mixture results findings using the
`plot_mixture_results` function. This will return a list of plots with
names corresponding to the interactions found.

``` r
mixture_plots <- plot_mixture_results(v_intxn_results = sim_results$`V-Specific Mix Results`)
mixture_plots$M1M2M3
```

![](README-plot%20sim%20mixture%20results-1.png)<!-- -->

This plot shows the ATE specific for each fold and for the weighted-mean
results over the fold with corresponding pooled variance. The rule is
the pooled rule which includes all observations that were indicated by
the fold specific rules.

``` r
marginal_results <- sim_results$`Pooled TMLE Marginal Results`
head(marginal_results) %>%
  kbl(caption = "Mixture Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Mixture Results
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Marginal ATE
</th>
<th style="text-align:right;">
Standard Error
</th>
<th style="text-align:right;">
Lower CI
</th>
<th style="text-align:right;">
Upper CI
</th>
<th style="text-align:right;">
P-value
</th>
<th style="text-align:right;">
P-value Adj
</th>
<th style="text-align:right;">
RMSE
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
M1\_2-M1\_1
</td>
<td style="text-align:right;">
-1.673242
</td>
<td style="text-align:right;">
0.1995122
</td>
<td style="text-align:right;">
-2.0642784
</td>
<td style="text-align:right;">
-1.282205
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.702335
</td>
</tr>
<tr>
<td style="text-align:left;">
M2\_2-M2\_1
</td>
<td style="text-align:right;">
1.500137
</td>
<td style="text-align:right;">
0.1672006
</td>
<td style="text-align:right;">
1.1724299
</td>
<td style="text-align:right;">
1.827844
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.667953
</td>
</tr>
<tr>
<td style="text-align:left;">
M3\_2-M3\_1
</td>
<td style="text-align:right;">
1.159470
</td>
<td style="text-align:right;">
0.1558593
</td>
<td style="text-align:right;">
0.8539916
</td>
<td style="text-align:right;">
1.464949
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.566421
</td>
</tr>
</tbody>
</table>

This plot shows the data-adaptively identified quantile comparisons for
each variable. Here `M1_2 - M1_1` shows the second quantile for variable
*M*1 minus the first quantile for *M*1. As expected, this difference is
positive and the other two are negative given how we simulated our data.

Similarly we can investigate and plot the v-fold specific estimates:

``` r
marginal_plots <- plot_marginal_results(v_marginal_results =  sim_results$`V-Specific Marg Results`, mix_comps = c(paste("M", seq(3), sep = "")))
marginal_plots$M1
```

![](README-plot%20sim%20marginal%20results-1.png)<!-- --> Same as the
mixtures plot, the marginal plot shows the ATE for an individual
variable in the mixture with corresponding ATE and variance estimates.
The top rule is the pooled rule for the reference category, the rule(s)
in the boxes are the pooled rules for each quantile that was found for
the variable of interest.

## Mixture and Marginal Results

Let’s first look at the mixture results for the model that had the
lowest RMSE:

``` r
mixture_results <- sim_results$`Mixture Results`
head(mixture_results) %>%
  kbl(caption = "Mixture Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Mixture Results
</caption>
<tbody>
<tr>
</tr>
</tbody>
</table>

In this table, Mixture ATE is the counterfactual mean difference if
everyone was exposed to this rule compared to if nobody was exposed to
this rule. Mixture interaction rule is the final rule created that
covers all individuals across the fold specific mixture rules. Coverage
is what proportion of individuals are covered by this rule and is an
indicator of rule stability.

These are the rules found for each individual variable in the vector of
exposures while controlling for other exposures and covariates.

``` r
marginal_results <- sim_results$`Marginal Results`
head(marginal_results) %>%
  kbl(caption = "Marginal Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Marginal Results
</caption>
<tbody>
<tr>
</tr>
</tbody>
</table>

Across the folds, the expected outcome given the cumulative sum of
marginal exposures is also estimated. That is, answering a question such
as “What is the exposure specific mean for each additional exposure
level.”

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/blind-contours/CVtreeMLE/issues).
Further details on filing issues are provided in our [contribution
guidelines](https://github.com/blind-contours/CVtreeMLE/blob/master/CONTRIBUTING.md).

------------------------------------------------------------------------

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/blind-contours/CVtreeMLE/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

------------------------------------------------------------------------

## Citation

After using the `CVtreeMLE` R package, please cite the following:

      @article{mccoyd2022CVtreeMLE-joss,
        author = {McCoy, David B; Hubbard, Alan; Van der Laan Mark},
        title = {{CVtreeMLE}: Efficient Estimation of Mixed Exposures using Data Adaptive Decision Trees and Cross-Validated Targeted Maximum Likelihood Estimation in {R}}
        year  = {2022},
        doi = {TBD},
        url = {TBD},
        journal = {Journal of Open Source Software},
        publisher = {The Open Journal}
      }

      @software{mccoyd2022CVtreeMLE-rpkg,
        author = {McCoy, David B; Hubbard, Alan; Van der Laan Mark},
        title = {{CVtreeMLE}: Efficient Estimation of Mixed Exposures using Data Adaptive Decision Trees and Cross-Validated Targeted Maximum Likelihood Estimation in {R}},
        year  = {2022},
        doi = {TBD},
        url = {https://CRAN.R-project.org/package=CVtreeMLE},
        note = {R package version 0.3.4}
      }

------------------------------------------------------------------------

## Related

-   [R/`sl3`](https://github.com/tlverse/sl3) - An R package providing
    implementation for Super Learner ensemble machine learning
    algorithms.

-   [R/`pre`](https://github.com/marjoleinF/pre) - An R package package
    for deriving prediction rule ensembles for binary, multinomial,
    (multivariate) continuous, count and survival responses.

-   [R/`partykit`](http://partykit.r-forge.r-project.org/partykit/) - A
    toolkit with infrastructure for representing, summarizing, and
    visualizing tree-structured regression and classification models.
    This unified infrastructure can be used for reading/coercing tree
    models from different sources (‘rpart,’ ‘RWeka,’ ‘PMML’) yielding
    objects that share functionality for print()/plot()/predict()
    methods.

-   [R/`SuperLearner`](https://github.com/ecpolley/SuperLearner) -
    Legacy R package providing implementation for Super Learner ensemble
    machine learning algorithms.

------------------------------------------------------------------------

## Funding

The development of this software was supported in part through grants
from the NIH-funded Biomedical Big Data Training Program at UC Berkeley
where I was a biomedical big data fellow.

------------------------------------------------------------------------

## License

© 2017-2021 [David B. McCoy](https://davidmccoy.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    Copyright (c) 2017-2021 David B. McCoy
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

------------------------------------------------------------------------

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-coyle2021sl3" class="csl-entry">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, Rachael V Phillips, and
Oleg Sofrygin. 2021. *<span class="nocase">sl3</span>: Modern Pipelines
for Machine Learning and Super Learning*.
<https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-Fokkema2020a" class="csl-entry">

Fokkema, Marjolein. 2020. “<span class="nocase">Fitting prediction rule
ensembles with R package pre</span>.” *Journal of Statistical Software*
92 (12). <https://doi.org/10.18637/jss.v092.i12>.

</div>

<div id="ref-Hubbard2016" class="csl-entry">

Hubbard, Alan E., Sara Kherad-Pajouh, and Mark J. Van Der Laan. 2016.
“<span class="nocase">Statistical Inference for Data Adaptive Target
Parameters</span>.” *International Journal of Biostatistics* 12 (1):
3–19. <https://doi.org/10.1515/ijb-2015-0013>.

</div>

<div id="ref-Zheng2010" class="csl-entry">

Zheng, Wenjing, and MJ van der Laan. 2010. “<span
class="nocase">Asymptotic theory for cross-validated targeted maximum
likelihood estimation</span>.” *U.C. Berkeley Division of Biostatistics
Working Paper Series*, no. 273.
<http://biostats.bepress.com/ucbbiostat/paper273/>.

</div>

</div>
