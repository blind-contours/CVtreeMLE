
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `R/CVtreeMLE` <img src="inst/figures/CVtreeMLE_sticker.png" height="300" align="right"/>

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
mixture, and the ![g](https://latex.codecogs.com/png.latex?g "g") and
![Q](https://latex.codecogs.com/png.latex?Q "Q") estimators needed for
the ATE. These rules and estimators created in training data are applied
to the validation data in order to calculate the final ATE target
parameter. In order to optimize the optimum bias-variance trade-off for
our causal parameter of interest we use cross-validated targeted minimum
loss based estimation (CV-TMLE). `CVtreeMLE` builds off of the CV-TMLE
general theorem of [cross-validated minimum loss based
estimation](https://biostats.bepress.com/cgi/viewcontent.cgi?article=1276&context=ucbbiostat)
(Zheng and Laan 2010) which allows the full utilization of loss based
super learning to obtain the initial estimators needed for our target
parameter without risk of overfitting. Thus, `CVtreeMLE` makes possible
the non-parametric estimation of the causal effects of a mixed exposure
that both results in interpretable results which are useful for public
policy and is asymptotically efficient.

`CVtreeMLE` integrates with the
[`sl3`package](https://github.com/tlverse/sl3) (Coyle et al. 2021) to
allow for ensemble machine learning to be leveraged in the estimation
procedure. `sl3` is used in the iterative backfitting procedure because
this step requires ensemble machine learning with an offset for the
decision tree predictions. In the
![Q](https://latex.codecogs.com/png.latex?Q "Q") and
![g](https://latex.codecogs.com/png.latex?g "g") mechanisms, `CVtreeMLE`
uses the legacy [`Super Learner`
package](https://github.com/tlverse/SuperLearner) (Coyle et al. 2021).
In the iterative backfitting procedure, for decision tree fitting on the
full mixture modeled together, the [`pre`
package](https://github.com/marjoleinF/pre)(Fokkema 2020) is used to fit
rule ensembles. In backfitting procedure to find thresholds in each
mixture component individually, the [`partykit`
package](http://partykit.r-forge.r-project.org/partykit/)\[partykit2015\].
In both instances,trees can be estimated with an offset for ensemble
machine learning predictions.

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
library(sl3)
library(kableExtra)
library(ggplot2)
library(jtools)


set.seed(429153)
```

Use the `simulate_mixture_cube` function to generate simulated data that
represents ground-truth. Here, we create three continuous mixture
variables, ![A](https://latex.codecogs.com/png.latex?A "A"), that are
correlated and baseline covariates,
![W](https://latex.codecogs.com/png.latex?W "W"), that are potential
confounders. Our outcome will be generated such that individuals with a
specific set of exposures have a different outcome compared to
individuals who are not exposed to this combination of exposure levels.

![](inst/figures/The_Cube.png) The above figure illustrates the data we
will generate using this function. Here, individuals exposed to
![M\_1](https://latex.codecogs.com/png.latex?M_1 "M_1") at values less
than 1.0, ![M\_2](https://latex.codecogs.com/png.latex?M_2 "M_2") at
levels more than 2.0, and
![M\_3](https://latex.codecogs.com/png.latex?M_3 "M_3") at levels at or
greater than 2.5 have an outcome of 6, compared to individuals not
exposed to this combination of thresholds who have an expected outcome
of 0 - thus our ATE is 6. Two covariates
![W](https://latex.codecogs.com/png.latex?W "W") confound this
relationship. Let’s simulate this scenario.

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
subspace_assoc_strength_betas <- c(0, 0, 0, 0, 0, 0, 6, 0) # index is the subspace to apply the outcome value to
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

In the current simulation, the threshold for
![M\_1](https://latex.codecogs.com/png.latex?M_1 "M_1") is 0.99,
![M\_2](https://latex.codecogs.com/png.latex?M_2 "M_2") is 2.0, and
![M\_3](https://latex.codecogs.com/png.latex?M_3 "M_3") is 2.5.
Therefore, if we were to specify,
`subspace_assoc_strength_betas <- c(0, 3, 0, 0, 0, 0, 0, 0)`, we would
simulate an outcome that is 3 where
![M\_1](https://latex.codecogs.com/png.latex?M_1 "M_1") is greater than
1 and ![M\_2](https://latex.codecogs.com/png.latex?M_2 "M_2") and
![M\_3](https://latex.codecogs.com/png.latex?M_3 "M_3") are less than
2.0 and 2.5 respectively. Now we can simulate the scenario shown in the
figure with:
`subspace_assoc_strength_betas <- c(0, 0, 0, 0, 0, 0, 6, 0)` with
additional confounding by
![W](https://latex.codecogs.com/png.latex?W "W") and random error.

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
backfitting procedure. These learners will fit
![Y\|W](https://latex.codecogs.com/png.latex?Y%7CW "Y|W") offset by
![Y\|A](https://latex.codecogs.com/png.latex?Y%7CA "Y|A") as we fit
decision trees to the exposure variables both jointly and individially.

``` r
lrnr_glm <- Lrnr_glm$new()
lrnr_bayesglm <- Lrnr_bayesglm$new()
lrnr_gam <- Lrnr_gam$new()
lrnr_lasso <- Lrnr_glmnet$new(alpha = 1)
lrnr_earth <- Lrnr_earth$new()
lrnr_ranger <- Lrnr_ranger$new()
lrnr_xgboost100 <- Lrnr_xgboost$new(nrounds = 100, early_stopping_rounds = 10)
lrnr_xgboost50 <- Lrnr_xgboost$new(nrounds = 50, early_stopping_rounds = 5)
lrnr_xgboost20 <- Lrnr_xgboost$new(nrounds = 20)
# put all the learners together (this is just one way to do it)
learners <- c(lrnr_glm, lrnr_bayesglm,
lrnr_gam, lrnr_ranger,
lrnr_xgboost100, lrnr_xgboost50, lrnr_xgboost20)

Q1_stack <- make_learner(Stack, learners)
```

This second stack of learners will be used in our
![Q](https://latex.codecogs.com/png.latex?Q "Q") and
![g](https://latex.codecogs.com/png.latex?g "g") mechanisms after we
identify rules using the first stack.

``` r
SL.library<- c('SL.randomForest',
               'SL.earth',
               "SL.glm",
               "SL.mean")
```

## Run `CVtreeMLE`

We will now pass the simulated data, learners, and variable names for
each node in
![O = W,A,Y](https://latex.codecogs.com/png.latex?O%20%3D%20W%2CA%2CY "O = W,A,Y")
to the `CVtreeMLE` function:

``` r
ptm <- proc.time()

sim_results <- CVtreeMLE(data = sim_data,
                                   W = c("W", "W2"),
                                   Y = "y",
                                   A = c(paste("M", seq(3), sep = "")),
                                   back_iter_SL = Q1_stack,
                                   SL.library = SL.library,
                                   n_folds = 5,
                                   family = "gaussian",
                                   H.AW_trunc_lvl = 10,
                                   parallel = TRUE,
                                   verbose = FALSE)

proc.time() - ptm
#>     user   system  elapsed 
#> 1157.040  106.988 3230.373
```

## Types of Models

`CVtreeMLE` fits four types of models:

1.  ![Y\|M\_i, W](https://latex.codecogs.com/png.latex?Y%7CM_i%2C%20W "Y|M_i, W")
    or the expected ![Y](https://latex.codecogs.com/png.latex?Y "Y")
    given ![M\_i](https://latex.codecogs.com/png.latex?M_i "M_i") and
    covariates ![W](https://latex.codecogs.com/png.latex?W "W") where
    other mixture components
    ![M\_{\\ne i}](https://latex.codecogs.com/png.latex?M_%7B%5Cne%20i%7D "M_{\ne i}")
    are controlled for - these are marginal rule models. That is,
    decision trees are fit to mixture component
    ![M\_i](https://latex.codecogs.com/png.latex?M_i "M_i") and Super
    Learner is fit to
    ![Y\|M\_{\\ne =i}, W](https://latex.codecogs.com/png.latex?Y%7CM_%7B%5Cne%20%3Di%7D%2C%20W "Y|M_{\ne =i}, W")
    in the iterative back-fitting procedure. In this way, we derive
    individual rules for each mixture compenent while controlling for
    other mixture components and W.

2.  ![Y\|M, W](https://latex.codecogs.com/png.latex?Y%7CM%2C%20W "Y|M, W")
    or the expected Y given
    ![M](https://latex.codecogs.com/png.latex?M "M") and covariates
    ![W](https://latex.codecogs.com/png.latex?W "W") where all mixture
    components are modeled collectively in ensemble partitioning while
    controlling for ![W](https://latex.codecogs.com/png.latex?W "W") -
    these are mixture rule models (multiple mixture variables included
    in a rule compared to 1. above). That is, decision trees are fit to
    the total mixture space and Super Learner is fit to
    ![Y\|W](https://latex.codecogs.com/png.latex?Y%7CW "Y|W") in the
    iterative back-fitting procedure. In this way, we derive rules for
    the total mixture while controlling for
    ![W](https://latex.codecogs.com/png.latex?W "W").

3.  The additive marginal model: this treats exposure as a cumulative
    sum of marginal rules found in the folds. That is, in each fold, the
    marginal fitting in model 1. is conducted, we simply sum up the
    rules found for each mixture component to derive an ordered factor
    variable that describes cumulative exposure. This is
    ![Y\| \\sum\_i^j A\_i, W](https://latex.codecogs.com/png.latex?Y%7C%20%5Csum_i%5Ej%20A_i%2C%20W "Y| \sum_i^j A_i, W")
    or the expected outcome given cumulative exposure while controlling
    for covariates.

4.  The non-additive marginal model: It could in fact be the case that
    there are interactions between the mixture, as represented as a
    vector of indicators, and covariates
    ![W](https://latex.codecogs.com/png.latex?W "W") or within the
    mixture itself. To capture this, we model
    ![Y\|M,W](https://latex.codecogs.com/png.latex?Y%7CM%2CW "Y|M,W")
    where now, each mixture variable
    ![M](https://latex.codecogs.com/png.latex?M "M") is represented as a
    binary indicator of the rule determined within the fold.

Of course, we want to only investigate statistical inference for our
target parameter for models that have the best fit. As such, we want to
review the RMSE for each of the models detailed above.

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
Model
</th>
<th style="text-align:right;">
RMSE
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
M1 &gt; 0.00859 & M1 &lt; 0.93418
</td>
<td style="text-align:right;">
1.9719464
</td>
</tr>
<tr>
<td style="text-align:left;">
M2 &gt; 2.99006 & M2 &lt; 3.99468
</td>
<td style="text-align:right;">
1.8367400
</td>
</tr>
<tr>
<td style="text-align:left;">
M3 &gt; 2.53437 & M3 &lt; 4.98025
</td>
<td style="text-align:right;">
1.8468551
</td>
</tr>
<tr>
<td style="text-align:left;">
M1 &gt; 0.021 & M1 &lt; 0.934 & M2 &gt; 2.034 & M2 &lt; 3.994 & M3 &gt;
2.534 & M3 &lt; 4.98
</td>
<td style="text-align:right;">
0.3791705
</td>
</tr>
<tr>
<td style="text-align:left;">
additive marginal model
</td>
<td style="text-align:right;">
0.9644628
</td>
</tr>
<tr>
<td style="text-align:left;">
non-additive marginal model
</td>
<td style="text-align:right;">
1.2904166
</td>
</tr>
</tbody>
</table>

In the above table the first three rows correspond to model type 1., or
marginal rules with respective RMSE. The fourth row corresponds to model
type 2. or the RMSE of a Super Learner fit with the exposure being a
mixture rule found when fitting decision trees on all the mixture
components simultaneously.

Lines five and six correspond to models 3. and 4. respectively.

As we can see, the model fit with the mixture rule as the *lowest RMSE*,
as we would expect given our simulated outcome was generated based on
this rule. We can also see that the *rule matches what we simulated*.

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
5.99961
</td>
<td style="text-align:right;">
0.0196168
</td>
<td style="text-align:right;">
5.961162
</td>
<td style="text-align:right;">
6.038058
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
0.3791705
</td>
<td style="text-align:left;">
M1 &gt; 0.021 & M1 &lt; 0.934 & M2 &gt; 2.034 & M2 &lt; 3.994 & M3 &gt;
2.534 & M3 &lt; 4.98
</td>
<td style="text-align:right;">
0.9666667
</td>
</tr>
</tbody>
</table>

In this table, Mixture ATE is the counterfactual mean difference if
everyone was exposed to this rule compared to if nobody was exposed to
this rule. Mixture interaction rule is the final rule created that
covers all individuals across the fold specific mixture rules. Coverage
is what proportion of individuals are covered by this rule and is an
indicator of rule stability. As we can see, `CVtreeMLE` identifies the
correct rule in the simulated data and estimates the correct ATE with
proper CI coverage.

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
<th style="text-align:left;">
Marginal Rules
</th>
<th style="text-align:right;">
Fraction Overlap
</th>
<th style="text-align:right;">
Min
</th>
<th style="text-align:right;">
Max
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
M1
</td>
<td style="text-align:right;">
1.331342
</td>
<td style="text-align:right;">
0.2460008
</td>
<td style="text-align:right;">
0.8491892
</td>
<td style="text-align:right;">
1.813495
</td>
<td style="text-align:right;">
1e-07
</td>
<td style="text-align:right;">
2e-07
</td>
<td style="text-align:right;">
1.971946
</td>
<td style="text-align:left;">
M1 &gt; 0.00859 & M1 &lt; 0.93418
</td>
<td style="text-align:right;">
0.5918919
</td>
<td style="text-align:right;">
0.0011497
</td>
<td style="text-align:right;">
2.988107
</td>
</tr>
<tr>
<td style="text-align:left;">
M2
</td>
<td style="text-align:right;">
1.336206
</td>
<td style="text-align:right;">
0.1797070
</td>
<td style="text-align:right;">
0.9839865
</td>
<td style="text-align:right;">
1.688425
</td>
<td style="text-align:right;">
0e+00
</td>
<td style="text-align:right;">
0e+00
</td>
<td style="text-align:right;">
1.836740
</td>
<td style="text-align:left;">
M2 &gt; 2.99006 & M2 &lt; 3.99468
</td>
<td style="text-align:right;">
0.4858300
</td>
<td style="text-align:right;">
0.0002522
</td>
<td style="text-align:right;">
3.994680
</td>
</tr>
<tr>
<td style="text-align:left;">
M3
</td>
<td style="text-align:right;">
1.491790
</td>
<td style="text-align:right;">
0.1734742
</td>
<td style="text-align:right;">
1.1517872
</td>
<td style="text-align:right;">
1.831793
</td>
<td style="text-align:right;">
0e+00
</td>
<td style="text-align:right;">
0e+00
</td>
<td style="text-align:right;">
1.846855
</td>
<td style="text-align:left;">
M3 &gt; 2.53437 & M3 &lt; 4.98025
</td>
<td style="text-align:right;">
0.9748954
</td>
<td style="text-align:right;">
0.0028738
</td>
<td style="text-align:right;">
4.980251
</td>
</tr>
</tbody>
</table>

Across the folds, the expected outcome given the cumulative sum of
marginal exposures is also estimated. That is, answering a question such
as “What is the exposure specific mean for each additional exposure
level.”

``` r
summary(sim_results$`Additive MSM`)
#> 
#> Call:
#> stats::glm(formula = QbarAW_additive_star ~ sum_marg_hits, data = mix_additive_data)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -2.33558  -0.48721   0.04473   0.45386   2.56758  
#> 
#> Coefficients:
#>                Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)    -0.07818    0.09966  -0.785 0.433110    
#> sum_marg_hits1  0.13036    0.11493   1.134 0.257223    
#> sum_marg_hits2  0.39839    0.11523   3.457 0.000592 ***
#> sum_marg_hits3  5.29362    0.14331  36.939  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.6257061)
#> 
#>     Null deviance: 1657.57  on 499  degrees of freedom
#> Residual deviance:  310.35  on 496  degrees of freedom
#> AIC: 1190.5
#> 
#> Number of Fisher Scoring iterations: 2
```

As we can see, the expected ATE when modeled as a cumulative exposure
does not match the simulations and has high RMSE, as we would expect
given the simulated data.

However, if this model did have the lowest RMSE we could investigate
further the cumulative impact of exposure to the marginal rules.

``` r
cumulative_sum_plot <- effect_plot(sim_results$`Additive MSM`, 
            pred = sum_marg_hits, 
            interval = TRUE, 
            y.label = "Expected Outcome",
            x.label = "Cumulative Exposure",
            cat.geom = "line",
            colors = "black")

cumulative_sum_plot
```

![](README-plot%20cumulative%20sum%20effects-1.png)<!-- -->

In the plot above, we see the expected outcome given exposure to none of
the rules found for each individual variable, exposure to any 1 rule for
![M\_1](https://latex.codecogs.com/png.latex?M_1 "M_1"),
![M\_2](https://latex.codecogs.com/png.latex?M_2 "M_2") or
![M\_3](https://latex.codecogs.com/png.latex?M_3 "M_3"), any two or all
three.

Lastly, it could be the case that model 4. has the lowest RMSE and we
want to investigate the respective marginal ATE, or combination of rule
exposures given that model.

The `fit_post_counterfactuals` function takes in the `CVtreeMLE` results
and uses the marginal combination data to calculate the ATE for new
counterfactuals using the fits found across the CV procedure. Below, we
run this to get ATE results if all individuals were exposed to the rules
found for each variable in the simulation compared to if none were
exposed when the marginal rules are modeled in a non-additive fashion.

``` r
post_fit_counterfactuals <- fit_post_counterfactuals(modeling_results = sim_results, 
                         target_mixtures = c("M1", "M2", "M3"), 
                         H.AW_trunc_lvl = 10, 
                         SL.library = SL.library,
                         p_adjust_n = 1)
#> [1] "Fitting SL to marginal rule for mixture"
#> Loading required package: nnls
#> [1] "Fitting SL to marginal rule for mixture"
#> [1] "Fitting SL to marginal rule for mixture"
#> [1] "Fitting SL to marginal rule for mixture"
#> [1] "Fitting SL to marginal rule for mixture"

post_fit_counterfactuals %>%
  kbl(caption = "Post Counterfactual Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Post Counterfactual Results
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
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
P-value adj
</th>
<th style="text-align:right;">
RMSE
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Mixture
</td>
<td style="text-align:right;">
3.619747
</td>
<td style="text-align:right;">
0.1614003
</td>
<td style="text-align:right;">
3.303408
</td>
<td style="text-align:right;">
3.936086
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.087848
</td>
</tr>
</tbody>
</table>

As we can see, the ATE from this model does not match the truth in
simulations and the RMSE is higher compared to the mixture rule fitting.
However, if the RMSE for this model was lowest, one could then
investigate the expected outcome under different combination of marginal
exposures.

## Vignette

For more details as to what’s under the hood in `CVtreeMLE` please see
the included vignette. There, additional applications are shown (on the
NIEHS mixtures workshop data) and results are compared to existing
mixture methods.

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
