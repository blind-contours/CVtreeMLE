
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `CVtreeMLE` \<img src=“man/figures/CVtreeMLE_sticker.png height=”300” align=“right”/\>

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
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070042.svg)](https://doi.org/10.5281/zenodo.4070042) -->
<!-- [![DOI](https://joss.theoj.org/papers/10.21105/joss.02447/status.svg)](https://doi.org/10.21105/joss.02447) -->
[![Codecov test
coverage](https://codecov.io/gh/blind-contours/CVtreeMLE/branch/main/graph/badge.svg)](https://app.codecov.io/gh/blind-contours/CVtreeMLE?branch=main)
<!-- badges: end -->

> Efficient Estimation of the Causal Effects of Joint Exposure using
> Data Adaptive Decision Trees and Cross-Validated Targeted Maximum
> Likelihood Estimation **Author:** David McCoy

------------------------------------------------------------------------

## What is `CVtreeMLE`?

To-date, there is no method which delivers robust statistical inference
based on leaves of a decision tree that also uses the full data.
Although decision trees are widely used algorithms because they provide
interpretable rule-based results, existing methods do not flexibly
adjust for covariates while fitting the tree to a specific set of
exposures. This is of particular interest for environmental
epidemiologists who may be interested in what parts of a pollution
exposure cause adverse health effects or to the medical analyst
interested in finding combinations of drugs which improve patient
outcomes. Current approaches use sample splitting techniques where the
rules found on a training set are applied to a test set to estimate the
conditional mean in each of tree leaves. This type of sample splitting
loses power. Furthermore, existing methods require the decision tree
algorithm to be chosen by the analyst which may not fit the data best,
introducing bias.

A gap exists in statistical estimation where the analyst wants to
understand what exposure levels, as part of a mixture, have the greatest
effect on some outcome while nonparametrically adjusting for baseline
covariates and data-adaptively identifying the best fitting tree applied
to the mixed exposure. This is the goal of `CVtreeMLE`.

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
data-adaptive parameters see (Hubbard, Kherad-Pajouh, and Van Der Laan
2016).

`CVtreeMLE` uses data-adaptive parameters by implementing V-fold
cross-validation (CV), that is, in 10-fold CV, the data is split 10
times (folds), where 90% of the data is used to determine rules in a
mixture, and the
![g](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g "g")
(probability of being exposed estimator, P(A\|W)) and
![Q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Q "Q")
(outcome estimator E(Y\|A,W)) estimators needed for the ATE. These rules
and estimators created in training data are applied to the validation
data in order to calculate the final ATE target parameter. This process
cycles 10 times and the parameters across the folds are pooled, utlizing
the full data. Because the partition nodes used in the tree may vary
across the folds, we provide a union rule and stability estimates which
informs how variable the rule is across the folds.

In order to optimize the optimum bias-variance trade-off for our causal
parameter of interest we use cross-validated targeted minimum loss based
estimation (CV-TMLE). `CVtreeMLE` builds off of the CV-TMLE general
theorem of cross-validated minimum loss based estimation Zheng and Laan
(2010) which allows the full utilization of loss based super learning to
obtain the initial estimators needed for our target parameter without
risk of overfitting. Thus, `CVtreeMLE` makes possible the nonparametric
estimation of the causal effects of a mixed exposure that both results
in interpretable results which are useful for public policy and are
asymptotically efficient.

`CVtreeMLE` integrates with the
[`sl3`package](https://github.com/tlverse/sl3) (Coyle et al. 2021) to
allow for ensemble machine learning to be leveraged in the estimation
procedure. `sl3` is used to create ensemble machine learning estimators
for the
![Q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Q "Q")
and
![g](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g "g")
mechanisms for the average treatment effect (ATE) target parameter.
`sl3` is also used in the iterative backfitting procedure. The iterative
backfitting algorithm fits E(Y\|A,W) = f(A) + h(W), where f(A) is a
Super Learner of decision trees and H(W) is an unrestricted Super
Learner. This semi-parametric additive model allows us to isolate
identification of nodes in A that best explain an outcome while flexibly
adjusting for covariates. In the iterative backfitting procedure, an
ensemble of decision trees are fit on the full mixture modeled together;
the [`pre` package](https://github.com/marjoleinF/pre)(Fokkema 2020) is
used to fit rule ensembles. In backfitting procedure to find thresholds
in each mixture component individually, a Super Learner of decision
trees generated from the [`partykit`
package](http://partykit.r-forge.r-project.org/partykit/)\[partykit2015\]
are created. In each case, the goal is to find the best fitting decision
tree from which we extract decision tree rules, we then calculate the
ATE for these rules.

Below we simulate data for a “mixture cube” and show basic results from
`CVtreeMLE` compared to ground-truth. For more details, please see the
vignette.

------------------------------------------------------------------------

## Installation

*Note:* Because `CVtreeMLE` package (currently) depends on `sl3`, an
enhancing dependency that allows ensemble machine learning to be used
for nuisance parameter estimation and `sl3` is not on CRAN the
`CVtreeMLE` package is not available on CRAN and must be downloaded
here.

For the latest features, install the most recent *stable version*  
of `CVtreeMLE`from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("blind-contours/CVtreeMLE@main")
```

`CVtreeMLE` uses newly added decision tree estimators added to `sl3` so
please download sl3 from:

``` r
remotes::install_github("tlverse/sl3@devel")
```

------------------------------------------------------------------------

## Example

First load the package and other packages needed

``` r
library(CVtreeMLE)
library(sl3)
library(kableExtra)
library(ggplot2)

set.seed(429153)
```

To illustrate how `CVtreeMLE` may be used to ascertain the effect of a
joint exposure, consider three exposures. We can represent these
exposures as a cube, this is depicted below:

![](man/figures/The_Cube.png)

This figure shows two regions we may be interested in within the mixture
space. On the left, the orange region is where
![M_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M_1 "M_1")
is less than 1,
![M_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M_2 "M_2")
is greater than 2 and
![M_3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M_3 "M_3")
is less than 2.5. If we were to generate this data, we’d expect
`CVtreeMLE` to find these thresholds in the mixture space. Likewise, in
blue this shows the region where
![M_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M_1 "M_1")
is less than 1,
![M_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M_2 "M_2")
is greater than 2 and
![M_3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;M_3 "M_3")
is greater than 2.5. We want to simulate an outcome that is highest in
this region and where the density of the region is based on covariates.

To do this, use the `simulate_mixture_cube` function to generate
simulated  
data that represents ground-truth. Here, we create three continuous
mixture variables,
![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A"),
that are correlated and baseline covariates,
![W](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;W "W"),
that are confounders. Our outcome will be generated such that
individuals with a specific set of exposures have a different outcome
compared to individuals who are not exposed to this combination of
exposure levels.

## Simulate Data

``` r
sim_data <- simulate_mixture_cube(n_obs = 800,
                                  splits = c(0.99, 2.0, 2.5),
                                  mins = c(0, 0, 0),
                                  maxs = c(3, 4, 5),
                                  subspace_assoc_strength_betas = c(0, 0, 0, 0,
                                                                    0, 0, 6, 0))
head(sim_data) %>%
  kbl(caption = "Simulated Data") %>%
  kable_classic(full_width = FALSE, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Simulated Data
</caption>
<thead>
<tr>
<th style="text-align:right;">
w
</th>
<th style="text-align:right;">
w2
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
0.0092552
</td>
<td style="text-align:right;">
-0.2031663
</td>
<td style="text-align:right;">
0.3790036
</td>
<td style="text-align:right;">
0.0344271
</td>
<td style="text-align:right;">
0.2936966
</td>
<td style="text-align:right;">
-0.1806314
</td>
</tr>
<tr>
<td style="text-align:right;">
0.0981536
</td>
<td style="text-align:right;">
0.2613579
</td>
<td style="text-align:right;">
1.3882993
</td>
<td style="text-align:right;">
2.3493205
</td>
<td style="text-align:right;">
3.8962661
</td>
<td style="text-align:right;">
0.3439479
</td>
</tr>
<tr>
<td style="text-align:right;">
-0.0950558
</td>
<td style="text-align:right;">
0.2606486
</td>
<td style="text-align:right;">
1.9012774
</td>
<td style="text-align:right;">
2.0433104
</td>
<td style="text-align:right;">
0.1834768
</td>
<td style="text-align:right;">
0.1821716
</td>
</tr>
<tr>
<td style="text-align:right;">
-0.0991683
</td>
<td style="text-align:right;">
-0.5569234
</td>
<td style="text-align:right;">
1.2717409
</td>
<td style="text-align:right;">
2.7877759
</td>
<td style="text-align:right;">
0.6990761
</td>
<td style="text-align:right;">
-0.6474042
</td>
</tr>
<tr>
<td style="text-align:right;">
0.1347455
</td>
<td style="text-align:right;">
0.3248600
</td>
<td style="text-align:right;">
1.6475103
</td>
<td style="text-align:right;">
3.3305123
</td>
<td style="text-align:right;">
4.1460804
</td>
<td style="text-align:right;">
0.4376962
</td>
</tr>
<tr>
<td style="text-align:right;">
-0.6767346
</td>
<td style="text-align:right;">
-0.5346316
</td>
<td style="text-align:right;">
0.2067554
</td>
<td style="text-align:right;">
2.0040681
</td>
<td style="text-align:right;">
0.3495084
</td>
<td style="text-align:right;">
-1.2134008
</td>
</tr>
</tbody>
</table>

Using the `simulate_mixture_cube` we generate 500 observations that are
exposed to three variables with min values being 0 for all and max
values being 3,4, and 5. In each variable we define split points 0.99,
2.0, and 2.5. Given the eight regions in the cube, the
`subspace_assoc_strength_betas` parameter is where we put the outcome in
a specific region.

The indices correspond to an area in the cube:

1.  All mixtures lower than specified thresholds
2.  M1 is higher but M2 and M3 are lower
3.  M2 is higher but M1 and M3 are lower
4.  M1 and M2 are higher and M3 is lower
5.  M3 is higher and M1 and M2 are lower
6.  M1 and M3 are higher and M2 is lower
7.  M2 and M3 are higher and M1 is lower
8.  All mixtures are higher than thresholds

So here - we put 6 at index 7 which means the outcome is 6 when M2 and
M3 are higher and M1 is lower than their respective split points. The
outcome is 0 in all other regions.

## Run `CVtreeMLE`

We will now pass the simulated data and variable names for each node in
![O = W,A,Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;O%20%3D%20W%2CA%2CY "O = W,A,Y")
to the `CVtreeMLE` function.

``` r
ptm <- proc.time()

sim_results <- CVtreeMLE(data = sim_data,
                         w = c("w", "w2"),
                         a = c(paste("M", seq(3), sep = "")),
                         y = "y",
                         n_folds = 5,
                         parallel_cv = TRUE,
                         seed = 2333,
                         parallel_type = "multi_session",
                         family = "gaussian",
                         num_cores = 2)

proc.time() - ptm
#>    user  system elapsed 
#>   83.73    4.45 1975.31
```

Note that above, there are default estimators for all parameters if they
are not passed to the function. Here we just use the out of the box
estimators that are defined in `utils_create_sls.R`. These estimators
are chosen to be both non-parametric but also not too computationally
demanding. Examples of estimators used by default are random forest,
xgboost, elastic net, and glms. Users can also pass in their own custom
stacks of learners. We also see here that, using 2 cores with these
learners on our simulated data with 500 observations and 6 variables,
our run time is 21 minutes. This can be greatly improved by increasing
the num_cores parameter.

## Results

We can look at the pooled TMLE results for this model. Let’s see if
`CVtreeMLE` identified the current rule in all our folds:

``` r
mixture_results <- sim_results$`Pooled TMLE Mixture Results`
mixture_results %>%
  dplyr::filter(Proportion_Folds == 1.0) %>%
  kbl(caption = "Mixture Results") %>%
  kable_classic(full_width = FALSE, html_font = "Cambria")
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
Union_Rule
</th>
<th style="text-align:right;">
Proportion_Folds
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
2.945
</td>
<td style="text-align:right;">
0.217
</td>
<td style="text-align:right;">
2.519
</td>
<td style="text-align:right;">
3.371
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
M1M2
</td>
<td style="text-align:right;">
1.552
</td>
<td style="text-align:left;">
M1 \> 0.002 & M1 \< 1.182 & M2 \> 1.946 & M2 \< 3.975
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:right;">
5.962
</td>
<td style="text-align:right;">
0.105
</td>
<td style="text-align:right;">
5.756
</td>
<td style="text-align:right;">
6.167
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
0.366
</td>
<td style="text-align:left;">
M1 \> 0.009 & M1 \< 0.989 & M2 \> 1.946 & M2 \< 3.975 & M3 \> 2.449 & M3
\< 4.99
</td>
<td style="text-align:right;">
1
</td>
</tr>
</tbody>
</table>

Above, the estimated mixture ATE for this rule is 5.96 (5.6 - 6.4),
which covers our true mixture ATE used to generate the data which was 6.
The estimated mixture ATE is interpreted as: the average counterfactual
mean outcome if all individuals were exposed to the rule shown in
`Union Rule` compared to if all individuals were unexposed is 5.96. That
is, those individuals who are exposed to this rule have an outcome that
is 5.96 higher compared to those that are not exposed to this rule. The
standard error, confidence intervals and p-values are derived from the
influence curve of this estimator.

We can also look at the v-fold specific results. This gives the analyst
the ability to investigate how stable the estimates and rules are. These
results are the same as standard sample splitting techniques and
therefore have proper variance estimates and p-values. Below we show the
v-fold specific interactions found with fold specific estimates of our
ATE target parameter and variance estimates from the fold specific IC.

``` r
mixture_v_results <- sim_results$`V-Specific Mix Results`
mixture_v_results$M1M2M3
#>     ate    se lower_ci upper_ci p_val p_val_adj  rmse
#> 1 5.006 0.424   4.1750   5.8380     0         0 0.834
#> 2 5.410 0.060   5.2920   5.5280     0         0 0.684
#> 3 5.704 0.072   5.5630   5.8440     0         0 0.156
#> 4 5.461 0.293   4.8860   6.0360     0         0 0.644
#> 5 5.552 0.032   5.4890   5.6150     0         0 0.444
#> 6 5.542 0.525   4.5137   6.5709     0         0 0.455
#>                                                                     mix_rule
#> 1                                      M1 <= 1.004 & M2 > 1.936 & M3 > 2.449
#> 2                                      M3 > 2.508 & M2 > 1.994 & M1 <= 0.989
#> 3                                      M1 <= 1.019 & M3 > 2.502 & M2 > 1.966
#> 4                                      M1 <= 0.995 & M2 > 2.006 & M3 > 2.341
#> 5                                      M1 <= 1.019 & M3 > 2.501 & M2 > 1.966
#> 6 M1 > 0.009 & M1 < 0.989 & M2 > 1.946 & M2 < 3.975 & M3 > 2.449 & M3 < 4.99
#>     fold variables
#> 1      1    M1M2M3
#> 2      2    M1M2M3
#> 3      3    M1M2M3
#> 4      4    M1M2M3
#> 5      5    M1M2M3
#> 6 Pooled    M1M2M3
```

In v-fold specific results we also give a pooled estimate. This is
different than the pooled TMLE estimate. Here we simply take the
weighted average of the fold specific ATEs and the harmonic mean of the
variances. This is similar to meta-analysis approaches.

We can plot our v-fold mixture results findings using the
`plot_mixture_results` function. This will return a list of plots with
names corresponding to the interactions found.

``` r
mixture_plots <- plot_mixture_results(v_intxn_results =
                                        sim_results$`V-Specific Mix Results`,
                                      hjust = 1.05)
mixture_plots$M1M2M3
```

![](man/figures/README-plot%20sim_mixture_results-1.png)<!-- --> This
plot shows the ATE specific for each fold and for the weighted-mean
results over the fold with corresponding pooled variance. The rule is
the union rule which includes all observations that were indicated by
the fold specific rules.

`CVtreeMLE` also data-adaptively identifies thresholds in the marginal
space. These around found here:

``` r
marginal_results <- sim_results$`Pooled TMLE Marginal Results`
marginal_results %>%
  dplyr::filter(`Proportion in Fold` == 1.0) %>%
  kbl(caption = "Marginal Results") %>%
  kable_classic(full_width = FALSE, html_font = "Cambria")
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
rmse
</th>
<th style="text-align:left;">
Variable Quantile
</th>
<th style="text-align:left;">
Marginal Rules
</th>
<th style="text-align:left;">
Min
</th>
<th style="text-align:left;">
Max
</th>
<th style="text-align:left;">
Proportion in Fold
</th>
<th style="text-align:left;">
var
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
M1_2-M1_1
</td>
<td style="text-align:right;">
-1.543
</td>
<td style="text-align:right;">
0.054
</td>
<td style="text-align:right;">
-1.648
</td>
<td style="text-align:right;">
-1.437
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.154
</td>
<td style="text-align:left;">
M1_2
</td>
<td style="text-align:left;">
M1 \> 0.9862 & M1 \< 2.99837
</td>
<td style="text-align:left;">
0.002
</td>
<td style="text-align:left;">
2.998
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
M1
</td>
</tr>
<tr>
<td style="text-align:left;">
M2_2-M2_1
</td>
<td style="text-align:right;">
1.447
</td>
<td style="text-align:right;">
0.101
</td>
<td style="text-align:right;">
1.250
</td>
<td style="text-align:right;">
1.644
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2.311
</td>
<td style="text-align:left;">
M2_2
</td>
<td style="text-align:left;">
M2 \> 1.96401 & M2 \< 3.97516
</td>
<td style="text-align:left;">
0.007
</td>
<td style="text-align:left;">
3.975
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
M2
</td>
</tr>
<tr>
<td style="text-align:left;">
M3_2-M3_1
</td>
<td style="text-align:right;">
1.528
</td>
<td style="text-align:right;">
0.128
</td>
<td style="text-align:right;">
1.277
</td>
<td style="text-align:right;">
1.780
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2.291
</td>
<td style="text-align:left;">
M3_2
</td>
<td style="text-align:left;">
M3 \> 2.5078 & M3 \< 4.99038
</td>
<td style="text-align:left;">
0.001
</td>
<td style="text-align:left;">
4.99
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
M3
</td>
</tr>
</tbody>
</table>

Additional details for this and other features are given in the
vignette.

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/blind-contours/CVtreeMLE/issues).
Further details on filing issues are provided in our [contribution
guidelines](https://github.com/blind-contours/%20CVtreeMLE/blob/master/CONTRIBUTING.md).

------------------------------------------------------------------------

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/blind-contours/%20CVtreeMLE/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

------------------------------------------------------------------------

## Citation

After using the `CVtreeMLE` R package, please cite the following:

(**article?**){mccoyd2022CVtreeMLE-joss, author = {McCoy, David B;
Hubbard, Alan; Van der Laan Mark}, title = {{CVtreeMLE}: Efficient
Estimation of Mixed Exposures using Data Adaptive Decision Trees and
Cross-Validated Targeted Maximum Likelihood Estimation in {R}} year =
{2022}, doi = {TBD}, url = {TBD}, journal = {Journal of Open Source
Software}, publisher = {The Open Journal} }

(**software?**){mccoyd2022CVtreeMLE-rpkg, author = {McCoy, David B;
Hubbard, Alan; Van der Laan Mark}, title = {{CVtreeMLE}: Efficient
Estimation of Mixed Exposures using Data Adaptive Decision Trees and
Cross-Validated Targeted Maximum Likelihood Estimation in {R}}, year =
{2022}, doi = {TBD}, url =
{<https://CRAN.R-project.org/package=CVtreeMLE>}, note = {R package
version 0.3.4} }

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
    models from different sources (‘rpart’, ‘RWeka’, ‘PMML’) yielding
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

© 2017-2021 David B. McCoy

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    Copyright (c) 2017-2023 David B. McCoy
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
