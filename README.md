
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `CVtreeMLE` <img src=“man/figures/CVtreeMLE_sticker.png" height=”300” align=“right”/>

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

> Discovery of Critical Thresholds in Mixed Exposures and Estimation of
> Policy Intervention Effects using Targeted Learning

**Author:** [David McCoy](https://davidmccoy.org)

------------------------------------------------------------------------

## What is `CVtreeMLE`?

People often encounter multiple simultaneous exposures (e.g. several
drugs or pollutants). Policymakers are interested in setting safe
limits, interdictions, or recommended dosage combinations based on a
combination of thresholds, one per exposure. Setting these thresholds is
difficult because all relevant interactions between exposures must be
accounted for. Previous statistical methods have used parametric
estimators which do not directly address the question of superadditive
or subadditive effects in a mixture, rely on unrealistic assumptions,
and do not result in a threshold based statistical quantity that is
directly relevant to policy regulators.

Here we present an estimator that a) identifies thresholds that
minimize/maximize the expected outcome marginalized over covariates and
other exposures; and which b) unbiasedly and efficiently estimates a
policy intervention which compares the expected outcome if everyone was
forced to these safe levels compared to the observed outcome under
observed exposure distribution. This is done by combining a custom
g-computation tree-based search algorithm with a targeted maximum
likelihood estimator using cross-validation.

This package takes in a mixed exposure, covariates, outcome, super
learner stacks of learners if determined (if not default are used),
number of folds, minimum observations in a region, if the desired region
is minimizer or maximizer and parallelization parameters.

The output are k-fold specific results for the region found in each fold
with valid inference, a pooled estimate of the overall oracle parameter
across all folds and pooled exposure sets if the region has some
inconsistency across the folds.

------------------------------------------------------------------------

## Installation

*Note:* Because `CVtreeMLE` package (currently) depends on `sl3` that
allows ensemble machine learning to be used for nuisance parameter
estimation and `sl3` is not on CRAN the `CVtreeMLE` package is not
available on CRAN and must be downloaded here.

There are many dependencies for `CVtreeMLE` so it’s easier to break up
installation of the various packages to ensure proper installation.

`CVtreeMLE` uses the `sl3` package to build ensemble machine learners
for each nuisance parameter. We have to install off the development
branch, first download these two packages for `sl3`

``` r
install.packages(c("ranger", "arm", "xgboost", "nnls"))
```

Now install `sl3` on devel:

``` r
remotes::install_github("tlverse/sl3@devel")
```

Make sure `sl3` installs correctly then install `CVtreeMLE`

``` r
remotes::install_github("blind-contours/CVtreeMLE@main")
```

`CVtreeMLE` has some other miscellaneous dependencies that are used in
the examples as well as in the plotting functions.

``` r
install.packages(c("kableExtra", "hrbrthemes", "viridis"))
```

------------------------------------------------------------------------

## Example

First load the package and other packages needed

``` r
library(CVtreeMLE)
library(sl3)
library(pre)
library(partykit)
library(kableExtra)
library(ggplot2)
library(here)
source(here("sandbox", "simulate_2d_data.R"))

set.seed(429153)
```

To illustrate how `CVtreeMLE` may be used to find and esitmate a region
that, if intervened on would lead to the biggest reduction in an outcome
we use simulated data that is described in our paper:

arXiv:2302.07976

Briefly, 2 discrete exposures are created based on a multinomial
regression from baseline covariates. We create an outcome with main
effect and squared interactions. One region has the minimum expected
outcome, when both exposures are equal to 1.

``` r

n <- 400
exposure_grid <- expand.grid(seq(5), seq(5))
labels <- apply(exposure_grid, 1, paste, collapse = " ")
exposure_grid <- cbind.data.frame(exposure_grid, labels)

# beta matrix
c_matrix <- matrix(c(0.1,0.7,0.8,0.9),
                   ncol  = 2,
                   nrow = 2)

# Generate simulated data -----------------
P_0_sim <- gen_covariates(n) %>% # gen covariates
  gen_multinom(exposure_grid = exposure_grid) %>% # assign obs to exposure cube
  assign_outcomes(exposure_grid = exposure_grid,
                  c_matrix = c_matrix)
#> Warning: `rbernoulli()` was deprecated in purrr 1.0.0.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.

intro_data <- P_0_sim$data
```

## Run `CVtreeMLE`

We will now pass the simulated data and variable names for each node in
O = W,A,Y to the `CVtreeMLE` function.

``` r
ptm <- proc.time()

sim_results <- CVtreeMLE(
  data = intro_data,
  w = c("age", "sex", "bmi"),
  a = c("region1", "region2"),
  y = "outcome_obs",
  n_folds = 5,
  parallel_cv = TRUE,
  seed = 2333,
  parallel_type = "multi_session",
  family = "continuous",
  num_cores = 2
)
#> [1] "Depth: 0 Parent mean: Inf Parent N: 320"
#> [1] "Best split found at depth: 0 Split on: region2 at point: 1"
#> [1] "Going deeper from depth: 0 with left rule: region2 <= 1 & "
#> [1] "Going deeper from depth: 0 with right rule: region2 > 1 & "
#> [1] "Depth: 1 Parent mean: 16.2745921683027 Parent N: 65"
#> [1] "Best split found at depth: 1 Split on: region1 at point: 3"
#> [1] "Going deeper from depth: 1 with left rule: region2 <= 1 & region1 <= 3 & "
#> [1] "Going deeper from depth: 1 with right rule: region2 <= 1 & region1 > 3 & "
#> [1] "Depth: 2 Parent mean: 12.6368174660957 Parent N: 28"
#> [1] "No best split found at depth: 2"
#> [1] "Depth: 2 Parent mean: 12.6368174660957 Parent N: 37"
#> [1] "No best split found at depth: 2"
#> [1] "Depth: 1 Parent mean: 16.2745921683027 Parent N: 255"
#> [1] "No best split found at depth: 1"
#> [1] "Depth: 0 Parent mean: Inf Parent N: 320"
#> [1] "Best split found at depth: 0 Split on: region2 at point: 1"
#> [1] "Going deeper from depth: 0 with left rule: region2 <= 1 & "
#> [1] "Going deeper from depth: 0 with right rule: region2 > 1 & "
#> [1] "Depth: 1 Parent mean: 16.6119824135153 Parent N: 61"
#> [1] "No best split found at depth: 1"
#> [1] "Depth: 1 Parent mean: 16.6119824135153 Parent N: 259"
#> [1] "No best split found at depth: 1"
#> [1] "Depth: 0 Parent mean: Inf Parent N: 320"
#> [1] "Best split found at depth: 0 Split on: region2 at point: 1"
#> [1] "Going deeper from depth: 0 with left rule: region2 <= 1 & "
#> [1] "Going deeper from depth: 0 with right rule: region2 > 1 & "
#> [1] "Depth: 1 Parent mean: 16.9582644829873 Parent N: 63"
#> [1] "No best split found at depth: 1"
#> [1] "Depth: 1 Parent mean: 16.9582644829873 Parent N: 257"
#> [1] "No best split found at depth: 1"
#> [1] "Depth: 0 Parent mean: Inf Parent N: 320"
#> [1] "Best split found at depth: 0 Split on: region2 at point: 1"
#> [1] "Going deeper from depth: 0 with left rule: region2 <= 1 & "
#> [1] "Going deeper from depth: 0 with right rule: region2 > 1 & "
#> [1] "Depth: 1 Parent mean: 16.3762713516552 Parent N: 66"
#> [1] "Best split found at depth: 1 Split on: region1 at point: 3"
#> [1] "Going deeper from depth: 1 with left rule: region2 <= 1 & region1 <= 3 & "
#> [1] "Going deeper from depth: 1 with right rule: region2 <= 1 & region1 > 3 & "
#> [1] "Depth: 2 Parent mean: 12.9263815293611 Parent N: 27"
#> [1] "No best split found at depth: 2"
#> [1] "Depth: 2 Parent mean: 12.9263815293611 Parent N: 39"
#> [1] "No best split found at depth: 2"
#> [1] "Depth: 1 Parent mean: 16.3762713516552 Parent N: 254"
#> [1] "No best split found at depth: 1"
#> [1] "Depth: 0 Parent mean: Inf Parent N: 320"
#> [1] "Best split found at depth: 0 Split on: region2 at point: 1"
#> [1] "Going deeper from depth: 0 with left rule: region2 <= 1 & "
#> [1] "Going deeper from depth: 0 with right rule: region2 > 1 & "
#> [1] "Depth: 1 Parent mean: 16.2722177201016 Parent N: 57"
#> [1] "Best split found at depth: 1 Split on: region1 at point: 3"
#> [1] "Going deeper from depth: 1 with left rule: region2 <= 1 & region1 <= 3 & "
#> [1] "Going deeper from depth: 1 with right rule: region2 <= 1 & region1 > 3 & "
#> [1] "Depth: 2 Parent mean: 12.5806656759671 Parent N: 25"
#> [1] "No best split found at depth: 2"
#> [1] "Depth: 2 Parent mean: 12.5806656759671 Parent N: 32"
#> [1] "No best split found at depth: 2"
#> [1] "Depth: 1 Parent mean: 16.2722177201016 Parent N: 263"
#> [1] "No best split found at depth: 1"

proc.time() - ptm
#>    user  system elapsed 
#>   4.785   0.604  82.810
```

We can look at the pooled TMLE results for this model. Let’s see if
`CVtreeMLE` identified the thresholds for the minimizing region:

``` r
mixture_results <- sim_results$`Pooled TMLE Mixture Results`
mixture_results %>%
  kbl(caption = "Oracle Results") %>%
  kable_classic(full_width = FALSE, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Oracle Results
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
Average_Rule
</th>
<th style="text-align:right;">
Proportion_Folds
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
-19.372
</td>
<td style="text-align:right;">
3.804
</td>
<td style="text-align:right;">
-26.827
</td>
<td style="text-align:right;">
-11.916
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1e-06
</td>
<td style="text-align:left;">
region2
</td>
<td style="text-align:right;">
11.370
</td>
<td style="text-align:left;">
region2 \<=1(1,1)
</td>
<td style="text-align:right;">
0.4
</td>
</tr>
<tr>
<td style="text-align:right;">
-1.230
</td>
<td style="text-align:right;">
0.192
</td>
<td style="text-align:right;">
-1.605
</td>
<td style="text-align:right;">
-0.855
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0e+00
</td>
<td style="text-align:left;">
region2-region1
</td>
<td style="text-align:right;">
0.238
</td>
<td style="text-align:left;">
region1 \<=3(3,3) & region2 \<=1(1,1)
</td>
<td style="text-align:right;">
0.6
</td>
</tr>
</tbody>
</table>

Above is the ARE for the region, region 1 = 1 and region 2 = 1. This was
found in all the folds so there is no variability in our oracle
estimate. The ARE is -0.2 which indicates that if we were to make a
policy that forced everyone into exposure levels of 1 for both exposures
the expected outcome would be 0.2 less than the current average.

We can also look at the v-fold specific results:

``` r
sim_results$`V-Specific Mix Results`
#>       are    se lower_ci upper_ci    p_val p_val_adj   rmse
#> 1  -0.920 0.313   -1.535   -0.306 0.003322  0.016611  0.236
#> 2 -20.565 5.532  -31.407   -9.723 0.000201  0.001006 11.648
#> 3 -18.192 5.193  -28.371   -8.013 0.000460  0.002302 11.067
#> 4  -0.388 0.327   -1.029    0.252 0.234842  1.000000  0.521
#> 5  -2.096 0.345   -2.772   -1.419 0.000000  0.000000  0.203
#>                      mix_rule fold       variables
#> 1 region2 <= 1 & region1 <= 3    1 region2-region1
#> 2                region2 <= 1    2         region2
#> 3                region2 <= 1    3         region2
#> 4 region2 <= 1 & region1 <= 3    4 region2-region1
#> 5 region2 <= 1 & region1 <= 3    5 region2-region1
```

Here we see the fold specific estimates which are used in the pooled
TMLE output we saw earlier.

Additional details for this and other features are given in the
vignette.

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/blind-contours/CVtreeMLE/issues).
Further details on filing issues are provided in our [contribution
guidelines](https://github.com/blind-contours/%20CVtreeMLE/blob/main/CONTRIBUTING.md).

------------------------------------------------------------------------

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/blind-contours/%20CVtreeMLE/blob/main/CONTRIBUTING.md)
prior to submitting a pull request.

------------------------------------------------------------------------

## Citation

After using the `CVtreeMLE` R package, please cite the following:

(**article?**){McCoy2023, doi = {10.21105/joss.04181}, url =
{<https://doi.org/10.21105/joss.04181>}, year = {2023}, publisher = {The
Open Journal}, volume = {8}, number = {82}, pages = {4181}, author =
{David McCoy and Alan Hubbard and Mark Van der Laan}, title =
{CVtreeMLE: Efficient Estimation of Mixed Exposures using Data Adaptive
Decision Trees and Cross-Validated Targeted Maximum Likelihood
Estimation in R}, journal = {Journal of Open Source Software} }

------------------------------------------------------------------------

## Related

- [R/`sl3`](https://github.com/tlverse/sl3) - An R package providing
  implementation for Super Learner ensemble machine learning algorithms.

- [R/`SuperLearner`](https://github.com/ecpolley/SuperLearner) - Legacy
  R package providing implementation for Super Learner ensemble machine
  learning algorithms.

------------------------------------------------------------------------

## Funding

The development of this software was supported in part through grants
from the NIH-funded Biomedical Big Data Training Program at UC Berkeley
where I was a biomedical big data fellow.

------------------------------------------------------------------------

## License

© 2017-2024 David B. McCoy

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    Copyright (c) 2017-2024 David B. McCoy
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
