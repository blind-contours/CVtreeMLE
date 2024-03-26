
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `CVtreeMLE` <img src="man/figures/CVtreeMLE_sticker.png" style="float:right; height:200px;">

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

This package operationalizes the methodology presented here:

<https://arxiv.org/abs/2302.07976>

People often encounter multiple simultaneous exposures (e.g. several
drugs or pollutants). Policymakers are interested in setting safe
limits, interdictions, or recommended dosage combinations based on a
combination of thresholds, one per exposure. Setting these thresholds is
difficult because all relevant interactions between exposures must be
accounted for. Previous statistical methods have used parametric
estimators which do not directly address the question of safe exposure
limits, rely on unrealistic assumptions, and do not result in a
threshold based statistical quantity that is directly relevant to policy
regulators.

Here we present an estimator that a) identifies thresholds that
minimize/maximize the expected outcome controlling for covariates and
other exposures; and which b) efficiently estimates a policy
intervention which compares the expected outcome if everyone was forced
to these safe levels compared to the observed outcome under observed
exposure distribution.

This is done by using cross-validation where in training folds of the
data, a custom g-computation tree-based search algorithm finds the
minimizing region, and an estimation sample is used to estimate the
policy intervention using targeted maximum likelihood estimation.

## Inputs and Outputs

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
library(kableExtra)
library(ggplot2)
seed <- 98484
set.seed(seed)
```

To illustrate how `CVtreeMLE` may be used to find and estimate a region
that, if intervened on would lead to the biggest reduction in an outcome
we use synthetic data from the National Institute of Environmental
Health:

## National Institute of Environmental Health Data

The 2015 NIEHS Mixtures Workshop was developed to determine if new
mixture methods detect ground-truth interactions built into the
simulated data. In this way we can simultaneously show `CVtreeMLE`
output, interpretation and validity.

For detailed information on this simulated data please see:

<https://github.com/niehs-prime/2015-NIEHS-MIxtures-Workshop>

``` r
niehs_data <- NIEHS_data_1

head(niehs_data) %>%
  kableExtra::kbl(caption = "NIEHS Data") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
NIEHS Data
</caption>
<thead>
<tr>
<th style="text-align:right;">
obs
</th>
<th style="text-align:right;">
Y
</th>
<th style="text-align:right;">
X1
</th>
<th style="text-align:right;">
X2
</th>
<th style="text-align:right;">
X3
</th>
<th style="text-align:right;">
X4
</th>
<th style="text-align:right;">
X5
</th>
<th style="text-align:right;">
X6
</th>
<th style="text-align:right;">
X7
</th>
<th style="text-align:right;">
Z
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
7.534686
</td>
<td style="text-align:right;">
0.4157066
</td>
<td style="text-align:right;">
0.5308077
</td>
<td style="text-align:right;">
0.2223965
</td>
<td style="text-align:right;">
1.1592634
</td>
<td style="text-align:right;">
2.4577556
</td>
<td style="text-align:right;">
0.9438601
</td>
<td style="text-align:right;">
1.8714406
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
19.611934
</td>
<td style="text-align:right;">
0.5293572
</td>
<td style="text-align:right;">
0.9339570
</td>
<td style="text-align:right;">
1.1210595
</td>
<td style="text-align:right;">
1.3350074
</td>
<td style="text-align:right;">
0.3096883
</td>
<td style="text-align:right;">
0.5190970
</td>
<td style="text-align:right;">
0.2418065
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
12.664050
</td>
<td style="text-align:right;">
0.4849759
</td>
<td style="text-align:right;">
0.7210988
</td>
<td style="text-align:right;">
0.4629027
</td>
<td style="text-align:right;">
1.0334138
</td>
<td style="text-align:right;">
0.9492810
</td>
<td style="text-align:right;">
0.3664090
</td>
<td style="text-align:right;">
0.3502445
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
15.600288
</td>
<td style="text-align:right;">
0.8275456
</td>
<td style="text-align:right;">
1.0457137
</td>
<td style="text-align:right;">
0.9699040
</td>
<td style="text-align:right;">
0.9045099
</td>
<td style="text-align:right;">
0.9107914
</td>
<td style="text-align:right;">
0.4299847
</td>
<td style="text-align:right;">
1.0007901
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
18.606498
</td>
<td style="text-align:right;">
0.5190363
</td>
<td style="text-align:right;">
0.7802400
</td>
<td style="text-align:right;">
0.6142188
</td>
<td style="text-align:right;">
0.3729743
</td>
<td style="text-align:right;">
0.5038126
</td>
<td style="text-align:right;">
0.3575472
</td>
<td style="text-align:right;">
0.5906156
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
18.525890
</td>
<td style="text-align:right;">
0.4009491
</td>
<td style="text-align:right;">
0.8639886
</td>
<td style="text-align:right;">
0.5501847
</td>
<td style="text-align:right;">
0.9011016
</td>
<td style="text-align:right;">
1.2907615
</td>
<td style="text-align:right;">
0.7990418
</td>
<td style="text-align:right;">
1.5097039
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

Briefly, this synthetic data can be considered the results of a
prospective cohort epidemiologic study. The outcome cannot cause the
exposures (as might occur in a cross-sectional study). Correlations
between exposure variables can be thought of as caused by common sources
or modes of exposure. The nuisance variable Z can be assumed to be a
potential confounder and not a collider. There are 7 exposures which
have a complicated dependency structure. $X_3$ and $X_6$ do not have an
impact on the outcome.

One issue is that many machine learning algorithms will fail given only
1 variable passed as a feature so let’s add some other covariates.

``` r
niehs_data$Z2 <- rbinom(nrow(niehs_data),
  size = 1,
  prob = 0.3
)

niehs_data$Z3 <- rbinom(nrow(niehs_data),
  size = 1,
  prob = 0.1
)
```

## Run `CVtreeMLE`

``` r
ptm <- proc.time()

niehs_results <- CVtreeMLE(
  data = as.data.frame(niehs_data),
  w = c("Z", "Z2", "Z3"),
  a = c(paste("X", seq(7), sep = "")),
  y = "Y",
  n_folds = 10,
  seed = seed,
  parallel_cv = TRUE,
  parallel = TRUE,
  family = "continuous",
  num_cores = 8,
  min_max = "min",
  min_obs = 25
)
proc.time() - ptm
#>    user  system elapsed 
#>  13.698   0.689 409.696
```

## Mixture Results

First let’s look at the k-fold specific estimates:

``` r
k_fold_results <- niehs_results$`V-Specific Mix Results`

k_fold_results %>%
  kableExtra::kbl(caption = "K-fold Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
K-fold Results
</caption>
<thead>
<tr>
<th style="text-align:right;">
are
</th>
<th style="text-align:right;">
se
</th>
<th style="text-align:right;">
lower_ci
</th>
<th style="text-align:right;">
upper_ci
</th>
<th style="text-align:right;">
p_val
</th>
<th style="text-align:right;">
p_val_adj
</th>
<th style="text-align:right;">
rmse
</th>
<th style="text-align:left;">
mix_rule
</th>
<th style="text-align:right;">
fold
</th>
<th style="text-align:left;">
variables
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
-0.036
</td>
<td style="text-align:right;">
11.477
</td>
<td style="text-align:right;">
-22.530
</td>
<td style="text-align:right;">
22.458
</td>
<td style="text-align:right;">
0.997515
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
2.689
</td>
<td style="text-align:left;">
X2 \<= 0.42
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X2
</td>
</tr>
<tr>
<td style="text-align:right;">
-0.173
</td>
<td style="text-align:right;">
7.816
</td>
<td style="text-align:right;">
-15.492
</td>
<td style="text-align:right;">
15.146
</td>
<td style="text-align:right;">
0.982327
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
2.759
</td>
<td style="text-align:left;">
X2 \<= 0.41
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
X2
</td>
</tr>
<tr>
<td style="text-align:right;">
0.486
</td>
<td style="text-align:right;">
15.208
</td>
<td style="text-align:right;">
-29.322
</td>
<td style="text-align:right;">
30.293
</td>
<td style="text-align:right;">
0.974526
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
2.906
</td>
<td style="text-align:left;">
X2 \<= 0.41
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
X2
</td>
</tr>
<tr>
<td style="text-align:right;">
0.969
</td>
<td style="text-align:right;">
15.328
</td>
<td style="text-align:right;">
-29.074
</td>
<td style="text-align:right;">
31.012
</td>
<td style="text-align:right;">
0.949580
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3.469
</td>
<td style="text-align:left;">
X2 \<= 0.39
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:left;">
X2
</td>
</tr>
<tr>
<td style="text-align:right;">
0.543
</td>
<td style="text-align:right;">
13.659
</td>
<td style="text-align:right;">
-26.229
</td>
<td style="text-align:right;">
27.314
</td>
<td style="text-align:right;">
0.968304
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3.358
</td>
<td style="text-align:left;">
X2 \<= 0.41
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:left;">
X2
</td>
</tr>
<tr>
<td style="text-align:right;">
0.537
</td>
<td style="text-align:right;">
15.434
</td>
<td style="text-align:right;">
-29.713
</td>
<td style="text-align:right;">
30.787
</td>
<td style="text-align:right;">
0.972228
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3.110
</td>
<td style="text-align:left;">
X2 \<= 0.42
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:left;">
X2
</td>
</tr>
<tr>
<td style="text-align:right;">
0.129
</td>
<td style="text-align:right;">
15.465
</td>
<td style="text-align:right;">
-30.182
</td>
<td style="text-align:right;">
30.439
</td>
<td style="text-align:right;">
0.993363
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3.191
</td>
<td style="text-align:left;">
X2 \<= 0.39
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
X2
</td>
</tr>
<tr>
<td style="text-align:right;">
-0.029
</td>
<td style="text-align:right;">
13.987
</td>
<td style="text-align:right;">
-27.443
</td>
<td style="text-align:right;">
27.386
</td>
<td style="text-align:right;">
0.998359
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
2.529
</td>
<td style="text-align:left;">
X2 \<= 0.39
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:left;">
X2
</td>
</tr>
<tr>
<td style="text-align:right;">
0.661
</td>
<td style="text-align:right;">
12.529
</td>
<td style="text-align:right;">
-23.895
</td>
<td style="text-align:right;">
25.217
</td>
<td style="text-align:right;">
0.957948
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
2.647
</td>
<td style="text-align:left;">
X2 \<= 0.39
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
X2
</td>
</tr>
<tr>
<td style="text-align:right;">
1.384
</td>
<td style="text-align:right;">
16.027
</td>
<td style="text-align:right;">
-30.029
</td>
<td style="text-align:right;">
32.797
</td>
<td style="text-align:right;">
0.931196
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3.422
</td>
<td style="text-align:left;">
X2 \<= 0.39
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
X2
</td>
</tr>
</tbody>
</table>

This indicates that the exposure X2 was found in every fold to have the
most minimizing impact on endocrine disruption if all individuals were
were forced to be exposed to levels less around 0.41. This resembles a
policy where, if everyone were still exposed to the other exposures but
we created a regulation that restricted individuals to only exposure of
X2 less than 0.41.

The pooled estimates, leveraging all the folds for our estimates oracle
target parameter looks like:

``` r
pooled_mixture_results <- niehs_results$`Oracle Region Results`

pooled_mixture_results %>%
  kableExtra::kbl(caption = "Oracle Mixture Results") %>%
  kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Oracle Mixture Results
</caption>
<thead>
<tr>
<th style="text-align:right;">
Region ARE
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
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.248
</td>
<td style="text-align:right;">
4.351
</td>
<td style="text-align:right;">
-8.28
</td>
<td style="text-align:right;">
8.777
</td>
<td style="text-align:right;">
0.954462
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

    @article{McCoy2023, 
    doi = {10.21105/joss.04181}, 
    url = {https://doi.org/10.21105/joss.04181}, 
    year = {2023}, publisher = {The Open Journal}, 
    volume = {8}, number = {82}, pages = {4181}, 
    author = {David McCoy and Alan Hubbard and Mark Van der Laan}, 
    title = {CVtreeMLE: Efficient Estimation of Mixed Exposures using Data Adaptive Decision Trees and Cross-Validated Targeted Maximum Likelihood Estimation in R}, 
    journal = {Journal of Open Source Software} }

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
