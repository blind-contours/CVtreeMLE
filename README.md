
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
for each nuisance parameter.

Install `sl3` on devel:

``` r
remotes::install_github("tlverse/sl3@devel")
```

Make sure `sl3` installs correctly then install `CVtreeMLE`

``` r
remotes::install_github("blind-contours/CVtreeMLE@main")
```

------------------------------------------------------------------------

## Example

First load the package and other packages needed

``` r
library(CVtreeMLE)
library(sl3)
library(dplyr)
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

# Convert continuous X variables to their corresponding deciles for example
niehs_data <- niehs_data %>%
  mutate(across(starts_with("X"), ~ ntile(., 10), .names = "decile_{col}"))

niehs_results <- CVtreeMLE(
  data = as.data.frame(niehs_data),
  w = c("Z", "Z2", "Z3"),
  a = c("decile_X1", "decile_X2", "decile_X3", "decile_X4", "decile_X5", "decile_X6", "decile_X7"),
  y = "Y",
  n_folds = 7,
  seed = seed,
  parallel_cv = TRUE,
  parallel = TRUE,
  family = "continuous",
  num_cores = 8,
  min_max = "min",
  max_depth = 2, 
  min_obs = 25
)
proc.time() - ptm
#>    user  system elapsed 
#>   7.184   0.816  91.497
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
-2.876
</td>
<td style="text-align:right;">
20.006
</td>
<td style="text-align:right;">
-42.087
</td>
<td style="text-align:right;">
36.335
</td>
<td style="text-align:right;">
0.885677
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
4.118
</td>
<td style="text-align:left;">
decile_X1 \<= 5 & decile_X7 \<= 3
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
decile_X1-decile_X7
</td>
</tr>
<tr>
<td style="text-align:right;">
-2.733
</td>
<td style="text-align:right;">
20.039
</td>
<td style="text-align:right;">
-42.009
</td>
<td style="text-align:right;">
36.544
</td>
<td style="text-align:right;">
0.891537
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
4.123
</td>
<td style="text-align:left;">
decile_X1 \<= 5 & decile_X7 \<= 4
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
decile_X1-decile_X7
</td>
</tr>
<tr>
<td style="text-align:right;">
-3.061
</td>
<td style="text-align:right;">
24.969
</td>
<td style="text-align:right;">
-51.999
</td>
<td style="text-align:right;">
45.878
</td>
<td style="text-align:right;">
0.902437
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5.026
</td>
<td style="text-align:left;">
decile_X1 \<= 5 & decile_X7 \<= 3
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
decile_X1-decile_X7
</td>
</tr>
<tr>
<td style="text-align:right;">
-2.196
</td>
<td style="text-align:right;">
18.755
</td>
<td style="text-align:right;">
-38.955
</td>
<td style="text-align:right;">
34.564
</td>
<td style="text-align:right;">
0.906806
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
4.452
</td>
<td style="text-align:left;">
decile_X1 \<= 5 & decile_X7 \<= 3
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:left;">
decile_X1-decile_X7
</td>
</tr>
<tr>
<td style="text-align:right;">
-3.257
</td>
<td style="text-align:right;">
21.951
</td>
<td style="text-align:right;">
-46.280
</td>
<td style="text-align:right;">
39.765
</td>
<td style="text-align:right;">
0.882028
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
4.489
</td>
<td style="text-align:left;">
decile_X1 \<= 5 & decile_X7 \<= 4
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:left;">
decile_X1-decile_X7
</td>
</tr>
<tr>
<td style="text-align:right;">
-3.371
</td>
<td style="text-align:right;">
16.739
</td>
<td style="text-align:right;">
-36.178
</td>
<td style="text-align:right;">
29.436
</td>
<td style="text-align:right;">
0.840393
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
4.000
</td>
<td style="text-align:left;">
decile_X1 \<= 5 & decile_X7 \<= 3
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:left;">
decile_X1-decile_X7
</td>
</tr>
<tr>
<td style="text-align:right;">
-3.441
</td>
<td style="text-align:right;">
21.135
</td>
<td style="text-align:right;">
-44.864
</td>
<td style="text-align:right;">
37.982
</td>
<td style="text-align:right;">
0.870658
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
4.491
</td>
<td style="text-align:left;">
decile_X1 \<= 5 & decile_X7 \<= 4
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
decile_X1-decile_X7
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
-3.12
</td>
<td style="text-align:right;">
7.776
</td>
<td style="text-align:right;">
-18.361
</td>
<td style="text-align:right;">
12.121
</td>
<td style="text-align:right;">
0.688247
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
guidelines](https://github.com/blind-contours/CVtreeMLE/blob/main/contributing.md)
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
