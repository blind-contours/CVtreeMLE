
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`CVtreeMLE`

<!-- badges: start -->

[![R-CMD-check](https://github.com/blind-contours/CVtreeMLE/workflows/R-CMD-check/badge.svg)](https://github.com/blind-contours/CVtreeMLE/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/blind-contours/CVtreeMLE/master.svg)](https://codecov.io/github/blind-contours/CVtreeMLE?branch=master)
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
Likelihood Estmation) R package is designed to provide facilities for
the construction of efficient estimators of an average treatment effect
(ATE) causal parameter defined as the counterfactual mean of an outcome
if all individuals were jointly exposed to a combination of exposure
levels in a mixed exposure compared to if all individuals were not
exposed. Here, a joint exposure is data-adaptively defined based on
decision trees applied to a set of exposure variables while flexibly
controlling for covariates non-parametrically. For more information on
data- adaptive parameters see (Dı́az and van der Laan 2012). In each
V-fold during cross- validation, `CVtreeMLE`: 1. Performs an iterative
backfitting of a Super Learner for Y\|W, the expected outcome given
covariates, and Y\|A, the expected outcome given the mixture variables.
Consider the estimator for Y\|W as h(x) and the estimator for Y\|A as
g(x). After initialization, h(x) is fit with an offset for predictions
from g(x). Likewise, g(x) is fit offset with predictions from h(x). This
procedure iterates until convergence, or when there is no change in fit
between either estimator.

2.  The same iterative backfitting is conducted but for Y\|M\_i and
    Y\|W,M\_ne\_i - or rather, iteratively backfitting decision trees to
    each mixture variable individually while controlling for covariates
    and other mixture variables in the complementary non-parametric
    Super Learner.

3.  If decision trees are found for the mixture measured as an
    interaction and/or any marginal mixture component which explain the
    outcome while flexibly controlling for covariates, the binary
    indicators are created which indicate observations that meet the
    respective rule(s) in the tree.

4.  For each decision tree, the propensity score is estimated (the
    probability of being exposed to the respective levels of the
    exposure(s)), or the g mechanism. Similarly, the Q mechanism, or the
    outcome estimator is created which estimates the expected outcome
    given exposure to the rule and covariates.

5.  Once these nuisance parameters are estimated, a pooled targeted
    maximum likelihood estimation is done to estimate the average
    treatment effect for the data-adaptively identified exposures
    (thresholds). This average treatment effect is our target parameter
    and is pooled over the validation data in each V-fold.

6.  Tables are given which include the ATE for each variable or sets of
    variables found across all the folds alongside variance estimates
    from the efficient influence curve. Average rules are created which
    are made from observations that meet each fold-specific rule.

`CVtreeMLE` integrates with the [`sl3`
package](https://github.com/tlverse/sl3) (**coyle2020sl3?**) to allow
for ensemble machine learning to be leveraged in the estimation
procedure. `sl3` is used in the iterative backfitting procedure because
this step requires ensemble machine learning with an offset for the
decision tree predictions. In the Q and g mechanisms, `CVtreeMLE` uses
the legacy [`Super Learner`
package](https://github.com/tlverse/SuperLearner) (**coyle2020sl3?**).
In the iterative backfitting procedure, for decision tree fitting on the
full mixture modeled together, the \[`pre` package(PRE PACKAGE
HERE)[citation](#citation) is used to fit rule ensembles. In backfitting
procedure to find thresholds in each mixture component individually, the
[`partykit` package](PARTYKIT%20PACKAGE%20HERE)[citation](#citation). In
both instances, trees can be estimated with an offset from ensemble
machine learning predictions.

For many practical applications (e.g., vaccine efficacy trials),
observed data is often subject to a two-phase sampling mechanism (i.e.,
through the use of a two-stage design). In such cases, efficient
estimators (of both varieties) must be augmented to construct unbiased
estimates of the population-level causal parameter. Rose and van der
Laan (2011) first introduced an augmentation procedure that relies on
introducing inverse probability of censoring (IPC) weights directly to
an appropriate loss function or to the efficient influence function
estimating equation. `CVtreeMLE` extends this approach to compute
IPC-weighted one-step and TML estimators of the counterfactual mean
outcome under a shift stochastic treatment regime. The package is
designed to implement the statistical methodology described in Hejazi et
al. (2020) and extensions thereof.

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
treatment, consider the following example:

``` r
library(CVtreeMLE)
library(sl3)
set.seed(429153)
# simulate simple data
n_obs <- 500


sim_data <- simulate_mixture_cube(
  n_obs = n_obs, ## number of observations
  splits = c(0.99, 2.0, 2.5), ## split points for each mixture
  mins = c(0, 0, 0), ## minimum values for each mixture
  maxs = c(3, 4, 5), ## maximum value for each mixture
  mu = c(0, 0, 0),
  sigma = matrix(c(1, 0.5, 0.8, 0.5, 1, 0.7, 0.8, 0.7, 1), nrow = 3, ncol = 3),
  w1_betas = c(0.0, 0.01, 0.03, 0.06, 0.1, 0.05, 0.2, 0.04), ## subspace probability relationship with covariate W1
  w2_betas = c(0.0, 0.04, 0.01, 0.07, 0.15, 0.1, 0.1, 0.04), ## subspace probability relationship with covariate W2
  mix_subspace_betas = c(0.00, 0.08, 0.05, 0.01, 0.05, 0.033, 0.07, 0.09), ## probability of mixture subspace (for multinomial outcome gen)
  subspace_assoc_strength_betas = c(0, 0, 0, 0, 0, 0, 6, 0), ## mixture subspace impact on outcome Y
  marginal_impact_betas = c(0, 0, 0), ## marginal impact of mixture component on Y
  eps_sd = 0.01, ## random error
  binary = FALSE
)

# make SL learner libraries for fitting g and Q
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

SL.library<- c('SL.randomForest',
               'SL.earth',
               "SL.glm",
               "SL.mean")

## run on simulated data
sim_results <- CVtreeMLE(data = sim_data,
                                   W = c("W", "W2"),
                                   Y = "y",
                                   A = c(paste("M", seq(3), sep = "")),
                                   back_iter_SL = Q1_stack,
                                   SL.library = SL.library,
                                   n_folds = 2,
                                   family = "gaussian",
                                   H.AW_trunc_lvl = 10,
                                   parallel = TRUE,
                                   verbose = FALSE)
#> [1] "iter:  1 SL:  0.27288655720532 ctree: 0.27288655720532 Diff:  0 Rule:"
#> [1] "iter:  2 SL:  0.272886304726474 ctree: 0.27288655720532 Diff:  2.52478846052284e-07 Rule:"
#> [1] "iter:  1 SL:  0.312034027924391 ctree: 0.312034027924391 Diff:  0 Rule:"
#> [1] "iter:  2 SL:  0.312034027924391 ctree: 0.312034027924391 Diff:  0 Rule:"
#> [1] "iter:  1 SL:  0.00328254089146472 ctree: 0.272450429463324 Diff:  0.26916788857186 Rule: M1 <= 0.698033232091907"                        
#> [2] "iter:  1 SL:  0.00328254089146472 ctree: 0.272450429463324 Diff:  0.26916788857186 Rule: M1 > 0.698033232091907 & M1 <= 1.92445913525497"
#> [3] "iter:  1 SL:  0.00328254089146472 ctree: 0.272450429463324 Diff:  0.26916788857186 Rule: M1 > 0.698033232091907 & M1 > 1.92445913525497" 
#> [1] "iter:  2 SL:  0.00366496758785031 ctree: 0.000901298530461516 Diff:  0.00276366905738879 Rule: M1 <= 0.698033232091907"
#> [2] "iter:  2 SL:  0.00366496758785031 ctree: 0.000901298530461516 Diff:  0.00276366905738879 Rule: M1 > 0.698033232091907" 
#> [1] "iter:  1 SL:  0.27224327763226 ctree: 0.272133220878407 Diff:  0.000110056753852394 Rule: M2 <= 2.00110418303427 & M2 <= 1.10148691524405"
#> [2] "iter:  1 SL:  0.27224327763226 ctree: 0.272133220878407 Diff:  0.000110056753852394 Rule: M2 <= 2.00110418303427 & M2 > 1.10148691524405" 
#> [3] "iter:  1 SL:  0.27224327763226 ctree: 0.272133220878407 Diff:  0.000110056753852394 Rule: M2 > 2.00110418303427"                          
#> [1] "iter:  2 SL:  0.271569590806612 ctree: 0.27224327763226 Diff:  0.000673686825648112 Rule: M2 <= 2.00110418303427 & M2 <= 1.10148691524405"
#> [2] "iter:  2 SL:  0.271569590806612 ctree: 0.27224327763226 Diff:  0.000673686825648112 Rule: M2 <= 2.00110418303427 & M2 > 1.10148691524405" 
#> [3] "iter:  2 SL:  0.271569590806612 ctree: 0.27224327763226 Diff:  0.000673686825648112 Rule: M2 > 2.00110418303427"                          
#> [1] "iter:  1 SL:  0.272804634550957 ctree: 0.272772560110105 Diff:  3.20744408526807e-05 Rule: M3 <= 2.53274239624712 & M3 <= 1.3571113984252"
#> [2] "iter:  1 SL:  0.272804634550957 ctree: 0.272772560110105 Diff:  3.20744408526807e-05 Rule: M3 <= 2.53274239624712 & M3 > 1.3571113984252" 
#> [3] "iter:  1 SL:  0.272804634550957 ctree: 0.272772560110105 Diff:  3.20744408526807e-05 Rule: M3 > 2.53274239624712"                         
#> [1] "iter:  2 SL:  0.272132675381959 ctree: 0.272804634550957 Diff:  0.000671959168998038 Rule: M3 <= 2.53274239624712 & M3 <= 1.21148852076342"
#> [2] "iter:  2 SL:  0.272132675381959 ctree: 0.272804634550957 Diff:  0.000671959168998038 Rule: M3 <= 2.53274239624712 & M3 > 1.21148852076342" 
#> [3] "iter:  2 SL:  0.272132675381959 ctree: 0.272804634550957 Diff:  0.000671959168998038 Rule: M3 > 2.53274239624712"                          
#> [1] "iter:  1 SL:  0.310644290661673 ctree: 0.310939596540529 Diff:  0.000295305878855401 Rule: M1 <= 0.934177674302425"                        
#> [2] "iter:  1 SL:  0.310644290661673 ctree: 0.310939596540529 Diff:  0.000295305878855401 Rule: M1 > 0.934177674302425 & M1 <= 1.87817514713179"
#> [3] "iter:  1 SL:  0.310644290661673 ctree: 0.310939596540529 Diff:  0.000295305878855401 Rule: M1 > 0.934177674302425 & M1 > 1.87817514713179" 
#> [1] "iter:  2 SL:  0.310112532942556 ctree: 0.310644290661673 Diff:  0.000531757719117454 Rule: M1 <= 0.934177674302425"
#> [2] "iter:  2 SL:  0.310112532942556 ctree: 0.310644290661673 Diff:  0.000531757719117454 Rule: M1 > 0.934177674302425" 
#> [1] "iter:  1 SL:  0.311455600437983 ctree: 0.311439709389857 Diff:  1.58910481259378e-05 Rule: M2 <= 2.00708027766665"
#> [2] "iter:  1 SL:  0.311455600437983 ctree: 0.311439709389857 Diff:  1.58910481259378e-05 Rule: M2 > 2.00708027766665" 
#> [1] "iter:  2 SL:  0.311452290431762 ctree: 0.311455600437983 Diff:  3.31000622177946e-06 Rule: M2 <= 2.00708027766665"                        
#> [2] "iter:  2 SL:  0.311452290431762 ctree: 0.311455600437983 Diff:  3.31000622177946e-06 Rule: M2 > 2.00708027766665 & M2 <= 2.92370001907995"
#> [3] "iter:  2 SL:  0.311452290431762 ctree: 0.311455600437983 Diff:  3.31000622177946e-06 Rule: M2 > 2.00708027766665 & M2 > 2.92370001907995" 
#> [1] "iter:  1 SL:  0.311794455780445 ctree: 0.311930500577001 Diff:  0.000136044796556101 Rule: M3 <= 2.47762214185431"
#> [2] "iter:  1 SL:  0.311794455780445 ctree: 0.311930500577001 Diff:  0.000136044796556101 Rule: M3 > 2.47762214185431" 
#> [1] "iter:  2 SL:  0.311261165216538 ctree: 0.311794455780445 Diff:  0.000533290563906885 Rule: M3 <= 2.5099802112645"                        
#> [2] "iter:  2 SL:  0.311261165216538 ctree: 0.311794455780445 Diff:  0.000533290563906885 Rule: M3 > 2.5099802112645 & M3 <= 3.74706103658575"
#> [3] "iter:  2 SL:  0.311261165216538 ctree: 0.311794455780445 Diff:  0.000533290563906885 Rule: M3 > 2.5099802112645 & M3 > 3.74706103658575"

sim_results$`Mixture Results`
#>   Mixture ATE Standard Error Lower CI Upper CI P-value P-value Adj   vars
#> 1    5.999865     0.01985289 5.960954 6.038776       0           0 M1M2M3
#>                                                    Mixture Interaction Rules
#> 1 M1 > 0.021 & M1 < 0.934 & M2 > 2.034 & M2 < 3.994 & M3 > 2.534 & M3 < 4.98
#>   Fraction Covered
#> 1        0.9666667
sim_results$`Marginal Results`
#>    Marginal ATE Standard Error  Lower CI Upper CI      P-value  P-value Adj
#> M1    1.3958243      0.1862784 1.0307254 1.760923 6.720602e-14 2.016181e-13
#> M2    1.1846510      0.2186746 0.7560567 1.613245 6.046683e-08 1.814005e-07
#> M3    0.9885222      0.2256091 0.5463365 1.430708 1.178267e-05 3.534802e-05
#>                 Marginal Rules Fraction Overlap          Min      Max
#> M1 M1 > 0.00859 & M1 < 0.69216        0.7945205 0.0011496947 2.988107
#> M2  M2 > 2.9311 & M2 < 3.99468        0.5141700 0.0002521947 3.994680
#> M3 M3 > 3.75935 & M3 < 4.98025        0.5364807 0.0028737772 4.980251
```

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

        @article{hejazi2020efficient,
          author = {Hejazi, Nima S and {van der Laan}, Mark J and Janes, Holly
            E and Gilbert, Peter B and Benkeser, David C},
          title = {Efficient nonparametric inference on the effects of
            stochastic interventions under two-phase sampling, with
            applications to vaccine efficacy trials},
          year = {2020},
          doi = {10.1111/biom.13375},
          url = {https://doi.org/10.1111/biom.13375},
          journal = {Biometrics},
          publisher = {Wiley Online Library}
        }

        @article{hejazi2020CVtreeMLE-joss,
          author = {Hejazi, Nima S and Benkeser, David C},
          title = {{CVtreeMLE}: Efficient estimation of the causal effects of
            stochastic interventions in {R}},
          year  = {2020},
          doi = {10.21105/joss.02447},
          url = {https://doi.org/10.21105/joss.02447},
          journal = {Journal of Open Source Software},
          publisher = {The Open Journal}
        }

        @software{hejazi2020CVtreeMLE-rpkg,
          author = {Hejazi, Nima S and Benkeser, David C},
          title = {{CVtreeMLE}: Efficient Estimation of the Causal Effects of
            Stochastic Interventions},
          year  = {2020},
          doi = {10.5281/zenodo.4070042},
          url = {https://CRAN.R-project.org/package=CVtreeMLE},
          note = {R package version 0.3.4}
        }

------------------------------------------------------------------------

## Related

-   [R/`tmle3shift`](https://github.com/tlverse/tmle3shift) - An R
    package providing an independent implementation of the same core
    routines for the TML estimation procedure and statistical
    methodology as is made available here, through reliance on a unified
    interface for Targeted Learning provided by the
    [`tmle3`](https://github.com/tlverse/tmle3) engine of the [`tlverse`
    ecosystem](https://github.com/tlverse).

-   [R/`medshift`](https://github.com/blind-contours/medshift) - An R
    package providing facilities to estimate the causal effect of
    stochastic treatment regimes in the mediation setting, including
    classical (IPW) and augmented double robust (one-step) estimators.
    This is an implementation of the methodology explored by Dı́az and
    Hejazi (2020).

-   [R/`haldensify`](https://github.com/blind-contours/haldensify) - A
    minimal package for estimating the conditional density treatment
    mechanism component of this parameter based on using the [highly
    adaptive lasso](https://github.com/tlverse/hal9001)
    (**coyle2020hal9001-rpkg?**; Hejazi, Coyle, and van der Laan 2020)
    in combination with a pooled hazard regression. This package
    implements a variant of the approach advocated by Dı́az and van der
    Laan (2011).

------------------------------------------------------------------------

## Funding

The development of this software was supported in part through grants
from the National Library of Medicine (award number [T32
LM012417](https://reporter.nih.gov/project-details/9248418)) and the
National Institute of Allergy and Infectious Diseases (award number [R01
AI074345](https://reporter.nih.gov/project-details/9926564)) of the
National Institutess of Health.

------------------------------------------------------------------------

## License

© 2017-2021 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    Copyright (c) 2017-2021 Nima S. Hejazi
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

<div id="ref-diaz2020causal" class="csl-entry">

Dı́az, Iván, and Nima S Hejazi. 2020. “Causal Mediation Analysis for
Stochastic Interventions.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 82 (3): 661–83.
<https://doi.org/10.1111/rssb.12362>.

</div>

<div id="ref-diaz2011super" class="csl-entry">

Dı́az, Iván, and Mark J van der Laan. 2011. “Super Learner Based
Conditional Density Estimation with Application to Marginal Structural
Models.” *The International Journal of Biostatistics* 7 (1): 1–20.

</div>

<div id="ref-diaz2012population" class="csl-entry">

———. 2012. “Population Intervention Causal Effects Based on Stochastic
Interventions.” *Biometrics* 68 (2): 541–49.

</div>

<div id="ref-hejazi2020hal9001-joss" class="csl-entry">

Hejazi, Nima S, Jeremy R Coyle, and Mark J van der Laan. 2020. “<span
class="nocase">hal9001</span>: Scalable Highly Adaptive Lasso Regression
in R.” *Journal of Open Source Software* 5 (53): 2526.
<https://doi.org/10.21105/joss.02526>.

</div>

<div id="ref-hejazi2020efficient" class="csl-entry">

Hejazi, Nima S, Mark J van der Laan, Holly E Janes, Peter B Gilbert, and
David C Benkeser. 2020. “Efficient Nonparametric Inference on the
Effects of Stochastic Interventions Under Two-Phase Sampling, with
Applications to Vaccine Efficacy Trials.” *Biometrics*.
<https://doi.org/10.1111/biom.13375>.

</div>

<div id="ref-rose2011targeted2sd" class="csl-entry">

Rose, Sherri, and Mark J van der Laan. 2011. “A Targeted Maximum
Likelihood Estimator for Two-Stage Designs.” *The International Journal
of Biostatistics* 7 (1): 1–21.

</div>

</div>
