% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_marg_nuisance_params.R
\name{est_marg_nuisance_params}
\alias{est_marg_nuisance_params}
\title{Estimate nuisance parameters for each marginal mixture component}
\usage{
est_marg_nuisance_params(
  at,
  av,
  w,
  aw_stack,
  family,
  a,
  no_marg_rules,
  marg_decisions,
  h_aw_trunc_lvl,
  parallel_cv,
  seed
)
}
\arguments{
\item{at}{Training data}

\item{av}{Validation data}

\item{w}{Vector of characters denoting covariates}

\item{aw_stack}{Super Learner library for fitting Q (outcome mechanism)
and g (treatment mechanism)}

\item{family}{Binomial or gaussian}

\item{a}{Vector of characters that denote the mixture components}

\item{no_marg_rules}{TRUE/FALSE if no marginal rules were
found across the folds}

\item{marg_decisions}{List of rules found within the fold for each
mixture component}

\item{h_aw_trunc_lvl}{Truncation level of the clever covariate
(induces more bias to reduce variance)}

\item{parallel_cv}{TRUE/FALSE if cv parallelization is used}

\item{seed}{Seed number}
}
\value{
Rules object. TODO: add more detail here.
}
\description{
For each marginal mixture component rule found, create a
g estimator for the probability of being exposed to the rule thresholds,
and a Q estimator for the outcome E(Y| A = a_mix, W).
Get estimates of g and Q using the validation data and
calculate the clever covariate used in the TMLE fluctuation step.
}
