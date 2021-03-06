% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_common_marginal_rules.R
\name{find_common_marginal_rules}
\alias{find_common_marginal_rules}
\title{Create a new rule based on observations that meet every rule across
the folds for each mixture.}
\usage{
find_common_marginal_rules(
  fold_rules,
  data,
  mix_comps,
  marginal_results,
  n_folds
)
}
\arguments{
\item{fold_rules}{List of rules found for each mixture component found
across the folds}

\item{data}{Full data which rules are evaluated}

\item{mix_comps}{Vector of mixture components}

\item{marginal_results}{Dataframe holding the results for each marginal
component rule}

\item{n_folds}{Total number of folds}
}
\value{
Rules object. TODO: add more detail here.

Dataframe with rules, threshold regions, proportion in folds and
min/max values
}
\description{
For each mixture component, a different rule could be found for
each fold. Therefore, it is necessary to create one rule for each mixture
component that can be interpreted
as a common rule across the folds. To do this, observations that meet all
rules for all folds are determined. Then a new rule is created for these
observations. A coverage metric is
calculated which is the fraction of observations in this common rule
compared to the sum of observations that met each rule. The coverage metric
can be thought of as an intersection over
union or the Jaccard coefficient for the rules.
}
