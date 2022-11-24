#' @title Calculate the ATE for each mixture rule
#'
#' @description Aggregate mixture rules found across the folds that have
#' the same variables. For each rule
#' extract the relevant nuisance parameter data calculated in the folds.
#' Given the validation data estimates across
#' the folds, for each tree do a TMLE update step to target the average
#' treatment effect. Update the initial counterfactuals,
#' calculate the influence curve and using the influence curve calculate
#' variance estimates and p-values.
#'
#' @param input_mix_rules List of dataframes of rules found for a mixture
#' across the folds
#' @param input_mix_data Nuisance parameter data for mixture rules found
#' @param y Outcome variable
#' across the folds
#' @param n_folds Number of folds used in cross-validation
#' @param no_mixture_rules TRUE/FALSE whether no mixture rules were found
#' across all the folds
#'
#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by
#'
#' @return A list with mixture analysis results which includes:
#' \itemize{
#'   \item \code{results}: A data frame with variable threshold combinations
#'   on the rows and ATE, variance and consistency estimates on the columns.
#'   \item \code{group_list}: A list of rule combinations found in the
#'   ensemble decision tree model grouped by variable sets in the rules and
#'   directions of the coefficient in the linear model. Also provided is the
#'   fold the rule was found and the RMSE of the model which used the rule.
#'   \item \code{mixture_data_list}: A list of data frames which houses the
#'   data for each mixture rule evaluated as an exposure, the baseline
#'   covariates, outcome, and nuisance parameter estimates for the respective
#'   rule.
#' }
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm
#' @importFrom stats qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export

calc_mixtures_ate <- function(input_mix_rules,
                              input_mix_data,
                              y,
                              n_folds,
                              no_mixture_rules) {
  fold_mix_rules <-
    data.table::rbindlist(unlist(input_mix_rules, recursive = FALSE))

  fold_mix_rules <-
    fold_mix_rules %>%
    dplyr::group_by(test, direction)

  groups <- fold_mix_rules %>%
    dplyr::group_by(test, direction)

  group_list <- dplyr::group_split(groups)

  mixture_results <-
    as.data.frame(matrix(
      data = NA,
      nrow = length(group_list),
      ncol = 8
    ))

  colnames(mixture_results) <-
    c(
      "Mixture ATE",
      "Standard Error",
      "Lower CI",
      "Upper CI",
      "P-value",
      "P-value Adj",
      "Vars",
      "RMSE"
    )

  mixture_data_list <- list()

  for (group in seq(group_list)) {
    intx_group <- group_list[[group]]
    intxn_rule_data_list <- list()
    vars <- unique(intx_group$test)

    for (i in seq(dim(intx_group)[1])) {
      intxn_rule <- intx_group[i, ]
      search_data <-
        as.data.frame(input_mix_rules[as.numeric(intxn_rule$fold)])
      srch_indx <- match(intxn_rule$rule, search_data[, 1])
      fold_data <- input_mix_data[[as.numeric(intxn_rule$fold)]]
      intx_rule_data <- fold_data[[1]][srch_indx]
      intxn_rule_data_list[[i]] <- intx_rule_data
    }

    mix_rule_data <- do.call(rbind, unlist(intxn_rule_data_list,
      recursive = FALSE
    ))

    flux_results <- fit_least_fav_submodel(
      h_aw = mix_rule_data$h_aw,
      y = y,
      data = mix_rule_data,
      qbar_aw = mix_rule_data$qbar_aw,
      qbar_1w = mix_rule_data$qbar_1w,
      qbar_0w = mix_rule_data$qbar_0w
    )

    mix_rule_data$qbar_aw_star <- flux_results$qbar_aw_star
    mix_rule_data$qbar_1w_star <- flux_results$qbar_1w_star
    mix_rule_data$qbar_0w_star <- flux_results$qbar_0w_star

    ate_results <- calc_ate_estimates(
      data = mix_rule_data,
      ate_var = "mix_ate",
      y = y,
      p_adjust_n = length(group_list)
    )

    ## calculate RMSE for Y| A = rule i, W
    sqrd_resids <- (mix_rule_data$qbar_aw_star - mix_rule_data[y])^2
    rmse <- sqrt(mean(sqrd_resids[, 1]))

    mixture_results$`Mixture ATE`[group] <- round(ate_results$ate, 3)
    mixture_results$`Standard Error`[group] <- round(ate_results$se, 3)
    mixture_results$`Lower CI`[group] <- round(ate_results$ci[1], 3)
    mixture_results$`Upper CI`[group] <- round(ate_results$ci[2], 3)
    mixture_results$`P-value`[group] <- round(ate_results$p_value, 6)
    mixture_results$`P-value Adj`[group] <- round(
      ate_results$adj_p_value,
      6
    )
    mixture_results$`Vars`[group] <- vars
    mixture_results$`RMSE`[group] <- round(rmse, 3)

    mixture_data_list[[group]] <- mix_rule_data
  }

  return(list(
    results = mixture_results,
    group_list = group_list,
    mixture_data_list = mixture_data_list
  ))
}
