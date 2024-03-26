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
#' across all the folds
#'
#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows
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
                              n_folds) {
  fold_mix_rules <-
    data.table::rbindlist(unlist(input_mix_rules, recursive = FALSE))

  fold_mix_data <- bind_rows(unlist(input_mix_data, recursive = FALSE))

  ## Calculate Pooled TMLE Across Max Region

  flux_results <- fit_least_fav_submodel(
    h_aw = fold_mix_data$h_aw,
    y = y,
    data = fold_mix_data,
    qbar_aw = fold_mix_data$qbar_aw,
    qbar_1w = fold_mix_data$qbar_1w
  )

  fold_mix_data$qbar_aw_star <- flux_results$qbar_aw_star
  fold_mix_data$qbar_1w_star <- flux_results$qbar_1w_star

  ate_results <- calc_ate_estimates(
    data = fold_mix_data,
    ate_var = "mix_ate",
    y = y,
    p_adjust_n = NULL
  )

  max_are_results <-
    as.data.frame(matrix(
      data = NA,
      nrow = 1,
      ncol = 5
    ))

  colnames(max_are_results) <-
    c(
      "Region ARE",
      "Standard Error",
      "Lower CI",
      "Upper CI",
      "P-value"
    )

  max_are_results$`Region ARE` <- round(ate_results$ate, 3)
  max_are_results$`Standard Error` <- round(ate_results$se, 3)
  max_are_results$`Lower CI` <- round(ate_results$ci[1], 3)
  max_are_results$`Upper CI` <- round(ate_results$ci[2], 3)
  max_are_results$`P-value` <- round(ate_results$p_value, 6)


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

  mixture_results_inv_variance <- mixture_results

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

    mix_rule_data <- bind_rows(intxn_rule_data_list)
    if (any(is.na(mix_rule_data$qbar_aw))) {
      next
    }

    flux_results <- fit_least_fav_submodel(
      h_aw = mix_rule_data$h_aw,
      y = y,
      data = mix_rule_data,
      qbar_aw = mix_rule_data$qbar_aw,
      qbar_1w = mix_rule_data$qbar_1w
    )

    mix_rule_data$qbar_aw_star <- flux_results$qbar_aw_star
    mix_rule_data$qbar_1w_star <- flux_results$qbar_1w_star

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

    if (nrow(intx_group) < n_folds) {
      n_1 <- nrow(intx_group)
      n_0 <- n_folds - n_1

      # Variance of the current pooled estimates
      var_pooled <- ate_results$se^2

      # Estimating the variance of the "null" estimate
      var_null <- var_pooled * n_1 / n_0
      sd_null <- sqrt(var_null)


      # Weights
      w_1 <- 1 / var_pooled
      w_0 <- 1 / var_null

      # Calculate new pooled estimate and its variance
      pooled_estimate <- ate_results$ate
      null_samples <- 0
      new_pooled_psi <- (w_1 * pooled_estimate + w_0 * null_samples) / (w_1 + w_0)
      var_new_pooled <- 1 / (w_1 + sum(w_0))

      # Calculate standard error and confidence intervals
      se_new_pooled <- sqrt(var_new_pooled)
      lower_ci <- new_pooled_psi - 1.96 * se_new_pooled
      upper_ci <- new_pooled_psi + 1.96 * se_new_pooled

      p_value_pooled <- 2 * stats::pnorm(abs(new_pooled_psi / se_new_pooled), lower.tail = F)

      mixture_results_inv_variance$`Mixture ATE`[group] <- round(new_pooled_psi, 3)
      mixture_results_inv_variance$`Standard Error`[group] <- round(se_new_pooled, 3)
      mixture_results_inv_variance$`Lower CI`[group] <- round(lower_ci, 3)
      mixture_results_inv_variance$`Upper CI`[group] <- round(upper_ci, 3)
      mixture_results_inv_variance$`P-value`[group] <- round(p_value_pooled, 6)
      mixture_results_inv_variance$`Vars`[group] <- vars
      mixture_results_inv_variance$`RMSE`[group] <- round(rmse, 3)
    }

    mixture_data_list[[group]] <- mix_rule_data
  }

  return(list(
    results = mixture_results,
    region_results = max_are_results,
    inv_var_results = mixture_results_inv_variance,
    group_list = group_list,
    mixture_data_list = mixture_data_list
  ))
}
