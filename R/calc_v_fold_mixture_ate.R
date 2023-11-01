#' @title Calculate the ATE for V-fold specific rules
#'
#' @description  For each rule in each fold extract the relevant nuisance
#' parameter data calculated in the folds. Given the validation data estimates
#' for each rule do a TMLE update step to target the average treatment effect.
#' Update the initial counterfactuals,calculate the influence curve
#' and using the influence curve calculate variance estimates and p-values.
#'
#' @param input_mix_rules List of dataframes of rules found for a mixture
#' across the folds
#' @param input_mix_data Nuisance parameter data for mixture rules found
#' across the folds
#' @param y Character indicating the outcome variable
#' @param n_folds Number of folds used in cross-validation
#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm
#' @importFrom stats qunif rnorm runif
#' @importFrom rlang :=
#' @return A data frame of mixture results for each variable level combination
#' for each fold
#'
#' @export

calc_v_fold_mixtures_ate <- function(input_mix_rules,
                                     input_mix_data,
                                     y,
                                     n_folds) {
  input_mix_rules <- unlist(input_mix_rules, recursive = FALSE)
  input_mix_data <- unlist(input_mix_data, recursive = FALSE)

  fold_mixture_results_list <- list()

  for (fold in seq(input_mix_data)) {
    fold_data <- input_mix_data[[fold]]
    fold_rules <- input_mix_rules[[fold]]

    rule_results_list <- list()

    for (rule in seq(fold_data)) {
      mix_data <- fold_data[[rule]]
      mix_rule_row <- fold_rules[rule, ]
      mix_rule <- mix_rule_row$description
      variables <- mix_rule_row$test

      no_rule_ind <- is.na(mix_rule_row$rule)

      if (no_rule_ind == TRUE) {
        break
      }

      flux_results <- fit_least_fav_submodel(
        h_aw = mix_data$h_aw,
        y = y,
        data = mix_data,
        qbar_aw = mix_data$qbar_aw,
        qbar_1w = mix_data$qbar_1w)

      mix_data$qbar_aw_star <- flux_results$qbar_aw_star
      mix_data$qbar_1w_star <- flux_results$qbar_1w_star
      # mix_data$qbar_0w_star <- flux_results$qbar_0w_star

      ate_results <- calc_ate_estimates(
        data = mix_data,
        ate_var = "mix_ate",
        y = y,
        p_adjust_n = length(input_mix_data),
        v_fold = TRUE
      )

      sqrd_resids <- (mix_data$qbar_aw_star - mix_data[y])^2
      rmse <- sqrt(mean(sqrd_resids[, 1], na.rm = TRUE))

      ate <- round(ate_results$ate, 3)
      se <- round(ate_results$se, 3)
      lower_ci <- round(ate_results$ci[1], 3)
      upper_ci <- round(ate_results$ci[2], 3)
      p_val <- round(ate_results$p_value, 6)
      p_val_adj <- round(ate_results$adj_p_value, 6)
      rmse <- round(rmse, 3)

      rule_results <- cbind.data.frame(
        ate, se, lower_ci, upper_ci,
        p_val, p_val_adj, rmse, mix_rule, fold, variables
      )

      rule_results_list[[rule]] <- rule_results
    }

    fold_results <- do.call(rbind, rule_results_list)
    fold_mixture_results_list[[fold]] <- fold_results
  }


  mixture_results <- do.call(rbind, fold_mixture_results_list)

  return(results = mixture_results)
}
