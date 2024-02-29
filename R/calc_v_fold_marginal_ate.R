#' @title Calculate the v-fold specifice ATE for each rule found for individual
#' mixture components.
#' @description Aggregate marginal rules for each mixture variable found across
#' the folds. For each rule
#' extract the relevant nuisance parameter data calculated in the folds. Given
#' the validation data estimates across
#' the folds, for each tree do a TMLE update step to target the average
#' treatment effect. Update the initial counterfactuals,
#' calculate the influence curve and using the influence curve calculate
#' variance estimates and p-values.
#'
#' @param marginal_data List of dataframes of nuisance parameter
#' data for each mixture
#' @param mix_comps Vector of characters indicating the mixture components
#' @param y Vector indicating the Y
#' @param n_folds Number of folds used in cross-validation
#' @param marginal_rules List of dataframes of marginal rules found across
#' the folds

#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif
#' @importFrom stats rnorm runif
#' @importFrom rlang :=
#' @return A list of the marginal results for each fold including:
#'  \itemize{
#'   \item \code{marginal_results}: A data frame with the data adpatively
#'   determined mixture component thresholds on the rows and ATE, variance,
#'   and RMSE estimates on the columns.
#'   \item \code{data}: A list of data frames for each mixture component
#'   threshold evaluated as the exposure, baseline covariates, outcome,
#'   nuisance parameter estimates, marginal ATE and the influence curve.
#'   }
#'
#' @export

calc_v_fold_marginal_ate <- function(marginal_data,
                                     mix_comps,
                                     marginal_rules,
                                     y,
                                     n_folds) {
  marg_data <- unlist(marginal_data, recursive = FALSE)
  marg_data <- unlist(marg_data, recursive = FALSE)

  marg_data <- marg_data[!sapply(marg_data, is.null)]

  marg_data <- marg_data[sapply(marg_data, FUN = function(i) {
    all(is.na(i$h_aw)) != TRUE
  })]

  marg_rules <- unlist(marginal_rules, recursive = FALSE)
  marg_rules <- marg_rules[!sapply(marg_rules, is.null)]
  marg_rules <- marg_rules[sapply(marg_rules, FUN = function(i) {
    dim(i)[1] != 0
  })]

  comp_labels <- list()
  rule_comp_cols <- list()
  rule_ref_cols <- list()

  index <- 0
  for (i in seq(marg_rules)) {
    fold_rules_groups <- marginal_group_split(marg_rules[[i]])

    for (j in seq(fold_rules_groups)) {
      rules <- fold_rules_groups[[j]]
      reference <- rules[rules$quantile == 1, ]
      reference_quant <- reference$var_quant_group
      reference_rule <- reference$rules
      comparisons <- rules[rules$quantile > 1, ]
      for (k in seq(nrow(comparisons))) {
        index <- index + 1
        comp_row <- comparisons[k, ]
        comp_label <- paste(comp_row$var_quant_group, reference_quant,
          sep = "-"
        )
        rule_comparison <- comp_row$rules
        rule_reference <- reference_rule
        comp_labels[[index]] <- comp_label
        rule_comp_cols[[index]] <- rule_comparison
        rule_ref_cols[[index]] <- rule_reference
      }
    }
  }

  rule_comparisons <- cbind(unlist(rule_ref_cols), unlist(rule_comp_cols))
  colnames(rule_comparisons) <- c("Reference_Rule", "Comparison_Rule")
  comp_labels <- unlist(comp_labels)

  marginal_results <-
    as.data.frame(matrix(
      data = NA,
      nrow = length(marg_data),
      ncol = 7
    ))

  colnames(marginal_results) <-
    c(
      "Marginal ATE",
      "Standard Error",
      "Lower CI",
      "Upper CI",
      "P-value",
      "P-value Adj",
      "RMSE"
    )

  marginal_results$Levels <- comp_labels

  updated_marginal_data <- list()

  for (i in seq(marg_data)) {
    marg_mix <- marg_data[[i]]

    if (all(is.na(marg_mix$h_aw)) == TRUE) {
      break
    }

    flux_results <- fit_least_fav_submodel(
      h_aw = marg_mix$h_aw,
      data = marg_mix,
      y = y,
      qbar_aw = marg_mix$qbar_aw,
      qbar_1w = marg_mix$qbar_1w
    )

    # ## back-scale Y
    # qbar_aw_star <- scale_to_original(scaled_vals = flux_results$qbar_aw_star,
    #                                  max_orig = max(marg_mix[y]),
    #                                  min_orig = min(marg_mix[y]))
    #
    # qbar_0w_star <- scale_to_original(scaled_vals = flux_results$qbar_0w_star,
    #                                  max_orig = max(marg_mix[y]),
    #                                  min_orig = min(marg_mix[y]))
    #
    # qbar_1w_star <- scale_to_original(scaled_vals = flux_results$qbar_1w_star,
    #                                  max_orig = max(marg_mix[y]),
    #                                  min_orig = min(marg_mix[y]))

    marg_mix$qbar_aw_star <- flux_results$qbar_aw_star
    marg_mix$qbar_0w_star <- flux_results$qbar_0w_star
    marg_mix$qbar_1w_star <- flux_results$qbar_1w_star

    ate_results <- calc_ate_estimates(
      data = marg_mix,
      ate_var = "marg_ate",
      y = y,
      p_adjust_n = length(marg_data)
    )

    sqrd_resids <- (marg_mix$qbar_aw_star - marg_mix[y])^2
    rmse <- sqrt(mean(sqrd_resids[, 1]))

    marginal_results$`Marginal ATE`[i] <- round(ate_results$ate, 3)
    marginal_results$`Standard Error`[i] <- round(ate_results$se, 3)
    marginal_results$`Lower CI`[i] <- round(ate_results$ci[1], 3)
    marginal_results$`Upper CI`[i] <- round(ate_results$ci[2], 3)
    marginal_results$`P-value`[i] <- round(ate_results$p_value, 6)
    marginal_results$`P-value Adj`[i] <- round(ate_results$adj_p_value, 3)
    marginal_results$RMSE[i] <- round(rmse, 3)

    updated_marginal_data[[i]] <- ate_results$data
  }

  marginal_results <- cbind(marginal_results, rule_comparisons)

  return(list(
    marginal_results = marginal_results,
    data = updated_marginal_data
  ))
}
