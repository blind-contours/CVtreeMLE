#' @title Calculate the v-fold specifice ATE for each rule found for individual mixture components.
#' @description Aggregate marginal rules for each mixture variable found across the folds. For each rule
#' extract the relevant nuisance parameter data calculated in the folds. Given the validation data estimates across
#' the folds, for each tree do a TMLE update step to target the average treatment effect. Update the initial counterfactuals,
#' calculate the influence curve and using the influence curve calculate variance estimates and p-values.
#'
#' @param marginal_data List of dataframes of nuisance parameter data for each mixture
#' @param mix_comps Vector of characters indicating the mixture components
#' @param Y Vector indicating the Y
#' @param n_folds Number of folds used in cross-validation
#' @param marginal_rules List of dataframes of marginal rules found across the folds

#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows

#' @return Rules object. TODO: add more detail here.
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export

calc_v_fold_marginal_ate <- function(marginal_data, mix_comps, marginal_rules, Y, n_folds) {
  marg_data <- list()

  x <- unlist(marginal_data, recursive = FALSE)
  x <- unlist(x, recursive = FALSE)
  x <- x[!sapply(x, is.null)]

  y <- unlist(marginal_rules, recursive = FALSE)
  y <- y[!sapply(y, is.null)]

  comp_labels <- list()
  rule_comp_cols <- list()
  rule_ref_cols <- list()

  index <- 0
  for (i in seq(y)) {

    fold_rules_groups <- marginal_group_split(y[[i]])

    for (j in seq(fold_rules_groups)) {
      # var_comps <- list()
      # rule_comps <- list()
      rules <- fold_rules_groups[[j]]
      reference <- rules[rules$quantile == 1, ]
      reference_quant <- reference$var_quant_group
      reference_rule <- reference$rules
      comparisons <- rules[rules$quantile > 1, ]
      for (k in seq(nrow(comparisons))) {
        index <- index + 1
        comp_row <- comparisons[k, ]
        comp_label <- paste(comp_row$var_quant_group, reference_quant, sep = "-")
        rule_comparison <- comp_row$rules
        rule_reference <- reference_rule

        # var_comps[k] <- comp_label
        comp_labels[[index]] <-  comp_label
        rule_comp_cols[[index]] <- rule_comparison
        rule_ref_cols[[index]] <- rule_reference
      }
    }
  }

  rule_comparisons <- cbind(unlist(rule_ref_cols), unlist(rule_comp_cols))
  colnames(rule_comparisons) <- c("Reference", "Comparison")
  comp_labels <- unlist(comp_labels)

  marginal_results <-
    as.data.frame(matrix(
      data = NA,
      nrow = length(x),
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

  marginal_results$comparison <- comp_labels

  updated_marginal_data <- list()

  ## find which marginal mixture rule was null across the folds to update bonferonni p-value

  no_null_indices <- list()


  for (i in seq(x)) {
    marg_mix <- x[[i]]

    flux_results <- fit_least_fav_submodel(H.AW = marg_mix$H.AW, data = marg_mix, QbarAW = marg_mix$QbarAW, Qbar1W = marg_mix$Qbar1W, Qbar0W = marg_mix$Qbar0W)

    ## back-scale Y
    QbarAW.star <- scale_to_original(scaled_vals = flux_results$QbarAW.star, max_orig = max(marg_mix[Y]), min_orig = min(marg_mix[Y]))
    Qbar0W.star <- scale_to_original(scaled_vals = flux_results$Qbar0W.star, max_orig = max(marg_mix[Y]), min_orig = min(marg_mix[Y]))
    Qbar1W.star <- scale_to_original(scaled_vals = flux_results$Qbar1W.star, max_orig = max(marg_mix[Y]), min_orig = min(marg_mix[Y]))

    marg_mix$QbarAW.star <- QbarAW.star
    marg_mix$Qbar0W.star <- Qbar0W.star
    marg_mix$Qbar1W.star <- Qbar1W.star

    ATE_results <- calc_ATE_estimates(data = marg_mix, ATE_var = "marg.ATE", outcome = Y, p_adjust_n = length(x))

    sqrd_resids <- (marg_mix$QbarAW.star - marg_mix[Y])^2
    RMSE <- sqrt(mean(sqrd_resids[, 1]))

    marginal_results$`Marginal ATE`[i] <- ATE_results$ATE
    marginal_results$`Standard Error`[i] <- ATE_results$SE
    marginal_results$`Lower CI`[i] <- ATE_results$CI[1]
    marginal_results$`Upper CI`[i] <- ATE_results$CI[2]
    marginal_results$`P-value`[i] <- ATE_results$`p-value`
    marginal_results$`P-value Adj`[i] <- ATE_results$`adj p-value`
    marginal_results$RMSE[i] <- RMSE

    updated_marginal_data[[i]] <- ATE_results$data
  }

  marginal_results <- cbind(marginal_results, rule_comparisons)

  return(list(marginal_results = marginal_results, data = updated_marginal_data))
}
