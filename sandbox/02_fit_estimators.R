
fit_estimators <- function(data,
                           covars,
                           exposures,
                           outcome,
                           seed,
                           P_0_data,
                           true_rule,
                           true_ate) {

  sim_results <- CVtreeMLE(data = data,
                           w = covars,
                           a = exposures,
                           y = outcome,
                           n_folds = 4,
                           num_cores = 20,
                           family = "continuous",
                           direction = "positive",
                           parallel = TRUE,
                           parallel_cv = FALSE,
                           seed = seed)

  max_ate_index <- which.max(abs(
    sim_results$`Pooled TMLE Mixture Results`$`Mixture ATE`))

  if (length(max_ate_index) == 1) {
    mixure_found_ind <- 1
  } else{
    mixure_found_ind <- 0
  }

  if (mixure_found_ind == 1) {
    mixture_results <- sim_results$`Pooled TMLE Mixture Results`[max_ate_index,]
    lower <- mixture_results$`Lower CI`
    upper <- mixture_results$`Upper CI`
    ate <- mixture_results$`Mixture ATE`
    mix_decisions <- mixture_results$Union_Rule

    DA_true_results <- calc_empir_truth(P_0_data, mix_decisions)
    DA_P0_truth <- DA_true_results[3]

    rule_binary <- data %>%
      dplyr::transmute(A = ifelse(eval(parse(text = mix_decisions)), 1, 0))

    rule_applied_P_0 <- P_0_data %>%
      dplyr::mutate(A = ifelse(eval(parse(text = mix_decisions)), 1, 0))

    true_rule_binary <- data %>%
      dplyr::mutate(A = ifelse(eval(parse(text = true_rule)), 1, 0))

    DA_rule_bias <- ate - DA_P0_truth
    bias <- true_ate - ate

    true_coverage = ifelse(
      (lower <= true_ate & true_ate <= upper), 1,0
    )

    da_covererage = ifelse(
      (lower <= DA_P0_truth & DA_P0_truth <= upper), 1,0
    )

    confusion_table <- table(true_rule_binary$A, rule_binary$A)

    true_pos <- confusion_table[2,2] / sum(true_rule_binary$A)
    true_neg <- confusion_table[1,1] / (length(true_rule_binary$A) -
                                          sum(true_rule_binary$A))
    false_pos <- confusion_table[1,2] / (length(true_rule_binary$A) -
                                           sum(true_rule_binary$A))
    false_neg <- confusion_table[2,1] / sum(true_rule_binary$A)

  }else{
    mixure_found_ind <- 0
    ate <- NULL
    lower <- NULL
    upper <- NULL
    true_pos <- NULL
    true_neg <- NULL
    false_pos <- NULL
    false_neg <- NULL
    DA_rule_bias <- NULL
    DA_P0_truth <- NULL
    mix_decisions <- NULL
    true_coverage <- NULL
    da_covererage <- NULL
    bias <- NULL
    true_rule <- NULL
  }

  sim_out <- list("Mix ind" = mixure_found_ind,
                  "ATE" = ate,
                  "Lower" = lower,
                  "Upper" = upper,
                  "DA Rule"= mix_decisions,
                  "True Pos" = true_pos,
                  "True Neg" = true_neg,
                  "False Pos" = false_pos,
                  "False Neg" = false_neg,
                  "DA Rule Bias" = DA_rule_bias,
                  "DA Truth" = DA_P0_truth,
                  "Bias" = bias,
                  "True Cov" = true_coverage,
                  "DA Cov" = da_covererage,
                  "True ATE" = true_ate,
                  "True Rule" = true_rule)
  return(sim_out)
}
