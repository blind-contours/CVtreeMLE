
fit_estimators <- function(data,
                           covars,
                           exposures,
                           outcome,
                           seed,
                           P_0_data,
                           true_rule,
                           true_ate,
                           true_region_mean,
                           true_region,
                           target_var_set,
                           exposure_dim,
                           cv,
                           region,
                           min_max) {

  sim_results <- CVtreeMLE(data = data,
                           w = covars,
                           a = exposures,
                           y = outcome,
                           n_folds = cv,
                           num_cores = 20,
                           family = "continuous",
                           parallel = TRUE,
                           parallel_cv = TRUE,
                           seed = seed,
                           region = region,
                           min_max = "max")

  tmle_pooled_mixture_results <- sim_results$`Pooled TMLE Mixture Results`[which.min(sim_results$`Pooled TMLE Mixture Results`$`Mixture ATE`), ]
  tmle_pooled_mixture_results <- na.omit(tmle_pooled_mixture_results)
  rule_mean_results <- sim_results$`Mixture Rules`

  ## tmle pooled results ------------------------
  # tmle_pooled_mixture_results <- sim_results$`Pooled TMLE Mixture Results`[max_ate_index,]
  tmle_pooled_lower <- tmle_pooled_mixture_results$`Lower CI`
  tmle_pooled_upper <- tmle_pooled_mixture_results$`Upper CI`
  tmle_pooled_ate <- tmle_pooled_mixture_results$`Mixture ATE`
  tmle_pooled_se <- tmle_pooled_mixture_results$`Standard Error`
  tmle_pooled_mix_decisions <- tmle_pooled_mixture_results$Average_Rule
  tmle_pooled_CI_range <- tmle_pooled_upper - tmle_pooled_lower

  ## tmle v-specific results ------------------------
  tmle_v_fold_mixture_results <- sim_results$`V-Specific Mix Results`
  tmle_v_fold_mixture_results <- na.omit(tmle_v_fold_mixture_results)
  # tmle_v_fold_mixture_results <- tmle_v_mixture_results[tmle_v_mixture_results$fold != "Pooled",]
  v_spec_mean_ate <- mean(tmle_v_fold_mixture_results$are)
  v_spec_mean_se <- mean(tmle_v_fold_mixture_results$se)
  v_spec_mean_lower <- mean(tmle_v_fold_mixture_results$lower_ci)
  v_spec_mean_upper <- mean(tmle_v_fold_mixture_results$upper_ci)
  v_spec_CI_range <- v_spec_mean_upper - v_spec_mean_lower


  ate_results <- list("tmle_pooled_ate" = tmle_pooled_ate,
                      "v_spec_mean_ate" = v_spec_mean_ate)

  se_results <- list("tmle_pooled_se" = tmle_pooled_se,
                     "v_spec_mean_se" = v_spec_mean_se)

  lower_ci_results <-  list("tmle_pooled_lower" = tmle_pooled_lower,
                            "v_spec_mean_lower" = v_spec_mean_lower)

  upper_ci_results <-  list("tmle_pooled_upper" = tmle_pooled_upper,
                            "v_spec_mean_upper" = v_spec_mean_upper)

  ## calc DA empirical truth ------------------------

  tmle_pooled_mix_decisions <- gsub("\\(.*?\\)", "", tmle_pooled_mix_decisions)


  da_p0_truth <- calc_empir_truth(data = P_0_data, rule = tmle_pooled_mix_decisions, exposure_dim = exposure_dim)


  tmle_pooled_da_bias <- da_p0_truth$pie - tmle_pooled_ate
  tmle_pooled_gt_bias <- true_ate - tmle_pooled_ate

  tmle_v_rule_spec_ates <- do.call(rbind,
                              lapply(X = tmle_v_fold_mixture_results$mix_rule,
                                     calc_empir_truth, data = P_0_data, exposure_dim = exposure_dim))
  tmle_v_rule_spec_ates <- as.data.frame(tmle_v_rule_spec_ates)

  v_spec_da_mean_bias <- mean(tmle_v_fold_mixture_results$are -
                             unlist(tmle_v_rule_spec_ates$pie))

  v_spec_gt_mean_bias <- mean(tmle_v_fold_mixture_results$are -
                               true_ate)

  bias_results <- list("tmle_pooled_da_bias" = tmle_pooled_da_bias,
                       "tmle_pooled_gt_bias" = tmle_pooled_gt_bias,
                       "v_spec_da_mean_bias" = v_spec_da_mean_bias,
                       "v_spec_gt_mean_bias" = v_spec_gt_mean_bias)

  ## calc DA empirical truth ------------------------

  rule_binary <- data %>%
    dplyr::transmute(A = ifelse(eval(parse(text = tmle_pooled_mix_decisions)), 1, 0))

  rule_applied_P_0 <- P_0_data %>%
    dplyr::mutate(A = ifelse(eval(parse(text = tmle_pooled_mix_decisions)), 1, 0))

  true_rule_binary <- data %>%
    dplyr::mutate(A = ifelse(eval(parse(text = true_rule)), 1, 0))


  confusion_table <- table(true_rule_binary$A, rule_binary$A)

  true_pos <- confusion_table[2,2] / sum(true_rule_binary$A)
  true_neg <- confusion_table[1,1] / (length(true_rule_binary$A) -
                                        sum(true_rule_binary$A))
  false_pos <- confusion_table[1,2] / (length(true_rule_binary$A) -
                                         sum(true_rule_binary$A))
  false_neg <- confusion_table[2,1] / sum(true_rule_binary$A)

  conf_table_results <- list("True Pos" = true_pos,
                             "True Neg" = true_neg,
                             "False Pos" = false_pos,
                             "False Neg" = false_neg)

  ## Coverage results --------------------------

  tmle_pooled_gt_coverage = ifelse(
    (tmle_pooled_lower <= true_ate & true_ate <= tmle_pooled_upper), 1,0
  )

  tmle_pooled_da_coverage = ifelse(
    (tmle_pooled_lower <= da_p0_truth$pie & da_p0_truth$pie <= tmle_pooled_upper), 1,0
  )

  v_spec_mean_da_cov <- mean(ifelse((tmle_v_fold_mixture_results$lower_ci
               <= tmle_v_rule_spec_ates$pie & tmle_v_rule_spec_ates$pie
               <= tmle_v_fold_mixture_results$upper_ci), 1, 0))

  v_spec_mean_gt_cov <- mean(ifelse((tmle_v_fold_mixture_results$lower_ci
                                     <= true_ate & true_ate
                                     <= tmle_v_fold_mixture_results$upper_ci), 1, 0))


  coverage_results <- list("tmle_pooled_gt_coverage" = tmle_pooled_gt_coverage,
                           "tmle_pooled_da_coverage" = tmle_pooled_da_coverage,
                           "v_spec_mean_da_cov" = v_spec_mean_da_cov,
                           "v_spec_mean_gt_cov" = v_spec_mean_gt_cov
                           )

  rule_min_bias <- mean(true_region_mean - rule_mean_results$coefficient)



sim_out <- c(ate_results,
             se_results,
             lower_ci_results,
             upper_ci_results,
             bias_results,
             conf_table_results,
             coverage_results,
             "true_ate" = true_ate,
             "true_pooled_da_ate" = da_p0_truth,
             "true_fold_da_ate" = mean(unlist(tmle_v_rule_spec_ates$pie)),
             "da rule" = tmle_pooled_mix_decisions,
             "rule_min_bias" = rule_min_bias,
             "pooled_tmle_CI_range" = tmle_pooled_CI_range,
             "v_spec_tmle_CI_range" = v_spec_CI_range)

  return(sim_out)
}
