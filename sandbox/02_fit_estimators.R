
fit_estimators <- function(data, true_rule) {

  sim_results <- CVtreeMLE(data = data,
                           w = c("w", "w2"),
                           a = c(paste("M", seq(3), sep = "")),
                           y = "y",
                           n_folds = 5,
                           num_cores = 20,
                           family = "gaussian",
                           parallel = TRUE,
                           parallel_cv = FALSE)

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

    rule_binary <- data %>%
      dplyr::transmute(A = ifelse(eval(parse(text = mix_decisions)), 1, 0))

    true_rule_binary <- data %>%
      dplyr::transmute(A = ifelse(eval(parse(text = true_rule)), 1, 0))

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
  }

  sim_out <- list("Mix ind" = mixure_found_ind, "ATE" = ate,
                  "Lower" = lower, "Upper" = upper,
                  "True Pos" = true_pos, "True Neg" = true_neg,
                  "False Pos" = false_pos, "False Neg" = false_neg)
  return(sim_out)
}
