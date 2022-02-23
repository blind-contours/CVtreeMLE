
fit_estimators <- function(data, true_rule) {
  lrnr_glm <- Lrnr_glm$new()
  ranger_lrnr <- Lrnr_ranger$new()
  lrnr_gam <- Lrnr_gam$new()

  # put all the learners together (this is just one way to do it)
  learners <- c(lrnr_glm, ranger_lrnr, lrnr_gam)

  Q1_stack <- make_learner(Stack, learners)

  lrnr_glmtree_001 <- Lrnr_glmtree$new(alpha = 0.5, maxdepth = 3)
  lrnr_glmtree_002 <- Lrnr_glmtree$new(alpha = 0.6,  maxdepth = 4)
  lrnr_glmtree_003 <- Lrnr_glmtree$new(alpha = 0.7, maxdepth = 2)
  lrnr_glmtree_004 <- Lrnr_glmtree$new(alpha = 0.8, maxdepth = 1)

  learners <- c( lrnr_glmtree_001, lrnr_glmtree_002, lrnr_glmtree_003, lrnr_glmtree_004)
  discrete_sl_metalrn <- Lrnr_cv_selector$new()

  tree_stack <- make_learner(Stack, learners)

  discrete_tree_sl <- Lrnr_sl$new(
    learners = tree_stack,
    metalearner = discrete_sl_metalrn
  )


  sim_results2 <- CVtreeMLE(data = data,
                           W = c("W", "W2"),
                           Y = "y",
                           A = c(paste("M", seq(3), sep = "")),
                           back_iter_SL = Q1_stack,
                           tree_SL = discrete_tree_sl,
                           n_folds = 5,
                           num_cores = 8,
                           family = "gaussian")

  max_ate_index <- which.max(abs(sim_results$`Pooled TMLE Mixture Results`$`Mixture ATE`))

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
    mix_decisions <- mixture_results$`Mixture Interaction Rules`

    rule_binary <- data %>%
      dplyr::transmute(A = ifelse(eval(parse(text = mix_decisions)), 1, 0))

    true_rule_binary <- data %>%
      dplyr::transmute(A = ifelse(eval(parse(text = true_rule)), 1, 0))

    confusion_table <- table(true_rule_binary$A, rule_binary$A)

    true_pos <- confusion_table[2,2] / sum(true_rule_binary$A)
    true_neg <- confusion_table[1,1] / (length(true_rule_binary$A) - sum(true_rule_binary$A))
    false_pos <- confusion_table[1,2] / (length(true_rule_binary$A) - sum(true_rule_binary$A))
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

  sim_out <- list("Mix ind" = mixure_found_ind, "ATE" = ate, "Lower" = lower, "Upper" = upper,
                  "True Pos" = true_pos, "True Neg" = true_neg,
                  "False Pos" = false_pos, "False Neg" = false_neg)
  return(sim_out)
}
