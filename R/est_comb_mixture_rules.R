est_comb_mixture_rules <- function(At, Av, W, Y, rules, no_rules, SL.library) {

  At_c <- At
  Av_c <- Av

  mix_interaction_train <- list()
  mix_interaction_valid <- list()

  rules <- as.data.frame(rules)

  if (dim(rules)[1] > 0 & no_rules == FALSE) {
    for (i in seq(dim(rules)[1])) {
      target_interaction_result <- rules[i, ]$description

      rule_name <- paste("mix_interaction", i, sep = "_")
      target_interaction_result <-
        gsub("-", " -", target_interaction_result)

      vec_train <- At_c %>%
        dplyr::mutate(!!(rule_name) := ifelse(eval(
          parse(text = target_interaction_result)
        ), 1, 0)) %>%
        dplyr::select(!!rule_name)

      vec_valid <- Av_c %>%
        dplyr::mutate(!!(rule_name) := ifelse(eval(
          parse(text = target_interaction_result)
        ), 1, 0)) %>%
        dplyr::select(!!rule_name)

      mix_interaction_train[[i]] <- vec_train
      mix_interaction_valid[[i]] <- vec_valid
    }

    mix_interactions_rules_train <-
      do.call(cbind, mix_interaction_train)

    At_mix_comb <-
      cbind(
        mix_interactions_rules_train,
        At_c[W]
      )

    mix_interactions_rules_valid <-
      do.call(cbind, mix_interaction_valid)

    Av_mix_comb <-
      cbind(
        mix_interactions_rules_valid,
        Av_c[W]
      )

    QbarAWSL_m <- SuperLearner(
      Y = At_c$y_scaled,
      X = At_mix_comb,
      SL.library = SL.library,
      family = "gaussian",
      verbose = FALSE
    )

    QbarAW <- bound_precision(predict(QbarAWSL_m, newdata = Av_mix_comb)$pred)

    Av_mix_comb$QbarAW_combo <- QbarAW
    Av_mix_comb$y_scaled <- Av_c$y_scaled
    Av_mix_comb$raw_outcome <- Av[, Y]
  } else {
    Av_mix_comb <- NA
    QbarAWSL_m <- NA
  }

  return(list("data" = Av_mix_comb, "learner" = QbarAWSL_m))
}
