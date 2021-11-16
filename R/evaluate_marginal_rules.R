evaluate_marginal_rules <- function(data, marg_decisions, mix_comps, marg_directions){

  marg_additive_data <- list()


  for (i in seq(marg_decisions)) {
    target_marg_rule <- marg_decisions[i][[1]]
    if (target_marg_rule != "1") {
      rule_name <- paste(mix_comps[i],"marg_rule",sep = "_")

      vec <- data %>%
        mutate(!!(rule_name) := ifelse(eval(parse(text = target_marg_rule)), 1, 0)) %>% dplyr::select(!!rule_name)

      if (marg_directions[[i]] == "negative") {
        vec <- 1- vec
      }

      marg_additive_data[[i]] <- vec
    }else{
      marg_additive_data[[i]] <- NA
    }
  }

  marg_rule_df <- do.call(cbind, marg_additive_data)
  marg_sums <- rowSums(marg_rule_df, na.rm = TRUE)

  data$sum_marg_hits <- marg_sums

  return(list(data = data, marginals = marg_rule_df))
}
