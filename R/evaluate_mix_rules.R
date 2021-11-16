evaluate_mixture_rules <- function(data, rules){

  rule_binary_match_list <- list()
  
  for (i in seq(dim(rules)[1])) {
    rule_t <- rules[i,]
    if (rule_t$description != "1"){
      mix_decisions <- rule_t$description
      rule_binary <- data %>%
        transmute(A = ifelse(eval(parse(text = mix_decisions)), 1, 0))
      # if (rule_t$direction == "negative") {
      #   rule_binary <- 1-rule_binary
      # }

      rule_binary_match_list[i] <- rule_binary
    }
  }

  rule_binary_match_df <- do.call(cbind, rule_binary_match_list)
  #rule_binary_match_sums <- rowSums(rule_binary_match_df)
  #rule_binary_match_sums_binary <- ifelse(rule_binary_match_sums == max(rule_binary_match_sums), 1, 0)
  colnames(rule_binary_match_df) <- rules$description

  binary_rule_matrix <- rule_binary_match_df

  # colnames(binary_rule_matrix)[dim(binary_rule_matrix)[2]] <- "all interactions"

  return(binary_rule_matrix)
}
