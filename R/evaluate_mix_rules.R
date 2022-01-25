#' Evaluate mixture rules found during the PRE process
#'
#' @param data Input data
#' @param rules Dataframe of the rules determined
#' @importFrom dplyr transmute
#' @return Rules object. TODO: add more detail here.

#'
#' @export

evaluate_mixture_rules <- function(data, rules) {
  rule_binary_match_list <- list()

  for (i in seq(dim(rules)[1])) {
    rule_t <- rules[i, ]
    if (rule_t$description != "1") {
      mix_decisions <- rule_t$description
      rule_binary <- data %>%
        dplyr::transmute(A = ifelse(eval(parse(text = mix_decisions)), 1, 0))

      rule_binary_match_list[i] <- rule_binary
    }
  }

  rule_binary_match_df <- base::do.call(cbind, rule_binary_match_list)

  colnames(rule_binary_match_df) <- rules$description

  binary_rule_matrix <- rule_binary_match_df

  return(binary_rule_matrix)
}
