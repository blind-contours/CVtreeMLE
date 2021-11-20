#' Evaluate mixture rules found during the rpart decision tree process
#'
#' @param data Input data
#' @param marg_decisions List of rules determined for each mixture component
#' @param mix_comps Vector of characters denoting mixture components
#' @param marg_directions Directionality of the marginal mixture rule

#' @importFrom dplyr transmute
#' @return Rules object. TODO: add more detail here.
#' @importFrom rlang :=
#'
#' @export

evaluate_marginal_rules <- function(data, marg_decisions, mix_comps, marg_directions) {
  marg_additive_data <- list()

  for (i in seq(marg_decisions)) {
    target_marg_rule <- marg_decisions[i][[1]]
    if (target_marg_rule != "1") {
      rule_name <- base::paste(mix_comps[i], "marg_rule", sep = "_")

      vec <- data %>%
        dplyr::mutate(!!(rule_name) := ifelse(eval(parse(text = target_marg_rule)), 1, 0)) %>%
        dplyr::select(!!rule_name)

      if (marg_directions[[i]] == "negative") {
        vec <- 1 - vec
      }

      marg_additive_data[[i]] <- vec
    } else {
      marg_additive_data[[i]] <- NA
    }
  }

  marg_rule_df <- do.call(cbind, marg_additive_data)
  marg_sums <- base::rowSums(marg_rule_df, na.rm = TRUE)

  data$sum_marg_hits <- marg_sums

  return(list(data = data, marginals = marg_rule_df))
}
