#' Evaluate mixture rules found during the rpart decision tree process
#'
#' @param data Input data
#' @param marg_decisions List of rules determined for each mixture component
#' @param mix_comps Vector of characters denoting mixture components

#' @importFrom dplyr transmute
#' @return Rules object. TODO: add more detail here.
#' @importFrom rlang :=
#'
#' @export

evaluate_marginal_rules <- function(data, marg_decisions, mix_comps) {
  marg_additive_data <- list()

  for (i in seq(dim(marg_decisions)[1])) {
    target_marg_rule <- marg_decisions[marg_decisions$target_m == mix_comps[i],]$rules
    if (target_marg_rule != "No Rules Found") {
      rule_name <- base::paste(mix_comps[i], "marg_rule", sep = "_")

      vec <- data %>%
        dplyr::mutate(!!(rule_name) := ifelse(eval(parse(text = target_marg_rule)), 1, 0)) %>%
        dplyr::select(!!rule_name)

      if (marg_decisions[marg_decisions$target_m == mix_comps[i],]$directions == "negative") {
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
