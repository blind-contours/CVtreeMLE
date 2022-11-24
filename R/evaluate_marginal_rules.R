#' Evaluate mixture rules found during the rpart decision tree process
#'
#' @param data Input data
#' @param marg_decisions List of rules determined for each mixture component
#' @param mix_comps Vector of characters denoting mixture components
#' @param no_marg_rules TRUE/FALSE if no marginal rules were found across the
#' folds

#' @importFrom dplyr transmute
#' @importFrom rlang :=
#' @return A list with marginal rules evaluated which includes:
#' \itemize{
#'   \item \code{data}: The fold raw data with the sum of evaluated rule
#'   binary vectors for each mixture component added.
#'   \item \code{marg_rule_df}: A binary matrix where each marginal rule has
#'   been evaluated to a binary vector.
#' }
#'
#' @export

evaluate_marginal_rules <- function(data,
                                    marg_decisions,
                                    no_marg_rules,
                                    mix_comps) {
  marg_additive_data <- list()

  if (no_marg_rules == FALSE) {
    for (i in seq(nrow(marg_decisions))) {
      target_marg_rule <- marg_decisions$rules[i]
      rule_name <- base::paste(marg_decisions$target_m[i], "marg_rule",
        marg_decisions$quantile[i],
        sep = "_"
      )

      vec <- data %>%
        dplyr::mutate(!!(rule_name) :=
          ifelse(eval(parse(text = target_marg_rule)), 1, 0)) %>%
        dplyr::select(!!rule_name)

      marg_additive_data[[i]] <- vec
    }

    marg_rule_df <- do.call(cbind, marg_additive_data)
    marg_sums <- base::rowSums(marg_rule_df, na.rm = TRUE)

    data$sum_marg_hits <- marg_sums
  } else {
    data$sum_marg_hits <- NA
    marg_rule_df <- NA
  }

  return(list(data = data, marginals = marg_rule_df))
}
