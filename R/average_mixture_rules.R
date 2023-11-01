#' @title Estimate the average rule. This is the rule that is the average
#' across thresholds in folds and gives the lower and upper bounds
#'
#' @description This function takes in a list of rules that are grouped by
#' variable sets. These rules for different variable sets may be slightly
#' different across the folds so we make a average rule for each variable set.
#' This entails creating a new rule that is is the average of cutpoints:
#'
#' @param group_list List of dataframes grouped by rules for variable sets
#' @param data Full data
#' @param mix_comps Mixture components of A
#' @param n_folds Number of folds used in cross-validation
#' @param no_mixture_rules TRUE/FALSE if no mixture rule was found
#' @param mixture_results data frame of results found for mixture rules
#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows

#' @return Rules object. TODO: add more detail here.
#' @importFrom stats as.formula glm p.adjust
#' @importFrom stats plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#' @return Data frame with mixture results including the union rule and
#' proportion found across the folds
#'
#' @export

average_mixture_rules <- function(group_list,
                                  data = data,
                                  mix_comps = mix_comps,
                                  n_folds,
                                  mixture_results) {

  average_rules <- list()
  fold_proportions <- list()

  for (i in seq_along(group_list)) {

    group <- as.data.frame(group_list[[i]])
    proportion_in_fold <- nrow(group) / n_folds
    fold_proportions[i] <- proportion_in_fold

    # If group contains only one rule, add it to the list and skip further processing
    if (nrow(group) == 1) {
      average_rules[[i]] <- group$description[1]
      next
    }

    vars <- mix_comps[mix_comps %in%
                        strsplit(group$description[1], split = " ")[[1]]]

    avg_rule_list <- list()

    for (var in vars) {
      all_min_values <- numeric()
      all_max_values <- numeric()

      for (rule in group$description) {

        pattern_for_min <- paste0(var, "\\s*>\\s*(\\d+\\.?\\d*)")
        pattern_for_max <- paste0(var, "\\s*<=\\s*(\\d+\\.?\\d*)")

        match_min <- regexec(pattern_for_min, rule)
        match_max <- regexec(pattern_for_max, rule)

        if (match_min[[1]][1] != -1) {
          min_value <- as.numeric(regmatches(rule, match_min)[[1]][2])
          all_min_values <- c(all_min_values, min_value)
        }

        if (match_max[[1]][1] != -1) {
          max_value <- as.numeric(regmatches(rule, match_max)[[1]][2])
          all_max_values <- c(all_max_values, max_value)
        }
      }

      if (is_empty(all_min_values)) {
        mean <- round(mean(all_max_values, na.rm = TRUE), 3)
        min <- round(min(all_max_values, na.rm = TRUE), 3)
        max <- round(max(all_max_values, na.rm = TRUE), 3)

        avg_rule_for_var <- paste0(var, " <=", mean, "(", min, ",", max, ")")

      }

      if (is_empty(all_max_values)) {
        mean <- round(mean(all_min_values, na.rm = TRUE), 3)
        min <- round(min(all_min_values, na.rm = TRUE), 3)
        max <- round(max(all_min_values, na.rm = TRUE), 3)

        avg_rule_for_var <- paste0(var, " >=", mean, "(", min, ",", max, ")")
      }


      if (avg_rule_for_var != "") {
        avg_rule_list <- c(avg_rule_list, avg_rule_for_var)
      }

    average_rules[[i]] <- paste(avg_rule_list, collapse = " & ")
  }

  }

  mixture_results$Average_Rule <- unlist(average_rules)
  mixture_results$Proportion_Folds <- unlist(fold_proportions)

  return(mixture_results)
}





