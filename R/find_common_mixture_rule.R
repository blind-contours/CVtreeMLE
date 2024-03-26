#' @title Estimate the union rule. This is the rule that covers all observations
#' across the folds.
#'
#' @description This function takes in a list of rules that are grouped by
#' variable sets. These rules for different variable sets may be slightly
#' different across the folds so we make a union rule for each variable set.
#' This entails creating a new rule that is essentially:
#' (rule fold 1) OR (rule fold 2) OR (rule fold 3) etc. We then evaluate this
#' rule on the input data to create a binary indicator of the union rule.
#' Then, for variables used in the rule, find the min in and max in regions
#' indicated by the rule for each variable. We then put these together with
#' AND statements to create the union rule.
#'
#' @param group_list List of dataframes grouped by rules for variable sets
#' @param data Full data
#' @param mix_comps Mixture components of A
#' @param n_folds Number of folds used in cross-validation
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

common_mixture_rules <- function(group_list,
                                 data = data,
                                 mix_comps = mix_comps,
                                 mixture_results,
                                 n_folds) {
  mixture_all_rules <- list()
  mixture_any_rules <- list()

  fractions <- list()
  fold_proportions <- list()
  mixture_data <- subset(data, select = mix_comps)

  total_rules <- list()

  for (i in seq(group_list)) {
    group <- group_list[[i]]
    group <- as.data.frame(group)

    vars <- mix_comps[mix_comps %in%
      strsplit(group$description[1], split = " ")[[1]]]

    intxn_rule <- paste(
      "(",
      paste(group$description, collapse = ")|("), ")"
    )

    intxn_data <- data %>%
      mutate("intxn_rule" := ifelse(eval(parse(text = intxn_rule)), 1, 0))

    new_rule <- list()

    for (var in vars) {
      var_min <-
        intxn_data %>%
        group_by(intxn_rule) %>%
        summarise(min = min(!!(as.name(var))))
      var_min <- subset(var_min, intxn_rule == 1, select = min)
      var_max <-
        intxn_data %>%
        group_by(intxn_rule) %>%
        summarise(max = max(!!(as.name(var))))
      var_max <- subset(var_max, intxn_rule == 1, select = max)

      augmented_rule <- paste(
        var, ">=", round(var_min, 3), "&", var,
        "<=", round(var_max, 3)
      )

      new_rule <- append(new_rule, augmented_rule)
    }

    interaction_rule <- paste(unlist(new_rule), collapse = " & ")

    rule_list <- group$description

    proportion_in_fold <- length(rule_list) / n_folds
    fold_rules_eval <- list()

    for (k in seq(rule_list)) {
      fold_rule <- rule_list[[k]]

      fold_rule <- mixture_data %>%
        transmute("fold_rule" := ifelse(eval(parse(text = fold_rule)), 1, 0))

      fold_rules_eval[[k]] <- fold_rule
    }

    fold_rules_df <- do.call(cbind, fold_rules_eval)
    colnames(fold_rules_df) <- paste("fold_", seq(nrow(group)))

    fold_rules_df$sum <- rowSums(fold_rules_df)

    fold_rules_df$all_folds <- as.numeric(fold_rules_df$sum == nrow(group))
    fold_rules_df <- cbind(fold_rules_df, mixture_data)

    total_count <- table(fold_rules_df$sum > 0)[[2]]

    if (dim(table(fold_rules_df$all_folds)) == 2) {
      count <- table(fold_rules_df$all_folds > 0)[[2]]
      fraction <- count / total_count
      new_rule <- list()

      for (var in vars) {
        var_min <-
          fold_rules_df %>%
          group_by(all_folds) %>%
          summarise(min = min(!!(as.name(var))))
        var_min <- subset(var_min, all_folds == 1, select = min)
        var_max <-
          fold_rules_df %>%
          group_by(all_folds) %>%
          summarise(max = max(!!(as.name(var))))
        var_max <- subset(var_max, all_folds == 1, select = max)

        augmented_rule <- paste(
          var, ">=", round(var_min, 3), "&",
          var, "<=", round(var_max, 3)
        )

        new_rule <- append(new_rule, augmented_rule)
      }

      rule <- paste(unlist(new_rule), collapse = " & ")

      mixture_any_rules[i] <- interaction_rule
      mixture_all_rules[i] <- rule

      fractions[i] <- fraction
    } else {
      mixture_all_rules[i] <- "No Obs Overlap across all folds"
      mixture_any_rules[i] <- interaction_rule
      fractions[i] <- NA
      total_rules[i] <- NA
    }

    fold_proportions[i] <- proportion_in_fold
  }


  mixture_results$Union_Rule <- unlist(mixture_any_rules)
  mixture_results$Proportion_Folds <- unlist(fold_proportions)

  return(mixture_results)
}
