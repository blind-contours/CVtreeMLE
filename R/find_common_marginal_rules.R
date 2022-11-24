#' @title Create a new rule based on observations that meet every rule across
#' the folds for each mixture.
#' @description For each mixture component, a different rule could be found for
#' each fold. Therefore, it is necessary to create one rule for each mixture
#' component that can be interpreted
#' as a common rule across the folds. To do this, observations that meet all
#' rules for all folds are determined. Then a new rule is created for these
#' observations. Specifically, we put OR statements between the rules
#' found across the folds then look at the min and max values in this new
#' region which encompasses all observations across the folds.
#'
#' @param fold_rules List of rules found for each mixture component found
#' across the folds
#' @param data Full data which rules are evaluated
#' @param mix_comps Vector of mixture components
#' @param marginal_results Data frame holding the results for each marginal
#' component rule
#' @param n_folds Total number of folds
#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis
#' @importFrom stats qnorm qunif rnorm runif
#' @importFrom rlang :=
#' @importFrom dplyr summarise group_by
#' @return Data frame with rules, threshold regions, proportion in folds and
#' min/max values
#'
#' @export
#'
find_common_marginal_rules <- function(fold_rules,
                                       data,
                                       mix_comps,
                                       marginal_results,
                                       n_folds) {
  fold_rules <- unlist(fold_rules, recursive = FALSE)
  fold_rules <- fold_rules[!sapply(fold_rules, is.null)]

  marg_rules <- list()
  total_rules <- list()
  mins <- list()
  maxs <- list()
  fold_proportions <- list()
  var_quantiles <- list()

  marginal_quantiles <- unique(do.call(rbind, fold_rules)$var_quant_group)

  for (i in seq(marginal_quantiles)) {
    var_quant <- marginal_quantiles[i]
    var_quantiles[i] <- var_quant
    mixture_data <- subset(data, select = mix_comps)

    rules <- list()
    for (rule in seq(length(fold_rules))) {
      temp <- fold_rules[[rule]]
      temp <- temp %>% filter(.data$var_quant_group == var_quant)
      rules[[rule]] <- temp
    }
    rules <- do.call(rbind, rules)
    var <- unique(rules$target_m)

    rules <- rules$rules

    proportion_in_fold <- length(rules) / n_folds
    fold_proportions[i] <- proportion_in_fold

    combined_rule <- rules

    total_rules[[i]] <- combined_rule

    marginal_rule <- paste("(", paste(rules, collapse = ")|("), ")")

    fold_rules_df <- data %>%
      mutate("all_folds" := ifelse(eval(parse(text = marginal_rule)), 1, 0))

    if (dim(table(fold_rules_df$all_folds)) == 2) {
      var_min_1 <-
        fold_rules_df %>%
        dplyr::group_by(all_folds) %>%
        summarise(min = min(!!(as.name(var))))
      var_min_1 <- subset(var_min_1, all_folds == 1, select = min)


      var_max_1 <-
        fold_rules_df %>%
        group_by(all_folds) %>%
        summarise(max = max(!!(as.name(var))))
      var_max_1 <- subset(var_max_1, all_folds == 1, select = max)

      rule <-
        paste(
          var,
          ">=",
          round(var_min_1, 5),
          "&",
          var,
          "<=",
          round(var_max_1, 5)
        )

      marg_rules[i] <- rule
      mins[i] <- round(min(mixture_data[var]), 3)
      maxs[i] <- round(max(mixture_data[var]), 3)
    } else {
      marg_rules[i] <- "No Overlapping"
      total_rules[i] <- NA
      mins[i] <- round(min(mixture_data[var]), 3)
      maxs[i] <- round(max(mixture_data[var]), 3)
    }
  }

  total_marginal_results <- as.data.frame(cbind(
    unlist(var_quantiles),
    unlist(marg_rules),
    unlist(mins),
    unlist(maxs),
    unlist(fold_proportions)
  ))

  colnames(total_marginal_results) <- c(
    "Variable Quantile",
    "Marginal Rules",
    "Min",
    "Max",
    "Proportion in Fold"
  )


  return(total_marginal_results)
}
