#' @title Create a new rule based on observations that meet every rule across the folds for each mixture.
#' @description For each mixture component, a different rule could be found for each fold. Therefore, it is necessary to create one rule for each mixture component that can be interpreted
#' as a common rule across the folds. To do this, observations that meet all rules for all folds are determined. Then a new rule is created for these observations. A coverage metric is
#' calculated which is the fraction of observations in this common rule compared to the sum of observations that met each rule. The coverage metric can be thought of as an intersection over
#' union or the Jaccard coefficient for the rules.
#'
#' @param fold_rules List of rules found for each mixture component found across the folds
#' @param data Full data which rules are evaluated
#' @param mix_comps Vector of mixture components
#' @param marginal_results Dataframe holding the results for each marginal component rule
#' @param n_folds Total number of folds

#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows

#' @return Rules object. TODO: add more detail here.
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#' @importFrom dplyr summarise group_by

#'
#' @export
#'
find_common_marginal_rules <- function(fold_rules, data, mix_comps, marginal_results, n_folds) {
  fold_rules <- unlist(fold_rules, recursive = FALSE)
  fold_rules <- fold_rules[!sapply(fold_rules, is.null)]

  ## Get final marginal rule for mixture 1 across folds
  marg_rules <- list()
  fractions <- list()
  total_rules <- list()
  mins <- list()
  maxs <- list()

  for (i in seq(nrow(fold_rules[[1]]))) {
    var <- fold_rules[[1]]$target_m[i]
    var_quant <- fold_rules[[1]]$var_quant_group[i]
    mixture_data <- subset(data, select = mix_comps)

    rules <- list()
    for (rule in seq(length(fold_rules))) {
      temp <- fold_rules[[rule]]
      temp <- temp %>% filter(.data$var_quant_group == var_quant)
      rules[[rule]] <- temp
    }

    rules <- do.call(rbind, rules)

    rules <- rules$rules

    combined_rule <- rules

    total_rules[[i]] <- combined_rule

    # combined_rule <- paste(unlist(combined_rule), collapse = " & ")
    fold_rules_eval <- list()

    for (k in seq(combined_rule)) {
      fold_rule <- combined_rule[[k]]

      fold_rule <- mixture_data %>%
        transmute("fold_rule" := ifelse(eval(parse(text = fold_rule)), 1, 0))

      fold_rules_eval[[k]] <- fold_rule
    }

    fold_rules_df <- do.call(cbind, fold_rules_eval)
    colnames(fold_rules_df) <- paste("fold_", seq(combined_rule))

    fold_rules_df$sum <- rowSums(fold_rules_df)

    fold_rules_df$all_folds <- as.numeric(fold_rules_df$sum == length(combined_rule))
    fold_rules_df <- cbind(fold_rules_df, mixture_data)

    if (dim(table(fold_rules_df$all_folds)) == 2) {
      count <- sum(as.numeric(fold_rules_df$sum == length(combined_rule)))
      total_count <- sum(as.numeric(fold_rules_df$sum > 0))
      fraction <- round(count / total_count, 3)

      var_min_1 <-
        fold_rules_df %>%
        dplyr::group_by(all_folds) %>%
        summarise(min = min(!!(as.name(var))))
      var_min_1 <- subset(var_min_1, all_folds == 1, select = min)

      var_min_0 <-
        fold_rules_df %>%
        group_by(all_folds) %>%
        summarise(min = min(!!(as.name(var))))
      var_min_0 <- subset(var_min_0, all_folds == 0, select = min)

      var_max_1 <-
        fold_rules_df %>%
        group_by(all_folds) %>%
        summarise(max = max(!!(as.name(var))))
      var_max_1 <- subset(var_max_1, all_folds == 1, select = max)

      var_max_0 <-
        fold_rules_df %>%
        group_by(all_folds) %>%
        summarise(max = max(!!(as.name(var))))
      var_max_0 <- subset(var_max_0, all_folds == 0, select = max)

      rule <-
        paste(
          var,
          ">",
          round(var_min_1, 5),
          "&",
          var,
          "<",
          round(var_max_1, 5)
        )

      marg_rules[i] <- rule
      fractions[i] <- fraction
      mins[i] <- round(min(mixture_data[var]), 3)
      maxs[i] <- round(max(mixture_data[var]), 3)
    } else {
      marg_rules[i] <- "No Overlapping"
      fractions[i] <- NA
      total_rules[i] <- NA
      mins[i] <- round(min(mixture_data[var]), 3)
      maxs[i] <- round(max(mixture_data[var]), 3)
    }
  }

  total_marginal_results <- as.data.frame(cbind(unlist(marg_rules), unlist(fractions), unlist(mins), unlist(maxs)))
  # marginal_results <-
  #   cbind(marginal_results, unlist(fractions))
  # marginal_results <-
  #   cbind(marginal_results, unlist(mins))
  # marginal_results <-
  #   cbind(marginal_results, unlist(maxs))


  colnames(total_marginal_results) <- c("Marginal Rules", "Fraction Overlap", "Min", "Max")


  return(total_marginal_results)
}
