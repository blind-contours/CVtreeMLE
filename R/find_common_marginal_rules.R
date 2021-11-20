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
#' @param fold_directions Direction of the rule, positive or negative - used to align the rules to create a common marginal rule
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
find_common_marginal_rules <- function(fold_rules, data, mix_comps, marginal_results, fold_directions, n_folds) {

  fold_rules <- unlist(fold_rules, recursive = FALSE)
  fold_rules <- fold_rules[!sapply(fold_rules, is.null)]

  fold_directions <- unlist(fold_directions, recursive = FALSE)
  fold_directions <- fold_directions[!sapply(fold_directions, is.null)]


  ## Get final marginal rule for mixture 1 across folds
  marg_rules <- list()
  fractions <- list()
  total_rules <- list()
  mins <- list()
  maxs <- list()

  for (i in seq(mix_comps)) {
    var <- mix_comps[i]
    mixture_data <- subset(data, select = mix_comps)
    rules <- sapply(fold_rules, "[", i)
    directions <- sapply(fold_directions, "[", i)

    if (any(unlist(rules) == "1") == FALSE) {
      combined_rule <- list()
      for (j in seq(rules)) {
        direction <- directions[[j]]
        rule <- rules[j]
        # rule <- gsub("-", " -", rule)

        if (direction == "negative") {
          rule_list <- c()
          neg_rules_flipped <- c()
          neg_rule <- mixture_data %>%
            transmute("pos_rule" := ifelse(eval(parse(text = rule)), 1, 0))
          pos_rule <- 1 - neg_rule
          target_neg_var_df <-
            cbind(subset(mixture_data, select = var), pos_rule)

          var_min_1 <-
            target_neg_var_df %>%
            dplyr::group_by(pos_rule) %>%
            dplyr::summarise(min = min(!!(as.name(var))))
          var_min_1 <- subset(var_min_1, pos_rule == 1, select = min)

          var_min_0 <-
            target_neg_var_df %>%
            dplyr::group_by(pos_rule) %>%
            dplyr::summarise(min = min(!!(as.name(var))))
          var_min_0 <- subset(var_min_0, pos_rule == 0, select = min)

          var_max_1 <-
            target_neg_var_df %>%
            dplyr::group_by(pos_rule) %>%
            dplyr::summarise(max = max(!!(as.name(var))))
          var_max_1 <- subset(var_max_1, pos_rule == 1, select = max)

          var_max_0 <-
            target_neg_var_df %>%
            dplyr::group_by(pos_rule) %>%
            dplyr::summarise(max = max(!!(as.name(var))))
          var_max_0 <- subset(var_max_0, pos_rule == 0, select = max)

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
        }

        combined_rule <- append(combined_rule, rule)
      }

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
      colnames(fold_rules_df) <- paste("fold_", seq(n_folds))

      fold_rules_df$sum <- rowSums(fold_rules_df)

      fold_rules_df$all_folds <- as.numeric(fold_rules_df$sum == n_folds)
      fold_rules_df <- cbind(fold_rules_df, mixture_data)

      if (dim(table(fold_rules_df$all_folds)) == 2) {
        count <- sum(as.numeric(fold_rules_df$sum == n_folds))
        total_count <- sum(as.numeric(fold_rules_df$sum > 0))
        fraction <- count / total_count

        var_min_1 <-
          fold_rules_df %>%
          group_by(all_folds) %>%
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
        mins[i] <- min(mixture_data[var])
        maxs[i] <- max(mixture_data[var])
      } else {
        marg_rules[i] <- "No Overlapping"
        fractions[i] <- NA
        total_rules[i] <- NA
        mins[i] <- min(mixture_data[var])
        maxs[i] <- max(mixture_data[var])
      }
    } else {
      marg_rules[i] <- NA
      fractions[i] <- NA
      total_rules[i] <- NA
      mins[i] <- min(mixture_data[var])
      maxs[i] <- max(mixture_data[var])
    }
  }

  marginal_results <-
    cbind(marginal_results, unlist(marg_rules))
  marginal_results <-
    cbind(marginal_results, unlist(fractions))
  marginal_results <-
    cbind(marginal_results, unlist(mins))
  marginal_results <-
    cbind(marginal_results, unlist(maxs))



  colnames(marginal_results)[7] <- "Marginal Rules"
  colnames(marginal_results)[8] <- "Fraction Overlap"
  colnames(marginal_results)[9] <- "Min"
  colnames(marginal_results)[10] <- "Max"

  return(marginal_results)
}
