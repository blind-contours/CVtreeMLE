#' @title Compute v-fold specific estimates and do a meta-analysis type pooled average ATE for mixture rules.
#' @description For each fold, estimates the ATE for a fold specific mixture rule. Also estimates a meta-analysis type average ATE and pooled variance. Creates a union rule that covers all the folds in the rules.
#'
#' @param v_fold_mixture_results List of dataframes for v-fold specific estimates of the fold-specific rule results
#' @param mix_comps Vector of characters indicating the mixture components
#' @param n_folds Number of folds used in cross-validation
#' @param data Input dataframe on which to evaluate the rules

#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows

#' @return Rules object. TODO: add more detail here.
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export


compute_meta_mix_results <- function(v_fold_mixture_results, mix_comps, n_folds, data) {
  v_fold_mixture_group <- v_fold_mixture_group_split(v_fold_mixture_results)
  # v_fold_mixture_group <- v_fold_mixture_results %>% dplyr::group_by(Variables)
  # v_fold_mixture_group <- dplyr::group_split(v_fold_mixture_group)

  v_fold_mixture_w_pooled <- list()
  intxn_names_list <- list()

  for (i in seq(v_fold_mixture_group)) {
    results_df <- v_fold_mixture_group[[i]]

    weighted_mean <- sum(results_df$`Mixture ATE` * (1 / results_df$`Standard Error`^2)) / sum((1 / results_df$`Standard Error`^2))
    weighted_RMSE <- sum(results_df$RMSE * (1 / results_df$`Standard Error`^2)) / sum((1 / results_df$`Standard Error`^2))

    pooled_se <- sqrt(1 / (1 / sum(results_df$`Standard Error`^2)))

    pooled_P_val <- round(2 * stats::pnorm(abs(weighted_mean / pooled_se), lower.tail = F), 5)

    pooled_CI <- c(
      round(weighted_mean + stats::qnorm(0.05 / 2, lower.tail = T) * pooled_se, 4),
      round(weighted_mean + stats::qnorm(0.05 / 2, lower.tail = F) * pooled_se, 4)
    )

    vars <- mix_comps[mix_comps %in% strsplit(results_df$`Mixture Interaction Rules`[1], split = " ")[[1]]]

    intxn_rule <- paste("(", paste(results_df$`Mixture Interaction Rules`, collapse = ")|("), ")")

    intxn_data <- data %>%
      dplyr::mutate("intxn_rule" := ifelse(eval(parse(text = intxn_rule)), 1, 0))

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

      augmented_rule <- paste(var, ">", round(var_min, 3), "&", var, "<", round(var_max, 3))

      new_rule <- append(new_rule, augmented_rule)
    }

    interaction_rule <- paste(unlist(new_rule), collapse = " & ")

    # split_rules <- strsplit(results_df$`Mixture Interaction Rules`, split = "&")
    #
    # average_rules <- list()

    # for (j in seq(vars)) {
    #   var <- vars[j]
    #   sign_list <- list()
    #   num_list <- list()
    #   for (k in seq(split_rules)) {
    #     rule <- split_rules[[k]]
    #     var_rule <- rule[(str_detect(rule, var))]
    #     sign_list[k] <- stringr::str_extract(var_rule, '>|>=|<|<=')
    #     num_list[k] <- round(as.numeric(stringr::str_extract(var_rule, ' [0-9]+([.][0-9]+)?')),2)
    #   }
    #
    #   rule_mean <- mean(unlist(num_list))
    #   rule_min <- min(unlist(num_list))
    #   rule_max <- max(unlist(num_list))
    #
    #   rule_range <- paste(rule_mean, " (", rule_min, "-", rule_max, ")", sep = "")
    #
    #   if (length(unique(unlist(sign_list))) == 1) {
    #     rule_sign <- unique(unlist(sign_list))
    #   }else{
    #     rule_sign <- "?"
    #   }
    #
    #   total_rule <- paste(var, rule_sign, rule_range)
    #   average_rules[j] <- total_rule
    # }
    #
    # interaction_rule <- paste(average_rules, collapse = " & ")

    average_results <- cbind(
      round(weighted_mean, 3),
      round(pooled_se, 3),
      pooled_CI[1], pooled_CI[2],
      pooled_P_val, pooled_P_val, round(weighted_RMSE, 3),
      interaction_rule,
      unique(results_df$Variables)
    )

    colnames(average_results) <- colnames(results_df)

    results <- as.data.frame(rbind(results_df, average_results))
    rownames(results) <- c(seq(nrow(results) - 1 ), "Pooled")
    intxn_names_list[[i]] <- unique(results$Variables)

    v_fold_mixture_w_pooled[[i]] <- results
  }

  names(v_fold_mixture_w_pooled) <- intxn_names_list


  return(v_fold_mixture_w_pooled)
}
