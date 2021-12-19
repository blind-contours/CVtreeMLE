find_common_mixture_rules <- function(group_list,
                                      data = data,
                                      mix_comps = mix_comps,
                                      mixture_results,
                                      n_folds) {
  if (anyNA(group_list) == FALSE) {
    ## Get final mixture rule across folds
    all_mixt_vals_contained_vector <- list()
    all_mixt_vals_contained_summary <- list()

    mixture_rules <- list()
    fractions <- list()
    mixture_data <- subset(data, select = mix_comps)

    all_rules <- list()
    total_rules <- list()

    for (i in seq(group_list)) {
      directions_in_set <- list()
      combined_rules <- list()

      group <- group_list[[i]]
      group <- as.data.frame(group)

      vars <- mix_comps[mix_comps %in% strsplit(group$description[1], split = " ")[[1]]]

      rule_list <- group$description

      # combined_rule <- paste(unlist(combined_rule), collapse = " & ")
      fold_rules_eval <- list()

      for (k in seq(rule_list)) {
        fold_rule <- rule_list[[k]]

        fold_rule <- mixture_data %>%
          transmute("fold_rule" := ifelse(eval(parse(text = fold_rule)), 1, 0))

        fold_rules_eval[[k]] <- fold_rule
      }

      fold_rules_df <- do.call(cbind, fold_rules_eval)
      colnames(fold_rules_df) <- paste("fold_", seq(n_folds))

      fold_rules_df$sum <- rowSums(fold_rules_df)

      fold_rules_df$all_folds <- as.numeric(fold_rules_df$sum == n_folds)
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

          augmented_rule <- paste(var, ">", round(var_min, 3), "&", var, "<", round(var_max, 3))

          new_rule <- append(new_rule, augmented_rule)
        }

        rule <- paste(unlist(new_rule), collapse = " & ")

        mixture_rules[i] <- rule
        fractions[i] <- fraction
      } else {
        mixture_rules[i] <- "No Obs Overlap across all folds"
        fractions[i] <- NA
        total_rules[i] <- NA
      }
    }


    mixture_results <-
      cbind(mixture_results, unlist(mixture_rules))

    colnames(mixture_results)[8] <- "Mixture Interaction Rules"

    mixture_results <-
      cbind(mixture_results, unlist(fractions))

    colnames(mixture_results)[9] <- "Fraction Covered"
  } else {
    mixture_results <-
      cbind(mixture_results, NA)

    colnames(mixture_results)[length(mixture_results)] <- "Mixture Interaction Rules"

    mixture_results <-
      cbind(mixture_results, NA)

    colnames(mixture_results)[length(mixture_results)] <- "Fraction Covered"
  }

  return(mixture_results)
}
