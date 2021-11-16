find_common_marginal_rules <- function(fold_rules, data = data, mix_comps = mix_comps, marginal_results = marginal_results, fold_directions, n_folds) {
  
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
            group_by(pos_rule) %>%
            summarise(min = min(!!(as.name(var))))
          var_min_1 <- subset(var_min_1, pos_rule == 1, select = min)

          var_min_0 <-
            target_neg_var_df %>%
            group_by(pos_rule) %>%
            summarise(min = min(!!(as.name(var))))
          var_min_0 <- subset(var_min_0, pos_rule == 0, select = min)

          var_max_1 <-
            target_neg_var_df %>%
            group_by(pos_rule) %>%
            summarise(max = max(!!(as.name(var))))
          var_max_1 <- subset(var_max_1, pos_rule == 1, select = max)

          var_max_0 <-
            target_neg_var_df %>%
            group_by(pos_rule) %>%
            summarise(max = max(!!(as.name(var))))
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

          # rule0 <-
          #   paste(var,
          #         ">",
          #         round(var_max_0, 5),
          #         "&",
          #         var,
          #         "<",
          #         round(var_max_1, 5))

          # rule <- paste("(",rule1,")", "|", "(", rule0,")")
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

        # rule0 <-
        #   paste(var,
        #         ">",
        #         round(var_max_1, 5),
        #         "&",
        #         var,
        #         "<",
        #         round(var_max_0, 5))

        # rule <- paste("(",rule1,")", "|", "(", rule0,")")

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



#         rule_df <-
#           do.call(rbind, str_split(neg_rules_flipped, " "))
#         rule_df <- as.data.frame(rule_df)
#         indexes <- which(!is.na(as.numeric(rule_df[1,])) == TRUE)
#
#         rule_df[indexes] <- sapply(rule_df[indexes], as.numeric)
#
#         rule_collection <- list()
#         for (k in seq(dim(rule_df)[2])) {
#           if (is.numeric(rule_df[, k]) == TRUE) {
#             target <- round(mean(rule_df[, k]), 3)
#           } else{
#             target <- names(which.max(table(rule_df[, k])))
#           }
#
#           rule_collection[k] <- target
#         }
#
#         rule_mods[[j]] <- paste(rule_collection, collapse = " ")
#       }
#
#       if (direction == "positive") {
#         if (grepl("&", rule)) {
#           rule <- gsub("&", " & ", rule)
#           rule_mods[[j]] <- rule
#         } else{
#           if (grepl("<", rule)) {
#             limit <- min(mixture_data[, var])
#             if (grepl("=", rule)) {
#               rule <- gsub("<=", " <= ", rule)
#               rule <- paste(var, ">=", round(limit, 3), "&", rule)
#               rule_mods[[j]] <- rule
#             } else{
#               rule <- gsub("<", " < ", rule)
#               rule <- paste(var, ">=", round(limit, 3), "&", rule)
#               rule_mods[[j]] <- rule
#             }
#           } else{
#             limit <- max(mixture_data[, var])
#             if (grepl("=", rule)) {
#               rule <- gsub(">=", " >= ", rule)
#               rule <- paste(rule, "&", var, "<=", round(limit, 3))
#               rule_mods[[j]] <- rule
#             } else{
#               rule <- gsub(">", " > ", rule)
#               rule <- paste(rule, "&", var, "<=", round(limit, 3))
#               rule_mods[[j]] <- rule
#             }
#
#           }
#
#         }
#       }
#     }
#
#
#     split_rules <- strsplit(unlist(rule_mods), "\\s+")
#     split_rules <-
#       lapply(split_rules, function(z) {
#         z[!is.na(z) & z != ""]
#       })
#     split_rules_df <- do.call(rbind, split_rules)
#     indexes <-
#       which(!is.na(as.numeric(split_rules_df[1,])) == TRUE) ## gives the NA introduced by coersion error so need to come back and fix this so user doesn't get messages
#
#     rules_cln_list <- list()
#     for (k in 1:length(indexes)) {
#       metrix_index <- indexes[k]
#       ave_metric <-
#         round(mean(as.numeric(split_rules_df[, metrix_index])), 3)
#       dir <-
#         names(which.max(table(split_rules_df[, metrix_index - 1])))
#       rule <- paste(var, dir, ave_metric)
#       rules_cln_list[k] <- rule
#     }
#
#     m_rule <- paste(rules_cln_list, collapse = " & ")
#     marg_ave_rules[[i]] <- m_rule
#   } else{
#     marg_ave_rules[[i]] <- NA
#   }
#
# }
#
# mix_comp_mins <- list()
# mix_comp_maxs <- list()
#
# for (i in seq(dim(mixture_data)[2])) {
#   data <- mixture_data[, i]
#   mix_comp_mins[i] <- min(data)
#   mix_comp_maxs[i] <- max(data)
# }
#
# names(mix_comp_mins) <- mix_comps
# names(mix_comp_maxs) <- mix_comps
#
# names(marg_ave_rules) <- mix_comps
# marginal_results <-
#   merge(marginal_results, t(as.data.frame(marg_ave_rules)), by = "row.names")
#
# colnames(marginal_results)[1] <- "Mixture Comps."
# rownames(marginal_results) <- marginal_results$`Mixture Comps.`
# marginal_results <- marginal_results[,-1]
#
# colnames(marginal_results)[length(marginal_results)] <-
#   "Marg. Rules"
#
#
# marginal_results <-
#   merge(marginal_results, t(as.data.frame(mix_comp_mins)), by = "row.names")
# colnames(marginal_results)[length(marginal_results)] <-
#   "min value"
# rownames(marginal_results) <- marginal_results$Row.names
# marginal_results <- marginal_results[,-1]
#
# marginal_results <-
#   merge(marginal_results, t(as.data.frame(mix_comp_maxs)), by = "row.names")
# rownames(marginal_results) <- marginal_results$Row.names
# marginal_results <- marginal_results[,-1]
# colnames(marginal_results)[length(marginal_results)] <-
#   "max value"
