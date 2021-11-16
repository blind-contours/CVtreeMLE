find_common_mixture_rules <- function(group_list, 
                                      data = data, 
                                      mix_comps = mix_comps, 
                                      mixture_results, 
                                      n_folds){
  
  if(anyNA(group_list) == FALSE){
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
      
      # for(rule_idx in 1:nrow(group)){
      #   rule_row <- group[rule_idx,]
      # 
      #   rule <- rule_row$description
      #   # direction <- rule_row$direction
      # 
      vars <- mix_comps[mix_comps %in% strsplit(group$description[1], split = " ")[[1]]]

        # if (direction == "negative") {
        #   mixture_data <- mixture_data %>%
        #     mutate("neg_rule" := ifelse(eval(parse(
        #       text = paste(rule_row$description, collapse = "")
        #     )), 1, 0))
        #
        #   mixture_data$pos_rule <- 1 - mixture_data$neg_rule
        #
        #   new_rule <- list()
        #
        #   for (var in vars) {
        #     if(length(vars) == 1){
        #     var_min_1 <-
        #       mixture_data %>%
        #       group_by(pos_rule) %>%
        #       summarise(min = min(!!(as.name(var))))
        #     var_min_1 <- subset(var_min_1, pos_rule == 1, select = min)
        #
        #     var_min_0 <-
        #       mixture_data %>%
        #       group_by(pos_rule) %>%
        #       summarise(min = min(!!(as.name(var))))
        #     var_min_0 <- subset(var_min_0, pos_rule == 0, select = min)
        #
        #     var_max_1 <-
        #       mixture_data %>%
        #       group_by(pos_rule) %>%
        #       summarise(max = max(!!(as.name(var))))
        #     var_max_1 <- subset(var_max_1, pos_rule == 1, select = max)
        #
        #     var_max_0 <-
        #       mixture_data %>%
        #       group_by(pos_rule) %>%
        #       summarise(max = max(!!(as.name(var))))
        #     var_max_0 <- subset(var_max_0, pos_rule == 0, select = max)
        #
        #     rule1 <-
        #       paste(
        #         var,
        #         ">",
        #         round(var_min_1, 5),
        #         "&",
        #         var,
        #         "<",
        #         round(var_min_0, 5)
        #       )
        #
        #     rule2 <-
        #       paste(var,
        #             ">",
        #             round(var_max_0, 5),
        #             "&",
        #             var,
        #             "<",
        #             round(var_max_1, 5))
        #
        #     augmented_rule <- paste("(",rule1,")", "|", "(", rule2,")")
        #
        #     new_rule <- append(new_rule, augmented_rule)
        #
        #     }else{
        #       var_min_1 <-
        #         mixture_data %>%
        #         group_by(pos_rule) %>%
        #         summarise(min = min(!!(as.name(var))))
        #       var_min_1 <- subset(var_min_1, pos_rule == 1, select = min)
        #
        #       var_max_1 <-
        #         mixture_data %>%
        #         group_by(pos_rule) %>%
        #         summarise(max = max(!!(as.name(var))))
        #       var_max_1 <- subset(var_max_1, pos_rule == 1, select = max)
        #
        #       augmented_rule <- paste("(", var, ">",round(var_min_1, 3), "&", var, "<", round(var_max_1,3), ")")
        #
        #       new_rule <- append(new_rule, augmented_rule)
        #       rule <- paste(unlist(new_rule), collapse = " & ")
        #       }
        #
        #   }
        #
        #   # rule <- paste(unlist(new_rule), collapse = " & ")
        # }
# 
#         combined_rules <- append(combined_rules, rule)
#       }
      
      # total_rules[[i]] <- combined_rules
      
      rule_list <- group$description
      
      # combined_rule <- paste(unlist(combined_rule), collapse = " & ")
      fold_rules_eval <- list()
      
      for(k in seq(rule_list)){
        
        fold_rule <- rule_list[[k]]
        
        fold_rule <- mixture_data %>%
          transmute("fold_rule" := ifelse(eval(parse(text = fold_rule)), 1, 0))
        
        fold_rules_eval[[k]] <- fold_rule
        
      }
      
      fold_rules_df <- do.call(cbind,fold_rules_eval)
      colnames(fold_rules_df) <- paste("fold_", seq(n_folds))
      
      fold_rules_df$sum <- rowSums(fold_rules_df)
      
      fold_rules_df$all_folds <- as.numeric(fold_rules_df$sum == n_folds)
      fold_rules_df <- cbind(fold_rules_df, mixture_data)
      
      total_count <-  table(fold_rules_df$sum > 0)[[2]]
      
      if(dim(table(fold_rules_df$all_folds)) == 2){
        
        count <- table(fold_rules_df$all_folds > 0)[[2]]
        fraction <- count/total_count
        new_rule <- list()
        
        for (var in vars) {
          var_min <-
            fold_rules_df %>% group_by(all_folds) %>% summarise(min = min(!!(as.name(var))))
          var_min <- subset(var_min, all_folds == 1, select = min)
          var_max <-
            fold_rules_df %>% group_by(all_folds) %>% summarise(max = max(!!(as.name(var))))
          var_max <- subset(var_max, all_folds == 1, select = max)
          
          augmented_rule <- paste(var, ">",round(var_min, 3), "&", var, "<", round(var_max,3))
          
          new_rule <- append(new_rule, augmented_rule)
          
        }
        
        rule <- paste(unlist(new_rule), collapse = " & ")
        
        mixture_rules[i] <- rule
        fractions[i] <- fraction
        
      }else{
        mixture_rules[i] <-  "No Obs Overlap across all folds"
        fractions[i] <-  NA
        total_rules[i] <-  NA
        
      }
      
    }
      
      
    # # combined_rule <- paste(unlist(combined_rules), collapse = " & ")
    #   
    #   mixture_data <- mixture_data %>%
    #     mutate("comb_rule" := ifelse(eval(parse(
    #       text = paste(combined_rule, collapse = "")
    #     )), 1, 0))
    #   
    #   new_rule <- list()
    #   for (var in vars) {
    #     var_min <-
    #       mixture_data %>% group_by(comb_rule) %>% summarise(min = min(!!(as.name(var))))
    #     var_min <- subset(var_min, comb_rule == 1, select = min)
    #     var_max <-
    #       mixture_data %>% group_by(comb_rule) %>% summarise(max = max(!!(as.name(var))))
    #     var_max <- subset(var_max, comb_rule == 1, select = max)
    #     
    #     augmented_rule <- paste(var, ">",round(var_min, 3), "&", var, "<", round(var_max,3))
    #     
    #     new_rule <- append(new_rule, augmented_rule)
    #     
    #   }
    #   
    #   rule <- paste(unlist(new_rule), collapse = " & ")
    #   mixture_rules[i] <- rule
    # 
    # 
    # }
    # 
    # 
    # 
          
          
      
    #   x <- strsplit(group$description, " ")
    #   ## find numerics
    #   y <- x[[1]]
    #   indexes <- which(!is.na(as.numeric(y)) == TRUE)
    #   x_df <- do.call(rbind, x)
    #   x_df <- as.data.frame(x_df)
    #   x_df[indexes] <- sapply(x_df[indexes], as.numeric)
    #   #means <- colMeans(x_df[indexes])
    #   vars <- unique(unlist(unique(x_df[indexes - 2])))
    #   rule_mod <- list()
    #   for (row_id in seq(dim(x_df)[1])) {
    #     pos_neg <- group[row_id, "direction"]
    #     row <- x_df[row_id,]
    #     
    # 
    #     #
    #     #   for (var in vars) {
    #     #     var_min <-
    #     #       target_neg_var_df %>% group_by(neg_rule) %>% summarise(min = min(!!(as.name(var))))
    #     #     var_min <- subset(var_min, neg_rule == 1, select = min)
    #     #     var_max <-
    #     #       target_neg_var_df %>% group_by(neg_rule) %>% summarise(max = max(!!(as.name(var))))
    #     #     var_max <- subset(var_max, neg_rule == 1, select = max)
    #     #
    #     #     neg_rules_flipped[var] <-
    #     #       paste(var,
    #     #             ">",
    #     #             round(var_min, 3),
    #     #             "&",
    #     #             var,
    #     #             "<",
    #     #             round(var_max, 3))
    #     #   }
    #     #
    #     #   rule <- paste(neg_rules_flipped, collapse = " & ")
    #     #
    #     # } else{
    #     pos_rule <- mixture_data %>%
    #       transmute("pos_rule" := ifelse(eval(parse(
    #         text = paste(row, collapse = "")
    #       )), 1, 0))
    #     
    #     target_pos_var_df <-
    #       cbind(subset(mixture_data, select = vars), pos_rule)
    #     pos_rules <- list()
    #     
    #     for (var in vars) {
    #       var_min <-
    #         target_pos_var_df %>% group_by(pos_rule) %>% summarise(min = min(!!(as.name(var))))
    #       var_min <- subset(var_min, pos_rule == 1, select = min)
    #       var_max <-
    #         target_pos_var_df %>% group_by(pos_rule) %>% summarise(max = max(!!(as.name(var))))
    #       var_max <- subset(var_max, pos_rule == 1, select = max)
    #       
    #       pos_rules[var] <-
    #         paste(var,
    #               ">",
    #               round(var_min, 3),
    #               "&",
    #               var,
    #               "<",
    #               round(var_max, 3))
    #     }
    #     
    #     rule <- paste(pos_rules, collapse = " & ")
    #     #rule <- paste("(", rule, ")")
    #     
    #     rule_mod[[row_id]] <- rule
    #   }
    #   
    #   rule_df <- do.call(rbind, str_split(rule_mod, " "))
    #   rule_df <- as.data.frame(rule_df)
    #   indexes <- which(!is.na(as.numeric(rule_df[1,])) == TRUE)
    #   
    #   rule_df[indexes] <- sapply(rule_df[indexes], as.numeric)
    #   
    #   rule_collection <- list()
    #   for (j in seq(dim(rule_df)[2])) {
    #     if (is.numeric(rule_df[, j]) == TRUE) {
    #       target <- round(mean(rule_df[, j]), 3)
    #     } else{
    #       target <- names(which.max(table(rule_df[, j])))
    #     }
    #     
    #     rule_collection[j] <- target
    #   }
    #   
    #   ave_rules_in_set[[i]] <-
    #     paste(rule_collection, collapse = " ")
    #   
    # }
    
    mixture_results <-
      cbind(mixture_results, unlist(mixture_rules))
    
    colnames(mixture_results)[7] <- "Mixture Interaction Rules"
    
    mixture_results <-
      cbind(mixture_results, unlist(fractions))
    
    colnames(mixture_results)[8] <- "Fraction Covered"
    
    
  }else{
    mixture_results <-
      cbind(mixture_results, NA)
    
    colnames(mixture_results)[length(mixture_results)] <- "Mixture Interaction Rules"
    
    mixture_results <-
      cbind(mixture_results, NA)
    
    colnames(mixture_results)[length(mixture_results)] <- "Fraction Covered"
  }
    
  
  return(mixture_results)
}
