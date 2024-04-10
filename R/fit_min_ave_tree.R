#' @title Fit minimum average tree
#'
#' @details Fits the min ave tree
#'
#' @param at Training dataframe
#' @param a Variable names in the mixture
#' @param w Variable names in the covariates
#' @param y Variable name for the outcome
#' @param fold Current fold in the cross-validation
#' @param min_max Min or Max oracle region to go for
#' @param min_obs Minimum number of observations needed to make a split
#' @param parallel_cv TRUE/FALSE indicator to parallelize cv
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @import dplyr
#' @importFrom stringr str_detect str_extract_all
#' @importFrom dplyr group_by filter top_n
#' @import ranger
#' @return A list of the mixture rule results within a fold including:
#'  \itemize{
#'   \item \code{rules}: A data frame with the data adpatively
#'   determined rules found in the \code{pre} model along with the coefficient,
#'   direction, fold, RMSE and other measures.
#'   \item \code{model}: The best fitting pre model found in the fold.
#'   }
#' @examples
#' data <- simulate_mixture_cube()
#' mix_comps <- c("M1", "M2", "M3")
#' W <- c("age", "sex", "bmi")
#' sls <- create_sls()
#' w_stack <- sls$W_stack
#' tree_stack <- sls$A_stack
#' example_output <- fit_pre_algorithm(
#'   at = data,
#'   a = mix_comps,
#'   w = W,
#'   y = "y",
#'   direction = "positive",
#'   w_stack = w_stack,
#'   fold = 1,
#'   max_iter = 1,
#'   verbose = FALSE,
#'   parallel = FALSE,
#'   seed = 6442
#' )
#' @export

fit_min_ave_tree_algorithm <- function(at,
                                       a,
                                       w,
                                       y,
                                       fold,
                                       parallel_cv,
                                       min_max,
                                       min_obs,
                                       max_depth = 2) {
  if (parallel_cv == TRUE) {
    future::plan(future::sequential, gc = TRUE)
  }

  # we need to impute any covariates before beginning the search algorithm
  at <- at %>%
    mutate(across(all_of(w), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))

  calculate_average <- function(region, w_names, var, outcome, depth) {

    if (depth == 0 ) {
      region_average <- mean(region[[outcome]])
      n <- nrow(region)
    }else{

      region_data <- region[, c(w_names, var, outcome)]

      outcome_model <- ranger(
          formula = as.formula(paste(outcome, "~.")),
          data = region_data, num.trees = 400
        )

        propensity_model <- ranger(
          formula = as.formula(paste(
            "split_indicator ~ ",
            paste(w_names, collapse = "+")
          )),
          data = data, num.trees = 400
        )

        propensity_predictions <- propensity_model$predictions
        propensity_predictions <- ifelse(propensity_predictions < 0.0001, 0.0001, propensity_predictions)
        propensity_predictions <- ifelse(propensity_predictions > 0.99, 0.99, propensity_predictions)

        split_1_haw <- calc_clever_covariate(ghat_1_w = 1/propensity_predictions,
                              data = data, exposure = "split_indicator",h_aw_trunc_lvl = 10)

        split_0_haw <- calc_clever_covariate(ghat_1_w = 1/1-propensity_predictions,
                                             data = data, exposure = "split_indicator",h_aw_trunc_lvl = 10)


        data_1 <- data
        data_1$split_indicator <- 1

        data_0 <- data
        data_0$split_indicator <- 0

        outcome_preds_aw <- predict(outcome_model, data = data)$predictions
        outcome_preds_1w <- predict(outcome_model, data = data_1)$predictions
        outcome_preds_0w <- predict(outcome_model, data = data_0)$predictions

        tmle_update_1 <- fit_least_fav_submodel(h_aw = split_1_haw, data = data, y = outcome, qbar_aw = outcome_preds_aw, qbar_1w = outcome_preds_1w)
        tmle_update_0 <- fit_least_fav_submodel(h_aw = split_0_haw, data = data, y = outcome, qbar_aw = outcome_preds_aw, qbar_1w = outcome_preds_0w)

        # Calculate the average predictions for left and right
        left_average <- mean(tmle_update_1$qbar_1w_star)
        right_average <- mean(tmle_update_0$qbar_1w_star)

      rf <- ranger(as.formula(paste(outcome, "~.")), data = region_data)

      region_average <- mean(rf$predictions)
    }
      return(list(average = region_average, n = n))

  }

  split_data <- function(data, split_var, split_point) {
    left <- subset(data, data[[split_var]] <= split_point)
    right <- subset(data, data[[split_var]] > split_point)
    return(list(left = left, right = right))
  }

  clean_rule <- function(rule) {
    rule <- gsub("&\\s*$", "", rule)
    return(rule)
  }

  rules_to_dataframe <- function(rules_list) {
    rules_df <- data.frame(
      RegionMean = numeric(),
      N = integer(),
      Rule = character(),
      Depth = integer(),
      stringsAsFactors = FALSE
    )

    for (rule in rules_list) {
      rules_df <- rbind(rules_df, data.frame(
        RegionMean = rule$RegionMean,
        N = rule$N,
        Rule = rule$Rule,
        Depth = rule$Depth,
        stringsAsFactors = FALSE
      ))
    }

    return(rules_df)
  }

  find_best_split <- function(data, split_variables, w_names,
                              outcome, min_obs = 10, parent_average,
                              min_max) {
    best_split <- NULL
    min_average <- Inf
    max_average <- -Inf

    for (var in split_variables) {
      unique_values <- sort(unique(data[[var]]))

      for (split_point in unique_values) {
        # Create a binary indicator for the split
        data$split_indicator <- ifelse(data[[var]] <= split_point, 1, 0)

        # Check for minimum observations in each split
        if (sum(data$split_indicator) < min_obs || sum(!data$split_indicator)
        < min_obs || length(unique(data[[outcome]])) == 1) {
          next
        }

        covars <- c(w_names,  a[a != var])

        outcome_model <- ranger(
          formula = as.formula(paste(
            outcome, "~ split_indicator +",
            paste(covars, collapse = "+")
          )),
          data = data, num.trees = 400
        )

        propensity_model <- ranger(
          formula = as.formula(paste(
            "split_indicator ~ ",
            paste(w_names, collapse = "+")
          )),
          data = data, num.trees = 400
        )

        propensity_predictions <- propensity_model$predictions
        propensity_predictions <- ifelse(propensity_predictions < 0.0001, 0.0001, propensity_predictions)
        propensity_predictions <- ifelse(propensity_predictions > 0.99, 0.99, propensity_predictions)

        split_1_haw <- calc_clever_covariate(ghat_1_w = 1/propensity_predictions,
                              data = data, exposure = "split_indicator",h_aw_trunc_lvl = 10)

        split_0_haw <- calc_clever_covariate(ghat_1_w = 1/1-propensity_predictions,
                                             data = data, exposure = "split_indicator",h_aw_trunc_lvl = 10)


        data_1 <- data
        data_1$split_indicator <- 1

        data_0 <- data
        data_0$split_indicator <- 0

        outcome_preds_aw <- predict(outcome_model, data = data)$predictions
        outcome_preds_1w <- predict(outcome_model, data = data_1)$predictions
        outcome_preds_0w <- predict(outcome_model, data = data_0)$predictions

        tmle_update_1 <- fit_least_fav_submodel(h_aw = split_1_haw, data = data, y = outcome, qbar_aw = outcome_preds_aw, qbar_1w = outcome_preds_1w)
        tmle_update_0 <- fit_least_fav_submodel(h_aw = split_0_haw, data = data, y = outcome, qbar_aw = outcome_preds_aw, qbar_1w = outcome_preds_0w)

        # Calculate the average predictions for left and right
        left_average <- mean(tmle_update_1$qbar_1w_star)
        right_average <- mean(tmle_update_0$qbar_1w_star)

        if (min_max == "min") {
          # Compare and update the best split
          if (left_average < min_average && left_average < parent_average) {
            min_average <- left_average
            best_split <- list(
              variable = var, point = split_point,
              side = "left", average = min_average
            )
          }

          if (right_average < min_average && right_average < parent_average) {
            min_average <- right_average
            best_split <- list(
              variable = var, point = split_point,
              side = "right", average = min_average
            )
          }
        } else {
          # Compare and update the best split
          if (left_average > max_average && left_average > parent_average) {
            max_average <- left_average
            best_split <- list(
              variable = var, point = split_point,
              side = "left", average = max_average
            )
          }

          if (right_average > max_average && right_average > parent_average) {
            max_average <- right_average
            best_split <- list(
              variable = var, point = split_point,
              side = "right", average = max_average
            )
          }
        }
      }
    }

    # Clean up the temporary variable
    data$split_indicator <- NULL

    if (!is.null(best_split)) {
      return(best_split)
    } else {
      return(NULL) # No split found that reduces the average below the parent's
    }
  }

  recursive_split_all_rules <- function(data,
                                        split_variables,
                                        w_names,
                                        depth = 0,
                                        max_depth = 3,
                                        outcome,
                                        path = "",
                                        min_obs = min_obs,
                                        min_max = "min",
                                        parent_average = NULL) {
    # Calculate the parent mean before attempting to find the best split
    if(depth == 0){
      parent_average <- mean(data[[outcome]])
      n <- nrow(data)
    }


    if (depth == max_depth || nrow(data) == 0) {
      current_rule <- list(
        "RegionMean" = parent_average, # Updated to use parent_stats
        "N" = nrow(data), # Updated to use parent_stats
        "Rule" = clean_rule(path),
        "Depth" = depth
      )
      return(list(current_rule))
    }

    # Find the best split using the parent mean calculated earlier
    best_split <- find_best_split(data, split_variables, w_names, outcome,
      min_obs,
      parent_average = parent_average,
      min_max = min_max
    )

    if (is.null(best_split)) { # If no best split found, return the current path
      current_rule <- list(
        "RegionMean" = parent_average,
        "N" = nrow(data),
        "Rule" = clean_rule(path),
        "Depth" = depth
      )
      return(list(current_rule))
    }

    # Perform the split
    splits <- split_data(data, best_split$variable, best_split$point)


    if (best_split$side == "left") {

      # Construct the rules for the left and right branches
      left_rule <- paste0(
        path, best_split$variable, " <= ",
        best_split$point, " & "
      )

      left_rules <- recursive_split_all_rules(
        data = splits$left,
        split_variables = split_variables,
        w_names = w_names,
        depth = depth + 1,
        max_depth = max_depth,
        outcome = outcome,
        path = left_rule,
        min_obs = min_obs,
        parent_average = best_split$average,
        min_max = min_max
      )

      right_rules <- NULL
    }else{
      right_rule <- paste0(
        path, best_split$variable, " > ",
        best_split$point, " & "
      )
      right_rules <- recursive_split_all_rules(
      data = splits$right,
      split_variables = split_variables,
      w_names = w_names,
      depth = depth + 1,
      max_depth = max_depth,
      outcome = outcome,
      path = right_rule,
      min_obs = min_obs,
      parent_average = best_split$average,
      min_max = min_max
      )
      left_rules <- NULL
      }




    # Combine the rules from the left and right branches and return them
    return(c(left_rules, right_rules))
  }

  at[, a] <- apply(at[, a], 2, round, 2)
  min_ave_tree_results <- recursive_split_all_rules(
    data = at,
    split_variables = a,
    w_names = w,
    max_depth = max_depth,
    outcome = y,
    min_max = min_max,
    min_obs = min_obs
  )

  tree <- rules_to_dataframe(min_ave_tree_results)

  if (min_max == "min") {
    region <- tree[which.min(tree$RegionMean), ]
  } else {
    region <- tree[which.max(tree$RegionMean), ]
  }

  rules <- data.frame(matrix(nrow = nrow(region), ncol = 7))
  colnames(rules) <- c(
    "rule", "coefficient", "description",
    "test", "direction", "fold", "RMSE"
  )

  rules$rule <- paste("rule", seq(nrow(region)), sep = "")
  rules$coefficient <- region$RegionMean
  rules$description <- region$Rule
  rules$effect_modifiers <- as.numeric(str_detect(
    region$Rule,
    paste(w, collapse = "|")
  ))
  rules$direction <- 1
  rules$RMSE <- NA
  rules$fold <- fold
  rules$test <- sort(paste(
    str_extract_all(
      region$Rule,
      paste(c(a, w), collapse = "|")
    )[[1]],
    collapse = "-"
  ))


  if (dim(rules)[1] == 0) {
    rules <- data.frame(matrix(nrow = 1, ncol = 7))
    colnames(rules) <- c(
      "rule", "coefficient", "description",
      "test", "direction", "fold", "RMSE"
    )
    rules$rule <- "None"
    rules$coefficient <- 0
    rules$description <- "No Rules Found"
    rules$test <- "No Rules Found"
    rules$effect_modifiers <- "None"
    rules$direction <- -1
    rules$fold <- fold
    rules$RMSE <- NA
  }

  return(list("rules" = rules, "model" = min_ave_tree_results))
}
