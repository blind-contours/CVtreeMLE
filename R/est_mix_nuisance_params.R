#' @title Estimate nuisance parameters for each mixture interaction identified
#'
#' @description For each mixture mixture interaction found, create a g estimator
#'  for the probability of being exposed to the rule thresholds,
#' and a Q estimator for the outcome E(Y| A = a_mix, W). Get estimates of g
#' and Q using the validation data and
#' calculate the clever covariate used in the TMLE fluctuation step.
#'
#' @param at Training data
#' @param av Validation data
#' @param w Vector of characters denoting covariates
#' @param y The outcome variable
#' @param no_mix_rules TRUE/FALSE indicator for if no mixture rules were found
#'
#' @param aw_stack Super Learner library for fitting Q (outcome mechanism) and
#' g (treatment mechanism)
#' @param family Binomial or continuous
#' @param rules Dataframe of rules found during the PRE fitting process
#' @param h_aw_trunc_lvl Truncation level of the clever covariate (induces more
#'  bias to reduce variance)
#' @param parallel_cv TRUE/FALSE if cv parallelization is used
#' @param seed Seed number
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by filter top_n
#' @importFrom stringr str_trim str_split
#' @return A list of dataframes where the nuisance parameters are added to
#' the raw data.

#'
#' @export

est_mix_nuisance_params <- function(at,
                                    av,
                                    w,
                                    a,
                                    y,
                                    aw_stack,
                                    family,
                                    rules,
                                    parallel_cv,
                                    seed,
                                    h_aw_trunc_lvl) {
  if (parallel_cv == TRUE) {
    future::plan(future::sequential, gc = TRUE)
  }

  identify_w_in_rule <- function(rule_description, w) {
    split_rule <- str_trim(unlist(str_split(rule_description, "&")))

    # Check which of the split elements contain variables from w
    variable_present <- sapply(split_rule, function(rule_part) {
      var_name <- unlist(str_split(rule_part, " "))[1]
      return(var_name %in% w)
    })

    # Return the elements that contain variables from w
    return(split_rule[variable_present])
  }

  set.seed(seed)

  mix_interaction_data <- list()

  for (interaction in seq(dim(rules)[1])) {
    rule_data <- rules[interaction, ]
    if (rule_data$effect_modifiers == 1) {
      w_rule <- paste(identify_w_in_rule(rule_data$description, w),
        collapse = " & "
      )
      a_rule <- paste(identify_w_in_rule(rule_data$description, a),
        collapse = " & "
      )
      # For 'at' dataset
      subset_at <- subset(at, eval(parse(text = w_rule)))

      # For 'av' dataset
      subset_av <- subset(av, eval(parse(text = w_rule)))

      interaction_rule_at <- subset_at %>%
        dplyr::transmute(A = ifelse(eval(parse(text = a_rule)), 1, 0))


      interaction_rule_av <- subset_av %>%
        dplyr::transmute(A = ifelse(eval(parse(text = a_rule)), 1, 0))

      at_mix <- subset_at
      av_mix <- subset_av
    } else {
      at_mix <- at
      av_mix <- av

      interaction_rule_at <- at_mix %>%
        dplyr::transmute(A = ifelse(eval(parse(text = rule_data$description)),
          1, 0
        ))

      interaction_rule_av <- av_mix %>%
        dplyr::transmute(A = ifelse(eval(parse(text = rule_data$description)),
          1, 0
        ))
    }

    if (dim(table(interaction_rule_at)) == 2) {
      at_mix$A_mix <- interaction_rule_at
      av_mix$A_mix <- interaction_rule_av

      exposures_used <- unique(unlist(strsplit(rule_data$test, split = "-")))
      other_a <- a[!a %in% exposures_used]

      task_at <- sl3::make_sl3_Task(
        data = at_mix,
        covariates = w,
        outcome = "A_mix",
        outcome_type = "binomial",
        folds = 10
      )

      task_av <- sl3::make_sl3_Task(
        data = av_mix,
        covariates = w,
        outcome = "A_mix",
        outcome_type = "binomial",
        folds = 10
      )

      cv_selector <- Lrnr_cv_selector$new(eval_function = loss_squared_error)
      sl <- Lrnr_sl$new(learners = aw_stack, metalearner = cv_selector)

      sl_fit <- suppressWarnings(sl$train(task_at))

      ghat_1w <- sl_fit$predict(task_av)

      h_aw <- calc_clever_covariate(
        ghat_1_w = ghat_1w,
        data = av_mix,
        exposure = "A_mix",
        h_aw_trunc_lvl = h_aw_trunc_lvl
      )

      task_at <- sl3::make_sl3_Task(
        data = at_mix,
        covariates = c(w, "A_mix", other_a),
        outcome = y,
        outcome_type = family,
        folds = 10
      )

      x_m1 <- av_mix
      x_m1$A_mix <- 1 # under exposure

      task_av <- sl3::make_sl3_Task(
        data = av_mix,
        covariates = c(w, "A_mix", other_a),
        outcome = y,
        outcome_type = family,
        folds = 10
      )

      task_av_1 <- sl3::make_sl3_Task(
        data = x_m1,
        covariates = c(w, "A_mix", other_a),
        outcome = y,
        outcome_type = family,
        folds = 10
      )

      sl_fit <- suppressWarnings(sl$train(task_at))

      qbar_aw <- sl_fit$predict(task_av)
      qbar_1w <- sl_fit$predict(task_av_1)

      ## add Qbar to the AV dataset
      av_mix$qbar_aw <- qbar_aw
      av_mix$qbar_1w <- qbar_1w

      av_mix$ghat_1w <- ghat_1w
      av_mix$h_aw <- h_aw

      mix_interaction_data[[interaction]] <- av_mix
    } else {
      av_mix$ghat_1w <- NA
      av_mix$h_aw <- NA
      av_mix$qbar_aw <- NA
      av_mix$qbar_1w <- NA
      av_mix$qbar_0w <- NA
      mix_interaction_data[[interaction]] <- av_mix
    }
  }


  return(list(data = mix_interaction_data))
}
