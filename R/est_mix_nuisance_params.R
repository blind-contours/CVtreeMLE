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
#' @param no_mix_rules TRUE/FALSE indicator for if no mixture rules were found
#'
#' @param aw_stack Super Learner library for fitting Q (outcome mechanism) and
#' g (treatment mechanism)
#' @param family Binomial or gaussian
#' @param rules Dataframe of rules found during the PRE fitting process
#' @param h_aw_trunc_lvl Truncation level of the clever covariate (induces more
#'  bias to reduce variance)
#' @param parallel_cv TRUE/FALSE if cv parallelization is used
#' @param seed Seed number
#' @import sl3
#' @importFrom pre pre
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by filter top_n
#' @return Rules object. TODO: add more detail here.

#'
#' @export

est_mix_nuisance_params <- function(at,
                                    av,
                                    w,
                                    no_mix_rules,
                                    aw_stack,
                                    family,
                                    rules,
                                    h_aw_trunc_lvl,
                                    parallel_cv,
                                    seed) {
  if (parallel_cv == TRUE) {
    future::plan(future::sequential, gc = TRUE)
  }

  set.seed(seed)

  at_mix <- at
  av_mix <- av

  if (no_mix_rules != TRUE) {
    at_rules_eval <-
      evaluate_mixture_rules(data = at_mix, rules = rules)
    av_rules_eval <-
      evaluate_mixture_rules(data = av_mix, rules = rules)

    mix_interaction_data <- list()

    for (interaction in seq(dim(at_rules_eval)[2])) {
      interaction_rule <- at_rules_eval[, interaction]

      if (dim(table(interaction_rule)) == 2) {
        at_mix$A_mix <- interaction_rule
        av_mix$A_mix <- av_rules_eval[, interaction]

        task_at <- sl3::make_sl3_Task(
          data = at_mix,
          covariates = w,
          outcome = "A_mix",
          outcome_type = "binomial",
          folds = 2
        )

        task_av <- sl3::make_sl3_Task(
          data = av_mix,
          covariates = w,
          outcome = "A_mix",
          outcome_type = "binomial"
        )

        discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new(
          sl3::loss_loglik_binomial)

        discrete_sl <- sl3::Lrnr_sl$new(
          learners = aw_stack,
          metalearner = discrete_sl_metalrn,
        )

        sl_fit <- suppressWarnings(discrete_sl$train(task_at))

        ghat_1w <- sl_fit$predict(task_av)

        h_aw <- calc_clever_covariate(ghat_1_w = ghat_1w,
                                      data = av_mix,
                                      exposure = "A_mix",
                                      h_aw_trunc_lvl = 10)

        task_at <- sl3::make_sl3_Task(
          data = at_mix,
          covariates = c(w, "A_mix"),
          outcome = "y_scaled",
          outcome_type = family,
          folds = 2
        )

        x_m1 <- x_m0 <- av_mix
        x_m1$A_mix <- 1 # under exposure
        x_m0$A_mix <- 0 # under control

        task_av <- sl3::make_sl3_Task(
          data = av_mix,
          covariates = c(w, "A_mix"),
          outcome = "y_scaled",
          outcome_type = family
        )

        task_av_1 <- sl3::make_sl3_Task(
          data = x_m1,
          covariates = c(w, "A_mix"),
          outcome = "y_scaled",
          outcome_type = family
        )

        task_av_0 <- sl3::make_sl3_Task(
          data = x_m0,
          covariates = c(w, "A_mix"),
          outcome = "y_scaled",
          outcome_type = family
        )

        discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new(
          sl3::loss_squared_error)

        discrete_sl <- sl3::Lrnr_sl$new(
          learners = aw_stack,
          metalearner = discrete_sl_metalrn,
        )

        sl_fit <- suppressWarnings(discrete_sl$train(task_at))

        sl_fit$predict(task_av)

        qbar_aw <- bound_precision(sl_fit$predict(task_av))
        qbar_1w <- bound_precision(sl_fit$predict(task_av_1))
        qbar_0w <- bound_precision(sl_fit$predict(task_av_0))

        ## add Qbar to the AV dataset
        av_mix$qbar_aw <- qbar_aw
        av_mix$qbar_1w <- qbar_1w
        av_mix$qbar_0w <- qbar_0w

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
  } else {
    mix_interaction_data <- list()

    av_mix$ghat_1w <- NA
    av_mix$h_aw <- NA
    av_mix$qbar_aw <- NA
    av_mix$qbar_1w <- NA
    av_mix$qbar_0w <- NA
    mix_interaction_data[[1]] <- av_mix
  }

  return(list(data = mix_interaction_data))
}
