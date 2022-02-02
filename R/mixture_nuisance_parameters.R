#' For each mixture rule found, create a g estimator for the probability of being exposed to the rule thresholds,
#' and a Q estimator for the outcome E(Y| A = a_mix, W). Get estimates of g and Q using the validation data and
#' calculate the clever covariate used in the TMLE fluctuation step.
#'
#' @param At Training data
#' @param Av Validation data
#' @param W Vector of characters denoting covariates
#' @param no_mix_rules TRUE/FALSE indicator for if no mixture rules were found
#'
#' @param Q1_stack Super Learner library for fitting Q (outcome mechanism) and g (treatment mechanism)
#' @param family Binomial or gaussian
#' @param rules Dataframe of rules found during the PRE fitting process
#' @param H.AW_trunc_lvl Truncation level of the clever covariate (induces more bias to reduce variance)
#'
#' @import sl3
#' @importFrom pre pre
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by filter top_n
#' @return Rules object. TODO: add more detail here.

#'
#' @export

est_mix_nuisance_params <- function(At, Av, W, no_mix_rules, Q1_stack, family, rules, H.AW_trunc_lvl) {
  future::plan(future::sequential, gc = TRUE)

  At_mix <- At
  Av_mix <- Av

  if (no_mix_rules != TRUE) {
    At_rules_eval <-
      evaluate_mixture_rules(data = At_mix, rules = rules)
    Av_rules_eval <-
      evaluate_mixture_rules(data = Av_mix, rules = rules)

    mix_interaction_data <- list()
    mix_directions <- list()

    for (interaction in seq(dim(At_rules_eval)[2])) {
      interaction_rule <- At_rules_eval[, interaction]

      if (dim(table(interaction_rule)) == 2) {
        At_mix$A_mix <- interaction_rule
        Av_mix$A_mix <- Av_rules_eval[, interaction]

        task_At <- sl3::make_sl3_Task(
          data = At_mix,
          covariates = W,
          outcome = "A_mix",
          outcome_type = "binomial",
          folds = 2
        )

        task_Av <- sl3::make_sl3_Task(
          data = Av_mix,
          covariates = W,
          outcome = "A_mix",
          outcome_type = "binomial",
          folds = 2
        )

        discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new()

        discrete_sl <- sl3::Lrnr_sl$new(
          learners = Q1_stack,
          metalearner = discrete_sl_metalrn,
        )

        sl_fit <- discrete_sl$train(task_At)

        gHat1W <- sl_fit$predict(task_Av)

        H.AW <- calc_clever_covariate(gHat1W = gHat1W, data = Av_mix, exposure = "A_mix", H.AW_trunc_lvl = 10)

        ## add treatment mechanism results to validation dataframe


        task_At <- sl3::make_sl3_Task(
          data = At_mix,
          covariates = c(W, "A_mix"),
          outcome = "y_scaled",
          outcome_type = family
        )

        X_m1 <- X_m0 <- Av_mix
        X_m1$A_mix <- 1 # under exposure
        X_m0$A_mix <- 0 # under control

        task_Av <- sl3::make_sl3_Task(
          data = Av_mix,
          covariates = c(W, "A_mix"),
          outcome = "y_scaled",
          outcome_type = family
        )

        task_Av_1 <- sl3::make_sl3_Task(
          data = X_m1,
          covariates = c(W, "A_mix"),
          outcome = "y_scaled",
          outcome_type = family
        )

        task_Av_0 <- sl3::make_sl3_Task(
          data = X_m0,
          covariates = c(W, "A_mix"),
          outcome = "y_scaled",
          outcome_type = family
        )

        sl_fit <- discrete_sl$train(task_At)

        sl_fit$predict(task_Av)

        QbarAW <- bound_precision(sl_fit$predict(task_Av))
        Qbar1W <- bound_precision(sl_fit$predict(task_Av_1))
        Qbar0W <- bound_precision(sl_fit$predict(task_Av_0))

        ## add Qbar to the AV dataset
        Av_mix$QbarAW <- QbarAW
        Av_mix$Qbar1W <- Qbar1W
        Av_mix$Qbar0W <- Qbar0W

        Av_mix$gHat1W <- gHat1W
        Av_mix$H.AW <- H.AW

        mix_interaction_data[[interaction]] <- Av_mix
      } else {
        Av_mix$gHat1W <- NA
        Av_mix$H.AW <- NA
        Av_mix$QbarAW <- NA
        Av_mix$Qbar1W <- NA
        Av_mix$Qbar0W <- NA
        mix_interaction_data[[interaction]] <- Av_mix
      }
    }
  } else {
    mix_interaction_data <- list()

    Av_mix$gHat1W <- NA
    Av_mix$H.AW <- NA
    Av_mix$QbarAW <- NA
    Av_mix$Qbar1W <- NA
    Av_mix$Qbar0W <- NA
    mix_interaction_data[[1]] <- Av_mix
  }

  return(list(data = mix_interaction_data))
}
