#' For each mixture rule found, create a g estimator for the probability of being exposed to the rule thresholds,
#' and a Q estimator for the outcome E(Y| A = a_mix, W). Get estimates of g and Q using the validation data and
#' calculate the clever covariate used in the TMLE fluctuation step.
#'
#' @param At Training data
#' @param Av Validation data
#' @param W Vector of characters denoting covariates
#' @param no_rules TRUE/FALSE indicator for if no mixture rules were found
#'
#' @param SL.library Super Learner library for fitting Q (outcome mechanism) and g (treatment mechanism)
#' @param family Binomial or gaussian
#' @param rules Dataframe of rules found during the PRE fitting process
#' @param H.AW_trunc_lvl Truncation level of the clever covariate (induces more bias to reduce variance)
#' @importFrom pre pre
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by filter top_n
#' @return Rules object. TODO: add more detail here.

#'
#' @export

est_mix_nuisance_params <- function(At, Av, W, no_rules, SL.library, family, rules, H.AW_trunc_lvl) {
  At_mix <- At
  Av_mix <- Av

  if (no_rules != TRUE) {
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

        ## select covariates
        X_Amix_T <- subset(At_mix,
          select = c(W)
        )

        X_Amix_V <- subset(Av_mix,
          select = c(W)
        )

        gHatSL <- SuperLearner::SuperLearner(
          Y = At_mix$A_mix,
          X = X_Amix_T,
          SL.library = SL.library,
          family = "binomial",
          verbose = FALSE
        )

        gHat1W <- predict(gHatSL, newdata = X_Amix_V)$pred

        H.AW <- calc_clever_covariate(gHat1W = gHat1W, data = Av_mix, exposure = "A_mix", H.AW_trunc_lvl = 10)

        ## add treatment mechanism results to validation dataframe

        X_train_mix <- At_mix[c("A_mix", W)]

        X_valid_mix <- Av_mix[c("A_mix", W)]

        ## QbarAW
        QbarAWSL_m <- SuperLearner::SuperLearner(
          Y = At_mix$y_scaled,
          X = X_train_mix,
          SL.library = SL.library,
          family = family,
          verbose = FALSE
        )

        X_m1 <- X_m0 <- X_train_mix
        X_m1$A_mix <- 1 # under exposure
        X_m0$A_mix <- 0 # under control

        QbarAW <- bound_precision(predict(QbarAWSL_m, newdata = X_train_mix)$pred)
        Qbar1W <- bound_precision(predict(QbarAWSL_m, newdata = X_m1)$pred)
        Qbar0W <- bound_precision(predict(QbarAWSL_m, newdata = X_m0)$pred)

        X_m1 <- X_m0 <- X_valid_mix
        X_m1$A_mix <- 1 # under exposure
        X_m0$A_mix <- 0 # under control

        QbarAW <- bound_precision(predict(QbarAWSL_m, newdata = X_valid_mix)$pred)
        Qbar1W <- bound_precision(predict(QbarAWSL_m, newdata = X_m1)$pred)
        Qbar0W <- bound_precision(predict(QbarAWSL_m, newdata = X_m0)$pred)

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
