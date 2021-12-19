#' For each marginal mixture component rule found, create a g estimator for the probability of being exposed to the rule thresholds,
#' and a Q estimator for the outcome E(Y| A = a_mix, W). Get estimates of g and Q using the validation data and
#' calculate the clever covariate used in the TMLE fluctuation step.
#'
#' @param At Training data
#' @param Av Validation data
#' @param W Vector of characters denoting covariates
#' @param SL.library Super Learner library for fitting Q (outcome mechanism) and g (treatment mechanism)
#' @param family Binomial or gaussian
#' @param A Vector of characters that denote the mixture components
#' @param marg_decisions List of rules found within the fold for each mixture component
#' @param H.AW_trunc_lvl Truncation level of the clever covariate (induces more bias to reduce variance)
#' @importFrom SuperLearner All
#' @importFrom magrittr %>%
#' @importFrom rlang :=

#' @importFrom dplyr group_by filter top_n
#' @return Rules object. TODO: add more detail here.

#'
#' @export

est_marg_nuisance_params <- function(At,
                                     Av,
                                     W,
                                     SL.library,
                                     family,
                                     A,
                                     marg_decisions,
                                     H.AW_trunc_lvl) {
  marginal_data <- list()
  marg_decisions$directions <- NA

  for (i in seq(A)) {
    At_c <- At
    Av_c <- Av

    target_m_rule <- marg_decisions[marg_decisions$target_m == A[i],]$rules

    if (target_m_rule != "No Rules Found") {
      rule_name <- paste(A[i], "marg_rule", sep = "_")

      At_c <- At_c %>%
        dplyr::mutate(!!(rule_name) := ifelse(eval(parse(text = target_m_rule)), 1, 0))

      Av_c <- Av_c %>%
        dplyr::mutate(!!(rule_name) := ifelse(eval(parse(text = target_m_rule)), 1, 0))

      ## covars
      X_Amarg_T <- At_c[W]

      X_Amarg_V <- Av_c[W]

      gHatSL <- SuperLearner::SuperLearner(
        Y = At_c[[rule_name]],
        X = X_Amarg_T,
        SL.library = SL.library,
        family = "binomial",
        verbose = FALSE
      )

      gHat1W <- bound_precision(predict(gHatSL, X_Amarg_V)$pred)

      H.AW <- calc_clever_covariate(gHat1W = gHat1W, data = Av_c, exposure = rule_name, H.AW_trunc_lvl = 10, type = "reg")

      ## add treatment mechanism results to Av_c dataframe
      Av_c$gHat1W <- gHat1W
      Av_c$H.AW <- H.AW

      X_train_mix <- At_c[c(rule_name, W)]
      X_valid_mix <- Av_c[c(rule_name, W)]

      ## QbarAW
      QbarAWSL_m <- SuperLearner::SuperLearner(
        Y = At_c$y_scaled,
        X = X_train_mix,
        SL.library = SL.library,
        family = family,
        verbose = FALSE
      )

      X_m1 <- X_m0 <- X_valid_mix
      X_m1[[rule_name]] <- 1 # under exposure
      X_m0[[rule_name]] <- 0 # under control

      QbarAW <- bound_precision(predict(QbarAWSL_m, newdata = X_valid_mix)$pred)
      Qbar1W <- bound_precision(predict(QbarAWSL_m, newdata = X_m1)$pred)
      Qbar0W <- bound_precision(predict(QbarAWSL_m, newdata = X_m0)$pred)

      flux_results <- fit_least_fav_submodel(H.AW, data = Av, QbarAW, Qbar1W, Qbar0W)

      QbarAW.star <- flux_results$QbarAW.star
      Qbar1W.star <- flux_results$Qbar1W.star
      Qbar0W.star <- flux_results$Qbar0W.star

      if (mean(Qbar1W.star - Qbar0W.star, na.rm = TRUE) > 0) {
        marg_decisions$directions[i] <- "positive"

        X_m1 <- X_m0 <- X_valid_mix
        X_m1[[rule_name]] <- 1 # under exposure
        X_m0[[rule_name]] <- 0 # under control

        QbarAW <- bound_precision(predict(QbarAWSL_m, newdata = X_valid_mix)$pred)
        Qbar1W <- bound_precision(predict(QbarAWSL_m, newdata = X_m1)$pred)
        Qbar0W <- bound_precision(predict(QbarAWSL_m, newdata = X_m0)$pred)

        Av_c$QbarAW <- QbarAW
        Av_c$Qbar1W <- Qbar1W
        Av_c$Qbar0W <- Qbar0W

        marginal_data[[i]] <- Av_c
      } else {
        marg_decisions$directions[i] <- "negative"

        At_c[, rule_name] <- 1 - At_c[, rule_name]

        ## covars
        X_Amarg_T <- At_c[W]
        X_Amarg_V <- Av_c[W]

        gHatSL <- SuperLearner::SuperLearner(
          Y = At_c[[rule_name]],
          X = X_Amarg_T,
          SL.library = SL.library,
          family = "binomial",
          verbose = FALSE
        )

        gHat1W <- predict(gHatSL, newdata = X_Amarg_V)$pred

        H.AW <- calc_clever_covariate(gHat1W, data = Av_c, exposure = rule_name, H.AW_trunc_lvl = 10, type = "reg")

        ## add treatment mechanism results to Av_c dataframe
        Av_c$gHat1W <- gHat1W
        Av_c$H.AW <- H.AW

        X_train_mix <-At_c[c(rule_name, W)]
        X_valid_mix <-Av_c[c(rule_name, W)]

        print(paste("Fitting SL of Y given W and rule for mixture", i))

        ## QbarAW
        QbarAWSL_m <- SuperLearner::SuperLearner(
          Y = At_c$y_scaled,
          X = X_train_mix,
          SL.library = SL.library,
          family = family,
          verbose = FALSE
        )

        X_m1 <- X_m0 <- X_valid_mix
        X_m1[[rule_name]] <- 1 # under exposure
        X_m0[[rule_name]] <- 0 # under control

        QbarAW <- bound_precision(predict(QbarAWSL_m, newdata = X_valid_mix)$pred)
        Qbar1W <- bound_precision(predict(QbarAWSL_m, newdata = X_m1)$pred)
        Qbar0W <- bound_precision(predict(QbarAWSL_m, newdata = X_m0)$pred)

        Av_c$QbarAW <- QbarAW
        Av_c$Qbar1W <- Qbar1W
        Av_c$Qbar0W <- Qbar0W

        marginal_data[[i]] <- Av_c
      }
    } else {
      print(
        paste(
          "No ATEs calculated in the validation for",
          A[i],
          "due to no rule found in training set for marginal impact"
        )
      )
      marg_decisions$directions[i] <- NA
    }
  }

  return(list(
    "data" = marginal_data,
    "marg_decisions" = marg_decisions
  ))
}
