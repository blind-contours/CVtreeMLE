#' For each marginal mixture component rule found, create a g estimator for the probability of being exposed to the rule thresholds,
#' and a Q estimator for the outcome E(Y| A = a_mix, W). Get estimates of g and Q using the validation data and
#' calculate the clever covariate used in the TMLE fluctuation step.
#'
#' @param At Training data
#' @param Av Validation data
#' @param covariates Vector of characters denoting covariates
#' @param SL.library Super Learner library for fitting Q (outcome mechanism) and g (treatment mechanism)
#' @param family Binomial or gaussian
#' @param mix_comps Vector of characters that denote the mixture components
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
                                     covariates,
                                     SL.library,
                                     family,
                                     mix_comps,
                                     marg_decisions,
                                     H.AW_trunc_lvl) {
  marginal_data <- list()
  marg_directions <- list()

  for (i in seq(mix_comps)) {
    At_c <- At
    Av_c <- Av

    target_m_rule <- marg_decisions[i][[1]]

    if (target_m_rule != "1") {
      rule_name <- paste(mix_comps[i], "marg_rule", sep = "_")

      At_c <- At_c %>%
        dplyr::mutate(!!(rule_name) := ifelse(eval(parse(text = target_m_rule)), 1, 0))

      Av_c <- Av_c %>%
        dplyr::mutate(!!(rule_name) := ifelse(eval(parse(text = target_m_rule)), 1, 0))

      ## covars
      X_Amarg_T <- subset(At_c,
        select = c(covariates)
      )

      X_Amarg_V <- subset(Av_c,
        select = c(covariates)
      )

      gHatSL <- SuperLearner::SuperLearner(
        Y = At_c[[rule_name]],
        X = X_Amarg_T,
        SL.library = SL.library,
        family = "binomial",
        verbose = FALSE
      )

      gHat1W <- predict(gHatSL, X_Amarg_V)$pred
      n <- length(gHat1W)
      gHat0W <- 1 - gHat1W

      gHatAW <- rep(NA, n)
      gHatAW[Av_c[rule_name] == 1] <- gHat1W[Av_c[rule_name] == 1]
      gHatAW[Av_c[rule_name] == 0] <- gHat0W[Av_c[rule_name] == 0]

      H.AW <-
        as.numeric(Av_c[[rule_name]] == 1) / gHat1W - as.numeric(Av_c[[rule_name]] ==
          0) / gHat0W

      H.AW <-
        ifelse(H.AW > H.AW_trunc_lvl, H.AW_trunc_lvl, H.AW)

      H.AW <-
        ifelse(H.AW < -H.AW_trunc_lvl, -H.AW_trunc_lvl, H.AW)

      ## add treatment mechanism results to Av_c dataframe
      Av_c$gHat1W <- gHat1W
      Av_c$H.AW <- H.AW

      X_train_mix <-
        subset(At_c, select = c(rule_name, covariates))

      X_valid_mix <-
        subset(Av_c, select = c(rule_name, covariates))

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

      QbarAW <- predict(QbarAWSL_m, newdata = X_valid_mix)$pred
      Qbar1W <- predict(QbarAWSL_m, newdata = X_m1)$pred
      Qbar0W <- predict(QbarAWSL_m, newdata = X_m0)$pred

      QbarAW[QbarAW >= 1.0] <- 0.999
      QbarAW[QbarAW <= 0] <- 0.001

      Qbar1W[Qbar1W >= 1.0] <- 0.999
      Qbar1W[Qbar1W <= 0] <- 0.001

      Qbar0W[Qbar0W >= 1.0] <- 0.999
      Qbar0W[Qbar0W <= 0] <- 0.001

      ## least optimal submodel
      logitUpdate <-
        stats::glm(
          Av_c$y_scaled ~ -1 + H.AW + offset(qlogis(QbarAW)),
          family = "quasibinomial",
          data = At
        )

      epsilon <- logitUpdate$coef
      QbarAW.star <- plogis(qlogis(QbarAW) + epsilon * H.AW)
      Qbar1W.star <- plogis(qlogis(Qbar1W) + epsilon * H.AW)
      Qbar0W.star <- plogis(qlogis(Qbar0W) + epsilon * H.AW)

      if (mean(Qbar1W.star - Qbar0W.star, na.rm = TRUE) > 0) {
        marg_directions[i] <- "positive"

        X_m1 <- X_m0 <- X_valid_mix
        X_m1[[rule_name]] <- 1 # under exposure
        X_m0[[rule_name]] <- 0 # under control

        QbarAW <- predict(QbarAWSL_m, newdata = X_valid_mix)$pred
        Qbar1W <- predict(QbarAWSL_m, newdata = X_m1)$pred
        Qbar0W <- predict(QbarAWSL_m, newdata = X_m0)$pred

        QbarAW[QbarAW >= 1.0] <- 0.999
        QbarAW[QbarAW <= 0] <- 0.001

        Qbar1W[Qbar1W >= 1.0] <- 0.999
        Qbar1W[Qbar1W <= 0] <- 0.001

        Qbar0W[Qbar0W >= 1.0] <- 0.999
        Qbar0W[Qbar0W <= 0] <- 0.001

        Av_c$QbarAW <- QbarAW
        Av_c$Qbar1W <- Qbar1W
        Av_c$Qbar0W <- Qbar0W

        marginal_data[[i]] <- Av_c
      } else {
        marg_directions[i] <- "negative"

        At_c[, rule_name] <- 1 - At_c[, rule_name]

        ## covars
        X_Amarg_T <- subset(At_c,
          select = c(covariates)
        )

        X_Amarg_V <- subset(Av_c,
          select = c(covariates)
        )

        gHatSL <- SuperLearner::SuperLearner(
          Y = At_c[[rule_name]],
          X = X_Amarg_T,
          SL.library = SL.library,
          family = "binomial",
          verbose = FALSE
        )

        gHat1W <- predict(gHatSL, newdata = X_Amarg_V)$pred
        gHat0W <- 1 - gHat1W

        gHatAW <- rep(NA, n)
        gHatAW[Av_c[[rule_name]] == 1] <-
          gHat1W[Av_c[[rule_name]] == 1]
        gHatAW[Av_c[[rule_name]] == 0] <-
          gHat0W[Av_c[[rule_name]] == 0] ## this gives an incorrect length error but the estimates all look fine so don't know what's up yet

        H.AW <-
          as.numeric(Av_c[[rule_name]] == 1) / gHat1W - as.numeric(Av_c[[rule_name]] ==
            0) / gHat0W

        H.AW <-
          ifelse(H.AW > H.AW_trunc_lvl, H.AW_trunc_lvl, H.AW)

        H.AW <-
          ifelse(H.AW < -H.AW_trunc_lvl, -H.AW_trunc_lvl, H.AW)

        ## add treatment mechanism results to Av_c dataframe
        Av_c$gHat1W <- gHat1W
        Av_c$H.AW <- H.AW

        X_train_mix <-
          subset(At_c, select = c(rule_name, covariates))
        X_valid_mix <-
          subset(Av_c, select = c(rule_name, covariates))

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

        QbarAW <- predict(QbarAWSL_m, newdata = X_valid_mix)$pred
        Qbar1W <- predict(QbarAWSL_m, newdata = X_m1)$pred
        Qbar0W <- predict(QbarAWSL_m, newdata = X_m0)$pred

        QbarAW[QbarAW >= 1.0] <- 0.999
        QbarAW[QbarAW <= 0] <- 0.001

        Qbar1W[Qbar1W >= 1.0] <- 0.999
        Qbar1W[Qbar1W <= 0] <- 0.001

        Qbar0W[Qbar0W >= 1.0] <- 0.999
        Qbar0W[Qbar0W <= 0] <- 0.001

        Av_c$QbarAW <- QbarAW
        Av_c$Qbar1W <- Qbar1W
        Av_c$Qbar0W <- Qbar0W

        marginal_data[[i]] <- Av_c
      }
    } else {
      print(
        paste(
          "No ATEs calculated in the validation for",
          mix_comps[i],
          "due to no rule found in training set for marginal impact"
        )
      )
      marg_directions[i] <- NA
    }
  }

  return(list(
    data = marginal_data,
    directions = marg_directions
  ))
}
