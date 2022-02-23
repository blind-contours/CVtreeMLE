#' Calculate the ATE and variance estimates
#'
#' @details Calculates the average treatment effect and variance estimates from the influence curve across folds

#' @param data Dataframe that includes treatement/exposure
#' @param ATE_var Variable that the ATE is assigned to in `data`
#' @param outcome Variable name in `data` for the raw outcome
#' @param p_adjust_n Number repeated measures to adjust p-value
#' @param v_fold TRUE/FALSE whether to calculate v-fold specific ATE estimates
#'
#' @return A \code{numeric} vector of the same length as \code{vals}, where
#'  the returned values are bounded to machine precision. This is intended to
#'  avoid numerical instability issues.
calc_ATE_estimates <- function(data, ATE_var, outcome, p_adjust_n, v_fold = FALSE) {
  # assertthat::assert_that(!(max(vals) > 1 | min(vals) < 0))
  data[ATE_var] <- data$Qbar1W.star - data$Qbar0W.star

  Thetas <-
    tapply(data[[ATE_var]], data$folds, mean, na.rm = TRUE)

  if (v_fold == TRUE) {
    data[, "Thetas"] <- Thetas
  } else {
    for (i in 1:length(Thetas)) {
      fold <- names(Thetas)[i]
      data[data$folds == fold, "Thetas"] <- Thetas[i][[1]]
    }
  }

  ICs <- base::by(data, data$folds, function(data) {
    result <- data["H.AW"] * (data[outcome] - data["QbarAW.star"]) + data["Qbar1W.star"] - data["Qbar0W.star"] - data["Thetas"]
    result
  })

  if (v_fold == TRUE) {
    data[, "IC"] <- ICs
  } else {
    for (i in seq(ICs)) {
      fold <- names(Thetas)[i]
      data[data$folds == fold, "IC"] <- ICs[fold]
    }
  }

  n <- dim(data)[1]
  varHat.IC <- stats::var(data$IC, na.rm = TRUE) / n
  se <- sqrt(varHat.IC)

  alpha <- 0.05

  Theta <- mean(Thetas)
  # obtain 95% two-sided confidence intervals:
  CI <- c(
    Theta + stats::qnorm(alpha / 2, lower.tail = T) * se,
    Theta + stats::qnorm(alpha / 2, lower.tail = F) * se
  )

  # p-value
  p.value <- 2 * stats::pnorm(abs(Theta / se), lower.tail = F)

  p.value.adjust <-
    stats::p.adjust(p.value, method = "bonferroni", n = p_adjust_n)

  return(list("ATE" = Theta, "SE" = se, "CI" = CI, "p-value" = p.value, "adj p-value" = p.value.adjust, "data" = data))
}

###############################################################################

#' Calculate the Clever Covariate for the TMLE step of the ATE
#'
#' @details Calculates the clever covariate or difference in inverse propensity scores for treatment used in the TMLE update step
#'
#' @param gHat1W \code{Numeric} vector of values in the unit interval that represent the predicted probabilities of being treated/exposed.
#' @param data Dataframe that includes treatement/exposure
#' @param exposure Variable name in `data` that the clever covariate is calculated for
#' @param H.AW_trunc_lvl Truncation level of the clever covariate
#' @param type Type of clever covariate to calculate, `reg` is standard difference in inverse weightings, `mod` is used for the additive ATE
#'
#' @return A \code{numeric} vector of the same length as \code{vals}, where
#'  the returned values are bounded to machine precision. This is intended to
#'  avoid numerical instability issues.
calc_clever_covariate <- function(gHat1W, data, exposure, H.AW_trunc_lvl, type = "reg") {
  # assertthat::assert_that(!(max(vals) > 1 | min(vals) < 0))
  if (type == "reg") {
    n <- length(gHat1W)
    gHat0W <- 1 - gHat1W

    gHatAW <- rep(NA, n)
    gHatAW[data[exposure] == 1] <- gHat1W[data[exposure] == 1]
    gHatAW[data[exposure] == 0] <- gHat0W[data[exposure] == 0]

    H.AW <-
      as.numeric(data[exposure] == 1) / gHat1W - as.numeric(data[exposure] == 0) /
        gHat0W

    H.AW <-
      ifelse(H.AW > H.AW_trunc_lvl, H.AW_trunc_lvl, H.AW)

    H.AW <-
      ifelse(H.AW < -H.AW_trunc_lvl, -H.AW_trunc_lvl, H.AW)
  } else {
    n <- length(gHat1W)

    gHatAW <- rep(NA, n)
    gHatAW[data[exposure] == 1] <- gHat1W[data[exposure] == 1]

    H.AW <-
      as.numeric(data[exposure] == 1) / gHat1W

    H.AW <-
      ifelse(H.AW > H.AW_trunc_lvl, H.AW_trunc_lvl, H.AW)

    H.AW <-
      ifelse(H.AW < -H.AW_trunc_lvl, -H.AW_trunc_lvl, H.AW)
  }
  return(H.AW)
}

###############################################################################

#' Least Favorable Submodel
#'
#' @details Does flucutation of initial outcome estimates using a least favorable submodel in the TMLE update step
#'
#' @param H.AW \code{numeric} vector of values for the clever covariate used in the fluctuation model
#' @param data dataframe that includes treatement/exposure
#' @param QbarAW Initial predictions for Y|A,W - A being observed A
#' @param Qbar1W Initial predictions for Y|A = 1,W - A being deterministically set to 1
#' @param Qbar0W Initial predictions for Y|A = 1,W - A being deterministically set to 0

#' @return A \code{numeric} vector of the same length as \code{vals}, where
#'  the returned values are bounded to machine precision. This is intended to
#'  avoid numerical instability issues.
fit_least_fav_submodel <- function(H.AW, data, QbarAW, Qbar1W, Qbar0W) {
  # assertthat::assert_that(!(max(vals) > 1 | min(vals) < 0))
  logitUpdate <-
    stats::glm(
      y_scaled ~ -1 + H.AW + offset(qlogis(bound_precision(QbarAW))),
      family = "quasibinomial",
      data = data
    )

  epsilon <- logitUpdate$coef
  QbarAW.star <- plogis(qlogis(bound_precision(QbarAW)) + epsilon * H.AW)
  Qbar1W.star <- plogis(qlogis(bound_precision(Qbar1W)) + epsilon * H.AW)
  Qbar0W.star <- plogis(qlogis(bound_precision(Qbar0W)) + epsilon * H.AW)

  return(list("QbarAW.star" = QbarAW.star, "Qbar1W.star" = Qbar1W.star, "Qbar0W.star" = Qbar0W.star))
}
