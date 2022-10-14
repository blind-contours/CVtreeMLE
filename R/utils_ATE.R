#' @title Calculate the ATE and variance estimates
#'
#' @details Calculates the average treatment effect and variance estimates
#' from the influence curve across folds

#' @param data Data frame that includes treatement/exposure
#' @param ate_var Variable that the ATE is assigned to in `data`
#' @param outcome Variable name in `data` for the raw outcome
#' @param p_adjust_n Number repeated measures to adjust p-value
#' @param v_fold TRUE/FALSE whether to calculate v-fold specific ATE estimates
#'
#' @return A \code{numeric} vector of the same length as \code{vals}, where
#'  the returned values are bounded to machine precision. This is intended to
#'  avoid numerical instability issues.
#' @export
calc_ate_estimates <- function(data,
                               ate_var,
                               outcome,
                               p_adjust_n,
                               v_fold = FALSE) {

  data[ate_var] <- data$qbar_1w_star - data$qbar_0w_star

  thetas <-
    tapply(data[[ate_var]], data$folds, mean, na.rm = TRUE)

  if (v_fold == TRUE) {
    data[, "thetas"] <- thetas
  } else {
    for (i in seq_along(thetas)) {
      fold <- names(thetas)[i]
      data[data$folds == fold, "thetas"] <- thetas[i][[1]]
    }
  }

  ics <- base::by(data, data$folds, function(data) {
    result <- data["h_aw"] * (data[outcome] - data["qbar_aw_star"]) +
      data["qbar_1w_star"] - data["qbar_0w_star"] - data["thetas"]
    result
  })

  if (v_fold == TRUE) {
    data[, "ic"] <- ics
  } else {
    for (i in seq(ics)) {
      fold <- names(thetas)[i]
      data[data$folds == fold, "ic"] <- ics[fold]
    }
  }

  n <- dim(data)[1]
  varhat_ic <- stats::var(data$ic, na.rm = TRUE) / n
  se <- sqrt(varhat_ic)

  alpha <- 0.05

  theta <- mean(thetas)
  # obtain 95% two-sided confidence intervals:
  ci <- c(
    theta + stats::qnorm(alpha / 2, lower.tail = TRUE) * se,
    theta + stats::qnorm(alpha / 2, lower.tail = FALSE) * se
  )

  # p-value
  p_value <- 2 * stats::pnorm(abs(theta / se), lower.tail = FALSE)

  p_value_adjust <-
    stats::p.adjust(p_value, method = "bonferroni", n = p_adjust_n)

  return(list("ate" = theta, "se" = se, "ci" = ci, "p_value" = p_value,
              "adj_p_value" = p_value_adjust, "data" = data))
}

###############################################################################

#' @title Calculate the Clever Covariate for the TMLE step of the ATE
#'
#' @details Calculates the clever covariate or difference in inverse propensity
#' scores for treatment used in the TMLE update step
#'
#' @param ghat_1_w \code{Numeric} vector of values in the unit interval that
#' represent the predicted probabilities of being treated/exposed.
#' @param data Dataframe that includes treatement/exposure
#' @param exposure Variable name in `data` that the clever covariate is
#' calculated for
#' @param h_aw_trunc_lvl Truncation level of the clever covariate
#' @param type Type of clever covariate to calculate, `reg` is standard
#' difference in inverse weightings, `mod` is used for the additive ATE
#'
#' @return A \code{numeric} vector of the same length as \code{vals}, where
#'  the returned values are bounded to machine precision. This is intended to
#'  avoid numerical instability issues.
#' @export
calc_clever_covariate <- function(ghat_1_w,
                                  data,
                                  exposure,
                                  h_aw_trunc_lvl,
                                  type = "reg") {
  if (type == "reg") {
    n <- length(ghat_1_w)

    ghat_1_w[ghat_1_w < 0.0001] <- 0.001
    ghat_1_w[ghat_1_w > 0.999] <- 0.99

    ghat_0_w <- 1 - ghat_1_w

    ghat_aw <- rep(NA, n)
    ghat_aw[data[exposure] == 1] <- ghat_1_w[data[exposure] == 1]
    ghat_aw[data[exposure] == 0] <- ghat_0_w[data[exposure] == 0]

    h_aw <-
      as.numeric(data[exposure] == 1) / ghat_1_w -
      as.numeric(data[exposure] == 0) / ghat_0_w

    h_aw <-
      ifelse(h_aw > h_aw_trunc_lvl, h_aw_trunc_lvl, h_aw)

    h_aw <-
      ifelse(h_aw < -h_aw_trunc_lvl, -h_aw_trunc_lvl, h_aw)
  } else {
    n <- length(ghat_1_w)

    ghat_aw <- rep(NA, n)
    ghat_aw[data[exposure] == 1] <- ghat_1_w[data[exposure] == 1]

    h_aw <-
      as.numeric(data[exposure] == 1) / ghat_1_w

    h_aw <-
      ifelse(h_aw > h_aw_trunc_lvl, h_aw_trunc_lvl, h_aw)

    h_aw <-
      ifelse(h_aw < -h_aw_trunc_lvl, -h_aw_trunc_lvl, h_aw)
  }
  return(h_aw)
}

###############################################################################

#' @title Least Favorable Submodel
#'
#' @details Does flucutation of initial outcome estimates using a least
#' favorable submodel in the TMLE update step
#'
#' @param h_aw \code{numeric} vector of values for the clever covariate
#' used in the fluctuation model
#' @param data dataframe that includes treatement/exposure
#' @param qbar_aw Initial predictions for Y|A,W - A being observed A
#' @param qbar_1w Initial predictions for Y|A = 1,W - A being deterministically
#' set to 1
#' @param qbar_0w Initial predictions for Y|A = 1,W - A being deterministically
#' set to 0

#' @return A \code{numeric} vector of the same length as \code{vals}, where
#'  the returned values are bounded to machine precision. This is intended to
#'  avoid numerical instability issues.
#' @export
fit_least_fav_submodel <- function(h_aw, data, qbar_aw, qbar_1w, qbar_0w) {

  data$y_scaled <- scale_to_unit(data[y])[[1]]

  logit_update <-
    stats::glm(
      y_scaled ~ -1 + h_aw + offset(qlogis(bound_precision(scale_to_unit(qbar_aw)))),
      family = "quasibinomial",
      data = data
    )

  epsilon <- logit_update$coef
  qbar_aw_star <- qbar_aw + epsilon * h_aw
  qbar_1w_star <- qbar_1w + epsilon * h_aw
  qbar_0w_star <- qbar_0w + epsilon * h_aw

  return(list("qbar_aw_star" = qbar_aw_star,
              "qbar_1w_star" = qbar_1w_star,
              "qbar_0w_star" = qbar_0w_star))
}
