#' @title Calculate the ATE and variance estimates
#'
#' @details Calculates the average treatment effect and variance estimates
#' from the influence curve across folds

#' @param data Data frame that includes treatement/exposure
#' @param ate_var Variable that the ATE is assigned to in `data`
#' @param y Variable name in `data` for the raw outcome
#' @param p_adjust_n Number repeated measures to adjust p-value
#' @param v_fold TRUE/FALSE whether to calculate v-fold specific ATE estimates
#' @param naive This is for a one step estimator so Psi is calculated without TMLE update
#'
#' @return A \code{numeric} vector of the same length as \code{vals}, where
#'  the returned values are bounded to machine precision. This is intended to
#'  avoid numerical instability issues.
#' @export
calc_ate_estimates <- function(data,
                               ate_var,
                               y,
                               p_adjust_n,
                               v_fold = FALSE,
                               naive = FALSE) {
  if (naive == TRUE) {
    data[, ate_var] <- data$qbar_1w - data[, y]
  } else {
    data[, ate_var] <- data$qbar_1w_star - data[, y]
  }


  if (v_fold == TRUE) {
    thetas <-
      tapply(data[[ate_var]], data$folds, mean, na.rm = TRUE)
    data[, "thetas"] <- thetas
  } else {
    data[, "thetas"] <- mean(data[, ate_var])
  }

  theta <- mean(data[, "thetas"])

  if (v_fold == TRUE) {
    ics <- base::by(data, data$folds, function(data) {
      result <- (data[, "h_aw"] * (data[, y] - data[, "qbar_aw_star"])) +
        (data[, "qbar_1w_star"] - mean(data[, "qbar_1w_star"])) # -
       #(data[, y] - mean(data[, y]))

    })
  } else {
    if (naive == TRUE) {
      ics <- (data[, "h_aw"] * (data[, y] - data[, "qbar_aw_star"])) +
        (data[, "qbar_1w_star"] - mean(data[, "qbar_1w_star"]))  #-
        #(data[, y] - mean(data[, y]))
    } else {
      ics <- (data[, "h_aw"] * (data[, y] - data[, "qbar_aw_star"])) +
        (data[, "qbar_1w_star"] - mean(data[, "qbar_1w_star"])) # -
        #(data[, y] - mean(data[, y]))
    }
  }


  data[, "ic"] <- ics

  n <- dim(data)[1]
  varhat_ic <- stats::var(data$ic, na.rm = TRUE) / n
  se <- sqrt(varhat_ic)

  alpha <- 0.05

  # obtain 95% two-sided confidence intervals:
  ci <- c(
    theta + stats::qnorm(alpha / 2, lower.tail = TRUE) * se,
    theta + stats::qnorm(alpha / 2, lower.tail = FALSE) * se
  )

  # p-value
  p_value <- 2 * stats::pnorm(abs(theta / se), lower.tail = FALSE)

  if (!is.null(p_adjust_n)) {
    p_value_adjust <-
      stats::p.adjust(p_value, method = "bonferroni", n = p_adjust_n)
  } else {
    p_value_adjust <- NULL
  }

  return(list(
    "ate" = theta, "se" = se, "ci" = ci, "p_value" = p_value,
    "adj_p_value" = p_value_adjust, "data" = data
  ))
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
    # n <- sum(data$A_mix) / nrow(data)

    h_aw <- 1 / ghat_1_w

    h_aw <-
      ifelse(h_aw > h_aw_trunc_lvl, h_aw_trunc_lvl, h_aw)

    h_aw <-
      ifelse(h_aw < -h_aw_trunc_lvl, -h_aw_trunc_lvl, h_aw)
  } else {
    n <- sum(data$A_mix) / nrow(data)

    h_aw <- 1 / ghat_1_w

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
#' @param y Outcome variable
#' used in the fluctuation model
#' @param data dataframe that includes treatement/exposure
#' @param qbar_aw Initial predictions for Y|A,W - A being observed A
#' @param qbar_1w Initial predictions for Y|A = 1,W - A being deterministically
#' set to 1

#' @return A \code{numeric} vector of the same length as \code{vals}, where
#'  the returned values are bounded to machine precision. This is intended to
#'  avoid numerical instability issues.
#' @export
fit_least_fav_submodel <- function(h_aw, data, y, qbar_aw, qbar_1w) {
  is_binary <- all(data[[y]] %in% c(0, 1))

  # Apply scaling for continuous outcomes
  if (!is_binary) {
    data$y_scaled <- scale_to_unit(data[[y]])
    qbar_aw_scaled <- scale_to_unit(qbar_aw)
    qbar_1w_scaled <- scale_to_unit(qbar_1w)
  } else {
    data$y_scaled <- data[[y]]
    qbar_aw_scaled <- qbar_aw
    qbar_1w_scaled <- qbar_1w
  }

  logit_update <- glm(
    y_scaled ~ -1 + h_aw,
    family = "quasibinomial",
    data = data,
    weights = qbar_aw_scaled
  )

  epsilon <- coef(logit_update)

  # Apply update step based on the type of outcome
  if (is_binary) {
    qbar_aw_star <- 1 / (1 + exp(-(log(qbar_aw / (1 - qbar_aw)) + epsilon * h_aw)))
    qbar_1w_star <- 1 / (1 + exp(-(log(qbar_1w / (1 - qbar_1w)) + epsilon * h_aw)))
  } else {
    qbar_aw_star <- qbar_aw + epsilon * h_aw
    qbar_1w_star <- qbar_1w + epsilon * h_aw
  }

  return(list(
    "qbar_aw_star" = qbar_aw_star,
    "qbar_1w_star" = qbar_1w_star
  ))
}
