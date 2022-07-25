#' @title Estimate nuisance parameters for each marginal mixture component
#'
#' @description For each marginal mixture component rule found, create a
#' g estimator for the probability of being exposed to the rule thresholds,
#' and a Q estimator for the outcome E(Y| A = a_mix, W).
#' Get estimates of g and Q using the validation data and
#' calculate the clever covariate used in the TMLE fluctuation step.
#'
#' @param at Training data
#' @param av Validation data
#' @param w Vector of characters denoting covariates
#' @param aw_stack Super Learner library for fitting Q (outcome mechanism)
#' and g (treatment mechanism)
#' @param family Binomial or gaussian
#' @param a Vector of characters that denote the mixture components
#' @param marg_decisions List of rules found within the fold for each
#' mixture component
#' @param h_aw_trunc_lvl Truncation level of the clever covariate
#' (induces more bias to reduce variance)
#' @param no_marg_rules TRUE/FALSE if no marginal rules were
#' found across the folds
#' @param parallel_cv TRUE/FALSE if cv parallelization is used
#' @param seed Seed number
#' @import sl3
#' @importFrom magrittr %>%
#' @importFrom rlang :=

#' @importFrom dplyr group_by filter top_n
#' @return Rules object. TODO: add more detail here.

#'
#' @export

est_marg_nuisance_params <- function(at,
                                     av,
                                     w,
                                     aw_stack,
                                     family,
                                     a,
                                     no_marg_rules,
                                     marg_decisions,
                                     h_aw_trunc_lvl,
                                     parallel_cv,
                                     seed) {
  if (parallel_cv == TRUE) {
    future::plan(future::sequential, gc = TRUE)
  }

  set.seed(seed)

  marginal_data <- list()

  marg_decisions_groups <- marginal_group_split(marg_decisions)

  if (no_marg_rules == FALSE) {
    for (i in seq(marg_decisions_groups)) {
      at_c <- at
      av_c <- av

      variable_decisions <- marg_decisions_groups[[i]]

      quant_one_row <- variable_decisions[variable_decisions$quantile == 1, ]
      quant_one_rule <- quant_one_row$rules

      at_c_ref_data <- at_c %>%
        dplyr::mutate("a" := ifelse(eval(parse(text = quant_one_rule)), 1, 0))

      av_c_ref_data <- av_c %>%
        dplyr::mutate("a" := ifelse(eval(parse(text = quant_one_rule)), 1, 0))

      at_c_ref_data <- at_c_ref_data[at_c_ref_data[, "a"] == 1, ]
      av_c_ref_data <- av_c_ref_data[av_c_ref_data[, "a"] == 1, ]

      at_c_ref_data[, "a"] <- 0
      av_c_ref_data[, "a"] <- 0

      quant_comparisons <- variable_decisions[variable_decisions$quantile > 1, ]

      in_group_marg_data <- list()

      for (j in seq(nrow(quant_comparisons))) {
        target_m_row <- quant_comparisons[j, ]

        at_c_comp_data <- at_c %>%
          dplyr::mutate("a" := ifelse(eval(parse(text =
                                                   target_m_row$rules)), 1, 0))

        av_c_comp_data <- av_c %>%
          dplyr::mutate("a" := ifelse(eval(parse(text =
                                                   target_m_row$rules)), 1, 0))

        at_c_comp_data <- at_c_comp_data[at_c_comp_data[, "a"] == 1, ]
        av_c_comp_data <- at_c_comp_data[at_c_comp_data[, "a"] == 1, ]

        at_data <- rbind(at_c_ref_data, at_c_comp_data)
        av_data <- rbind(av_c_ref_data, av_c_comp_data)

        task_at <- sl3::make_sl3_Task(
          data = at_data,
          covariates = w,
          outcome = "a",
          outcome_type = "binomial",
          folds = 2
        )

        task_av <- sl3::make_sl3_Task(
          data = av_data,
          covariates = w,
          outcome = "a",
          outcome_type = "binomial"
        )

        discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new(
          sl3::loss_loglik_binomial)

        discrete_sl <- sl3::Lrnr_sl$new(
          learners = aw_stack,
          metalearner = discrete_sl_metalrn,
        )

        sl_fit <- suppressWarnings(discrete_sl$train(task_at))

        ghat_1w <- bound_precision(sl_fit$predict(task_av))

        h_aw <- calc_clever_covariate(ghat_1_w = ghat_1w,
                                      data = av_data,
                                      exposure = "a",
                                      h_aw_trunc_lvl = 10,
                                      type = "reg")

        av_data$ghat_1w <- ghat_1w
        av_data$h_aw <- h_aw

        task_at <- sl3::make_sl3_Task(
          data = at_data,
          covariates = c(w, "a"),
          outcome = "y_scaled",
          outcome_type = family,
          folds = 2
        )

        x_m1 <- x_m0 <- av_data
        x_m1$a <- 1 # under exposure
        x_m0$a <- 0 # under control

        task_av <- sl3::make_sl3_Task(
          data = av_data,
          covariates = c(w, "a"),
          outcome = "y_scaled",
          outcome_type = family
        )

        task_av_1 <- sl3::make_sl3_Task(
          data = x_m1,
          covariates = c(w, "a"),
          outcome = "y_scaled",
          outcome_type = family
        )

        task_av_0 <- sl3::make_sl3_Task(
          data = x_m0,
          covariates = c(w, "a"),
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

        q_bar_aw <- bound_precision(sl_fit$predict(task_av))
        q_bar_1w <- bound_precision(sl_fit$predict(task_av_1))
        q_bar_0w <- bound_precision(sl_fit$predict(task_av_0))

        av_data$qbar_aw <- q_bar_aw
        av_data$qbar_1w <- q_bar_1w
        av_data$qbar_0w <- q_bar_0w

        in_group_marg_data[[j]] <- av_data
      }
      marginal_data[[i]] <- in_group_marg_data
    }

    marginal_data <- unlist(marginal_data, recursive = FALSE, use.names = FALSE)
  } else {
    av$ghat_1w <- NA
    av$h_aw <- NA
    av$qbar_aw <- NA
    av$qbar_1w <- NA
    av$qbar_0w <- NA
    marginal_data[[1]] <- av
  }

  return(marginal_data)
}
