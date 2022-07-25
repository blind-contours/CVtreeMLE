#' @title Estimate the expected outcome for the combination of marginal
#' thresholds identified in the fold.
#' @description Estimate the expected outcome given exposure to the combination
#' of marginal exposures. This is different
#' compared to the cumulative sum; whereas with the cumulative sum, the exposure
#'  is the additive effect of each
#' marginal rule found in the fold, here each marginal rule is included in a
#' Super Learner as a binary vector and
#' therefore this can pick up possible nonlinearity between the combination of
#' binary exposures.
#'
#' @param at Training data
#' @param av Validation data
#' @param y Outcome variable
#' @param w Vector of characters denoting covariates
#' @param marg_rule_train Data frame of binary vectors for marginal rules
#' identified in the training fold
#' @param marg_rule_valid Data frame of binary vectors for marginal rules
#' identified in the validation fold
#' @param aw_stack Super Learner library for fitting Q (outcome mechanism) and
#' g (treatment mechanism)
#' @param family Outcome type family
#' @param no_marg_rules TRUE/FALSE if no marginal rules were found across all
#' @param seed Seed number
#' folds
#' @param parallel_cv TRUE/FALSE if parallel CV is used
#' @import sl3
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#' @importFrom dplyr group_by filter top_n
#' @return Rules object. TODO: add more detail here.
#'
#' @export

est_comb_exposure <- function(at,
                              av,
                              y = "y_scaled",
                              w,
                              marg_rule_train,
                              marg_rule_valid,
                              no_marg_rules,
                              aw_stack,
                              family,
                              parallel_cv,
                              seed) {
  if (parallel_cv == TRUE) {
    future::plan(future::sequential, gc = TRUE)
  }

  set.seed(seed)

  if (no_marg_rules == FALSE) {
    at_mc <- at
    av_mc <- av

    at_marg_comb <-
      cbind(marg_rule_train, at_mc[w], at_mc["y_scaled"])

    av_marg_comb <-
      cbind(marg_rule_valid, av_mc[w], av_mc["y_scaled"])

    task_at <- sl3::make_sl3_Task(
      data = at_marg_comb,
      covariates = c(colnames(marg_rule_train), w),
      outcome = "y_scaled",
      outcome_type = family
    )

    task_av <- sl3::make_sl3_Task(
      data = av_marg_comb,
      covariates = c(colnames(marg_rule_valid), w),
      outcome = "y_scaled",
      outcome_type = family
    )

    discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new(sl3::loss_squared_error)

    discrete_sl <- sl3::Lrnr_sl$new(
      learners = aw_stack,
      metalearner = discrete_sl_metalrn,
    )

    sl_fit <- suppressWarnings(discrete_sl$train(task_at))

    qbar_aw <- bound_precision(sl_fit$predict(task_av))
    qbar_aw <- scale_to_original(scaled_vals = qbar_aw,
                                 max_orig = max(at_mc[y]),
                                min_orig = min(at_mc[y]))

    av_marg_comb$qbar_aw_combo <- qbar_aw
    av_marg_comb$y_scaled <- av_mc$y_scaled
    av_marg_comb$raw_outcome <- av[, y]
  } else {
    av_marg_comb <- NA
    sl_fit <- NA
  }

  return(list("data" = av_marg_comb, "learner" = sl_fit))
}
