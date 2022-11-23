#' @title Estimate the NDE nuisance parameters for each mixture
#'  interaction identified
#'
#' @description For each mixture mixture interaction found, create a g estimator
#' for the probability of being exposed to the rule thresholds,
#'a Q estimator for the outcome E(Y| A = a_mix, Z, W), an e estimator for
#' e(A|Z,W) and a Phi estimator that regresses the outcome counterfactuals on
#' the covariates of A =0.  Get estimates of g, e,phi, and Q using the
#' validation data and calculate the clever covariate used in the
#' TMLE fluctuation step.
#'
#' @param at Training data
#' @param av Validation data
#' @param w Vector of characters denoting covariates
#' @param z Vector of characters denoting mediators
#' @param y The outcome variable
#' @param no_mix_rules TRUE/FALSE indicator for if no mixture rules were found
#'
#' @param aw_stack Super Learner library for fitting Q (outcome mechanism) and
#' g (treatment mechanism)
#' @param psi_z_stack Super Learner library for fitting Q_diff (difference
#' in conditional outcomes for A=1 and A = 0, and Z = z, W = w) onto W of
#' controls observations
#' @param family Binomial or continuous
#' @param rules Dataframe of rules found during the PRE fitting process
#' @param h_aw_trunc_lvl Truncation level of the clever covariate (induces more
#'  bias to reduce variance)
#' @param parallel_cv TRUE/FALSE if cv parallelization is used
#' @param seed Seed number
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by filter top_n
#' @return A list of dataframes where the nuisance parameters are added to
#' the raw data.

#'
#' @export

est_mix_nuisance_params_nde <- function(at,
                                        av,
                                        w,
                                        z,
                                        y,
                                        no_mix_rules,
                                        aw_stack,
                                        psi_z_stack,
                                        family,
                                        rules,
                                        h_aw_trunc_lvl = 10,
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

        # get g mech estimates ---------------------------

        task_at_g <- sl3::make_sl3_Task(
          data = at_mix,
          covariates = w,
          outcome = "A_mix",
          outcome_type = "binomial",
          folds = 2
        )

        task_av_g <- sl3::make_sl3_Task(
          data = av_mix,
          covariates = w,
          outcome = "A_mix",
          outcome_type = "binomial"
        )

        sl <- sl3::Lrnr_sl$new(learners = aw_stack,
                               metalearner = sl3::Lrnr_nnls$new())

        sl_fit_g <- suppressWarnings(sl$train(task_at_g))

        at_g1_est <- sl_fit_g$predict(task_at_g)
        at_g0_est <- 1 - at_g1_est

        av_g1_est <- sl_fit_g$predict(task_av_g)
        av_g0_est <- 1 - av_g1_est

        # get e mech estimates ---------------------------

        task_at_e <- sl3::make_sl3_Task(
          data = at_mix,
          covariates = c(w,z),
          outcome = "A_mix",
          outcome_type = "binomial",
          folds = 2
        )

        task_av_e <- sl3::make_sl3_Task(
          data = av_mix,
          covariates = c(w,z),
          outcome = "A_mix",
          outcome_type = "binomial"
        )

        sl_fit_e <- suppressWarnings(sl$train(task_at_e))

        at_e_est <- sl_fit_e$predict(task_at_e)

        # calculate training clever covariate  ---------------------------

        at_treatment_indicator <- at_mix$A_mix
        at_control_indicator <- 1 - at_mix$A_mix

        at_e1_est <- at_e_est
        at_e0_est <- 1 - at_e1_est

        at_hy <- (at_treatment_indicator / at_g1_est) * ((at_e0_est/at_g0_est)*(at_g1_est/at_e1_est)) - (at_control_indicator / at_g0_est)
        at_hz <- 1 / at_g0_est

        at_mix$hy <- at_hy
        at_mix$hz <- at_hz

        # calculate validation clever covariate  ---------------------------

        av_e_est <- sl_fit_e$predict(task_av_e)

        av_treatment_indicator <- av_mix$A_mix
        av_control_indicator <- 1 - av_mix$A_mix

        av_e1_est <- av_e_est
        av_e0_est <- 1 - av_e1_est

        av_hy <- (av_treatment_indicator / av_g1_est) * ((av_e0_est/av_g0_est)*(av_g1_est/av_e1_est)) - (av_control_indicator / av_g0_est)
        av_hz <- 1 / av_g0_est

        av_mix$hy <- av_hy
        av_mix$hz <- av_hz

        # calculate Qy Initial ---------------------------

        at_1w <- at_0w <- at_mix
        at_1w$A_mix <- 1 # under exposure
        at_0w$A_mix <- 0 # under control

        task_at <- sl3::make_sl3_Task(
          data = at_mix,
          covariates = c(w, z, "A_mix"),
          outcome = y,
          outcome_type = family,
          folds = 2
        )

        task_at_1 <- sl3::make_sl3_Task(
          data = at_1w,
          covariates = c(w, z, "A_mix"),
          outcome = y,
          outcome_type = family,
          folds = 2
        )

        task_at_0 <- sl3::make_sl3_Task(
          data = at_0w,
          covariates = c(w, z, "A_mix"),
          outcome = y,
          outcome_type = family,
          folds = 2
        )

        av_1w <- av_0w <- av_mix

        av_1w$A_mix <- 1 # under exposure
        av_0w$A_mix <- 0 # under control

        task_av <- sl3::make_sl3_Task(
          data = av_mix,
          covariates = c(w, z, "A_mix"),
          outcome = y,
          outcome_type = family
        )

        task_av_1 <- sl3::make_sl3_Task(
          data = av_1w,
          covariates = c(w, z, "A_mix"),
          outcome = y,
          outcome_type = family
        )

        task_av_0 <- sl3::make_sl3_Task(
          data = av_0w,
          covariates = c(w, z, "A_mix"),
          outcome = y,
          outcome_type = family
        )

        sl_fit_y <- suppressWarnings(sl$train(task_at))

        # calculate Qz Initial ---------------------------

        at_mix$qbar_azw <- sl_fit_y$predict(task_at)
        at_mix$qbar_1zw <- sl_fit_y$predict(task_at_1)
        at_mix$qbar_0zw <- sl_fit_y$predict(task_at_0)

        at_Q_y_tmle_udpates <- fit_least_fav_submodel(h_aw = at_hy,
                                                   data = at_mix,
                                                   y = y,
                                                   qbar_aw = at_mix$qbar_azw,
                                                   qbar_1w = at_mix$qbar_1zw,
                                                   qbar_0w = at_mix$qbar_0zw)


        at_q_diff <- at_Q_y_tmle_udpates$qbar_1w_star - at_Q_y_tmle_udpates$qbar_0w_star
        at_mix$q_diff <- at_q_diff

        av_mix$qbar_azw <- sl_fit_y$predict(task_av)
        av_mix$qbar_1zw <- sl_fit_y$predict(task_av_1)
        av_mix$qbar_0zw <- sl_fit_y$predict(task_av_0)

        av_Q_y_tmle_udpates <- fit_least_fav_submodel(h_aw = av_hy,
                                                      data = av_mix,
                                                      y = y,
                                                      qbar_aw = av_mix$qbar_azw,
                                                      qbar_1w = av_mix$qbar_1zw,
                                                      qbar_0w = av_mix$qbar_0zw)


        av_q_diff <- av_Q_y_tmle_udpates$qbar_1w_star - av_Q_y_tmle_udpates$qbar_0w_star
        av_mix$q_diff <- av_q_diff

        av_mix$qbar_azw_star <- av_Q_y_tmle_udpates$qbar_aw_star
        av_mix$qbar_1zw_star <- av_Q_y_tmle_udpates$qbar_1w_star
        av_mix$qbar_0zw_star <- av_Q_y_tmle_udpates$qbar_0w_star

        at_control_w_data <- at_mix %>% dplyr::filter(A_mix == 0)
        # av_control_w_data <- av_mix %>% filter(A_mix == 0)

        task_at <- sl3::make_sl3_Task(
          data = at_control_w_data,
          covariates = w,
          outcome = "q_diff",
          outcome_type = family,
          folds = 2
        )
        task_av <- sl3::make_sl3_Task(
          data = av_mix,
          covariates = w,
          outcome = "q_diff",
          outcome_type = family,
          folds = 2
        )

        sl_fit_diff_w <- suppressWarnings(sl$train(task_at))
        qdiff_w <- sl_fit_diff_w$predict(task_av)

        av_mix$qdiff_w <- qdiff_w

        logit_update <-
          stats::glm(
            q_diff ~ -1 + hz + offset(qdiff_w),
            data = av_mix
          )

        epsilon <- logit_update$coef
        qbar_diff_star <- qdiff_w + epsilon * av_mix$hz

        av_mix$psi_z <- qbar_diff_star
        # compute individual scores for DY, DA, DW
        d_y <- av_hy * (av_mix[,y] - av_Q_y_tmle_udpates$qbar_aw_star)
        d_z <- av_control_indicator * av_hz * (av_Q_y_tmle_udpates$qbar_1w_star - av_Q_y_tmle_udpates$qbar_0w_star - qbar_diff_star)
        d_w <- qbar_diff_star

        # parameter and influence function
        theta <- mean(qbar_diff_star)
        eif <- d_y + d_z + d_w - theta

        av_mix$theta <- theta
        av_mix$eif <- eif

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
