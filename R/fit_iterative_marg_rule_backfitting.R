#' @title Iteratively back-fit a Super Learner on marginal mixture components
#'  and covariates
#'
#' @details Fit the semi-parametric additive model E(Y) = f(A) + h(W) where
#' f(A) is a Super Learner of decision trees applied to each mixture component
#' and h(W) is a Super Learner applied to the covariates. Each estimator is fit
#' offset by the predictions of the other until convergence where convergence
#' is essentially no difference between the model fits. If a partitioning set
#' is found in f(A) return the rules which are the data-adaptively identified
#' thresholds for the mixture component that maximize the between group
#' difference while controlling for covariates.
#'
#' @param mix_comps A vector of characters indicating variables for the
#' mixture components
#' @param at Training data
#' @param w A vector of characters indicating variables that are covariates
#' @param w_stack Stack of algorithms made in SL 3 used in ensemble machine
#' learning to fit Y|W
#' @param tree_stack Stack of algorithms made in SL for the decision tree
#' estimation
#' @param fold Current fold in the cross-validation
#' @param max_iter Max number of iterations of iterative backfitting algorithm
#' @param verbose Run in verbose setting
#' @param parallel_cv Parallelize the cross-validation (TRUE/FALSE)
#' @param seed Numeric, seed number for consistent results
#' @return A matrix with the rules for each mixture component
#' @import partykit
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect
#' @importFrom stats na.omit
#' @importFrom purrr discard
#'
#' @examples
#'data <- simulate_mixture_cube()
#'data$y_scaled <- data$y
#'mix_comps <- c("M1", "M2", "M3")
#'W <- c("w", "w2")
#' sls <- create_sls()
#' w_stack <- sls$W_stack
#' tree_stack <- sls$A_stack
#' example_output <- fit_marg_rule_backfitting(mix_comps = mix_comps,
#'                                                      at = data,
#'                                                      w = W,
#'                                                      w_stack = w_stack,
#'                                                      tree_stack = tree_stack,
#'                                                      fold = 1,
#'                                                      max_iter = 1,
#'                                                      verbose = FALSE,
#'                                                      parallel_cv = FALSE,
#'                                                      seed = set.seed())
#'
#' @export


fit_marg_rule_backfitting <- function(mix_comps,
                                                at,
                                                w,
                                                w_stack,
                                                tree_stack,
                                                fold,
                                                max_iter,
                                                verbose,
                                                parallel_cv,
                                                seed) {
  if (parallel_cv == TRUE) {
    future::plan(future::sequential, gc = TRUE)
  }

  set.seed(seed)

  marg_decisions <- list()

  at$Qbar_ne_M_W_initial <- 0
  at$Qbar_M_W_initial <- 0
  at$Qbar_ne_M_W_now <- 0
  at$Qbar_M_W_now <- 0

  for (i in seq(mix_comps)) {
    target_m <- mix_comps[i]
    covars_m <- c(mix_comps[-i], w)

    task <- sl3::make_sl3_Task(
      data = at,
      covariates = covars_m,
      outcome = "y_scaled",
      outcome_type = "continuous"
    )

    discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new(sl3::loss_squared_error)

    discrete_sl <- sl3::Lrnr_sl$new(
      learners = w_stack,
      metalearner = discrete_sl_metalrn,
    )

    sl_fit <- suppressWarnings(discrete_sl$train(task))

    qbar_ne_m_w_initial <- sl_fit$predict()

    at[, "Qbar_ne_M_W_initial"] <- qbar_ne_m_w_initial

    task <- sl3::make_sl3_Task(
      data = at,
      covariates = target_m,
      outcome = "y_scaled",
      outcome_type = "continuous"
    )

    tree_sl <- sl3::Lrnr_sl$new(
      learners = tree_stack,
      metalearner = discrete_sl_metalrn,
    )

    glmtree_fit <- tree_sl$train(task)

    qbar_m_w_initial <- glmtree_fit$predict()

    at[, "Qbar_M_W_initial"] <- qbar_m_w_initial

    iter <- 0
    stop <- FALSE

    at_no_offset <- data.table::copy(at)
    at_no_offset$Qbar_M_W_initial <- 0
    at_no_offset$Qbar_ne_M_W_initial <- 0

    while (stop == FALSE) {
      iter <- iter + 1

      task_offset <- sl3::sl3_Task$new(
        data = at,
        covariates = covars_m,
        outcome = "y_scaled",
        outcome_type = "continuous",
        offset = "Qbar_M_W_initial"
      )

      task_no_offset <- sl3::sl3_Task$new(
        data = at_no_offset,
        covariates = covars_m,
        outcome = "y_scaled",
        outcome_type = "continuous",
        offset = "Qbar_M_W_initial"
      )

      sl_fit_backfit_offset <- discrete_sl$train(task_offset)
      sl_fit_backfit_no_offset <- discrete_sl$train(task_no_offset)

      preds_offset <- sl_fit_backfit_offset$predict()
      preds_no_offset <- sl_fit_backfit_no_offset$predict()

      at[, "Qbar_ne_M_W_now"] <- preds_no_offset

      task <- sl3::make_sl3_Task(
        data = at,
        covariates = target_m,
        outcome = "y_scaled",
        outcome_type = "continuous",
        offset = "Qbar_ne_M_W_initial"
      )

      task_no_offset <- sl3::make_sl3_Task(
        data = at_no_offset,
        covariates = target_m,
        outcome = "y_scaled",
        outcome_type = "continuous",
        offset = "Qbar_ne_M_W_initial"
      )

      glmtree_fit_offset <- tree_sl$train(task)

      glmtree_model_preds_offset <- glmtree_fit_offset$predict(task)
      glmtree_model_preds_no_offset <- glmtree_fit_offset$predict(
        task_no_offset)

      at[, "Qbar_M_W_now"] <- glmtree_model_preds_no_offset

      curr_diff <- abs(glmtree_model_preds_offset - preds_offset)

      qbar_ne_m_w_initial <- at$Qbar_ne_M_W_now
      at$Qbar_ne_M_W_initial <- qbar_ne_m_w_initial
      at$Qbar_M_W_initial <- at$Qbar_M_W_now

      selected_learner <- glmtree_fit_offset$learner_fits[[
        which(glmtree_fit_offset$coefficients == 1)]]

      if (verbose) {
        if (iter == 1) {
          print(paste(
            "Fold: ", fold, "|",
            "Process: ", target_m, "Marginal Decision Backfitting", "|",
            "Iteration: ", iter, "|",
            "Delta: ", "None", "|",
            "Diff: ", mean(curr_diff), "|",
            "Rules:", list_rules_party(selected_learner$fit_object)
          ))
        } else {
          print(paste(
            "Fold: ", fold, "|",
            "Process: ", target_m, "Marginal Decision Backfitting", "|",
            "Iteration: ", iter, "|",
            "Delta: ", mean(curr_diff - prev_diff), "|",
            "Diff: ", mean(curr_diff), "|",
            "Rules:", list_rules_party(selected_learner$fit_object)
          ))
        }
      }

      if (iter == 1) {
        stop <- FALSE
        prev_diff <- curr_diff
      } else if (abs(mean(curr_diff - prev_diff)) <= 0.001) {
        stop <- TRUE
      } else if (iter >= max_iter) {
        stop <- TRUE
      } else {
        stop <- FALSE
        prev_diff <- curr_diff
      }
    }

    rules <- list_rules_party(selected_learner$fit_object)
    quantile <- seq(length(rules))


    if (length(rules) == 1) {
      if (rules == "") {
        rules <- "No Rules Found"
      }
    }
    rules <- as.data.frame(cbind(rules, fold, target_m, quantile))

    backfit_resids <- (at$y_scaled - glmtree_model_preds_offset)^2
    backfit_rmse <- sqrt(mean(backfit_resids))

    rules$RMSE <- backfit_rmse


    marg_decisions[[i]] <- rules
  }

  marg_decisions <- do.call(rbind, marg_decisions)

  return(marg_decisions)
}
