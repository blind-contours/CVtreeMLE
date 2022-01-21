#' Iteratively Backfit a Super Learner, h(x) = Y|M_neX,W, and Indiviidual Decision Algorithms, g(x), for each Y|M_x until Convergence.
#'
#' @details Same procedure as `fit_iterative_mix_rule_backfitting` but the Super Learner controls for W and other mixture variables not equal to the mixture variable of interest, `glmtree` is used to find partitions
#' instead of `pre`.
#'
#' @param mix_comps A vector of characters indicating variables for the mixture components
#' @param At Training data
#' @param W A vector of characters indicating variables that are covariates#'
#' @param Q1_stack Stack of algorithms made in SL 3 used in ensemble machine learning to fit Y|W
#' @param tree_SL Stack of algorithms made in SL for the decision tree estimation
#' @param fold Current fold in the cross-validation
#' @param max_iter Max number of iterations of iterative backfitting algorithm
#' @param verbose Run in verbose setting
#' @return Rules object. TODO: add more detail here.
#' @import partykit
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect
#' @importFrom stats na.omit
#' @importFrom purrr discard
#'
#' @export

fit_iterative_marg_rule_backfitting <- function(mix_comps,
                                                At,
                                                W,
                                                Q1_stack,
                                                tree_SL,
                                                fold,
                                                max_iter,
                                                verbose) {

  future::plan(future::sequential, gc = TRUE)

  n <- dim(At)[1]
  marg_decisions <- list()

  At$Qbar_ne_M_W_initial <- 0
  At$Qbar_M_W_initial <- 0
  At$Qbar_ne_M_W_now <- 0
  At$Qbar_M_W_now <- 0

  for (i in seq(mix_comps)) {
    target_m <- mix_comps[i]
    covars_m <- c(mix_comps[-i], W)

    # Y <- At %>% dplyr::select(target_m)
    # X_m_part <- At %>% dplyr::select(covariates, covars_m)

    # binary_mix_var <- all(na.omit(Y[[1]]) %in% 0:1)

    task <- sl3::make_sl3_Task(
      data = At,
      covariates = covars_m,
      outcome = "y_scaled",
      outcome_type = "continuous"
    )

    discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new()

    discrete_sl <- sl3::Lrnr_sl$new(
      learners = Q1_stack,
      metalearner = discrete_sl_metalrn,
    )
    # delayed_sl_fit <- delayed_learner_train(discrete_sl, task)

    sl_fit <- discrete_sl$train(task)

    Qbar_ne_M_W_initial <- sl_fit$predict()

    At[,"Qbar_ne_M_W_initial"] <- Qbar_ne_M_W_initial

    task <- sl3::make_sl3_Task(
      data = At,
      covariates = target_m,
      outcome = "y_scaled",
      outcome_type = "continuous"
    )

    # formula <-
      # as.formula(paste("y_scaled", "~", target_m))

    ctree_fit <- tree_SL$train(task)
    # ctree_fit <- partykit::glmtree(formula,
    #                      data = At,
    #                      alpha = alpha,
    #                      prune = "AIC",
    #                      minsize = minsize,
    #                      maxdepth = max_depth)

    Qbar_M_W_initial <- ctree_fit$predict()

    At[,"Qbar_M_W_initial"] <- Qbar_M_W_initial

    iter <- 0
    stop <- FALSE

    At_no_offset <- data.table::copy(At)
    At_no_offset$Qbar_M_W_initial <- 0
    At_no_offset$Qbar_ne_M_W_initial <- 0

    while (stop == FALSE) {
      iter <- iter + 1

      task_offset <- sl3::sl3_Task$new(
        data = At,
        covariates = covars_m,
        outcome = "y_scaled",
        outcome_type = "continuous",
        offset = "Qbar_M_W_initial"
      )

      task_no_offset <- sl3::sl3_Task$new(
        data = At_no_offset,
        covariates = covars_m,
        outcome = "y_scaled",
        outcome_type = "continuous",
        offset = "Qbar_M_W_initial"
      )

      sl_fit_backfit_offset <- discrete_sl$train(task_offset)
      sl_fit_backfit_no_offset <- discrete_sl$train(task_no_offset)


      # preds_offset <- sl_fit_backfit$predict()
      preds_offset <- sl_fit_backfit_offset$predict()
      preds_no_offset <- sl_fit_backfit_no_offset$predict()

      At[,"Qbar_ne_M_W_now"] <- preds_no_offset

      task <- sl3::make_sl3_Task(
        data = At,
        covariates = target_m,
        outcome = "y_scaled",
        outcome_type = "continuous",
        offset = "Qbar_ne_M_W_initial"
      )

      task_no_offset <- sl3::make_sl3_Task(
        data = At_no_offset,
        covariates = target_m,
        outcome = "y_scaled",
        outcome_type = "continuous",
        offset = "Qbar_ne_M_W_initial"
      )

      ctree_fit <- tree_SL$train(task)

      glmtree_model_preds_offset <- ctree_fit$predict(task)
      At[, "Qbar_M_W_now"] <- ctree_fit$predict(task_no_offset)

      # delta_h <- mean(abs(At$Qbar_ne_M_W_initial - At$Qbar_ne_M_W_now))
      # delta_g <- mean(abs(At$Qbar_M_W_initial - At$Qbar_M_W_now))

      curr_diff <- abs(glmtree_model_preds_offset - preds_offset)

      Qbar_ne_M_W_initial <-  At$Qbar_ne_M_W_now
      At$Qbar_ne_M_W_initial <- Qbar_ne_M_W_initial
      At$Qbar_M_W_initial <- At$Qbar_M_W_now

      selected_learner <- ctree_fit$learner_fits[[which(ctree_fit$coefficients == 1)]]

      if (verbose){
        if (iter == 1) {
          print(paste("Fold: ", fold, "|",
                      "Process: ", target_m,  "Marginal Decision Backfitting", "|",
                      "Iteration: ", iter, "|",
                      "Delta: ", "None", "|",
                      "Diff: ", mean(curr_diff), "|",
                      "Rules:", list.rules.party(selected_learner$fit_object)))
        }else{
          print(paste("Fold: ", fold, "|",
                      "Process: ", target_m, "Marginal Decision Backfitting", "|",
                      "Iteration: ", iter, "|",
                      "Delta: ", mean(curr_diff - prev_diff), "|",
                      "Diff: ", mean(curr_diff), "|",
                      "Rules:", list.rules.party(selected_learner$fit_object)))
        }
      }

      if (iter == 1) {
        stop <- FALSE
        prev_diff <- curr_diff
      } else if (abs(mean(curr_diff - prev_diff)) <= 0.001) {
        stop <- TRUE
      } else if (iter >= max_iter){
        stop <- TRUE
      } else {
        stop <- FALSE
        prev_diff <- curr_diff
      }
    }

    rules <- list.rules.party(selected_learner$fit_object)
    quantile <- seq(length(rules))


    if (length(rules) == 1) {
      if(rules == "") {
        rules <- "No Rules Found"
      }
    }
    rules <- as.data.frame(cbind(rules, fold, target_m, quantile))

    backfit_resids <- (At$y_scaled - glmtree_model_preds_offset)^2
    backfit_RMSE <- sqrt(mean(backfit_resids))

    rules$RMSE <- backfit_RMSE


    marg_decisions[[i]] <- rules

  }

  marg_decisions <- do.call(rbind, marg_decisions)

  return(marg_decisions)
}

