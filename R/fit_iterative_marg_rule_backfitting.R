#' Iteratively Backfit a Super Learner, h(x) = Y|M_neX,W, and Indiviidual Decision Algorithms, g(x), for each Y|M_x until Convergence.
#'
#' @details Same procedure as `fit_iterative_mix_rule_backfitting` but the Super Learner controls for W and other mixture variables not equal to the mixture variable of interest, `glmtree` is used to find partitions
#' instead of `pre`.
#'
#' @param mix_comps A vector of characters indicating variables for the mixture components
#' @param At Training data
#' @param W A vector of characters indicating variables that are covariates#'
#' @param Q1_stack Stack of algorithms made in SL 3 used in ensemble machine learning to fit Y|W
#' @param fold Current fold in the cross-validation
#' @param verbose Run in verbose setting
#' @return Rules object. TODO: add more detail here.
#' @import partykit
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect
#' @importFrom stats na.omit
#' @importFrom purrr discard
#' @importFrom rpart.plot rpart.rules
#'
#' @export

fit_iterative_marg_rule_backfitting <- function(mix_comps,
                                                At,
                                                W,
                                                Q1_stack,
                                                fold,
                                                verbose) {


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

    formula <-
      as.formula(paste("y_scaled", "~", target_m))

    ctree_fit <- partykit::glmtree(formula,
                         data = At,
                         alpha = 0.9,
                         prune = "AIC",
                         minsize = n/4,
                         maxdepth = 4)

    Qbar_M_W_initial <- predict(ctree_fit, newdata = At)

    At[,"Qbar_M_W_initial"] <- Qbar_M_W_initial

    iter <- 0
    stop <- FALSE

    At_no_offset <- data.table::copy(At)
    At_no_offset$Qbar_M_W_initial <- 0
    At_no_offset$Qbar_ne_M_W_initial <- 0

    while (stop == FALSE) {
      iter <- iter + 1

      task <- sl3::make_sl3_Task(
        data = At,
        covariates = covars_m,
        outcome = "y_scaled",
        outcome_type = "continuous",
        offset = "Qbar_M_W_initial"
      )

      sl_fit_backfit <- discrete_sl$train(task)

      sl_fit_backfit_no_offset <- sl3::sl3_Task$new(
        data = At_no_offset,
        covariates = covars_m,
        outcome = "y_scaled",
        outcome_type = "continuous",
        offset = "Qbar_M_W_initial"
      )
      # preds_offset <- sl_fit_backfit$predict()
      preds_no_offset <- sl_fit_backfit$predict(sl_fit_backfit_no_offset)


      At[,"Qbar_ne_M_W_now"] <- preds_no_offset

      ctree_fit <- partykit::glmtree(formula,
                                     data = At,
                                     offset = Qbar_ne_M_W_initial,
                                     alpha = 0.9,
                                     prune = "AIC",
                                     minsize = n/4)

      At[, "Qbar_M_W_now"] <- predict(ctree_fit, newdata = At_no_offset)

      delta_h <- mean(abs(At$Qbar_ne_M_W_initial - At$Qbar_ne_M_W_now))
      delta_g <- mean(abs(At$Qbar_M_W_initial - At$Qbar_M_W_now))

      diff <- abs(delta_g - delta_h)

      Qbar_ne_M_W_initial <-  At$Qbar_ne_M_W_now
      At$Qbar_ne_M_W_initial <- Qbar_ne_M_W_initial
      At$Qbar_M_W_initial <- At$Qbar_M_W_now

      if (verbose){
        print(paste("iter: ", iter, "SL: ", delta_h, "ctree:", delta_g, "Diff: ", diff, "Rule:",  list.rules.party(ctree_fit)))
      }

      if (iter == 1) {
        stop <- FALSE
      }else if(diff <= 0.01) {
        stop <- TRUE
      }else{
        # prev_diff <- diff
        stop <- FALSE
      }
    }

    rules <- list.rules.party(ctree_fit)
    rules <- rules[length(rules)]

    if (rules == "") {
      rules <- "No Rules Found"
    }
    rules <- cbind(rules, fold, target_m)

    marg_decisions[[i]] <- rules

  }

  marg_decisions <- do.call(rbind, marg_decisions)

  return(marg_decisions)
}

