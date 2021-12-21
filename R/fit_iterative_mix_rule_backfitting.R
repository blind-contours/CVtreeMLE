#' @title Iteratively Backfit a Super Learner, h(x) = Y|W, and an Ensemble Decision Algorithm, g(x), Y|M_x until Convergence.
#'
#' @details Performs an iterative backfitting algorithm to flexibly adjust for covariates W while finding the best fitting set of mixture rules to partition the space in M.
#'
#' @param At Training dataframe
#' @param A Variable names in the mixture
#' @param W Variable names in the covariates
#' @param Y Variable name for the outcome
#' @param Q1_stack Stack of algorithms made in SL 3 used in ensemble machine learning to fit Y|W
#' @param fold Current fold in the cross-validation
#' @import sl3
#' @importFrom pre pre maxdepth_sampler
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by filter top_n
#' @return Rules object. TODO: add more detail here.

#' @export

fit_iterative_mix_rule_backfitting <- function(At, A, W, Y, Q1_stack, fold) {
  pre_boot_list <- list()

  task <- sl3::make_sl3_Task(
    data = At, covariates = W,
    outcome = "y_scaled", outcome_type = "continuous"
  )

  discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new()

  discrete_sl <- sl3::Lrnr_sl$new(
    learners = Q1_stack,
    metalearner = discrete_sl_metalrn,
  )
  # delayed_sl_fit <- delayed_learner_train(discrete_sl, task)

  sl_fit <- discrete_sl$train(task)

  ## calculate remaining variance unexplained by W in residuals
  QbarW_initial <- sl_fit$predict()

  At <- cbind(At, QbarW_initial)

  formula <-
    as.formula(paste("y_scaled", "~", paste(A, collapse = "+")))

  pre_model_t0 <- pre::pre(formula,
                           data = At,
                           family = "gaussian",
                           use.grad = FALSE,
                           tree.unbiased = TRUE,
                           removecomplements = TRUE,
                           removeduplicates = TRUE,
                           maxdepth = pre::maxdepth_sampler(),
                           sampfrac = min(1, (11 * sqrt(dim(At)[1]) + 1) / dim(At)[1]),
                           nfolds = 10
  )

  QbarAW_initial <- predict(pre_model_t0)

  At <- cbind(At, QbarAW_initial)

  iter <- 0
  stop <- FALSE

  At_no_offset <- data.table::copy(At)
  At_no_offset$QbarW_initial <- 0
  At_no_offset$QbarAW_initial <- 0

  while (stop == FALSE) {
    iter <- iter + 1

    task <- sl3::make_sl3_Task(
      data = At, covariates = W,
      outcome = "y_scaled",
      outcome_type = "continuous",
      offset = "QbarAW_initial"
    )

    sl_fit_backfit <- discrete_sl$train(task)

    sl_fit_backfit_no_offset <- sl3::sl3_Task$new(
      data = At_no_offset,
      covariates = W,
      outcome = "y_scaled",
      outcome_type = "continuous",
      offset = "QbarAW_initial"
    )
    # preds_offset <- sl_fit_backfit$predict()
    preds_no_offset <- sl_fit_backfit$predict(sl_fit_backfit_no_offset)

    At$QbarW_now <- preds_no_offset

    pre_model <- pre::pre(formula,
                          data = At,
                          family = "gaussian",
                          use.grad = FALSE,
                          tree.unbiased = TRUE,
                          removecomplements = TRUE,
                          removeduplicates = TRUE,
                          maxdepth = pre::maxdepth_sampler(),
                          sampfrac = min(1, (11 * sqrt(dim(At)[1]) + 1) / dim(At)[1]),
                          nfolds = 10,
                          offset = At$QbarW_initial
    )

    At$QbarAW_now <- predict(pre_model, newoffset = At_no_offset$QbarW_initial)

    pre_model_coefs <- stats::coef(pre_model, penalty.par.val = "lambda.min")

    pre_coefs_no_zero <- base::subset(pre_model_coefs, coefficient != 0)
    pre_coefs_no_zero$test <- apply(pre_coefs_no_zero, 1, function(x) pull_out_rule_vars(x, A))

    rules <- pre_coefs_no_zero
    rules$boot_num <- iter
    pre_boot_list[[iter]] <- rules

    delta_h <- mean(abs(At$QbarW_initial - At$QbarW_now))
    delta_g <- mean(abs(At$QbarAW_initial - At$QbarAW_now))

    diff <- abs(delta_g - delta_h)

    print(paste("iter: ", iter, "SL: ", delta_h, "ctree:", delta_g, "Diff: ", diff, "Rule:"))


    At$QbarW_initial <- At$QbarW_now
    At$QbarAW_initial <- At$QbarAW_now

    if (iter == 1) {
      stop <- FALSE
    } else if (diff < 0.01) {
      stop <- TRUE
    } else {
      stop <- FALSE
    }
  }

  pre_boot_df <- do.call(rbind, pre_boot_list)
  pre_boot_df$direction <- ifelse(pre_boot_df$coefficient > 0, 1, 0)

  rules <- pre_boot_df %>%
    dplyr::group_by(test, direction) %>%
    dplyr::filter(dplyr::n() >= iter/2)

  rules <- rules[!is.na(rules$test), ]
  rules <- rules[!rules$test == 0, ]

  rules <- rules %>%
    dplyr::group_by(test, direction) %>%
    dplyr::top_n(1, abs(coefficient))

  rules$fold <- fold

  if (dim(rules)[1] == 0) {
    rules <- data.frame(matrix(nrow = 1, ncol = 7))
    colnames(rules) <- c("rule", "coefficient", "description", "test", "boot_num", "direction", "fold")
    rules$rule <- "None"
    rules$coefficient <- 0
    rules$description <- "No Rules Found"
    rules$test <- "No Rules Found"
    rules$boot_num <- 0
    rules$direction <- -1
    rules$fold <- fold
  }

  return(rules)
}
