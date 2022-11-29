#' @title Iteratively Backfit a Super Learner, h(x) = Y|W, and an Ensemble
#'  Decision Algorithm, g(x), Y|M_x until Convergence.
#'
#' @details Performs an iterative backfitting algorithm to flexibly adjust
#' for covariates W while finding the best fitting set of mixture rules
#' to partition the space in M.
#'
#' @param at Training dataframe
#' @param a Variable names in the mixture
#' @param w Variable names in the covariates
#' @param y Variable name for the outcome
#' @param direction Positive/negative - max or min coefficient to keep in
#' the ensemble
#' @param w_stack Stack of algorithms made in SL 3 used in ensemble machine
#' learning to fit Y|W
#' @param fold Current fold in the cross-validation
#' @param max_iter Max number of iterations of iterative backfitting algorithm
#' @param verbose Run in verbose setting
#' @param parallel_cv TRUE/FALSE indicator to parallelize cv
#' @param seed Seed number for consistent results
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr group_by filter top_n
#' @return A list of the mixture rule results within a fold including:
#'  \itemize{
#'   \item \code{rules}: A data frame with the data adpatively
#'   determined rules found in the \code{pre} model along with the coefficient,
#'   direction, fold, RMSE and other measures.
#'   \item \code{model}: The best fitting pre model found in the fold.
#'   }
#' @examples
#' data <- simulate_mixture_cube()
#' mix_comps <- c("M1", "M2", "M3")
#' W <- c("age", "sex", "bmi")
#' sls <- create_sls()
#' w_stack <- sls$W_stack
#' tree_stack <- sls$A_stack
#' example_output <- fit_mix_rule_backfitting(
#'   at = data,
#'   a = mix_comps,
#'   w = W,
#'   y = "y",
#'   direction = "positive",
#'   w_stack = w_stack,
#'   fold = 1,
#'   max_iter = 1,
#'   verbose = FALSE,
#'   parallel = FALSE,
#'   seed = 6442
#' )
#' @export

fit_mix_rule_backfitting <- function(at,
                                     a,
                                     w,
                                     y,
                                     direction,
                                     w_stack,
                                     fold,
                                     max_iter,
                                     verbose,
                                     parallel_cv,
                                     seed) {
  if (parallel_cv == TRUE) {
    future::plan(future::sequential, gc = TRUE)
  }

  set.seed(seed)

  pre_boot_list <- list()

  task <- sl3::make_sl3_Task(
    data = at,
    covariates = w,
    outcome = y,
    outcome_type = "continuous"
  )

  discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new(sl3::loss_squared_error)

  discrete_sl <- sl3::Lrnr_sl$new(
    learners = w_stack,
    metalearner = discrete_sl_metalrn,
  )

  sl_fit <- discrete_sl$train(task)

  qbar_w_initial <- sl_fit$predict()

  at <- cbind(at, qbar_w_initial)

  formula <-
    as.formula(paste(y, "~", paste(a, collapse = "+")))

  pre_model_t0 <- pre::pre(formula,
    data = at,
    family = "gaussian",
    use.grad = FALSE,
    tree.unbiased = TRUE,
    removecomplements = TRUE,
    removeduplicates = TRUE,
    maxdepth = pre::maxdepth_sampler(),
    sampfrac = min(1, (11 * sqrt(dim(at)[1]) + 1) /
      dim(at)[1]),
    nfolds = 10,
    par.final = FALSE,
    par.init = FALSE
  )

  qbar_aw_initial <- predict(pre_model_t0)

  at <- cbind(at, qbar_aw_initial)

  iter <- 0
  stop <- FALSE

  at_no_offset <- data.table::copy(at)
  at_no_offset$qbar_w_initial <- 0
  at_no_offset$qbar_aw_initial <- 0

  while (stop == FALSE) {
    iter <- iter + 1

    task <- sl3::make_sl3_Task(
      data = at,
      covariates = w,
      outcome = y,
      outcome_type = "continuous",
      offset = "qbar_aw_initial"
    )

    sl_fit_backfit <- discrete_sl$train(task)

    sl_fit_backfit_no_offset <- sl3::sl3_Task$new(
      data = at_no_offset,
      covariates = w,
      outcome = y,
      outcome_type = "continuous",
      offset = "qbar_aw_initial"
    )
    preds_no_offset <- sl_fit_backfit$predict(sl_fit_backfit_no_offset)
    preds_offset <- sl_fit_backfit$predict(task)

    at$QbarW_now <- preds_no_offset

    pre_model <- pre::pre(formula,
      data = at,
      family = "gaussian",
      use.grad = FALSE,
      tree.unbiased = TRUE,
      removecomplements = TRUE,
      removeduplicates = TRUE,
      maxdepth = pre::maxdepth_sampler(),
      sampfrac = min(1, (11 * sqrt(dim(at)[1]) + 1) /
        dim(at)[1]),
      nfolds = 10,
      offset = at$qbar_w_initial,
      par.final = FALSE,
      par.init = FALSE
    )

    pre_model_preds_offset <- predict(pre_model,
      newoffset = at$qbar_w_initial
    )

    at$qbar_aw_now <- predict(pre_model,
      newoffset = at_no_offset$qbar_w_initial
    )

    pre_model_coefs <- stats::coef(pre_model, penalty.par.val = "lambda.min")

    pre_coefs_no_zero <- base::subset(pre_model_coefs, coefficient != 0)
    pre_coefs_no_zero$test <- apply(pre_coefs_no_zero, 1, function(x) {
      pull_out_rule_vars(x, a)
    })

    rules <- pre_coefs_no_zero
    rules$boot_num <- iter
    pre_boot_list[[iter]] <- rules

    curr_diff <- abs(pre_model_preds_offset - preds_offset)

    if (verbose) {
      if (iter == 1) {
        print(paste(
          "Fold: ", fold, "|",
          "Process: ", "Mixture Decision Backfitting", "|",
          "Iteration: ", iter, "|",
          "Delta: ", "None", "|",
          "Diff: ", mean(curr_diff), "|",
          "Rules:", dim(pre_coefs_no_zero)[1]
        ))
      } else {
        print(paste(
          "Fold: ", fold, "|",
          "Process: ", "Mixture Decision Backfitting", "|",
          "Iteration: ", iter, "|",
          "Delta: ", mean(curr_diff - prev_diff), "|",
          "Diff: ", mean(curr_diff), "|",
          "Rules:", dim(pre_coefs_no_zero)[1]
        ))
      }
    }

    at$qbar_w_initial <- at$qbar_w_now
    at$qbar_aw_initial <- at$qbar_aw_now

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

  pre_boot_df <- do.call(rbind, pre_boot_list)
  pre_boot_df$direction <- ifelse(pre_boot_df$coefficient > 0, 1, 0)

  rules <- pre_boot_df %>% dplyr::filter(.data$boot_num == iter)
  rules <- rules[!is.na(rules$test), ]
  # rules <- rules[!rules$test == 0, ]

  if (direction == "positive") {
    rules <- rules %>%
      dplyr::group_by(test) %>%
      dplyr::slice_max(n = 1, coefficient)
  } else {
    rules <- rules %>%
      dplyr::group_by(test) %>%
      dplyr::slice_min(n = 1, coefficient)
  }

  rules$fold <- fold

  backfit_resids <- (at[, y] - pre_model_preds_offset)^2
  backfit_rmse <- sqrt(mean(backfit_resids))

  rules$RMSE <- backfit_rmse

  if (dim(rules)[1] == 0) {
    rules <- data.frame(matrix(nrow = 1, ncol = 7))
    colnames(rules) <- c(
      "rule", "coefficient", "description",
      "test", "boot_num", "direction", "fold"
    )
    rules$rule <- "None"
    rules$coefficient <- 0
    rules$description <- "No Rules Found"
    rules$test <- "No Rules Found"
    rules$boot_num <- 0
    rules$direction <- -1
    rules$fold <- fold
    rules$RMSE <- NA
  }

  return(list("rules" = rules, "model" = pre_model))
}
