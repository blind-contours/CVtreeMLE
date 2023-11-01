#' @title Just Fit Ensemble Decision Trees Including Covariates
#'
#' @details Does not fit iteratively, includes W in rules if found as well.
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
#' @importFrom stringr str_detect
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
#' example_output <- fit_pre_algorithm(
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

fit_pre_algorithm <- function(at,
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

  vars_to_keep <- c(a, w, y)
  at <- at[, vars_to_keep]

  formula <-
    as.formula(paste(y, "~."))

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
    par.final = FALSE,
    par.init = FALSE
  )


  pre_model_coefs <- stats::coef(pre_model, penalty.par.val = "lambda.min")

  pre_coefs_no_zero <- base::subset(pre_model_coefs, coefficient != 0)
  sets <- c(a, w)
  pre_coefs_no_zero$test <- apply(pre_coefs_no_zero, 1, function(x) {
    pull_out_rule_vars(x, sets)
  })

  effect_modifiers <- apply(pre_coefs_no_zero, 1, function(x) {
    any(str_detect(x[[4]], a)) & any(str_detect(x[[4]], w))
  })

  targets <- apply(pre_coefs_no_zero, 1, function(x) {
    any(str_detect(x[[4]], a)) & any(str_detect(x[[4]], w)) | any(str_detect(x[[4]], a))
  })

  pre_coefs_no_zero$effect_modifiers <- as.numeric(effect_modifiers)

  data_sub <- pre_coefs_no_zero[targets, ]

  rules <- data_sub


  rules$direction <- ifelse(rules$coefficient > 0, 1, 0)

  rules <- rules[!is.na(rules$test), ]

  if (direction == "positive") {
    rules <- rules %>%
      dplyr::group_by(test) %>%
      dplyr::slice_max(n = 1, coefficient)
  } else {
    rules <- rules %>%
      dplyr::group_by(test) %>%
      dplyr::slice_min(n = 1, coefficient)
  }
  # rules <- rules[str_detect(rules$rule, "rule"), ]
  rules$fold <- fold
  rules$RMSE <- NA
  rules <- rules[str_detect(rules$rule, "rule"), ]



  if (dim(rules)[1] == 0) {
    rules <- data.frame(matrix(nrow = 1, ncol = 7))
    colnames(rules) <- c(
      "rule", "coefficient", "description",
      "test", "direction", "fold", "RMSE"
    )
    rules$rule <- "None"
    rules$coefficient <- 0
    rules$description <- "No Rules Found"
    rules$test <- "No Rules Found"
    rules$effect_modifiers <- "None"
    rules$direction <- -1
    rules$fold <- fold
    rules$RMSE <- NA
  }

  return(list("rules" = rules, "model" = pre_model))
}
