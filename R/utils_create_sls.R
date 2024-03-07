#' @title Create default Super Learner estimators for the data adaptive
#' and nuisance parameters used in `CVtreeMLE`
#' @description If Super Learners are not passed to the stack arguments
#' this function is used to create some default ensemble machine
#' learning estimators for each parameter. The default estimators are
#' fast but also flexible.
#'
#' @return List of ensemble estimators
#' @import sl3
#' @export

create_sls <- function() {
  lrnr_ranger <- sl3::Lrnr_ranger$new()
  lrnr_glm <- sl3::make_learner(sl3::Lrnr_glm)
  lrnr_elasticnet <- sl3::make_learner(sl3::Lrnr_glmnet, alpha = .5)

  grid_params <- list(
    max_depth = c(3, 5, 8),
    eta = c(0.001, 0.1, 0.3),
    nrounds = c(100, 200, 500)
  )

  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)

  xgb_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
    do.call(sl3::Lrnr_xgboost$new, as.list(tuning_params))
  })

  learners <- c(
    lrnr_ranger, lrnr_elasticnet, lrnr_glm, xgb_learners[[4]], xgb_learners[[15]], xgb_learners[[20]], xgb_learners[[27]]
  )

  w_stack <- sl3::make_learner(sl3::Stack, learners)
  aw_stack <- sl3::make_learner(sl3::Stack, learners)


  return(list(
    "W_stack" = w_stack,
    "AW_stack" = aw_stack
  ))
}
