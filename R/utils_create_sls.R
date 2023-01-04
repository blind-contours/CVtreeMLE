#' @title Create default Super Learner estimators for the data adaptive
#' and nuisance parameters used in `CVtreeMLE`
#' @description If Super Learners are not passed to the stack arguments
#' this function is used to create some default ensemble machine
#' learning estimators for each parameter. The default estimators are
#' fast but also flexible.
#'
#' @return List of ensemble estimators
#' @export

create_sls <- function() {
  lrnr_ranger <- sl3::Lrnr_ranger$new()
  lrnr_glm <- sl3::make_learner(sl3::Lrnr_glm)
  lrnr_elasticnet <- sl3::make_learner(sl3::Lrnr_glmnet, alpha = .5)

  grid_params <- list(
    max_depth = c(3, 5, 8),
    eta = c(0.001, 0.1, 0.3),
    nrounds = 100
  )

  grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)

  xgb_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
    do.call(sl3::Lrnr_xgboost$new, as.list(tuning_params))
  })

  learners <- c(
    lrnr_glm,
    lrnr_ranger, lrnr_elasticnet,
    xgb_learners
  )

  w_stack <- sl3::make_learner(sl3::Stack, learners)
  aw_stack <- sl3::make_learner(sl3::Stack, learners)

  lrnr_glmtree_001 <- sl3::Lrnr_glmtree$new(
    alpha = 0.05,
    maxdepth = 2,
    bonferroni = TRUE,
    minsize = 50
  )

  lrnr_glmtree_002 <- sl3::Lrnr_glmtree$new(
    alpha = 0.05,
    maxdepth = 3,
    bonferroni = TRUE,
    minsize = 50
  )

  lrnr_glmtree_003 <- sl3::Lrnr_glmtree$new(
    alpha = 0.1,
    maxdepth = 2,
    bonferroni = TRUE,
    minsize = 50
  )

  lrnr_glmtree_004 <- sl3::Lrnr_glmtree$new(
    alpha = 0.1,
    maxdepth = 1,
    bonferroni = TRUE,
    minsize = 50
  )

  lrnr_glmtree_005 <- sl3::Lrnr_glmtree$new(
    alpha = 0.1,
    maxdepth = 2,
    bonferroni = TRUE,
    minsize = 50
  )

  lrnr_glmtree_006 <- sl3::Lrnr_glmtree$new(
    alpha = 0.1,
    maxdepth = 3,
    bonferroni = FALSE,
    minsize = 50
  )

  lrnr_glmtree_007 <- sl3::Lrnr_glmtree$new(
    alpha = 0.1,
    maxdepth = 2,
    bonferroni = FALSE,
    minsize = 50
  )

  learners <- c(
    lrnr_glmtree_001, lrnr_glmtree_002,
    lrnr_glmtree_003, lrnr_glmtree_004,
    lrnr_glmtree_005, lrnr_glmtree_006,
    lrnr_glmtree_007
  )

  tree_stack <- sl3::make_learner(sl3::Stack, learners)

  return(list(
    "W_stack" = w_stack,
    "AW_stack" = aw_stack,
    "A_stack" = tree_stack
  ))
}
