#' @title Create default Super Learner estimators for the data adaptive
#' and nuisance parameters used in `CVtreeMLE`
#' @description If Super Learners are not passed to the stack arguments
#' this function is used to create some default ensemble machine
#' learning estimators for each parameter. The default estimators are
#' fast but also flexible.
#'
#' @import sl3
#' @return List of ensemble estimators
#' @export

create_sls <- function() {

  mean_lrnr <- Lrnr_mean$new()
  lrnr_ranger <- Lrnr_ranger$new()
  lrnr_glm <- make_learner(Lrnr_glm)
  lrnr_elasticnet <- make_learner(Lrnr_glmnet, alpha = .5)
  lrnr_xgboost <- Lrnr_xgboost$new()

  learners <- c(mean_lrnr, lrnr_glm,
                lrnr_ranger, lrnr_elasticnet,
                lrnr_xgboost)

  w_stack <- make_learner(Stack, learners)
  aw_stack <- make_learner(Stack, learners)

  lrnr_glmtree_001 <- Lrnr_glmtree$new(alpha = 0.1,
                                       maxdepth = 3,
                                       bonferroni = TRUE,
                                       minsize = 20)

  lrnr_glmtree_002 <- Lrnr_glmtree$new(alpha = 0.2,
                                       maxdepth = 4,
                                       bonferroni = TRUE,
                                       minsize = 30)

  lrnr_glmtree_003 <- Lrnr_glmtree$new(alpha = 0.3,
                                       maxdepth = 2,
                                       bonferroni = TRUE,
                                       minsize = 30)

  lrnr_glmtree_004 <- Lrnr_glmtree$new(alpha = 0.5,
                                       maxdepth = 1,
                                       bonferroni = TRUE,
                                       minsize = 20)

  lrnr_glmtree_005 <- Lrnr_glmtree$new(alpha = 0.5,
                                       maxdepth = 4,
                                       bonferroni = TRUE,
                                       minsize = 20)

  lrnr_glmtree_006 <- Lrnr_glmtree$new(alpha = 0.1,
                                       maxdepth = 4,
                                       bonferroni = FALSE,
                                       minsize = 20)

  lrnr_glmtree_007 <- Lrnr_glmtree$new(alpha = 0.1,
                                       maxdepth = 5,
                                       bonferroni = FALSE,
                                       minsize = 20)

  learners <- c(lrnr_glmtree_001, lrnr_glmtree_002,
                lrnr_glmtree_003, lrnr_glmtree_004,
                lrnr_glmtree_005, lrnr_glmtree_006,
                lrnr_glmtree_007)

  tree_stack <- make_learner(Stack, learners)

  return(list("W_stack" = w_stack,
              "AW_stack" = aw_stack,
              "A_stack" = tree_stack))
}
