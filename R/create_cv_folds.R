#' @title Stratified CV to insure balance (by one grouping variable, Y)
#' @description Creates a dummy variable that partitions the data into v
#' equal sized groups for v-fold CV.
#' @param v number of folds
#' @param y Outcome variable. If binary will be used for stratification.
#' @param verbose If T will display extra output.
#'
#' @return Vector of fold assignments.
#'
#' @importFrom cvTools cvFolds
#' @export
create_cv_folds <- function(v, y, verbose = FALSE) {
  ys <- unique(y)
  nys <- length(ys)
  nn <- length(y)
  # Binary outcome so we can do stratified fold generation.
  if (nys == 2) {
    out <- rep(NA, nn)
    for (i in 1:nys) {
      # Record how many observations have this Y value.
      n <- sum(y == ys[i])
      folds <- cvTools::cvFolds(n, K = v,  R = 1, type = "random")$which

      out[y == ys[i]] <- folds
    }
    if (verbose) {
      cat("Cross-validation fold breakdown:\n")
      print(table(y, "Fold" = out, useNA = "ifany"))
    }
  } else {
    xx <- cvTools::cvFolds(nn, K = v,  R = 1, type = "random")$which
    out <- xx
  }
  return(out)
}
