est_comb_exposure <- function(At, Av, Y, W, marg_rule_train, marg_rule_valid, SL.library) {

  if (all(is.na(marg_rule_train)) == FALSE) {
    At_mc <- At
    Av_mc <- Av

    marg_rule_df <-
      marg_rule_train[, colSums(is.na(marg_rule_train)) < nrow(marg_rule_train)]

    At_marg_comb <-
      cbind(marg_rule_df, subset(At_mc, select = W))

    At_marg_comb <- At_marg_comb[, colSums(is.na(At_marg_comb)) < nrow(At_marg_comb)]

    Av_marg_comb <-
      cbind(marg_rule_valid, subset(Av_mc, select = W))

    Av_marg_comb <- Av_marg_comb[, colSums(is.na(Av_marg_comb)) < nrow(Av_marg_comb)]

    QbarAWSL_m <- SuperLearner::SuperLearner(
      Y = At_mc$y_scaled,
      X = At_marg_comb,
      SL.library = SL.library,
      family = "gaussian",
      verbose = FALSE
    )

    QbarAW <- bound_precision(predict(QbarAWSL_m, newdata = Av_marg_comb)$pred)

    Av_marg_comb$QbarAW_combo <- QbarAW
    Av_marg_comb$y_scaled <- Av_mc$y_scaled
    Av_marg_comb$raw_outcome <- Av[, Y]
  }else{
    Av_marg_comb <- NA
    QbarAWSL_m <- NA

  }

  return(list("data" = Av_marg_comb, "learner" = QbarAWSL_m))
}
