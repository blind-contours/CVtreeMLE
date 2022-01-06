#' @title Estimate the expected outcome for the combination of marginal thresholds identified in the fold.
#' @description Estimate the expected outcome given exposure to the combination of marginal exposures. This is different
#' compared to the cumulative sum; whereas with the cumulative sum, the exposure is the additive effect of each
#' marginal rule found in the fold, here each marginal rule is included in a Super Learner as a binary vector and
#' therefore this can pick up possible nonlinearity between the combination of binary exposures.
#'
#' @param At Training data
#' @param Av Validation data
#' @param Y Outcome variable
#' @param W Vector of characters denoting covariates
#' @param marg_rule_train Dataframe of binary vectors for marginal rules identified in the training fold
#' @param marg_rule_valid Dataframe of binary vectors for marginal rules identified in the validation fold
#' @param SL.library Super Learner library for fitting Q (outcome mechanism) and g (treatment mechanism)
#'
#' @import SuperLearner
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#' @importFrom dplyr group_by filter top_n
#' @return Rules object. TODO: add more detail here.
#'
#' @export

est_comb_exposure <- function(At, Av, Y, W, marg_rule_train, marg_rule_valid, SL.library) {

  if (all(is.na(marg_rule_train)) == FALSE) {
    At_mc <- At
    Av_mc <- Av

    marg_rule_df <-
      marg_rule_train[, colSums(is.na(marg_rule_train)) < nrow(marg_rule_train)]

    At_marg_comb <-
      cbind(marg_rule_df, At_mc[W])

    At_marg_comb <- At_marg_comb[, colSums(is.na(At_marg_comb)) < nrow(At_marg_comb)]

    Av_marg_comb <-
      cbind(marg_rule_valid, Av_mc[W])

    Av_marg_comb <- Av_marg_comb[, colSums(is.na(Av_marg_comb)) < nrow(Av_marg_comb)]

    QbarAWSL_m <- SuperLearner(
      Y = At_mc$y_scaled,
      X = At_marg_comb,
      SL.library = SL.library,
      family = "gaussian",
      verbose = FALSE
    )

    QbarAW <- bound_precision(predict(QbarAWSL_m, newdata = Av_marg_comb)$pred)
    QbarAW <- scale_to_original(scaled_vals = QbarAW, max_orig = max(At_mc[Y]), min_orig = min(At_mc[Y]))

    Av_marg_comb$QbarAW_combo <- QbarAW
    Av_marg_comb$y_scaled <- Av_mc$y_scaled
    Av_marg_comb$raw_outcome <- Av[, Y]
  }else{
    Av_marg_comb <- NA
    QbarAWSL_m <- NA

  }

  return(list("data" = Av_marg_comb, "learner" = QbarAWSL_m))
}
