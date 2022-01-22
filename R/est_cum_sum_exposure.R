#' #' @title Calculate the cumulative sum exposure based on marginals found in the fold
#' #' @description For each marginal mixture component rule found, create a g estimator for the probability of being exposed to the rule thresholds,
#' #' and a Q estimator for the outcome E(Y| A = a_mix, W). Get estimates of g and Q using the validation data and
#' #' calculate the clever covariate used in the TMLE fluctuation step.
#' #'
#' #' @param At_c Training data
#' #' @param Av_c Validation data
#' #' @param W Vector of characters denoting covariates
#' #' @param Q1_stack Super Learner library for fitting Q (outcome mechanism) and g (treatment mechanism)
#' #' @param H.AW_trunc_lvl Truncation level of the clever covariate (induces more bias to reduce variance)
#' #' @param no_marg_rules TRUE/FALSE if no marginal rules were found across the folds
#' #'
#' #' @importFrom magrittr %>%
#' #' @importFrom rlang :=
#' #' @importFrom dplyr group_by filter top_n
#' #' @return Rules object. TODO: add more detail here.
#' #'
#' #' @export
#'
#'
#' est_cum_sum_exposure <- function(At_c, Av_c, W, Q1_stack, no_marg_rules, H.AW_trunc_lvl) {
#'   if (no_marg_rules == FALSE) {
#'     H.AW.list <- list()
#'
#'     for (i in rev(seq(table(At_c$sum_marg_hits)))) {
#'       level <- table(At_c$sum_marg_hits)[i]
#'       if (level[[1]] < 10 & as.numeric(names(level)) != min(At_c$sum_marg_hits)) {
#'         At_c$sum_marg_hits[At_c$sum_marg_hits == as.numeric(names(level))] <- as.numeric(names(level)) - 1
#'         Av_c$sum_marg_hits[Av_c$sum_marg_hits == as.numeric(names(level))] <- as.numeric(names(level)) - 1
#'       } else if (level[[1]] < 10 & as.numeric(names(level)) == min(At_c$sum_marg_hits)) {
#'         At_c$sum_marg_hits[At_c$sum_marg_hits == as.numeric(names(level))] <- 0
#'         Av_c$sum_marg_hits[Av_c$sum_marg_hits == as.numeric(names(level))] <- 0
#'       }
#'     }
#'
#'
#'     for (i in sort(unique(At_c$sum_marg_hits[At_c$sum_marg_hits != 0]))) {
#'       target.lvl <- sort(unique(At_c$sum_marg_hits))[sort(unique(At_c$sum_marg_hits)) == i]
#'       X_train_covars <- At_c[W]
#'       X_valid_covars <- Av_c[W]
#'
#'       At_c$binarized_cat <-
#'         as.numeric(At_c$sum_marg_hits == target.lvl)
#'
#'       Av_c$binarized_cat <-
#'         as.numeric(Av_c$sum_marg_hits == target.lvl)
#'
#'
#'
#'       gHat1W <- predict(gHatSL, newdata = X_valid_covars)$pred
#'
#'       H.AW <- calc_clever_covariate(gHat1W, data = Av_c, exposure = "binarized_cat", H.AW_trunc_lvl = 10, type = "mod")
#'
#'       H.AW.list[[i]] <- H.AW
#'     }
#'
#'     HA.W.by.lvl <- do.call(cbind, H.AW.list)
#'     HA.W.lvl.sums <- rowSums(HA.W.by.lvl, na.rm = TRUE)
#'
#'     X_train_mix <- At_c[, c("sum_marg_hits", W)]
#'     X_valid_mix <- Av_c[, c("sum_marg_hits", W)]
#'
#'
#'
#'     QbarAW <- bound_precision(predict(QbarAWSL_m, newdata = X_valid_mix)$pred)
#'
#'     Av_c$QbarAW_additive <- QbarAW
#'     Av_c$HAW_additive <- HA.W.lvl.sums
#'   } else {
#'     Av_c$QbarAW_additive <- NA
#'     Av_c$HAW_additive <- NA
#'   }
#'
#'   return(Av_c)
#' }
