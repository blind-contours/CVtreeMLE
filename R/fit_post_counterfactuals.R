#' Given a vector of target mixture components, load a \code{CVtreeMLE} results object, get the lists of rules, nuisance parameters, and learners used. Apply counterfactuals to the target
#' variables of interest and return a dataframe with the ATE, variance estimates and p-values for these counterfactuals.
#'
#' Currently uses PRE Predictive Rule Ensembles package but could be generalized
#'
#' @param modeling_results A \code{CVtreeMLE} results object
#' @param target_mixtures Vector of characters indicating which mixture compnents to calculate a new counterfactual parameter on given rules found across the folds for each marginal
#' mixture component
#' @param H.AW_trunc_lvl Truncation level of the clever covariate - reduces variance but increases bias of the ATE
#' @param SL.library Library of algorithms used to fit the counterfactuals
#' @param p_adjust_n p-value adjustment for multiple comparisons

#' @importFrom stats glm pnorm qnorm
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by filter top_n
#' @return Rules object. TODO: add more detail here.

#'
#' @export

fit_post_counterfactuals <- function(modeling_results, target_mixtures, H.AW_trunc_lvl, SL.library, p_adjust_n) {
  ## TODO: Fix whatever the fuck is going on here with the NULL compounding in lists from foreach

  marg_comb_fold_data <- modeling_results$`Fold Results Marg Comb Data`
  marg_comb_sls <- modeling_results$`Fold Super Learners Marg Comb`

  marg_comb_sls <- unlist(marg_comb_sls,recursive=FALSE, use.names = FALSE)
  marg_comb_sls <- marg_comb_sls[!sapply(marg_comb_sls,is.null)]

  marg_comb_fold_data <- groupby_fold(marg_comb_fold_data)

  marg_comb_fold_data <- dplyr::group_split(marg_comb_fold_data)

  updated_resuls_list <- list()

  for (i in seq(marg_comb_sls)) {
    sl <- marg_comb_sls[[i]]
    combo_data <- marg_comb_fold_data[[i]]
    combo_data <- combo_data[-1]
    n <- dim(combo_data)[1]
    folds <- i

    matches <- unique(grep(paste(target_mixtures, collapse = "|"),
                           names(combo_data),
                           value = TRUE
    ))


    combo_data_1 <- combo_data_0 <- combo_data
    combo_data_1[, matches] <- 1 # under exposure
    combo_data_0[, matches] <- 0 # under control

    combo_data_1 <- combo_data_1 %>% dplyr::select(-QbarAW_combo, -y_scaled, -raw_outcome)
    combo_data_0 <- combo_data_0 %>% dplyr::select(-QbarAW_combo, -y_scaled, -raw_outcome)

    Qbar1W <- predict(sl, newdata = combo_data_1)$pred
    Qbar0W <- predict(sl, newdata = combo_data_0)$pred

    Qbar1W[Qbar1W >= 1.0] <- 0.999
    Qbar1W[Qbar1W <= 0] <- 0.001

    Qbar0W[Qbar0W >= 1.0] <- 0.999
    Qbar0W[Qbar0W <= 0] <- 0.001

    combo_data$Qbar1W <- Qbar1W
    combo_data$Qbar0W <- Qbar0W

    #######################################################################################
    ############# GET PROPENSITY SCORE FOR RULE THAT MAKES THE COUNTEFACTUALS #############
    #######################################################################################

    # TODO: add try catch if no observations with both marginal exposure rules found
    combo_data$comb_rule <- rowSums(combo_data[, matches])
    combo_data$comb_rule <- ifelse(combo_data$comb_rule == max(combo_data$comb_rule), 1, 0)

    # Confirm that some observations exist who were exposed to both rules in the fold
    if (max(combo_data$comb_rule) == 0 ) {
      stop("No observations were exposed to all marginal rules")
    }

    X_Amix_V <- combo_data %>% dplyr::select(-matches, -QbarAW_combo, -y_scaled, -Qbar1W, -Qbar0W, -comb_rule, -raw_outcome)

    print("Fitting SL to marginal rule for mixture")

    gHatSL <- SuperLearner(
      Y = combo_data$comb_rule,
      X = X_Amix_V,
      SL.library = SL.library,
      family = "binomial",
      verbose = FALSE
    )

    gHat1W <- gHatSL$SL.predict
    gHat0W <- 1 - gHat1W

    gHatAW <- rep(NA, n)
    gHatAW[combo_data$comb_rule == 1] <- gHat1W[combo_data$comb_rule == 1]
    gHatAW[combo_data$comb_rule == 0] <- gHat0W[combo_data$comb_rule == 0]

    H.AW <- as.numeric(combo_data$comb_rule == 1) / gHat1W - as.numeric(combo_data$comb_rule == 0) / gHat0W

    H.AW <-
      ifelse(H.AW > H.AW_trunc_lvl, H.AW_trunc_lvl, H.AW)

    H.AW <-
      ifelse(H.AW < -H.AW_trunc_lvl, -H.AW_trunc_lvl, H.AW)


    ## add treatment mechanism results to validation dataframe
    combo_data$gHat1W <- gHat1W
    combo_data$H.AW <- H.AW
    combo_data$folds <- i


    updated_resuls_list[[i]] <- combo_data
  }

  updated_results_counterfactuals <- do.call(rbind, updated_resuls_list)

  flux_model_results <- fit_least_fav_submodel(H.AW = updated_results_counterfactuals$H.AW,
                         data = updated_results_counterfactuals,
                         QbarAW = scale_to_unit(updated_results_counterfactuals$QbarAW_combo),
                         Qbar1W = updated_results_counterfactuals$Qbar1W,
                         Qbar0W = updated_results_counterfactuals$Qbar0W)

  ## back-scale Y
  QbarAW.star <- scale_to_original(scaled_vals = flux_model_results$QbarAW.star,
                                   max_orig = max(updated_results_counterfactuals$raw_outcome),
                                   min(updated_results_counterfactuals$raw_outcome))

  Qbar1W.star <- scale_to_original(scaled_vals = flux_model_results$Qbar1W.star,
                                   max_orig = max(updated_results_counterfactuals$raw_outcome),
                                   min(updated_results_counterfactuals$raw_outcome))

  Qbar0W.star <- scale_to_original(scaled_vals = flux_model_results$Qbar0W.star,
                                   max_orig = max(updated_results_counterfactuals$raw_outcome),
                                   min(updated_results_counterfactuals$raw_outcome))

  updated_results_counterfactuals$QbarAW.star <- QbarAW.star
  updated_results_counterfactuals$Qbar0W.star <- Qbar0W.star
  updated_results_counterfactuals$Qbar1W.star <- Qbar1W.star

  ATE_results <- calc_ATE_estimates(data = updated_results_counterfactuals, ATE_var = "mix.ATE", outcome = "raw_outcome", p_adjust_n = p_adjust_n)


  ## TODO: Fix this RMSE output to reflect full combo data not the counterfactuals

  ## calculate RMSE
  RMSE <- (updated_results_counterfactuals$QbarAW.star - updated_results_counterfactuals$raw_outcome)^2
  non_additive_RMSE <- sqrt(mean(RMSE))

  mixture_results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
  colnames(mixture_results) <- c("Mixture ATE", "Standard Error", "Lower CI", "Upper CI", "P-value", "P-value adj", "RMSE")
  rownames(mixture_results) <- "Mixture"

  mixture_results$`Mixture ATE` <- ATE_results$ATE
  mixture_results$`Standard Error` <- ATE_results$SE
  mixture_results$`Lower CI` <- ATE_results$CI[1]
  mixture_results$`Upper CI` <- ATE_results$CI[2]
  mixture_results$`P-value` <- ATE_results$`p-value`
  mixture_results$`P-value adj` <- ATE_results$`adj p-value`
  mixture_results$`RMSE` <- non_additive_RMSE

  return(mixture_results)
}
