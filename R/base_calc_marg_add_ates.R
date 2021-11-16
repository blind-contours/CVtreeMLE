fit_post_counterfactuals <- function(modeling_results, target_mixtures, H.AW_trunc_lvl) {
  ## TODO: Fix whatever the fuck is going on here with the NULL compounding in lists from foreah
  
  marg_comb_fold_data <- modeling_results$`Fold Results Marg Comb Data`
  marg_comb_sls <- modeling_results$`Fold Super Learners Marg Comb`
  
  # marg_comb_fold_data <- unlist(marg_comb_fold_data,recursive=FALSE, use.names = FALSE)
  marg_comb_sls <- unlist(marg_comb_sls,recursive=FALSE, use.names = FALSE)
  
  # marg_comb_fold_data <- marg_comb_fold_data[!sapply(marg_comb_fold_data,is.null)]
  marg_comb_sls <- marg_comb_sls[!sapply(marg_comb_sls,is.null)]
  

  updated_resuls_list <- list()

  for (i in seq(marg_comb_sls)) {
    sl <- marg_comb_sls[[i]]
    combo_data <- marg_comb_fold_data[[i]]
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

    gHatSL <- suppressWarnings(SuperLearner(
      Y = combo_data$comb_rule,
      X = X_Amix_V,
      SL.library = SL.library,
      family = "binomial",
      verbose = FALSE
    ))

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

  ## least optimal submodel
  logitUpdate <- glm(y_scaled ~ -1 + H.AW + offset(qlogis(QbarAW_combo)),
    family = "quasibinomial", data = updated_results_counterfactuals
  )

  epsilon <- logitUpdate$coef

  QbarAW.star <- plogis(qlogis(updated_results_counterfactuals$QbarAW_combo) + epsilon * updated_results_counterfactuals$H.AW)
  Qbar1W.star <- plogis(qlogis(updated_results_counterfactuals$Qbar1W) + epsilon * updated_results_counterfactuals$H.AW)
  Qbar0W.star <- plogis(qlogis(updated_results_counterfactuals$Qbar0W) + epsilon * updated_results_counterfactuals$H.AW)

  ## back-scale Y
  QbarAW.star <- QbarAW.star * (max(updated_results_counterfactuals$raw_outcome) - min(updated_results_counterfactuals$raw_outcome)) + min(updated_results_counterfactuals$raw_outcome)
  Qbar0W.star <- Qbar0W.star * (max(updated_results_counterfactuals$raw_outcome) - min(updated_results_counterfactuals$raw_outcome)) + min(updated_results_counterfactuals$raw_outcome)
  Qbar1W.star <- Qbar1W.star * (max(updated_results_counterfactuals$raw_outcome) - min(updated_results_counterfactuals$raw_outcome)) + min(updated_results_counterfactuals$raw_outcome)

  updated_results_counterfactuals$QbarAW.star <- QbarAW.star
  updated_results_counterfactuals$Qbar0W.star <- Qbar0W.star
  updated_results_counterfactuals$Qbar1W.star <- Qbar1W.star

  updated_results_counterfactuals$mix.ATE <- updated_results_counterfactuals$Qbar1W.star - updated_results_counterfactuals$Qbar0W.star

  Thetas <- tapply(updated_results_counterfactuals$mix.ATE, updated_results_counterfactuals$folds, mean, na.rm = TRUE)

  for (i in seq(Thetas)) {
    updated_results_counterfactuals[updated_results_counterfactuals$folds == i, "Thetas"] <- Thetas[i][[1]]
  }

  ICs <- base::by(updated_results_counterfactuals, updated_results_counterfactuals$folds, function(updated_results_counterfactuals) {
    result <- with(updated_results_counterfactuals, H.AW * (updated_results_counterfactuals$raw_outcome - QbarAW.star) + Qbar1W.star - Qbar0W.star - Thetas)
    result
  })

  for (i in seq(ICs)) {
    updated_results_counterfactuals[updated_results_counterfactuals$folds == i, "IC"] <- ICs[i]
  }

  n <- dim(updated_results_counterfactuals)[1]
  varHat.IC <- var(updated_results_counterfactuals$IC, na.rm = TRUE) / n
  se <- sqrt(varHat.IC)

  alpha <- 0.05

  Theta <- mean(Thetas)
  # obtain 95% two-sided confidence intervals:
  CI <- c(
    Theta + qnorm(alpha / 2, lower.tail = T) * se,
    Theta + qnorm(alpha / 2, lower.tail = F) * se
  )

  # p-value
  p.value <- 2 * pnorm(abs(Theta / se), lower.tail = F)

  
  ## TODO: Fix this RMSE output to reflect full combo data not the counterfactuals
  
  ## calculate RMSE
  RMSE <- sqrt((updated_results_counterfactuals$QbarAW.star - updated_results_counterfactuals$raw_outcome)^2 / n)
  RMSE <- mean(RMSE[, 1])


  mixture_results <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 5))
  colnames(mixture_results) <- c("Mixture ATE", "Standard Error", "Lower CI", "Upper CI", "P-value")
  rownames(mixture_results) <- "Mixture"

  mixture_results$`Mixture ATE` <- Theta
  mixture_results$`Standard Error` <- se
  mixture_results$`Lower CI` <- CI[1]
  mixture_results$`Upper CI` <- CI[2]
  mixture_results$`P-value` <- p.value
  mixture_results$`Model RMSE` <- RMSE

  return(mixture_results)
}
