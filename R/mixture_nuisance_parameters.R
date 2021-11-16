mixture_nuisance_parameters <- function(At, Av, covariates, no_rules, 
                                        SL.library, family, rules, fold_k, mix_fold_data, 
                                        fold_results_mix_rules, H.AW_trunc_lvl){
  
  At_mix <- At
  Av_mix <- Av
  
  if (no_rules != TRUE) {
    At_rules_eval <-
      evaluate_mixture_rules(data = At_mix, rules = rules)
    Av_rules_eval <-
      evaluate_mixture_rules(data = Av_mix, rules = rules)
    
    mix_interaction_data <- list()
    mix_directions <- list()
    
    for (interaction in seq(dim(At_rules_eval)[2])) {
      interaction_rule <- At_rules_eval[, interaction]
      
      if (dim(table(interaction_rule)) == 2) {
        At_mix$A_mix <- interaction_rule
        Av_mix$A_mix <- Av_rules_eval[, interaction]
        
        ## select covariates
        X_Amix_T <- subset(At_mix,
                           select = c(covariates))
        
        X_Amix_V <- subset(Av_mix,
                           select = c(covariates))
        
        print("Fitting SL to marginal rule for mixture")
        
        gHatSL <- SuperLearner(
          Y = At_mix$A_mix,
          X = X_Amix_T,
          SL.library = SL.library,
          family = "binomial",
          verbose = FALSE
        )
        
        gHat1W <- predict(gHatSL, newdata = X_Amix_V)$pred
        n <- length(gHat1W)
        gHat0W <- 1 - gHat1W
        
        gHatAW <- rep(NA, n)
        gHatAW[Av_mix$A_mix == 1] <- gHat1W[Av_mix$A_mix == 1]
        gHatAW[Av_mix$A_mix == 0] <- gHat0W[Av_mix$A_mix == 0]
        
        H.AW <-
          as.numeric(Av_mix$A_mix == 1) / gHat1W - as.numeric(Av_mix$A_mix == 0) /
          gHat0W
        
        H.AW <-
          ifelse(H.AW > H.AW_trunc_lvl, H.AW_trunc_lvl, H.AW)
        
        H.AW <-
          ifelse(H.AW < -H.AW_trunc_lvl, -H.AW_trunc_lvl, H.AW)
        
        ## add treatment mechanism results to validation dataframe
        
        
        X_train_mix <-
          subset(At_mix, select = c("A_mix", covariates))
        
        X_valid_mix <-
          subset(Av_mix, select = c("A_mix", covariates))
        
        print(paste(
          "Fitting SL Y given W and marginal rule for mixture",
          interaction
        ))
        
        ##QbarAW
        QbarAWSL_m <- SuperLearner(
          Y = At_mix$y_scaled,
          X = X_train_mix,
          SL.library = SL.library,
          family = family,
          verbose = FALSE
        )
        
        X_m1 <- X_m0 <- X_train_mix
        X_m1$A_mix <- 1 # under exposure
        X_m0$A_mix <- 0 # under control
        
        QbarAW <- predict(QbarAWSL_m, newdata = X_train_mix)$pred
        Qbar1W <- predict(QbarAWSL_m, newdata = X_m1)$pred
        Qbar0W <- predict(QbarAWSL_m, newdata = X_m0)$pred
        
        QbarAW[QbarAW >= 1.0] <- 0.999
        QbarAW[QbarAW <= 0] <- 0.001
        
        Qbar1W[Qbar1W >= 1.0] <- 0.999
        Qbar1W[Qbar1W <= 0] <- 0.001
        
        Qbar0W[Qbar0W >= 1.0] <- 0.999
        Qbar0W[Qbar0W <= 0] <- 0.001
        
        # if(mean(Qbar1W - Qbar0W) > 0){
        #   mix_directions[[interaction]] <- "positive"
          
        X_m1 <- X_m0 <- X_valid_mix
        X_m1$A_mix <- 1 # under exposure
        X_m0$A_mix <- 0 # under control
        
        QbarAW <- predict(QbarAWSL_m, newdata = X_valid_mix)$pred
        Qbar1W <- predict(QbarAWSL_m, newdata = X_m1)$pred
        Qbar0W <- predict(QbarAWSL_m, newdata = X_m0)$pred
        
        QbarAW[QbarAW >= 1.0] <- 0.999
        QbarAW[QbarAW <= 0] <- 0.001
        
        Qbar1W[Qbar1W >= 1.0] <- 0.999
        Qbar1W[Qbar1W <= 0] <- 0.001
        
        Qbar0W[Qbar0W >= 1.0] <- 0.999
        Qbar0W[Qbar0W <= 0] <- 0.001
        
        ## add Qbar to the AV dataset
        Av_mix$QbarAW <- QbarAW
        Av_mix$Qbar1W <- Qbar1W
        Av_mix$Qbar0W <- Qbar0W
        
        Av_mix$gHat1W <- gHat1W
        Av_mix$H.AW <- H.AW
        
        mix_interaction_data[[interaction]] <- Av_mix
          
        # }else{
        #   mix_directions[[interaction]] <- "negative"
        #   
        #   At_mix$A_mix <- abs(interaction_rule - 1)
        #   Av_mix$A_mix <- abs(Av_rules_eval[, interaction]-1)
        #   
        #   ## select covariates
        #   X_Amix_T <- subset(At_mix,
        #                      select = c(covariates))
        #   
        #   X_Amix_V <- subset(Av_mix,
        #                      select = c(covariates))
        #   
        #   print("Fitting SL to marginal rule for mixture")
        #   
        #   gHatSL <- SuperLearner(
        #     Y = At_mix$A_mix,
        #     X = X_Amix_T,
        #     SL.library = SL.library,
        #     family = "binomial",
        #     verbose = FALSE
        #   )
        #   
        #   gHat1W <- predict(gHatSL, newdata = X_Amix_V)$pred
        #   gHat0W <- 1 - gHat1W
        #   
        #   gHatAW <- rep(NA, n)
        #   gHatAW[Av_mix$A_mix == 1] <- gHat1W[Av_mix$A_mix == 1]
        #   gHatAW[Av_mix$A_mix == 0] <- gHat0W[Av_mix$A_mix == 0]
        #   
        #   H.AW <-
        #     as.numeric(Av_mix$A_mix == 1) / gHat1W - as.numeric(Av_mix$A_mix == 0) /
        #     gHat0W
        #   
        #   H.AW <-
        #     ifelse(H.AW > H.AW_trunc_lvl, H.AW_trunc_lvl, H.AW)
        #   
        #   H.AW <-
        #     ifelse(H.AW < -H.AW_trunc_lvl, -H.AW_trunc_lvl, H.AW)
        #   
        #   ## add treatment mechanism results to validation dataframe
        #   
        #   
        #   X_train_mix <-
        #     subset(At_mix, select = c("A_mix", covariates))
        #   
        #   X_valid_mix <-
        #     subset(Av_mix, select = c("A_mix", covariates))
        #   
        #   print(paste(
        #     "Fitting SL Y given W and marginal rule for mixture",
        #     interaction
        #   ))
        #   
        #   ##QbarAW
        #   QbarAWSL_m <- SuperLearner(
        #     Y = At_mix$y_scaled,
        #     X = X_train_mix,
        #     SL.library = SL.library,
        #     family = family,
        #     verbose = FALSE
        #   )
        #   
        #   X_m1 <- X_m0 <- X_valid_mix
        #   X_m1$A_mix <- 1 # under exposure
        #   X_m0$A_mix <- 0 # under control
        #   
        #   QbarAW <- predict(QbarAWSL_m, newdata = X_valid_mix)$pred
        #   Qbar1W <- predict(QbarAWSL_m, newdata = X_m1)$pred
        #   Qbar0W <- predict(QbarAWSL_m, newdata = X_m0)$pred
        #   
        #   QbarAW[QbarAW >= 1.0] <- 0.999
        #   QbarAW[QbarAW <= 0] <- 0.001
        #   
        #   Qbar1W[Qbar1W >= 1.0] <- 0.999
        #   Qbar1W[Qbar1W <= 0] <- 0.001
        #   
        #   Qbar0W[Qbar0W >= 1.0] <- 0.999
        #   Qbar0W[Qbar0W <= 0] <- 0.001
        #   
        #   ## add Qbar to the AV dataset
        #   Av_mix$QbarAW <- QbarAW
        #   Av_mix$Qbar1W <- Qbar1W
        #   Av_mix$Qbar0W <- Qbar0W
        #   
        #   Av_mix$gHat1W <- gHat1W
        #   Av_mix$H.AW <- H.AW
        #   
        #   mix_interaction_data[[interaction]] <- Av_mix
        # }
        
        
      } else{
        print(
          "Only one A level resulted when applying training rules to validation set, no estimates made"
        )
        Av_mix$gHat1W <- NA
        Av_mix$H.AW <- NA
        Av_mix$QbarAW <- NA
        Av_mix$Qbar1W <- NA
        Av_mix$Qbar0W <- NA
        mix_interaction_data[[interaction]] <- Av_mix
        
        mix_directions <- NA
      }
    }
  } else{
    print(
      "No ATEs calculated in the validation set for the mixture rule due to no rule found in training set for mixture"
    )
    
    mix_interaction_data <- list()
    
    Av_mix$gHat1W <- NA
    Av_mix$H.AW <- NA
    Av_mix$QbarAW <- NA
    Av_mix$Qbar1W <- NA
    Av_mix$Qbar0W <- NA
    mix_interaction_data[[1]] <- Av_mix
    
    mix_directions <- NA
  }

  
 return(list(data =  mix_interaction_data)) 
}
