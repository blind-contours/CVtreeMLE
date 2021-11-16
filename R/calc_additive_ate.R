calc_additive_ate <- function(additive_data, outcome){
  
  ## TODO: figure out why when using doParallel sometimes there is a NULL list appended to the list of lists
  additive_data <- unlist(additive_data,recursive=FALSE)
  additive_data <- additive_data[!sapply(additive_data,is.null)]
  
  if (length(additive_data) == n_folds) {
    mix_additive_data = do.call(rbind, lapply(1:length(additive_data), function(i) {
      fold = additive_data[[i]]
      fold
    }))
    
    ## least optimal submodel
    logitUpdate <-
      glm(y_scaled ~ -1 + HAW_additive + offset(qlogis(QbarAW_additive)) ,
          family = 'quasibinomial',
          data = mix_additive_data)
    
    epsilon <- logitUpdate$coef
    
    QbarAW_additive_star <-
      plogis(
        qlogis(mix_additive_data$QbarAW_additive) + epsilon * mix_additive_data$HAW_additive
      )
    
    mix_additive_data$QbarAW_additive_star <- QbarAW_additive_star
   
    mix_additive_data$QbarAW_additive <-
      with(mix_additive_data,
           QbarAW_additive * (max(mix_additive_data[outcome]) - min(mix_additive_data[outcome])) + min(mix_additive_data[outcome]))
    
    ## Calculate Additive RMSE for non-updated predictions
    additive.MSM.RMSE <-
      sqrt((mix_additive_data$QbarAW_additive - mix_additive_data[outcome]) ^
             2 / length(mix_additive_data[outcome])
      )
    additive.MSM.RMSE <- mean(additive.MSM.RMSE[, 1])
    
    
    ## Calculate Additive RMSE for updated predictions
    mix_additive_data$QbarAW_additive_star <-
      with(mix_additive_data,
           QbarAW_additive_star * (max(mix_additive_data[outcome]) - min(mix_additive_data[outcome])) + min(mix_additive_data[outcome]))
    
    updated_additive.MSM.RMSE <-
      sqrt((
        mix_additive_data$QbarAW_additive_star - mix_additive_data[outcome]
      ) ^ 2 / length(mix_additive_data[outcome])
      )
    updated_additive.MSM.RMSE <-
      mean(updated_additive.MSM.RMSE[, 1])
    
    #############################################################################################################
    ################################ RUN MSM THROUGH THE CUMULATIVE SUM EXPECTATIONS ############################
    #############################################################################################################
    mix_additive_data$sum_marg_hits <- as.factor(mix_additive_data$sum_marg_hits)
    MSM_results <- glm(QbarAW_additive_star~sum_marg_hits, data = mix_additive_data)
    
  } else{
    print('No additive results found')
    additive.MSM.RMSE <- NA
    updated_additive.MSM.RMSE <- NA
    MSM_results <- NA
    
  }
  
  return(list(data = mix_additive_data, 
              RMSE = additive.MSM.RMSE, 
              RMSE_star = updated_additive.MSM.RMSE,
              MSM = MSM_results))
  
}
