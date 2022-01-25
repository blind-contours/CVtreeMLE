#' Evaluate mixture rules found during the rpart decision tree process
#'
#' @param additive_data Data for the cumulative exposure across the folds for each fold specific set of rules
#' @param n_folds Number of folds used in the cross-validation
#' @param Y Character indicating the Y of interest
#' @importFrom dplyr transmute
#' @return Rules object. TODO: add more detail here.
#' @importFrom rlang :=
#'
#' @export

calc_additive_ate <- function(additive_data, Y, n_folds) {

  ## TODO: figure out why when using doParallel sometimes there is a NULL list appended to the list of lists
  additive_data <- unlist(additive_data, recursive = FALSE)
  additive_data <- additive_data[!sapply(additive_data, is.null)]

  mix_additive_data <- do.call(rbind, lapply(1:length(additive_data), function(i) {
    fold <- additive_data[[i]]
    fold
  }))

  if (any(is.na(mix_additive_data$QbarAW_additive)) == FALSE) {

    ## least optimal submodel
    logitUpdate <-
      stats::glm(y_scaled ~ -1 + HAW_additive + offset(qlogis(bound_precision(QbarAW_additive))),
        family = "quasibinomial",
        data = mix_additive_data
      )

    epsilon <- logitUpdate$coef

    QbarAW_additive_star <-
      stats::plogis(
        stats::qlogis(bound_precision(mix_additive_data$QbarAW_additive)) + epsilon * mix_additive_data$HAW_additive
      )

    mix_additive_data$QbarAW_additive_star <- QbarAW_additive_star

    mix_additive_data$QbarAW_additive <- scale_to_original(mix_additive_data$QbarAW_additive, max(mix_additive_data[Y]), min(mix_additive_data[Y]))
    mix_additive_data$QbarAW_additive_star <- scale_to_original(mix_additive_data$QbarAW_additive_star, max(mix_additive_data[Y]), min(mix_additive_data[Y]))


    ## Calculate Additive RMSE for non-updated predictions

    sqrd_resids <- (mix_additive_data$QbarAW_additive_star - mix_additive_data[Y])^2
    cum_sum_RMSE <- sqrt(mean(sqrd_resids[, 1]))


    # cum_sum_RMSE <-
    #   sqrt((mix_additive_data$QbarAW_additive - mix_additive_data[Y]) ^
    #          2 / length(mix_additive_data[Y])
    #   )
    # additive.MSM.RMSE <- sum(additive.MSM.RMSE[, 1])



    #
    #     updated_additive.MSM.RMSE <-
    #       sqrt((mix_additive_data$QbarAW_additive_star - mix_additive_data[Y])^2/ length(mix_additive_data[Y]))
    #     updated_additive.MSM.RMSE <-
    #       sum(updated_additive.MSM.RMSE[, 1])


    ################################ RUN MSM THROUGH THE CUMULATIVE SUM EXPECTATIONS ############################

    mix_additive_data$sum_marg_hits <- as.factor(mix_additive_data$sum_marg_hits)
    MSM_results <- stats::glm(QbarAW_additive_star ~ sum_marg_hits, data = mix_additive_data)
  } else {
    print("No additive results found")
    additive.MSM.RMSE <- NA
    updated_additive.MSM.RMSE <- NA
    MSM_results <- NA
  }

  return(list(
    data = mix_additive_data,
    RMSE_star = cum_sum_RMSE,
    MSM = MSM_results
  ))
}
