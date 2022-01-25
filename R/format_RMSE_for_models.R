#' @title Creates the RMSE output
#' @description Create a simple data frame with the marginal, mixture, additive marginal and non-additive marginal models and their root mean squared error.
#'
#' @param marginal_results Marginal results dataframe
#' @param mixture_results Mixture results dataframe
#' @param additive_RMSE_star Additive marginal model RMSE
#' @param marg_combo_data Non-additive marginal dataframe
#'
#' @export
#'
format_RMSE_for_models <- function(marginal_results, mixture_results, additive_RMSE_star, marg_combo_data) {
  marginal_rules_RMSE <- cbind("Marginal Model", cbind(rownames(marginal_results), marginal_results$RMSE))
  mixture_rules_RMSE <- cbind("Mixture Model", cbind(mixture_results$`Mixture Interaction Rules`, mixture_results$RMSE))
  # additive_RMSE <- cbind("Additive Model", cbind("additive marginal model", additive_RMSE_star))
  non_additive_sqr_resids <- (marg_combo_data$QbarAW_combo - marg_combo_data$raw_outcome)^2
  non_additive_RMSE <- sqrt(mean(non_additive_sqr_resids))
  non_additive_RMSE <- cbind("Combination Model", cbind("non-additive marginal model", non_additive_RMSE))

  models_RMSE <- as.data.frame(rbind(marginal_rules_RMSE, mixture_rules_RMSE, non_additive_RMSE))

  colnames(models_RMSE) <- c("Model Type", "Rule", "RMSE")
  models_RMSE$RMSE <- as.numeric(models_RMSE$RMSE)

  return(models_RMSE)
}
