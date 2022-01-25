#' Aggregate mixture rules found across the folds that have the same variables and directions. For each rule
#' extract the relevant nuisance parameter data calculated in the folds. Given the validation data estimates across
#' the folds, for each tree do a TMLE update step to target the average treatment effect. Update the initial counterfactuals,
#' calculate the influence curve and using the influence curve calculate variance estimates and p-values.
#'
#' @param input_mix_rules List of dataframes of rules found for a mixture across the folds
#' @param input_mix_data Nuisance parameter data for mixture rules found across the folds
#' @param outcome Character indicating the outcome variable
#' @param n_folds Number of folds used in cross-validation
#' @param no_mixture_rules TRUE/FALSE if no mixture rules were found across the folds


#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by

#' @return Rules object. TODO: add more detail here.
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export

calc_v_fold_mixtures_ate <- function(input_mix_rules, input_mix_data, outcome, n_folds, no_mixture_rules) {
  input_mix_rules <- unlist(input_mix_rules, recursive = FALSE, use.names = FALSE)
  input_mix_rules <- input_mix_rules[!sapply(input_mix_rules, is.null)]

  input_mix_data <- unlist(input_mix_data, recursive = FALSE, use.names = FALSE)
  input_mix_data <- unlist(input_mix_data, recursive = FALSE, use.names = FALSE)
  input_mix_data <- input_mix_data[!sapply(input_mix_data, is.null)]


  if (no_mixture_rules == FALSE) {
    mixture_results <-
      as.data.frame(matrix(
        data = NA,
        nrow = length(input_mix_data),
        ncol = 9
      ))

    colnames(mixture_results) <-
      c(
        "Mixture ATE",
        "Standard Error",
        "Lower CI",
        "Upper CI",
        "P-value",
        "P-value Adj",
        "RMSE",
        "Mixture Interaction Rules",
        "Variables"
      )

    mixture_data_list <- list()

    for (group in seq(input_mix_data)) {
      mix_data <- as.data.frame(input_mix_data[[group]])

      flux_results <- fit_least_fav_submodel(mix_data$H.AW, mix_data, mix_data$QbarAW, mix_data$Qbar1W, mix_data$Qbar0W)
      QbarAW.star <- flux_results$QbarAW.star
      Qbar1W.star <- flux_results$Qbar1W.star
      Qbar0W.star <- flux_results$Qbar0W.star

      ## back-scale Y
      mix_data$QbarAW.star <- scale_to_original(scaled_vals = QbarAW.star, max_orig = max(mix_data[outcome]), min(mix_data[outcome]))
      mix_data$Qbar0W.star <- scale_to_original(scaled_vals = Qbar0W.star, max_orig = max(mix_data[outcome]), min(mix_data[outcome]))
      mix_data$Qbar1W.star <- scale_to_original(scaled_vals = Qbar1W.star, max_orig = max(mix_data[outcome]), min(mix_data[outcome]))

      ATE_results <- calc_ATE_estimates(data = mix_data, ATE_var = "mix.ATE", outcome = outcome, p_adjust_n = length(input_mix_data), v_fold = TRUE)

      ## TODO:
      ## Figure out why the scale_to_original function still outputs NA when the qnorm should be scaled to 0-1.

      ## calculate RMSE for Y| A = rule i, W
      sqrd_resids <- (mix_data$QbarAW.star - mix_data[outcome])^2
      RMSE <- sqrt(mean(sqrd_resids[, 1], na.rm = TRUE))

      mixture_results$`Mixture ATE`[[group]] <- round(ATE_results$ATE, 3)
      mixture_results$`Standard Error`[group] <- round(ATE_results$SE, 3)
      mixture_results$`Lower CI`[group] <- round(ATE_results$CI[1], 3)
      mixture_results$`Upper CI`[group] <- round(ATE_results$CI[2], 3)
      mixture_results$`P-value`[group] <- round(ATE_results$`p-value`, 6)
      mixture_results$`P-value Adj`[group] <- round(ATE_results$`adj p-value`, 6)
      mixture_results$`RMSE`[group] <- round(RMSE, 3)

      mixture_data_list[[group]] <- mix_data
    }

    mixture_results[, "Mixture Interaction Rules"] <- do.call(rbind, input_mix_rules)$description
    mixture_results[, "Variables"] <- do.call(rbind, input_mix_rules)$test
  } else {
    mixture_results <-
      as.data.frame(matrix(
        data = NA,
        nrow = 1,
        ncol = 9
      ))

    colnames(mixture_results) <-
      c(
        "Mixture ATE",
        "Standard Error",
        "Lower CI",
        "Upper CI",
        "P-value",
        "P-value Adj",
        "RMSE",
        "Mixture Interaction Rules",
        "Variables"
      )


    mixture_data_list <- NA
  }
  return(list(
    results = mixture_results,
    mixture_data_list = mixture_data_list
  ))
}
