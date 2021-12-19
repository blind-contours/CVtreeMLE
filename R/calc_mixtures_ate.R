#' Aggregate mixture rules found across the folds that have the same variables and directions. For each rule
#' extract the relevant nuisance parameter data calculated in the folds. Given the validation data estimates across
#' the folds, for each tree do a TMLE update step to target the average treatment effect. Update the initial counterfactuals,
#' calculate the influence curve and using the influence curve calculate variance estimates and p-values.
#'
#' @param input_mix_rules List of dataframes of rules found for a mixture across the folds
#' @param input_mix_data Nuisance parameter data for mixture rules found across the folds
#' @param outcome Character indicating the outcome variable
#' @param n_folds Number of folds used in cross-validation


#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by

#' @return Rules object. TODO: add more detail here.
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export

calc_mixtures_ate <- function(input_mix_rules, input_mix_data, outcome, n_folds){

  input_mix_rules <- unlist(input_mix_rules,recursive=FALSE, use.names = FALSE)
  input_mix_rules <- input_mix_rules[!sapply(input_mix_rules,is.null)]

  input_mix_data <- unlist(input_mix_data,recursive=FALSE, use.names = FALSE)
  input_mix_data <- input_mix_data[!sapply(input_mix_data,is.null)]


  if(anyNA(unlist(input_mix_rules)) == FALSE){

    fold_mix_rules <-
      data.table::rbindlist(input_mix_rules)
    fold_mix_rules <-
      fold_mix_rules %>% dplyr::group_by(test, direction) %>% dplyr::filter(n() >= n_folds)

    groups <- fold_mix_rules %>%
      dplyr::group_by(test, direction)

    group_list <- dplyr::group_split(groups)

    mixture_results <-
      as.data.frame(matrix(
        data = NA,
        nrow = length(group_list),
        ncol = 5
      ))

    colnames(mixture_results) <-
      c("Mixture ATE",
        "Standard Error",
        "Lower CI",
        "Upper CI",
        "P-value")

    for (group in seq(group_list)) {
      intx_group <- group_list[[group]]
      intxn_rule_data_list <- list()
      vars <- intx_group$test[1]

      for (i in seq(dim(intx_group)[1])) {
        intxn_rule <- intx_group[i,]
        search_data <-
          as.data.frame(input_mix_rules[as.numeric(intxn_rule$fold)])
        srch_indx <- match(intxn_rule$rule , search_data$rule)
        fold_data <- input_mix_data[[as.numeric(intxn_rule$fold)]]
        intx_rule_data <- fold_data[[srch_indx]]
        intxn_rule_data_list[[i]] <- intx_rule_data
      }

      # Extract the results from each CV-TMLE fold and rbind into a single dataframe.
      mix_data = do.call(rbind, intxn_rule_data_list)

      flux_results <- fit_least_fav_submodel(mix_data$H.AW, mix_data, mix_data$QbarAW, mix_data$Qbar1W, mix_data$Qbar0W)
      QbarAW.star <- flux_results$QbarAW.star
      Qbar1W.star <- flux_results$Qbar1W.star
      Qbar0W.star <- flux_results$Qbar0W.star

      ## back-scale Y
      mix_data$QbarAW.star <- scale_to_original(scaled_vals = QbarAW.star, max_orig = max(mix_data[outcome]), min(mix_data[outcome]))
      mix_data$Qbar0W.star <- scale_to_original(scaled_vals = Qbar0W.star, max_orig = max(mix_data[outcome]), min(mix_data[outcome]))
      mix_data$Qbar1W.star <- scale_to_original(scaled_vals = Qbar1W.star, max_orig = max(mix_data[outcome]), min(mix_data[outcome]))

      ATE_results <- calc_ATE_estimates(data = mix_data, ATE_var = "mix.ATE", outcome = outcome, p_adjust_n = length(group_list))

      ## calculate RMSE
      RMSE <-
        sqrt((mix_data$QbarAW.star - mix_data[outcome]) ^ 2 / length(mix_data[outcome]))
      RMSE <- mean(RMSE[, 1])

      mixture_results$`Mixture ATE`[group] <- ATE_results$ATE
      mixture_results$`Standard Error`[group] <- ATE_results$SE
      mixture_results$`Lower CI`[group] <- ATE_results$CI[1]
      mixture_results$`Upper CI`[group] <-ATE_results$CI[2]
      mixture_results$`P-value`[group] <- round(ATE_results$`p-value`,6)
      mixture_results$`P-value Adj`[group] <- round(ATE_results$`adj p-value`, 6)
      mixture_results$`vars`[group] <- vars

    }
  }else{
    mixture_results <-
      as.data.frame(matrix(
        data = NA,
        nrow = 1,
        ncol = 5
      ))

    colnames(mixture_results) <-
      c("Mixture ATE",
        "Standard Error",
        "Lower CI",
        "Upper CI",
        "P-value")

    mixture_results$`Mixture ATE` <- NA
    mixture_results$`Standard Error` <- NA
    mixture_results$`Lower CI` <- NA
    mixture_results$`Upper CI` <- NA
    mixture_results$`P-value` <- NA
    mixture_results$`P-value Adj` <- NA

    group_list <- NA

    RMSE <- NA

  }
  return(list(results = mixture_results, group_list = group_list, mix_RMSE = RMSE))
}
