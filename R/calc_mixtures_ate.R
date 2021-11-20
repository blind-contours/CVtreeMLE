#' Aggregate mixture rules found across the folds that have the same variables and directions. For each rule
#' extract the relevant nuisance parameter data calculated in the folds. Given the validation data estimates across
#' the folds, for each tree do a TMLE update step to target the average treatment effect. Update the initial counterfactuals,
#' calculate the influence curve and using the influence curve calculate variance estimates and p-values.
#'
#' @param input_mix_rules List of dataframes of rules found for a mixture across the folds
#' @param input_mix_data Nuisance parameter data for mixture rules found across the folds
#' @param outcome Character indicating the outcome variable

#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by

#' @return Rules object. TODO: add more detail here.
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export

calc_mixtures_ate <- function(input_mix_rules, input_mix_data, outcome){

  input_mix_rules <- unlist(input_mix_rules,recursive=FALSE, use.names = FALSE)
  input_mix_rules <- input_mix_rules[!sapply(input_mix_rules,is.null)]

  input_mix_data <- unlist(input_mix_data,recursive=FALSE, use.names = FALSE)
  input_mix_data <- input_mix_data[!sapply(input_mix_data,is.null)]


  if(anyNA(unlist(input_mix_rules)) == FALSE){

    fold_mix_rules <-
      data.table::rbindlist(input_mix_rules, idcol = "fold", fill=TRUE)
    fold_mix_rules <-
      fold_mix_rules %>% group_by(test, direction) %>% filter(n() >= n_folds)

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

      ## least optimal submodel
      logitUpdate <-
        glm(y_scaled ~ -1 + H.AW + offset(qlogis(QbarAW)) ,
            family = 'quasibinomial',
            data = mix_data)

      epsilon <- logitUpdate$coef

      QbarAW.star <-
        stats::plogis(stats::qlogis(mix_data$QbarAW) + epsilon * mix_data$H.AW)
      Qbar1W.star <-
        stats::plogis(stats::qlogis(mix_data$Qbar1W) + epsilon * mix_data$H.AW)
      Qbar0W.star <-
        stats::plogis(stats::qlogis(mix_data$Qbar0W) + epsilon * mix_data$H.AW)

      ## back-scale Y
      QbarAW.star <-
        QbarAW.star * (max(mix_data[outcome]) - min(mix_data[outcome])) + min(mix_data[outcome])
      Qbar0W.star <-
        Qbar0W.star * (max(mix_data[outcome]) - min(mix_data[outcome])) + min(mix_data[outcome])
      Qbar1W.star <-
        Qbar1W.star * (max(mix_data[outcome]) - min(mix_data[outcome])) + min(mix_data[outcome])

      mix_data$QbarAW.star <- QbarAW.star
      mix_data$Qbar0W.star <- Qbar0W.star
      mix_data$Qbar1W.star <- Qbar1W.star

      mix_data$mix.ATE <- mix_data$Qbar1W.star - mix_data$Qbar0W.star

      Thetas <-
        tapply(mix_data$mix.ATE, mix_data$folds, mean, na.rm = TRUE)

      for (i in seq(Thetas)) {
        mix_data[mix_data$folds == i, "Thetas"] <- Thetas[i][[1]]
      }

      ICs = base::by(mix_data, mix_data$folds, function(mix_data) {
        result = with(mix_data,
                      H.AW * (mix_data[outcome] - QbarAW.star) + Qbar1W.star - Qbar0W.star - Thetas)
        result
      })

      for (i in seq(ICs)) {
        mix_data[mix_data$folds == i, "IC"] <- ICs[i]
      }

      n <- dim(mix_data)[1]
      varHat.IC <- stats::var(mix_data$IC, na.rm = TRUE) / n
      se <- sqrt(varHat.IC)

      alpha <- 0.05

      Theta <- mean(Thetas)
      # obtain 95% two-sided confidence intervals:
      CI <- c(Theta + stats::qnorm(alpha / 2, lower.tail = T) * se,
              Theta + stats::qnorm(alpha / 2, lower.tail = F) * se)

      # p-value
      p.value <- 2 * stats::pnorm(abs(Theta / se), lower.tail = F)

      p.value.adjust <-
        stats::p.adjust(p.value, method = "bonferroni", n = length(group_list))

      ## calculate RMSE
      RMSE <-
        sqrt((mix_data$QbarAW.star - mix_data[outcome]) ^ 2 / length(mix_data[outcome]))
      RMSE <- mean(RMSE[, 1])

      mixture_results$`Mixture ATE`[group] <- Theta
      mixture_results$`Standard Error`[group] <- se
      mixture_results$`Lower CI`[group] <- CI[1]
      mixture_results$`Upper CI`[group] <- CI[2]
      mixture_results$`P-value`[group] <- round(p.value,6)
      mixture_results$`P-value Adj`[group] <- round(p.value.adjust, 6)
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
