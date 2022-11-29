#' @title Calculate the NDE for each mixture rule
#'
#' @description Aggregate mixture rules found across the folds that have
#' the same variables. For each rule
#' extract the relevant nuisance parameter data calculated in the folds.
#' Given the validation data estimates across
#' the folds, for each tree do a TMLE update step to target the natural
#' direct effect. Update the initial counterfactuals,
#' calculate the influence curve and using the influence curve calculate
#' variance estimates and p-values.
#'
#' @param input_mix_rules List of dataframes of rules found for a mixture
#' across the folds
#' @param input_mix_data Nuisance parameter data for mixture rules found
#' @param y Outcome variable
#' across the folds
#' @param n_folds Number of folds used in cross-validation
#' @param no_mixture_rules TRUE/FALSE whether no mixture rules were found
#' across all the folds
#'
#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by
#'
#' @return A list with mixture analysis results which includes:
#' \itemize{
#'   \item \code{results}: A data frame with variable threshold combinations
#'   on the rows and ATE, variance and consistency estimates on the columns.
#'   \item \code{group_list}: A list of rule combinations found in the
#'   ensemble decision tree model grouped by variable sets in the rules and
#'   directions of the coefficient in the linear model. Also provided is the
#'   fold the rule was found and the RMSE of the model which used the rule.
#'   \item \code{mixture_data_list}: A list of data frames which houses the
#'   data for each mixture rule evaluated as an exposure, the baseline
#'   covariates, outcome, and nuisance parameter estimates for the respective
#'   rule.
#' }
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm
#' @importFrom stats qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export

calc_mixtures_nde <- function(input_mix_rules,
                              input_mix_data,
                              y,
                              n_folds,
                              no_mixture_rules) {
  fold_mix_rules <-
    data.table::rbindlist(unlist(input_mix_rules, recursive = FALSE))

  fold_mix_rules <-
    fold_mix_rules %>%
    dplyr::group_by(test, direction)

  groups <- fold_mix_rules %>%
    dplyr::group_by(test, direction)

  group_list <- dplyr::group_split(groups)

  nde_mixture_results <-
    as.data.frame(matrix(
      data = NA,
      nrow = length(group_list),
      ncol = 8
    ))

  colnames(nde_mixture_results) <-
    c(
      "Mixture NDE",
      "Standard Error",
      "Lower CI",
      "Upper CI",
      "P-value",
      "P-value Adj",
      "Vars",
      "RMSE"
    )

  nde_mixture_data_list <- list()

  for (group in seq(group_list)) {
    intx_group <- group_list[[group]]
    intxn_rule_data_list <- list()
    vars <- unique(intx_group$test)

    for (i in seq(dim(intx_group)[1])) {
      intxn_rule <- intx_group[i, ]
      search_data <-
        as.data.frame(input_mix_rules[as.numeric(intxn_rule$fold)])
      srch_indx <- match(intxn_rule$rule, search_data[, 1])
      fold_data <- input_mix_data[[as.numeric(intxn_rule$fold)]]
      intx_rule_data <- fold_data[[1]][srch_indx]
      intxn_rule_data_list[[i]] <- intx_rule_data
    }

    mix_rule_data <- do.call(rbind, unlist(intxn_rule_data_list,
      recursive = FALSE
    ))
    logit_update <-
      stats::glm(
        q_diff ~ -1 + hz + offset(qdiff_w),
        data = mix_rule_data
      )

    epsilon <- logit_update$coef
    qbar_diff_star <- mix_rule_data$qdiff_w + epsilon * mix_rule_data$hz

    av_treatment_indicator <- mix_rule_data$A_mix
    av_control_indicator <- 1 - mix_rule_data$A_mix

    # compute individual scores for DY, DA, DW
    d_y <- mix_rule_data$hy * (mix_rule_data[, y] - mix_rule_data$qbar_azw_star)
    d_z <- av_control_indicator * mix_rule_data$hz * (mix_rule_data$qbar_1zw_star - mix_rule_data$qbar_0zw_star - qbar_diff_star)
    d_w <- qbar_diff_star

    # parameter and influence function
    theta <- mean(qbar_diff_star)
    eif <- d_y + d_z + d_w - theta

    n <- dim(mix_rule_data)[1]
    varhat_ic <- stats::var(eif, na.rm = TRUE) / n
    se <- sqrt(varhat_ic)

    alpha <- 0.05

    # obtain 95% two-sided confidence intervals:
    ci <- c(
      theta + stats::qnorm(alpha / 2, lower.tail = TRUE) * se,
      theta + stats::qnorm(alpha / 2, lower.tail = FALSE) * se
    )

    # p-value
    p_value <- 2 * stats::pnorm(abs(theta / se), lower.tail = FALSE)

    p_adjust_n <- length(group_list)

    p_value_adjust <-
      stats::p.adjust(p_value, method = "bonferroni", n = p_adjust_n)


    ## calculate RMSE for Y| A = rule i, W
    sqrd_resids <- (mix_rule_data$qbar_azw_star - mix_rule_data[y])^2
    rmse <- sqrt(mean(sqrd_resids[, 1]))

    nde_mixture_results$`Mixture NDE`[group] <- round(theta, 3)
    nde_mixture_results$`Standard Error`[group] <- round(se, 3)
    nde_mixture_results$`Lower CI`[group] <- round(ci[1], 3)
    nde_mixture_results$`Upper CI`[group] <- round(ci[2], 3)
    nde_mixture_results$`P-value`[group] <- round(p_value, 6)
    nde_mixture_results$`P-value Adj`[group] <- round(
      p_value_adjust,
      6
    )
    nde_mixture_results$`Vars`[group] <- vars
    nde_mixture_results$`RMSE`[group] <- round(rmse, 3)

    nde_mixture_data_list[[group]] <- mix_rule_data
  }

  return(list(
    results = nde_mixture_results,
    group_list = group_list,
    mixture_data_list = nde_mixture_data_list
  ))
}
