#' @title Calculate the ATE for each rule found for individual mixture
#' components.
#' @description Aggregate marginal rules for each mixture variable found
#' across the folds. For each rule
#' extract the relevant nuisance parameter data calculated in the folds. Given
#' the validation data estimates across
#' the folds, for each tree do a TMLE update step to target the average
#'  treatment effect. Update the initial counterfactuals,
#' calculate the influence curve and using the influence curve calculate
#' variance estimates and p-values.
#'
#' @param marginal_data List of dataframes of nuisance parameter data for
#' each mixture
#' @param mix_comps Vector of characters indicating the mixture components
#' @param marginal_rules List of dataframes for rules found across the folds
#' @param y Vector indicating the Y
#' @param n_folds Number of folds used in cross-validation
#'
#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows
#'
#' @return A list of the marginal results including:
#' ' \itemize{
#'   \item \code{marginal_results}: A data frame with the data adpatively
#'   determined mixture component thresholds on the rows and ATE, variance,
#'   and RMSE estimates on the columns.
#'   \item \code{data}: A list of data frames for each mixture component
#'   threshold evaluated as the exposure, baseline covariates, outcome,
#'   nuisance parameter estimates, marginal ATE and the influence curve.
#'   }
#' @importFrom stats as.formula glm p.adjust
#' @importFrom stats plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export

calc_marginal_ate <- function(marginal_data,
                              mix_comps,
                              marginal_rules,
                              y,
                              n_folds) {
  marg_data <- list()

  m_data <- unlist(marginal_data, recursive = FALSE)
  m_rules <- unlist(marginal_rules, recursive = FALSE)

  input <- do.call(rbind, m_rules)

  fold_rules_groups <- marginal_group_split(input)

  comp_labels <- list()

  for (i in seq(fold_rules_groups)) {
    var_comps <- list()
    rules <- fold_rules_groups[[i]]
    reference <- rules[rules$quantile == 1, ]
    reference_quant <- unique(reference$var_quant_group)
    comparisons <- rules[rules$quantile > 1, ]
    comparisons <- comparisons[!duplicated(comparisons$var_quant_group), ]

    for (j in seq(nrow(comparisons))) {
      comp_row <- comparisons[j, ]
      comp_label <- paste(comp_row$var_quant_group, reference_quant, sep = "-")
      var_comps[[j]] <- comp_label
    }

    comp_labels[[i]] <- var_comps
  }

  comp_labels <- unlist(comp_labels)

  all_fold_var_quants <- unique(unlist(unique(sapply(m_data, names))))

  all_fold_var_quants <- all_fold_var_quants[!is.na(all_fold_var_quants)]

  for (i in seq(all_fold_var_quants)) {
    var_quant <- all_fold_var_quants[i]
    marg_data[[var_quant]] <-
      dplyr::bind_rows(sapply(m_data, "[", var_quant))
  }

  marginal_results <-
    as.data.frame(matrix(
      data = NA,
      nrow = length(marg_data),
      ncol = 7
    ))

  colnames(marginal_results) <-
    c(
      "Marginal ATE",
      "Standard Error",
      "Lower CI",
      "Upper CI",
      "P-value",
      "P-value Adj",
      "rmse"
    )

  rownames(marginal_results) <- comp_labels

  updated_marginal_data <- list()

  no_null_indices <- list()

  for (i in seq(length(marg_data))) {
    check <- marg_data[[i]]
    if (anyNA(check$QbarAW) != TRUE) {
      no_null_indices <- append(no_null_indices, i)
    } else {
      no_null_indices <- no_null_indices
    }
  }

  no_null_indices <- unlist(no_null_indices)

  for (i in seq(marg_data)) {
    marg_mix <- marg_data[[i]]
    flux_results <- fit_least_fav_submodel(h_aw = marg_mix$h_aw,
                                           data = marg_mix,
                                           qbar_aw = marg_mix$qbar_aw,
                                           qbar_1w = marg_mix$qbar_1w,
                                           qbar_0w = marg_mix$qbar_0w)

    ## back-scale Y
    qbar_aw_star <- scale_to_original(scaled_vals = flux_results$qbar_aw_star,
                                     max_orig = max(marg_mix[y]),
                                     min_orig = min(marg_mix[y]))

    qbar_0w_star <- scale_to_original(scaled_vals = flux_results$qbar_0w_star,
                                     max_orig = max(marg_mix[y]),
                                     min_orig = min(marg_mix[y]))

    qbar_1w_star <- scale_to_original(scaled_vals = flux_results$qbar_1w_star,
                                     max_orig = max(marg_mix[y]),
                                     min_orig = min(marg_mix[y]))

    marg_mix$qbar_aw_star <- qbar_aw_star
    marg_mix$qbar_0w_star <- qbar_0w_star
    marg_mix$qbar_1w_star <- qbar_1w_star

    ate_results <- calc_ate_estimates(data = marg_mix,
                                      ate_var = "marg_ATE",
                                      outcome = y,
                                      p_adjust_n = length(no_null_indices))

    sqrd_resids <- (marg_mix$qbar_aw_star - marg_mix[y])^2
    rmse <- sqrt(mean(sqrd_resids[, 1]))

    marginal_results$`Marginal ATE`[i] <- round(ate_results$ate, 3)
    marginal_results$`Standard Error`[i] <- round(ate_results$se, 3)
    marginal_results$`Lower CI`[i] <- round(ate_results$ci[1], 3)
    marginal_results$`Upper CI`[i] <- round(ate_results$ci[2], 3)
    marginal_results$`P-value`[i] <- round(ate_results$p_value, 6)
    marginal_results$`P-value Adj`[i] <- round(ate_results$adj_p_value, 6)
    marginal_results$rmse[i] <- round(rmse, 3)

    updated_marginal_data[[i]] <- ate_results$data
  }

  return(list(marginal_results = marginal_results,
              data = updated_marginal_data))
}
