#' @title Compute v-fold specific estimates and do a meta-analysis type pooled
#' average ATE for each individual marginal rule.
#' @description For each fold, estimates the ATE for a fold specific mixture
#' rule. Also estimates a meta-analysis type average ATE and pooled variance.
#' Creates a union rule that covers all the folds in the rules.
#'
#' @param v_fold_marginal_results List of dataframes for v-fold specific
#' estimates of the fold-specific rule results for marginals
#' @param mix_comps Vector of characters indicating the mixture components
#' @param n_folds Number of folds used in cross-validation
#' @param data Input dataframe on which to evaluate the rules

#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows

#' @return Rules object. TODO: add more detail here.
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm
#' @importFrom stats qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export

compute_meta_marg_results <- function(v_fold_marginal_results,
                                      data,
                                      mix_comps,
                                      n_folds) {
  v_fold_marg_group <- v_fold_marginal_qgroup_split(v_fold_marginal_results)

  v_fold_marginal_w_pooled <- list()
  list_names <- list()

  for (i in seq(v_fold_marg_group)) {
    results_df <- v_fold_marg_group[[i]]

    var <- mix_comps[mix_comps %in% strsplit(results_df$Reference_Rule[1],
      split = " "
    )[[1]]]


    weighted_mean <- sum(results_df$`Marginal ATE` *
      (1 / results_df$`Standard Error`^2)) /
      sum((1 / results_df$`Standard Error`^2))

    weighted_rmse <- sum(results_df$RMSE *
      (1 / results_df$`Standard Error`^2)) /
      sum((1 / results_df$`Standard Error`^2))

    pooled_se <- sqrt(1 / (1 / sum(results_df$`Standard Error`^2)))

    pooled_p_val <- round(2 * stats::pnorm(
      abs(weighted_mean /
        pooled_se),
      lower.tail = FALSE
    ), 5)

    pooled_ci <- c(
      round(weighted_mean + stats::qnorm(
        0.05 /
          2,
        lower.tail = TRUE
      ) * pooled_se, 4),
      round(weighted_mean + stats::qnorm(
        0.05 /
          2,
        lower.tail = FALSE
      ) * pooled_se, 4)
    )

    ref_rule <- paste(results_df$Reference_Rule, collapse = "|")
    comp_rule <- paste(results_df$Comparison_Rule, collapse = "|")

    ref_data <- data %>%
      mutate("ref_rule" := ifelse(eval(parse(text = ref_rule)), 1, 0))

    comp_data <- data %>%
      mutate("comp_rule" := ifelse(eval(parse(text = comp_rule)), 1, 0))

    ## ref rule

    ref_min <-
      ref_data %>%
      dplyr::group_by(ref_rule) %>%
      summarise(min = min(!!(as.name(var))))

    ref_min <- subset(ref_min, ref_rule == 1, select = min)

    ref_max <-
      ref_data %>%
      group_by(ref_rule) %>%
      summarise(max = max(!!(as.name(var))))

    ref_max <- subset(ref_max, ref_rule == 1, select = max)

    ref_rule <- paste(
      var, " >= ", round(ref_min[[1]], 3), "&", var, "<=",
      round(ref_max[[1]], 3)
    )

    ## comp rule

    comp_min <-
      comp_data %>%
      dplyr::group_by(comp_rule) %>%
      summarise(min = min(!!(as.name(var))))

    comp_min <- subset(comp_min, comp_rule == 1, select = min)

    comp_max <-
      comp_data %>%
      group_by(comp_rule) %>%
      summarise(max = max(!!(as.name(var))))

    comp_max <- subset(comp_max, comp_rule == 1, select = max)

    comp_rule <- paste(
      var, " >= ", round(comp_min[[1]], 3), "&", var, "<=",
      round(comp_max[[1]], 3)
    )

    average_results <- cbind(
      round(weighted_mean, 3),
      round(pooled_se, 3),
      round(pooled_ci[1], 3),
      round(pooled_ci[2], 3),
      round(pooled_p_val, 6),
      round(pooled_p_val, 6),
      round(weighted_rmse, 3),
      unique(results_df$Levels),
      ref_rule,
      comp_rule
    )

    colnames(average_results) <- colnames(results_df)

    list_names[[i]] <- unique(results_df$Levels)

    results <- as.data.frame(rbind(results_df, average_results))
    results$fold <- c(seq(nrow(results) - 1), "Pooled")

    v_fold_marginal_w_pooled[[i]] <- results
  }

  names(v_fold_marginal_w_pooled) <- list_names


  return(v_fold_marginal_w_pooled)
}
