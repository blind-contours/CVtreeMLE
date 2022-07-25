#' @title Compute v-fold specific estimates and do a meta-analysis type
#' pooled average ATE for mixture rules.
#' @description For each fold, estimates the ATE for a fold specific
#' mixture rule. Also estimates a meta-analysis type average ATE and pooled
#' variance. Creates a union rule that covers all the folds in the rules.
#'
#' @param v_fold_mixture_results List of dataframes for v-fold specific
#' estimates of the fold-specific rule results
#' @param mix_comps Vector of characters indicating the mixture components
#' @param n_folds Number of folds used in cross-validation
#' @param data Input dataframe on which to evaluate the rules

#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows

#' @return Rules object. TODO: add more detail here.
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis
#' @importFrom stats qnorm qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export


meta_mix_results <- function(v_fold_mixture_results,
                             mix_comps,
                             n_folds,
                             data) {

  v_fold_mixture_group <- v_fold_mixture_group_split(v_fold_mixture_results)

  v_fold_mixture_w_pooled <- list()
  intxn_names_list <- list()

  for (i in seq(v_fold_mixture_group)) {
    results_df <- v_fold_mixture_group[[i]]

    weighted_mean <- sum(results_df$ate *
                           (1 / results_df$se^2)) / sum((1 / results_df$se^2))

    weighted_rmse <- sum(results_df$rmse *
                           (1 / results_df$se^2)) / sum((1 / results_df$se^2))

    pooled_se <- sqrt(1 / (1 / sum(results_df$se^2)))

    pooled_p_val <- round(2 *
                            stats::pnorm(abs(weighted_mean / pooled_se),
                                         lower.tail = FALSE), 5)

    pooled_ci <- c(
      round(weighted_mean + stats::qnorm(0.05 / 2, lower.tail = TRUE) *
              pooled_se, 4),
      round(weighted_mean + stats::qnorm(0.05 / 2, lower.tail = FALSE) *
              pooled_se, 4)
    )

    vars <- mix_comps[mix_comps %in% strsplit(results_df$mix_rule[1],
                                              split = " ")[[1]]]

    intxn_rule <- paste("(", paste(results_df$mix_rule, collapse = ")|("), ")")

    intxn_data <- data %>%
      dplyr::mutate("intxn_rule" := ifelse(eval(parse(text = intxn_rule)),
                                           1, 0))

    new_rule <- list()

    for (var in vars) {
      var_min <-
        intxn_data %>%
        group_by(intxn_rule) %>%
        summarise(min = min(!!(as.name(var))))
      var_min <- subset(var_min, intxn_rule == 1, select = min)
      var_max <-
        intxn_data %>%
        group_by(intxn_rule) %>%
        summarise(max = max(!!(as.name(var))))
      var_max <- subset(var_max, intxn_rule == 1, select = max)

      augmented_rule <- paste(var, ">", round(var_min, 3), "&", var, "<",
                              round(var_max, 3))

      new_rule <- append(new_rule, augmented_rule)
    }

    interaction_rule <- paste(unlist(new_rule), collapse = " & ")

    average_results <- cbind.data.frame(
      round(weighted_mean, 3),
      round(pooled_se, 3),
      pooled_ci[1], pooled_ci[2],
      pooled_p_val, pooled_p_val, round(weighted_rmse, 3),
      interaction_rule,
      "Pooled",
      unique(results_df$variables)
    )

    colnames(average_results) <- colnames(results_df)

    results <- as.data.frame(rbind.data.frame(results_df, average_results))
    intxn_names_list[[i]] <- unique(results$variables)

    v_fold_mixture_w_pooled[[i]] <- results
  }

  names(v_fold_mixture_w_pooled) <- intxn_names_list


  return(v_fold_mixture_w_pooled)
}
