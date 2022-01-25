#' @title Calculate the ATE for each rule found for individual mixture components.
#' @description Aggregate marginal rules for each mixture variable found across the folds. For each rule
#' extract the relevant nuisance parameter data calculated in the folds. Given the validation data estimates across
#' the folds, for each tree do a TMLE update step to target the average treatment effect. Update the initial counterfactuals,
#' calculate the influence curve and using the influence curve calculate variance estimates and p-values.
#'
#' @param marginal_data List of dataframes of nuisance parameter data for each mixture
#' @param mix_comps Vector of characters indicating the mixture components
#' @param marginal_rules List of dataframes for rules found across the folds
#' @param Y Vector indicating the Y
#' @param n_folds Number of folds used in cross-validation

#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows

#' @return Rules object. TODO: add more detail here.
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export

calc_marginal_ate <- function(marginal_data, mix_comps, marginal_rules, Y, n_folds) {
  marg_data <- list()

  x <- unlist(marginal_data, recursive = FALSE)
  x <- x[!sapply(x, is.null)]

  y <- unlist(marginal_rules, recursive = FALSE)
  y <- y[!sapply(y, is.null)]

  fold_rules_groups <- marginal_group_split(y[[1]])

  # fold_rules_groups <- y[[1]] %>% dplyr::group_by(target_m)
  # fold_rules_groups <- dplyr::group_split(fold_rules_groups)

  comp_labels <- list()

  for (i in seq(fold_rules_groups)) {
    var_comps <- list()
    rules <- fold_rules_groups[[i]]
    reference <- rules[rules$quantile == 1, ]
    reference_quant <- reference$var_quant_group
    comparisons <- rules[rules$quantile > 1, ]
    for (j in seq(nrow(comparisons))) {
      comp_row <- comparisons[j, ]
      comp_label <- paste(comp_row$var_quant_group, reference_quant, sep = "-")
      var_comps[[j]] <- comp_label
    }

    comp_labels[[i]] <- var_comps
  }

  comp_labels <- unlist(comp_labels)


  for (i in seq(length(x[[1]]))) {
    # for(j in seq(marginal_data)){
    marg_data[[i]] <-
      dplyr::bind_rows(sapply(x, "[", i))
  }

  marginal_results <-
    as.data.frame(matrix(
      data = NA,
      nrow = length(x[[1]]),
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
      "RMSE"
    )

  target_quantile <- as.data.frame(marginal_rules[[1]])

  rownames(marginal_results) <- comp_labels

  updated_marginal_data <- list()

  ## find which marginal mixture rule was null across the folds to update bonferonni p-value

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
    if (i %in% no_null_indices == TRUE) {
      marg_mix <- marg_data[[i]]

      if (length(unique(marg_mix$folds)) == n_folds) {
        flux_results <- fit_least_fav_submodel(H.AW = marg_mix$H.AW, data = marg_mix, QbarAW = marg_mix$QbarAW, Qbar1W = marg_mix$Qbar1W, Qbar0W = marg_mix$Qbar0W)

        ## back-scale Y
        QbarAW.star <- scale_to_original(scaled_vals = flux_results$QbarAW.star, max_orig = max(marg_mix[Y]), min_orig = min(marg_mix[Y]))
        Qbar0W.star <- scale_to_original(scaled_vals = flux_results$Qbar0W.star, max_orig = max(marg_mix[Y]), min_orig = min(marg_mix[Y]))
        Qbar1W.star <- scale_to_original(scaled_vals = flux_results$Qbar1W.star, max_orig = max(marg_mix[Y]), min_orig = min(marg_mix[Y]))

        marg_mix$QbarAW.star <- QbarAW.star
        marg_mix$Qbar0W.star <- Qbar0W.star
        marg_mix$Qbar1W.star <- Qbar1W.star

        ATE_results <- calc_ATE_estimates(data = marg_mix, ATE_var = "marg.ATE", outcome = Y, p_adjust_n = length(no_null_indices))

        sqrd_resids <- (marg_mix$QbarAW.star - marg_mix[Y])^2
        RMSE <- sqrt(mean(sqrd_resids[, 1]))

        marginal_results$`Marginal ATE`[i] <- round(ATE_results$ATE, 3)
        marginal_results$`Standard Error`[i] <- round(ATE_results$SE, 3)
        marginal_results$`Lower CI`[i] <- round(ATE_results$CI[1], 3)
        marginal_results$`Upper CI`[i] <- round(ATE_results$CI[2], 3)
        marginal_results$`P-value`[i] <- round(ATE_results$`p-value`, 6)
        marginal_results$`P-value Adj`[i] <- round(ATE_results$`adj p-value`, 6)
        marginal_results$RMSE[i] <- round(RMSE, 3)

        updated_marginal_data[[i]] <- ATE_results$data
      } else {
        print(
          paste(
            "A marginal rule for mixture component",
            i,
            "was not found for all folds, dropping this mixture from results due to inconsistency"
          )
        )
        marginal_results$`Marginal ATE`[i] <- NA
        marginal_results$`Standard Error`[i] <- NA
        marginal_results$`Lower CI`[i] <- NA
        marginal_results$`Upper CI`[i] <- NA
        marginal_results$`P-value`[i] <- NA
        marginal_results$`P-value Adj`[i] <- NA
        marginal_results$RMSE[i] <- NA

        updated_marginal_data[[i]] <- NA
      }
    } else {
      print(
        paste(
          "A marginal rule for mixture component",
          i,
          "was not found for one folds, dropping this mixture from results due to inconsistency"
        )
      )
      marginal_results$`Marginal ATE`[i] <- NA
      marginal_results$`Standard Error`[i] <- NA
      marginal_results$`Lower CI`[i] <- NA
      marginal_results$`Upper CI`[i] <- NA
      marginal_results$`P-value`[i] <- NA
      marginal_results$`P-value Adj`[i] <- NA
      marginal_results$RMSE[i] <- NA

      updated_marginal_data[[i]] <- NA
    }
  }
  return(list(marginal_results = marginal_results, data = updated_marginal_data))
}
