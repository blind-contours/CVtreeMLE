#' @title Calculate the ATE for each rule found for individual mixture components.
#' @description Aggregate marginal rules for each mixture variable found across the folds. For each rule
#' extract the relevant nuisance parameter data calculated in the folds. Given the validation data estimates across
#' the folds, for each tree do a TMLE update step to target the average treatment effect. Update the initial counterfactuals,
#' calculate the influence curve and using the influence curve calculate variance estimates and p-values.
#'
#' @param marginal_data List of dataframes of nuisance parameter data for each mixture
#' @param mix_comps Vector of characters indicating the mixture components
#' @param Y Vector indicating the Y
#' @param n_folds Number of folds used in cross-validation

#' @importFrom data.table rbindlist
#' @importFrom dplyr group_by bind_rows

#' @return Rules object. TODO: add more detail here.
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#'
#' @export

calc_marginal_ate <- function(marginal_data, mix_comps, Y, n_folds){

  marg_data <- list()

  x <- unlist(marginal_data,recursive=FALSE)
  x <- x[!sapply(x,is.null)]

  for (i in seq(mix_comps)) {
    # for(j in seq(marginal_data)){
    marg_data[[i]] <-
      dplyr::bind_rows(sapply(x, "[", i))
    }

  marginal_results <-
    as.data.frame(matrix(
      data = NA,
      nrow = length(mix_comps),
      ncol = 7
    ))

  colnames(marginal_results) <-
    c("Marginal ATE",
      "Standard Error",
      "Lower CI",
      "Upper CI",
      "P-value",
      "P-value Adj",
      "RMSE")

  rownames(marginal_results) <- mix_comps

  updated_marginal_data <- list()

  ## find which marginal mixture rule was null across the folds to update bonferonni p-value

  no_null_indices <- list()

  for (i in seq(length(marg_data))) {
    check <- marg_data[i]
    if (length(check[[1]]) != 0) {
      no_null_indices <- append(no_null_indices, i)
    } else{
      no_null_indices <- no_null_indices
    }
  }

  no_null_indices <- unlist(no_null_indices)

  for (i in seq(marg_data)) {
    if(i %in% no_null_indices == TRUE){
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
        RMSE <- sqrt(mean(sqrd_resids[,1]))

        marginal_results$`Marginal ATE`[i] <- ATE_results$ATE
        marginal_results$`Standard Error`[i] <- ATE_results$SE
        marginal_results$`Lower CI`[i] <- ATE_results$CI[1]
        marginal_results$`Upper CI`[i] <- ATE_results$CI[2]
        marginal_results$`P-value`[i] <- ATE_results$`p-value`
        marginal_results$`P-value Adj`[i] <- ATE_results$`adj p-value`
        marginal_results$RMSE[i] <- RMSE

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
