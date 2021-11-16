flatten_tmle3mixrule_results <- function(results){
  num_folds <- length(results$fold_mix_ates)
  num_marginals <- length(results$fold_marg_ates[[1]])
  marg_columns <- num_marginals*4

  marg_results_df <- as.data.frame(matrix(data = NA, nrow = num_folds, ncol = marg_columns))
  mix_results_df <- as.data.frame(matrix(data = NA, nrow = length(results$fold_mix_ates), ncol = 4))

  ## mixture data
  mix_ATEs <- unlist(results$fold_mix_ates)
  mix_SEs <- unlist(results$fold_mix_se)
  mix_CIs <- unlist(results$fold_mix_CIs)
  mix_CIs <- matrix(mix_CIs, nrow = num_folds, byrow = TRUE)
  mix_pvals <- unlist(results$fold_mix_pval)


  mix_results_df <- cbind(mix_ATEs, mix_SEs, mix_CIs, mix_pvals)
  if (is.na(mix_CIs) == FALSE) {
    colnames(mix_results_df) <- c("Mixture ATE", "Mixture SE", "Mixture Low 95% CI", "Mixture High 95% CI", "Mixture ATE p-value")
  }
  fold_mix_averages <- colMeans(mix_results_df)

  mix_results <- rbind(mix_results_df,fold_mix_averages)
  rownames(mix_results) <- c(paste("Fold", seq(num_folds), sep = " "), "Fold Average")

  mix_results <- as.data.frame(mix_results)
  ###### Marginals Data Constr ###########
  even_indexes <- seq(2,num_marginals*2,2)
  odd_indexes <- seq(1,num_marginals*2,2)

  ## ATEs
  marg_ATEs <- unlist(results$fold_marg_ates)
  marg_ATEs <- matrix(marg_ATEs, nrow = num_folds, byrow = TRUE)
  colnames(marg_ATEs) <- c(paste("Marg ATE for Mix Comp", seq(num_marginals), sep = " "))
  ## SEs
  marg_SEs <- unlist(results$fold_marg_se)
  marg_SEs <- matrix(marg_SEs, nrow = num_folds, byrow = TRUE)
  colnames(marg_SEs) <- c(paste("Marg SE for Mix Comp", seq(num_marginals), sep = " "))

  marg_CIs <- unlist(results$fold_marg_CI)
  marg_CIs <- matrix(marg_CIs, nrow = num_folds, byrow = TRUE)
  marg_CIs <- as.data.frame(marg_CIs)

  colnames(marg_CIs)[odd_indexes] <- c(paste("Marg CI Low for Mix Comp", seq(num_marginals), sep = " "))
  colnames(marg_CIs)[even_indexes] <- c(paste("Marg CI High for Mix Comp", seq(num_marginals), sep = " "))

  ## pvals
  marg_pvals <- unlist(results$fold_marg_pval)
  marg_pvals <- matrix(marg_pvals, nrow = num_folds, byrow = TRUE)
  colnames(marg_pvals) <- c(paste("Marg p-val for Mix Comp", seq(num_marginals), sep = " "))

  ## make additive data
  additive_marginal_ATEs <- rowSums(marg_ATEs)
  additive_low_CI <- rowSums(marg_CIs[,odd_indexes])
  additive_high_CI <- rowSums(marg_CIs[,even_indexes])
  additive_results_df <- cbind(additive_marginal_ATEs, additive_low_CI, additive_high_CI)
  fold_additive_averages <- colMeans(additive_results_df)
  additive_results <- rbind(additive_results_df, fold_additive_averages )
  rownames(additive_results) <- c(paste("Fold", seq(num_folds), sep = " "), "Fold Average")
  additive_results <- as.data.frame(additive_results)

  marg_results_df <- cbind(marg_ATEs, marg_SEs, marg_CIs, marg_pvals)

  fold_marg_averages <- colMeans(marg_results_df)
  marg_results <- rbind(marg_results_df,fold_marg_averages)
  rownames(marg_results) <- c(paste("Fold", seq(num_folds), sep = " "), "Fold Average")

  marg_results <- as.data.frame(marg_results)

  ## put together the additive and mixture results
  mix_additive <- cbind(mix_results, additive_results)
  if (is.na(mix_additive[1,1]) == FALSE) {
    mix_additive$Synergy_ATE <- mix_additive$`Mixture ATE` - mix_additive$additive_marginal_ATEs
  }


  # Compile results
  results_list = list("Mixture-Additive Results" = mix_additive, "Marginal Results" = marg_results)

  return(results_list)

}
