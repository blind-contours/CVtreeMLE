
create_add_marg_ATE_table <- function(updated_marginal_data, A, Y, n, n_folds) {
  additive_results <-
    as.data.frame(matrix(
      data = NA,
      nrow = 1,
      ncol = 5
    ))
  colnames(additive_results) <-
    c(
      "Exp. Additive ATE",
      "Standard Error",
      "Lower CI",
      "Upper CI",
      "P-value"
    )

  if (all(is.na(updated_marginal_data)) == FALSE) {

    ## identify the mixture components which had at least one fold where no rule was identified to associate with y - is therefore NA or there is no association with y
    dropped_mixed_marginals <- A[is.na(updated_marginal_data)]
    updated_marginal_data_RM_na <-
      updated_marginal_data[!is.na(updated_marginal_data)]

    ## for those mixtures that are not NA, get the marginal ATEs for each rule and take the rowSums for the additive ATE calculation
    ind.ATEs <- purrr::map(updated_marginal_data_RM_na, "marg.ATE")
    marginal.ATEs <- do.call(cbind, ind.ATEs)
    marginal.ATEs <- rowSums(marginal.ATEs, na.rm = TRUE)

    ## same thing but now get the IC for each ATE given each rule, bind them together (across the rules) and take sum across the rows to get the sum IC

    ind.ICs <- purrr::map(updated_marginal_data_RM_na, "IC")
    marginal.ICs <- do.call(cbind, ind.ICs)
    additive.IC <- rowSums(marginal.ICs, na.rm = TRUE)

    ave.additive.Marg <- mean(marginal.ATEs)
    varHat.IC <- stats::var(additive.IC, na.rm = TRUE) / n
    se <- sqrt(varHat.IC)

    alpha <- 0.05

    Theta <- ave.additive.Marg
    # obtain 95% two-sided confidence intervals:
    CI <- c(
      Theta + stats::qnorm(alpha / 2, lower.tail = T) * se,
      Theta + stats::qnorm(alpha / 2, lower.tail = F) * se
    )

    # p-value
    p.value <- round(2 * stats::pnorm(abs(Theta / se), lower.tail = F), 4)

    additive_results$`Exp. Additive ATE` <- Theta
    additive_results$`Standard Error` <- se
    additive_results$`Lower CI` <- CI[1]
    additive_results$`Upper CI` <- CI[2]
    additive_results$`P-value` <- p.value
  } else {
    additive_results$`Exp. Additive ATE` <- NA
    additive_results$`Standard Error` <- NA
    additive_results$`Lower CI` <- NA
    additive_results$`Upper CI` <- NA
    additive_results$`P-value` <- NA
  }

  return(additive_results)
}
