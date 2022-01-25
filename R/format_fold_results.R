
format_fold_results <- function(marginal_rules,
                                mix_rules,
                                mix_combo_data,
                                marg_combo_data) {


  # TODO: turn this into a function instead of repeat coding

  ##########################################################################################################
  ########################################## MARGINALS  ####################################################
  ##########################################################################################################
  marginal_rules <- unlist(marginal_rules, recursive = FALSE, use.names = FALSE)
  marginal_rules <- marginal_rules[!sapply(marginal_rules, is.null)]
  marginal_rules <- do.call(rbind, marginal_rules)



  ##########################################################################################################
  ########################################## MIXTURES  #####################################################
  ##########################################################################################################
  mix_rules <- unlist(mix_rules, recursive = FALSE, use.names = FALSE)
  mix_rules <- mix_rules[!sapply(mix_rules, is.null)]
  mix_rules <-
    data.table::rbindlist(mix_rules)

  ############################################################################################################
  ########################################## COMBO DATA  #####################################################
  ############################################################################################################
  mix_combo_data <- unlist(mix_combo_data, recursive = FALSE, use.names = FALSE)
  mix_combo_data <- mix_combo_data[!sapply(mix_combo_data, is.null)]

  if (all(is.na(mix_combo_data)) == FALSE) {
    mix_combo_data <- data.table::rbindlist(mix_combo_data, idcol = "fold", fill = TRUE)
  }



  marg_combo_data <- unlist(marg_combo_data, recursive = FALSE, use.names = FALSE)
  marg_combo_data <- marg_combo_data[!sapply(marg_combo_data, is.null)]
  marg_combo_data <- data.table::rbindlist(marg_combo_data, idcol = "fold", fill = TRUE)

  return(list(
    "marginal_rules" = marginal_rules,
    "mix_rules" = mix_rules,
    "mix_combo_data" = mix_combo_data,
    "marg_combo_data" = marg_combo_data
  ))
}
