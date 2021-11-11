

CVtreeMLE <- function(W,
                      A,
                      Y,
                      data,
                      SL.library,
                      n_folds,
                      binary,
                      H.AW_trunc_lvl,
                      rule_boot_n,
                      n_stable_tree,
                      minbucket,
                      verbose,
                      parallel) {
  if (binary == TRUE) {
    family <- "binomial"
    ## create the CV folds
    data[, "y_scaled"] <- data[[Y]]
    data$folds <- create_cv_folds(n_folds, data$y_scaled)
  } else {
    family <- "gaussian"
    data[, "y_scaled"] <-
      ((data[[Y]] - min(data[[Y]])) / (max(data[[Y]]) - min(data[[Y]])))
    ## create the CV folds
    data$folds <- create_cv_folds(n_folds, data$y_scaled)
  }

  if (parallel == TRUE) {
    library(parallel, quietly = TRUE)
    library(doParallel, quietly = TRUE)
    library(foreach, quietly = TRUE)
    num_cores <- detectCores()
    doParallel::registerDoParallel(num_cores)
    print(paste("Parallelizing over", num_cores, "Cores"))
  }

  data <- as.data.frame(data)
  n <- dim(data)[1]

  fold_results_mix_rules <- fold_results_marg_rules <- fold_results_marg_directions <- list()
  mix_fold_data <- mix_fold_directions <- list()

  fold_results_marginal_data <- fold_results_marginal_additive_data <- list()

  fold_results_marginal_combo_data <- fold_results_marginal_Sls <- list()

  fold_results_mix_combo_data <- fold_results_mix_Sls <- list()

  results <- foreach(fold_k = unique(data$folds), .combine = "c", .multicombine = TRUE) %dopar% {
    print(paste("Starting Fold", fold_k))
    ## get training and validation samples
    At <- data[data$folds != fold_k, ]
    Av <- data[data$folds == fold_k, ]


    ## subset At to only Ws for SL Y~W
    X <- subset(At,
      select = c(W)
    )

    if (verbose) cat("Fitting Y|W using SL \n")

    ## train Y~W learner
    QbarWSL <- suppressWarnings(SuperLearner(
      Y = At$y_scaled,
      X = X,
      SL.library = SL.library,
      family = family,
      verbose = FALSE
    ))

    ## calculate remaining variance unexplained by W in residuals
    QbarW <- predict(QbarWSL, newdata = X)$pred
    QbarW_res <-
      QbarW - At$y_scaled ## include variance detection as an option - most people would want residual

    total.res <- QbarW_res
    At <- cbind(At, total.res)

    # res.diff <- 1

    if (verbose) cat("Fitting (Y_hat-Y)|M using Ensemble Decision Trees \n")

    formula <-
      as.formula(paste("total.res~", paste(A, collapse = "+")))

    ##############################################################################################################################################
    ########### Bootstrap (with replacement) the training data and fit pre rule_boot_n times #####################################################
    ########### Group by rules with same variables that are used at all bootstraps and use the rule with highest coefficient #####################
    ##############################################################################################################################################

    rules <-
      boot_ensemble_trees(At, formula, rule_boot_n, n_stable_tree, A)

    if (dim(rules)[1] == 0) {
      no_rules <- TRUE
    } else if (rules$rule == "(Intercept)") {
      no_rules <- TRUE
    } else {
      no_rules <- FALSE
    }


    if (no_rules == TRUE) {
      print(
        "Ensemble decision trees found no mixture rules that explain (Y_hat-Y)|W"
      )

      no_rule_df <- data.frame(
        rule = NA,
        coefficient = NA,
        description = NA,
        test = NA,
        boot_num = NA
      )


      fold_results_mix_rules[[fold_k]] <- no_rule_df
    } else {
      fold_results_mix_rules[[fold_k]] <- rules
    }

    #########################################################################################################
    ################## BEGIN FITTING DECISION TREES TO EACH MARGINAL COMPONENT OF THE MIXTURE ###############
    #########################################################################################################

    marg_decisions <- fit_marginal_decision_trees(
      A,
      At,
      W,
      SL.library,
      family,
      minbucket
    )

    fold_results_marg_rules[[fold_k]] <- marg_decisions

    #######################################################################################################
    ####### =################### NUISANCE PARAMETER FITTING ON TRAINING ####################################
    #######################################################################################################

    if (verbose) {
      cat(
        "Training learners for nuisance parameters using training data and rules identified in training"
      )
    }

    mix_nuisance_params <- mixture_nuisance_parameters(
      At, Av, W, no_rules,
      SL.library, family, rules, fold_k,
      mix_fold_data, fold_results_mix_rules, H.AW_trunc_lvl
    )

    mix_interaction_data <- mix_nuisance_params$data
    # mix_directions <- mix_nuisance_params$directions

    mix_fold_data[[fold_k]] <- mix_interaction_data
    # fold_results_mix_rules[[fold_k]]$direction <- suppressWarnings(unlist(mix_directions))

    ###################################################################################################
    ################ MARGINAL ESTIMATES FROM RULE IN TRAINING SET #####################################
    ###################################################################################################

    marg_nuisance_params <- marginal_nuisance_parameters(
      At, Av, W, SL.library, family,
      A, marg_decisions, fold_k, fold_results_marginal_data,
      fold_results_marg_directions, H.AW_trunc_lvl
    )

    fold_results_marginal_data[[fold_k]] <- marg_nuisance_params$data
    fold_results_marg_directions[[fold_k]] <- marg_nuisance_params$directions

    marginal_data <- marg_nuisance_params$data
    marg_directions <- marg_nuisance_params$directions


    ################################################################################################################
    ################## CALCULATE THE EXPECTED OUTCOME FOR EACH ADDED MARGINAL MIXTURE RULE EXPOSURE ################
    ################################################################################################################

    At_c <- At
    Av_c <- Av

    At_c <-
      evaluate_marginal_rules(
        data = At_c,
        marg_decisions = marg_decisions,
        mix_comps = mix_comps,
        marg_directions = marg_directions
      )
    Av_c <-
      evaluate_marginal_rules(
        data = Av_c,
        marg_decisions = marg_decisions,
        mix_comps = mix_comps,
        marg_directions = marg_directions
      )

    marg_rule_train <- At_c$marginals
    marg_rule_valid <- Av_c$marginals

    At_c <- At_c$data
    Av_c <- Av_c$data


    H.AW.list <- list()


    for (i in rev(seq(table(At_c$sum_marg_hits)))) {
      level <- table(At_c$sum_marg_hits)[i]
      if (level[[1]] < 10 & as.numeric(names(level)) != min(At_c$sum_marg_hits)) {
        At_c$sum_marg_hits[At_c$sum_marg_hits == as.numeric(names(level))] <- as.numeric(names(level)) - 1
        Av_c$sum_marg_hits[Av_c$sum_marg_hits == as.numeric(names(level))] <- as.numeric(names(level)) - 1
      } else if (level[[1]] < 10 & as.numeric(names(level)) == min(At_c$sum_marg_hits)) {
        At_c$sum_marg_hits[At_c$sum_marg_hits == as.numeric(names(level))] <- 0
        Av_c$sum_marg_hits[Av_c$sum_marg_hits == as.numeric(names(level))] <- 0
      }
    }


    for (i in sort(unique(At_c$sum_marg_hits[At_c$sum_marg_hits != 0]))) {
      target.lvl <- sort(unique(At_c$sum_marg_hits))[sort(unique(At_c$sum_marg_hits)) == i]
      X_train_covars <- At_c %>% dplyr::select(W)
      X_valid_covars <- Av_c %>% dplyr::select(W)

      At_c$binarized_cat <-
        as.numeric(At_c$sum_marg_hits == target.lvl)

      Av_c$binarized_cat <-
        as.numeric(Av_c$sum_marg_hits == target.lvl)


      gHatSL <- suppressWarnings(SuperLearner(
        Y = At_c$binarized_cat,
        X = X_train_covars,
        SL.library = SL.library,
        family = "binomial",
        verbose = FALSE
      ))

      gHat1W <- predict(gHatSL, newdata = X_valid_covars)$pred
      gHat0W <- 1 - gHat1W

      gHatAW <- rep(NA, dim(Av_c)[1])

      gHatAW[Av_c$binarized_cat == 1] <-
        gHat1W[Av_c$binarized_cat == 1]

      gHatAW[Av_c$binarized_cat == 0] <-
        gHat0W[Av_c$binarized_cat == 0]

      H.AW <-
        as.numeric(Av_c$binarized_cat == 1) / gHat1W #- as.numeric(Av_c$binarized_cat==0)/gHat0W

      H.AW <-
        ifelse(H.AW > H.AW_trunc_lvl, H.AW_trunc_lvl, H.AW)

      H.AW <-
        ifelse(H.AW < -H.AW_trunc_lvl, -H.AW_trunc_lvl, H.AW)

      H.AW.list[[i]] <- H.AW
    }

    HA.W.by.lvl <- do.call(cbind, H.AW.list)
    HA.W.lvl.sums <- rowSums(HA.W.by.lvl, na.rm = TRUE)


    X_train_mix <-
      subset(At_c, select = c("sum_marg_hits", W))
    X_valid_mix <-
      subset(Av_c, select = c("sum_marg_hits", W))

    print("Fitting SL of Y given additive marginal rule fits")

    ## QbarAW
    QbarAWSL_m <- suppressWarnings(SuperLearner(
      Y = At_c$y_scaled,
      X = X_train_mix,
      SL.library = SL.library,
      family = family,
      verbose = FALSE
    ))

    QbarAW <- predict(QbarAWSL_m, newdata = X_valid_mix)$pred

    QbarAW[QbarAW >= 1.0] <- 0.999
    QbarAW[QbarAW <= 0] <- 0.001

    Av_c$QbarAW_additive <- QbarAW
    Av_c$HAW_additive <- HA.W.lvl.sums

    fold_results_marginal_additive_data[[fold_k]] <- Av_c

    #########################################################################################################################
    ################## CALCULATE THE EXPECTED OUTCOME FOR EACH COMBINATION OF MARGINAL MIXTURE RULE EXPOSURE ################
    #########################################################################################################################
    At_mc <- At
    Av_mc <- Av

    marg_rule_df <-
      marg_rule_train[, colSums(is.na(marg_rule_train)) < nrow(marg_rule_train)]

    At_marg_comb <-
      cbind(marg_rule_df, subset(At_mc, select = W))

    At_marg_comb <- At_marg_comb[, colSums(is.na(At_marg_comb)) < nrow(At_marg_comb)]

    Av_marg_comb <-
      cbind(marg_rule_valid, subset(Av_mc, select = W))

    Av_marg_comb <- Av_marg_comb[, colSums(is.na(Av_marg_comb)) < nrow(Av_marg_comb)]

    QbarAWSL_m <- SuperLearner(
      Y = At_mc$y_scaled,
      X = At_marg_comb,
      SL.library = SL.library,
      family = family,
      verbose = FALSE
    )

    QbarAW <- predict(QbarAWSL_m, newdata = Av_marg_comb)$pred

    QbarAW[QbarAW >= 1.0] <- 0.999
    QbarAW[QbarAW <= 0] <- 0.001

    Av_marg_comb$QbarAW_combo <- QbarAW
    Av_marg_comb$y_scaled <- Av_mc$y_scaled
    Av_marg_comb$raw_outcome <- Av[, Y]

    fold_results_marginal_combo_data[[fold_k]] <- Av_marg_comb
    fold_results_marginal_Sls[[fold_k]] <- QbarAWSL_m


    #########################################################################################################################
    ################## CALCULATE THE EXPECTED OUTCOME FOR EACH COMBINATION OF THE MIXTURE INTERACTION RULES  ################
    #########################################################################################################################

    At_c <- At
    Av_c <- Av

    mix_interaction_train <- list()
    mix_interaction_valid <- list()

    rules <- as.data.frame(rules)

    if (dim(rules)[1] > 0 & no_rules == FALSE) {
      for (i in seq(dim(rules)[1])) {
        target_interaction_result <- rules[i, ]$description

        rule_name <- paste("mix_interaction", i, sep = "_")
        target_interaction_result <-
          gsub("-", " -", target_interaction_result)

        vec_train <- At_c %>%
          mutate(!!(rule_name) := ifelse(eval(
            parse(text = target_interaction_result)
          ), 1, 0)) %>%
          dplyr::select(!!rule_name)

        vec_valid <- Av_c %>%
          mutate(!!(rule_name) := ifelse(eval(
            parse(text = target_interaction_result)
          ), 1, 0)) %>%
          dplyr::select(!!rule_name)

        # if (rules[i,]$direction == "negative") {
        #   vec_train <- 1 - vec_train
        #   vec_valid <- 1 - vec_valid
        # }
        mix_interaction_train[[i]] <- vec_train
        mix_interaction_valid[[i]] <- vec_valid
      }

      mix_interactions_rules_train <-
        do.call(cbind, mix_interaction_train)
      At_mix_comb <-
        cbind(
          mix_interactions_rules_train,
          subset(At_c, select = W)
        )

      mix_interactions_rules_valid <-
        do.call(cbind, mix_interaction_valid)

      Av_mix_comb <-
        cbind(
          mix_interactions_rules_valid,
          subset(Av_c, select = W)
        )

      QbarAWSL_m <- SuperLearner(
        Y = At_c$y_scaled,
        X = At_mix_comb,
        SL.library = SL.library,
        family = family,
        verbose = FALSE
      )

      QbarAW <- predict(QbarAWSL_m, newdata = Av_mix_comb)$pred

      QbarAW[QbarAW >= 1.0] <- 0.999
      QbarAW[QbarAW <= 0] <- 0.001

      Av_mix_comb$QbarAW_combo <- QbarAW
      Av_mix_comb$y_scaled <- Av_c$y_scaled
      Av_mix_comb$raw_outcome <- Av[, Y]

      fold_results_mix_combo_data[[fold_k]] <- Av_mix_comb
      fold_results_mix_Sls[[fold_k]] <- QbarAWSL_m
    } else {
      fold_results_mix_combo_data[[fold_k]] <- NA
      fold_results_mix_Sls[[fold_k]] <- NA
    }

    # results_list <- list(
    #   setNames(fold_results_mix_rules, paste("Mix Rules", fold_k)),
    #   setNames(mix_fold_data, paste("Mix Data", fold_k)),
    #   setNames(mix_fold_directions, paste("Mix Directions", fold_k)),
    #   setNames(fold_results_marg_rules, paste("Marginal Rules", fold_k)),
    #   setNames(fold_results_marg_directions, paste("Marginal Directions", fold_k)),
    #   setNames(fold_results_marginal_data, paste("Marginal Data", fold_k)),
    #   setNames(fold_results_marginal_additive_data, paste("Additive Data", fold_k)),
    #   setNames(fold_results_marginal_combo_data, paste("Marginal Combination Data", fold_k)),
    #   setNames(fold_results_marginal_Sls, paste("Marginal SLs", fold_k)),
    #   setNames(fold_results_mix_combo_data, paste("Mixture Combination Data", fold_k)),
    #   setNames(fold_results_mix_Sls, paste("Mixture SLs", fold_k))
    # )

    results_list <- list(
      fold_results_mix_rules,
      mix_fold_data,
      mix_fold_directions,
      fold_results_marg_rules,
      fold_results_marg_directions,
      fold_results_marginal_data,
      fold_results_marginal_additive_data,
      fold_results_marginal_combo_data,
      fold_results_marginal_Sls,
      fold_results_mix_combo_data,
      fold_results_mix_Sls
    )

    names(results_list) <- c(
      paste("mix rules", fold_k),
      paste("mix data", fold_k),
      paste("mix directions", fold_k),
      paste("marg rules", fold_k),
      paste("marg directions", fold_k),
      paste("marg data", fold_k),
      paste("additive data", fold_k),
      paste("marg combo data", fold_k),
      paste("marginal SLs", fold_k),
      paste("mix combo data", fold_k),
      paste("mix SLs", fold_k)
    )




    results_list
  }

  ########################################################################################################
  ##################################### AGGREGATE DATA ACROSS FOLDS ######################################
  ########################################################################################################

  marginal_data <- results[grepl("marg data", names(results))]
  marginal_rules <- results[grepl("marg rules", names(results))]
  marginal_directions <- results[grepl("marg directions", names(results))]

  mix_data <- results[grepl("mix data", names(results))]
  mix_rules <- results[grepl("mix rules", names(results))]
  mix_directions <- results[grepl("mix directions", names(results))]

  additive_data <- results[grepl("additive data", names(results))]
  marg_combo_data <- results[grepl("marg combo data", names(results))]

  mix_combo_data <- results[grepl("mix combo data", names(results))]

  marg_Sls <- results[grepl("marginal SLs", names(results))]
  mix_Sls <- results[grepl("mix SLs", names(results))]



  ########################################################################################################
  ######################### POOLED TMLE ATE FOR MIXTURES FOUND ACROSS FOLDS.##############################
  ########################################################################################################

  mixture_results <- calc_mixtures_ate(
    input_mix_rules = mix_rules,
    input_mix_data = mix_data,
    outcome = outcome
  )

  group_list <- mixture_results$group_list
  mix_RMSE <- mixture_results$mix_RMSE
  mixture_results <- mixture_results$results

  ########################################################################################################
  ###################### POOLED TMLE ATE FOR ADDITIVE MARGINAL FOUND ACROSS FOLDS ########################
  ########################################################################################################

  additive_results <- calc_additive_ate(additive_data, outcome)

  additive_RMSE <- additive_results$RMSE
  additive_RMSE_star <- additive_results$RMSE_star
  additive_MSM <- additive_results$MSM

  ########################################################################################################
  ###################### POOLED TMLE ATE FOR MARGINAL RULES FOUND ACROSS FOLDS ###########################
  ########################################################################################################

  marginal_results <- calc_marginal_ate(marginal_data, A, outcome)

  updated_marginal_data <- marginal_results$data
  marginal_results <- marginal_results$marginal_results


  #################################################
  ######### Additive Results Config. ##############
  #################################################


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
  varHat.IC <- var(additive.IC, na.rm = TRUE) / n
  se <- sqrt(varHat.IC)

  alpha <- 0.05

  Theta <- ave.additive.Marg
  # obtain 95% two-sided confidence intervals:
  CI <- c(
    Theta + qnorm(alpha / 2, lower.tail = T) * se,
    Theta + qnorm(alpha / 2, lower.tail = F) * se
  )

  # p-value
  p.value <- round(2 * pnorm(abs(Theta / se), lower.tail = F), 4)

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
  additive_results$`Exp. Additive ATE` <- Theta
  additive_results$`Standard Error` <- se
  additive_results$`Lower CI` <- CI[1]
  additive_results$`Upper CI` <- CI[2]
  additive_results$`P-value` <- p.value

  ########################################################
  ######## Calculate Mixture Rules Over Folds ############
  ########################################################

  mixture_results <- find_common_mixture_rules(group_list,
    data = data,
    mix_comps = mix_comps,
    mixture_results = mixture_results,
    n_folds = n_folds
  )



  marginal_results <- find_common_marginal_rules(
    fold_rules = marginal_rules,
    data = data,
    mix_comps = mix_comps,
    marginal_results = marginal_results,
    fold_directions = marginal_directions,
    n_folds = n_folds
  )

  ##########################################################################################################
  ################################ FORMAT OUTPUTS FOR OUTPUT  ##############################################
  ##########################################################################################################

  # TODO: turn this into a function instead of repeat coding

  ##########################################################################################################
  ########################################## MARGINALS  ####################################################
  ##########################################################################################################
  marginal_rules <- unlist(marginal_rules, recursive = FALSE, use.names = FALSE)
  marginal_rules <- marginal_rules[!sapply(marginal_rules, is.null)]
  marginal_rules <- do.call(rbind, marginal_rules)


  marginal_directions <- unlist(marginal_directions, recursive = FALSE, use.names = FALSE)
  marginal_directions <- marginal_directions[!sapply(marginal_directions, is.null)]
  marginal_directions <- do.call(rbind, marginal_directions)

  marginal_rules_directions <- rbind(marginal_rules, marginal_directions)


  ##########################################################################################################
  ########################################## MIXTURES  #####################################################
  ##########################################################################################################
  mix_rules <- unlist(mix_rules, recursive = FALSE, use.names = FALSE)
  mix_rules <- mix_rules[!sapply(mix_rules, is.null)]
  mix_rules <-
    rbindlist(mix_rules, idcol = "fold", fill = TRUE)

  ############################################################################################################
  ########################################## COMBO DATA  #####################################################
  ############################################################################################################
  mix_combo_data <- unlist(mix_combo_data, recursive = FALSE, use.names = FALSE)
  mix_combo_data <- mix_combo_data[!sapply(mix_combo_data, is.null)]

  marg_combo_data <- unlist(marg_combo_data, recursive = FALSE, use.names = FALSE)
  marg_combo_data <- marg_combo_data[!sapply(marg_combo_data, is.null)]


  results_list <- list(
    "Marginal Results" = marginal_results,
    "Mixture Results" = mixture_results,
    "Additive Results" = additive_results,
    "Marginal Rules" = marginal_rules_directions,
    "Mixture Rules" = mix_rules,
    "Additive Model RMSE no star" = additive_RMSE,
    "Additive Model RMSE star" = additive_RMSE_star,
    "Additive MSM" = additive_MSM,
    "Mixture Model RMSE" = mix_RMSE,
    "Fold Results Mix Comb Data" = mix_combo_data,
    "Fold Results Marg Comb Data" = marg_combo_data,
    "Fold Super Learners Mix Comb" = mix_Sls,
    "Fold Super Learners Marg Comb" = marg_Sls
  )


  return(results_list)
}
