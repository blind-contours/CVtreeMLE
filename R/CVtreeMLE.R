
#' @title Estimate the joint impact of a mixed exposure using cross-validated targeted learning and decision trees
#'
#' @param data Data frame of (W,A,Y) variables of interest.
#' @param W A vector of characters indicating which variables in the data to use as covariates.
#' @param A A vector of characters indicating which variables in the data to use as exposures.
#' @param Y A character indicating which variable in the data to use as the outcome.
#' @param back_iter_SL Stack of estimators used in the SL during the iterative backfitting for `Y|W`, this should be an SL3 object
#' @param SL.library Library of algorithms used by Super Learner for estimating both the propensity of exposure and the outcome.
#' @param n_folds Number of cross-validation folds.
#' @param family Family ('binomial' or 'gaussian').
#' @param H.AW_trunc_lvl Truncation level for the clever covariate.
#' @param verbose If true, creates a sink file where detailed information and diagnostics of the run process are given.
#' @param parallel Use parallel processing if a backend is registered; enabled by default.
#'
#' @return Results object. TODO: add more detail here.
#'
#' @importFrom foreach %dopar% %do%
#' @importFrom magrittr %>%
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom SuperLearner All
#' @importFrom rlang :=
#' @importFrom dplyr filter
#' @importFrom data.table rbindlist
#'
#' @section Authors:
#' David McCoy, University of California, Berkeley
#'
#' @section References:
#' Benjamini, Y., & Hochberg, Y. (1995). \emph{Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing}. Journal of the
#' royal statistical society. Series B (Methodological), 289-300.
#'
#' Gruber, S., & van der Laan, M. J. (2012). \emph{tmle: An R Package for
#' Targeted Maximum Likelihood Estimation}. Journal of Statistical Software,
#' 51(i13).
#'
#' Hubbard, A. E., Kherad-Pajouh, S., & van der Laan, M. J. (2016).
#' \emph{Statistical Inference for Data Adaptive Target Parameters}. The
#' international journal of biostatistics, 12(1), 3-19.
#'
#' Hubbard, A., Munoz, I. D., Decker, A., Holcomb, J. B., Schreiber, M. A.,
#' Bulger, E. M., ... & Rahbar, M. H. (2013). \emph{Time-Dependent Prediction
#' and Evaluation of Variable Importance Using SuperLearning in High Dimensional
#' Clinical Data}. The journal of trauma and acute care surgery, 75(1 0 1), S53.
#'
#' Hubbard, A. E., & van der Laan, M. J. (2016). \emph{Mining with inference:
#' data-adaptive target parameters (pp. 439-452)}. In P. Buhlmann et al. (Ed.),
#' \emph{Handbook of Big Data}. CRC Press, Taylor & Francis Group, LLC: Boca
#' Raton, FL.
#'
#' van der Laan, M. J. (2006). \emph{Statistical inference for variable
#' importance}. The International Journal of Biostatistics, 2(1).
#'
#' van der Laan, M. J., & Pollard, K. S. (2003). \emph{A new algorithm for
#' hybrid hierarchical clustering with visualization and the bootstrap}. Journal
#' of Statistical Planning and Inference, 117(2), 275-303.
#'
#' van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). \emph{Super
#' learner}. Statistical applications in genetics and molecular biology, 6(1).
#'
#' van der Laan, M. J., & Rose, S. (2011). \emph{Targeted learning: causal
#' inference for observational and experimental data}. Springer Science &
#' Business Media.
#'
#' @examples
#' @export

CVtreeMLE <- function(W,
                      A,
                      Y,
                      data,
                      back_iter_SL,
                      SL.library,
                      n_folds,
                      family,
                      H.AW_trunc_lvl,
                      parallel,
                      verbose) {

  ######################
  # Argument checks.

  # Confirm that data has at least two columns.
  if (length(A) < 2L) {
    stop("Vector of mixture variable names must have at least two variable characters")
  }

  if (any(sapply(data[,A], is.factor))) {
    print("Factor variable detected in exposures, converting to numeric")
    data[,A] <- sapply(data[,A], as.numeric)
  }

  if (any(sapply(data[,W], is.factor))) {
    print("Factor variable detected in the covariate sapce, converting to numeric to avoid issues in different factors between training and validation folds")
    data[,W] <- sapply(data[,W], as.numeric)
  }

  # Ensure that Y is numeric; e.g. can't be a factor.
  stopifnot(class(data[,Y]) %in% c("numeric", "integer"))

  if (family == "binomial" &&
      (min(data[,Y], na.rm = TRUE) < 0 || max(data[,Y], na.rm = TRUE) > 1)) {
    stop("With binomial family Y must be bounded by [0, 1]. Specify family=\"gaussian\" otherwise.")
  }

  if (!family %in% c("binomial", "gaussian")) {
    stop('Family must be either "binomial" or "gaussian".')
  }

  ######################
  # Check NA values, impute with mean and create indicator variable included in W

  impute_results <- impute_NA_vals(data)

  data <- impute_results$data
  W <- c(W, impute_results$`impute cols`)


  if (family == "binomial") {
    ## create the CV folds
    data[, "y_scaled"] <- data[Y]
    data$folds <- create_cv_folds(n_folds, data$y_scaled)
  } else {
    data[, "y_scaled"] <- scale_to_unit(data[Y])
    ## create the CV folds
    data$folds <- create_cv_folds(n_folds, data$y_scaled)
  }

  if (parallel == TRUE) {
    num_cores <- parallel::detectCores()
    doParallel::registerDoParallel(num_cores)
    future::plan(future::multisession)
  }

  data <- as.data.frame(data)
  n <- dim(data)[1]

  fold_results_mix_rules <- fold_results_marg_rules <- fold_results_marg_directions <- list()
  mix_fold_data <- mix_fold_directions <- list()

  fold_results_marginal_data <- fold_results_marginal_additive_data <- list()

  fold_results_marginal_combo_data <- fold_results_marginal_Sls <- list()

  fold_results_mix_combo_data <- fold_results_mix_Sls <- list()

  ###################### Fit Iterative Backwards Rule Ensembles to the Mixture Components Across Fold Data ########################################


  fold_mixture_rules <- foreach::foreach(fold_k = unique(data$folds), .combine = "rbind", .multicombine = TRUE) %do% {

    ## get training and validation samples

    At <- data[data$folds != fold_k, ]
    Av <- data[data$folds == fold_k, ]

    rules <-
      fit_iterative_mix_rule_backfitting(At = At, A = A, W = W, Y = Y, Q1_stack = back_iter_SL, fold = fold_k, verbose)

    rules
  }

  # only keep rules with variables and directions that are found in all the folds

  fold_mixture_rules <- fold_mixture_rules %>%
    dplyr::group_by(test, direction) %>%
    dplyr::filter(dplyr::n() >=  n_folds)

  ###################### Fit Iterative Backwards Rule Ensembles to the Marginal Mixture Components Across Fold Data ########################################


  fold_marginal_rules <- foreach::foreach(fold_k = unique(data$folds), .combine = "rbind", .multicombine = TRUE) %do% {
    ## get training and validation samples

    At <- data[data$folds != fold_k, ]
    Av <- data[data$folds == fold_k, ]


    marg_decisions <- fit_iterative_marg_rule_backfitting(
      mix_comps = A,
      At = At,
      W = W,
      Q1_stack = back_iter_SL,
      fold = fold_k,
      verbose)

    marg_decisions
  }

  fold_marginal_rules <- as.data.frame(fold_marginal_rules)
  fold_marginal_rules$fold <- as.numeric(fold_marginal_rules$fold)

  results <- foreach::foreach(fold_k = unique(data$folds), .combine = "c", .multicombine = TRUE) %dopar% {

    At <- data[data$folds != fold_k, ]
    Av <- data[data$folds == fold_k, ]

    rules <- filter_rules(fold_mixture_rules, fold_k = fold_k)
    marg_decisions <- filter_rules(fold_marginal_rules, fold_k = fold_k)

    if (dim(rules)[1] == 1) {
      if (rules$description == "No Rules Found") {
        no_rules <- TRUE
      }else{
        no_rules <- FALSE
      }
    }else{
      no_rules <- FALSE
    }

    fold_results_mix_rules[[fold_k]] <- rules


    ################## BEGIN FITTING DECISION TREES TO EACH MARGINAL COMPONENT OF THE MIXTURE ###############



    ######################### NUISANCE PARAMETER FITTING ON TRAINING ####################################

    mix_nuisance_params <- est_mix_nuisance_params(
      At,
      Av,
      W,
      no_rules,
      SL.library,
      family,
      rules,
      H.AW_trunc_lvl
    )

    mix_interaction_data <- mix_nuisance_params$data

    mix_fold_data[[fold_k]] <- mix_interaction_data

    ################ MARGINAL ESTIMATES FROM RULE IN TRAINING SET #####################################

    marg_nuisance_params <- est_marg_nuisance_params(At,
                                                     Av,
                                                     W,
                                                     SL.library,
                                                     family,
                                                     A,
                                                     marg_decisions,
                                                     H.AW_trunc_lvl
    )

    marg_decisions <- marg_nuisance_params$marg_decisions

    fold_results_marg_rules[[fold_k]] <- marg_decisions


    fold_results_marginal_data[[fold_k]] <- marg_nuisance_params$data

    marginal_data <- marg_nuisance_params$data

    ################## CALCULATE THE EXPECTED OUTCOME FOR EACH ADDED MARGINAL MIXTURE RULE EXPOSURE ################

    At_c <- At
    Av_c <- Av

    At_c <-
      evaluate_marginal_rules(
        data = At_c,
        marg_decisions = marg_decisions,
        mix_comps = A)
    Av_c <-
      evaluate_marginal_rules(
        data = Av_c,
        marg_decisions = marg_decisions,
        mix_comps = A
      )

    marg_rule_train <- At_c$marginals
    marg_rule_valid <- Av_c$marginals

    At_c <- At_c$data
    Av_c <- Av_c$data

    cum_sum_results <- est_cum_sum_exposure(At_c, Av_c, W, SL.library, H.AW_trunc_lvl)

    fold_results_marginal_additive_data[[fold_k]] <- cum_sum_results

    ################## CALCULATE THE EXPECTED OUTCOME FOR EACH COMBINATION OF MARGINAL MIXTURE RULE EXPOSURE ################

    comb_results <- est_comb_exposure(At, Av, Y, W, marg_rule_train, marg_rule_valid, SL.library)

    fold_results_marginal_combo_data[[fold_k]] <- comb_results$data
    fold_results_marginal_Sls[[fold_k]] <- comb_results$learner

    ################## CALCULATE THE EXPECTED OUTCOME FOR EACH COMBINATION OF THE MIXTURE INTERACTION RULES  ################

    mixed_comb_results <- est_comb_mixture_rules(At, Av, W, Y, rules, no_rules, SL.library)

    fold_results_mix_combo_data[[fold_k]] <- mixed_comb_results$data
    fold_results_mix_Sls[[fold_k]] <- mixed_comb_results$learner

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

  ##################################### AGGREGATE DATA ACROSS FOLDS ######################################

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


  ######################### POOLED TMLE ATE FOR MIXTURES FOUND ACROSS FOLDS.##############################

  mixture_results <- calc_mixtures_ate(
    input_mix_rules = mix_rules,
    input_mix_data = mix_data,
    outcome = Y,
    n_folds = n_folds
  )

  group_list <- mixture_results$group_list
  mixture_results <- mixture_results$results

  ########################################################################################################
  ###################### POOLED TMLE ATE FOR ADDITIVE MARGINAL FOUND ACROSS FOLDS ########################
  ########################################################################################################

  additive_results <- calc_additive_ate(additive_data, Y, n_folds)

  additive_RMSE_star <- additive_results$RMSE_star
  additive_MSM <- additive_results$MSM

  ########################################################################################################
  ###################### POOLED TMLE ATE FOR MARGINAL RULES FOUND ACROSS FOLDS ###########################
  ########################################################################################################

  marginal_results <- calc_marginal_ate(marginal_data, mix_comps = A, Y, n_folds)

  updated_marginal_data <- marginal_results$data
  marginal_results <- marginal_results$marginal_results


  #################################################
  ######### Additive Results Config. ##############
  #################################################

  additive_results <- create_add_marg_ATE_table(updated_marginal_data, A, Y, n, n_folds)


  ########################################################
  ######## Calculate Mixture Rules Over Folds ############
  ########################################################

  mixture_results <- find_common_mixture_rules(group_list,
    data = data,
    mix_comps = A,
    mixture_results = mixture_results,
    n_folds = n_folds
  )


  marginal_results <- find_common_marginal_rules(
    fold_rules = marginal_rules,
    data = data,
    mix_comps = A,
    marginal_results = marginal_results,
    n_folds = n_folds
  )

  ##########################################################################################################
  ################################ FORMAT OUTPUTS FOR OUTPUT  ##############################################
  ##########################################################################################################

  formatted_results <-format_fold_results(marginal_rules,
                                          mix_rules,
                                          mix_combo_data,
                                          marg_combo_data)

  marginal_rules <- formatted_results$marginal_rules
  mixture_rules <- formatted_results$mix_rules
  mix_combo_data <- formatted_results$mix_combo_data
  marg_combo_data <- formatted_results$marg_combo_data

  RMSE_df <- format_RMSE_for_models(marginal_results, mixture_results, additive_RMSE_star, marg_combo_data)

  results_list <- list(
    "Model RMSEs" = RMSE_df,
    "Marginal Results" = marginal_results,
    "Mixture Results" = mixture_results,
    "Additive Results" = additive_results,
    "Additive MSM" = additive_MSM,
    "Mixture Rules" = mixture_rules,
    "Marginal Rules" = marginal_rules,
    "Fold Results Mix Comb Data" = mix_combo_data,
    "Fold Results Marg Comb Data" = marg_combo_data,
    "Fold Super Learners Mix Comb" = mix_Sls,
    "Fold Super Learners Marg Comb" = marg_Sls
  )


  return(results_list)
}
