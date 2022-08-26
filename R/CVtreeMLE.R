#' @title Fit ensemble decision trees to a vector of exposures and use targeted
#' maximum likelihood estimation to determine the average treatment effect
#' in each leaf of best fitting tree
#'
#' @description Fit ensemble decision trees on a mixed exposure while
#' controlling for covariates using iterative backfitting of two Super Learners.
#' If partitioning nodes are identified, use these partitions as a rule-based
#' exposure. The CV-TMLE framework is used to create training and estimation
#' samples. Trees are fit to the training and the average treatment effect (ATE)
#' of the rule-based exposure is estimated in the validation folds. Any type of
#' mixed exposure (continuous, binary, multinomial) is accepted. The ATE for
#' multiple mixture components (interactions) are given as well as marginal
#' effects if data-adaptively identified.
#'
#' @param data Data frame of (W,A,Y) variables of interest.
#' @param w A character vector indicating which variables in the data to use
#'  as baseline covariates.
#' @param a A character vector indicating which variables in the data to use
#'  as exposures.
#' @param y A character indicating which variable in the data to use as the
#' outcome.
#' @param a_stack Stack of estimators used in the Super Learner during
#' the iterative backfitting for `Y|A`, this should be an SL3 object.
#' If not provided, \code{utils_create_sls} is used to create default
#' decision tree estimators used in the ensemble.
#' @param w_stack Stack of estimators used in the Super Learner during the
#' iterative backfitting for `Y|W`, this should be an SL3 stack. If not
#' provided, \code{utils_create_sls} is used to create default
#' estimators used in the ensemble.
#' @param aw_stack Stack of estimators used in the Super Learner for the Q and g
#' mechanisms. If not
#' provided, \code{utils_create_sls} is used to create default
#' estimators used in the ensemble.
#' @param n_folds Number of cross-validation folds.
#' @param seed Pass in a seed number for consistency of results. If not provided
#' a default seed is generated.
#' @param family Family ('binomial' or 'gaussian').
#' @param h_aw_trunc_lvl Truncation level for the clever covariate.
#' @param max_iter Max number of iterations of iterative backfitting algorithm.
#' Default is `5`.
#' @param verbose If true, a message will be printed indicating what process
#' is being started by CVtreeMLE. During the iterative backfitting procedure
#' the iteration in each fold along with the difference in model fit and
#' rules determined will be printed. After the iterative backfitting procedures
#' a table of rules determined in each fold will be printed.
#' @param parallel Use parallel processing if a backend is registered; enabled
#' by default.
#' @param parallel_cv Use parallel processing on CV procedure vs. parallel
#' processing on Super Learner model fitting
#' @param parallel_type default is `multi_session`, if parallel is true
#' which type of parallelization to do `multi_session` or `multicore`
#' @param num_cores If using parallel, the number of cores to parallelize over
#'
#' @details The function performs the following functions.
#'  \enumerate{
#'  \item Imputes missing values with the mean and creates dummy indicator
#'  variables for imputed variables.
#'  \item Separate out covariates into factors and continuous (ordered).
#'  \item Create a variable which indicates the fold number assigned to each
#'  observation.
#'  \item Fit iterative backfitting algorithm onto the mixed exposure
#'  which applies ensemble decision trees to the mixed exposure and an
#'  unrestricted Super Learner on the covariates. Algorithms are fit, offset by
#'  their compliment until there is virtually no difference between the model
#'  fits. Extract partition nodes found for the mixture. This is done on each
#'  training fold data.
#'  \item Fit iterative backfitting algorithm onto each individual mixture
#'  component which applies ensemble decision trees to the mixed exposure and an
#'  unrestricted Super Learner on the covariates. Algorithms are fit, offset by
#'  their compliment until there is virtually no difference between the model
#'  fits. Extract partition nodes found for the mixture. This is done on each
#'  training fold data.
#'  \item Estimate nuisance parameters (Q and g estimates) for mixture
#'  interaction rule
#'  \item Estimate nuisance parameters (Q and g estimates) for marginal
#'  rules
#'  \item Estimate the Q outcome mechanism over all the marginal rules for
#'  later user input for targeted ATE for different marginal combinations based
#'  on data-adaptively identified thresholds.
#'  \item Use the mixture rules and data and do a TMLE fluctuation step to
#'  target the ATE for the given rule across all the folds. Calculate
#'  proportion of folds the rule is found.
#'  \item Use the marginal rules and data and do a TMLE fluctuation step to
#'  target the ATE for the given rule across all the folds. Calculate
#'  proportion of folds the rule is found.
#'  \item Calculate V-fold specific TMLE estimates of the rules.
#'  \item For the mixture rules, calculate a union rule or the rule that covers
#'  all the observations across the folds that the respective variable set in
#'  the rule.
#'  \item For the marginal rules, calculate a union rule or the rule that covers
#'  all the observations across the folds that the respective variable set in
#'  the rule.
#' }
#'
#' @importFrom foreach %dopar% %do%
#' @importFrom magrittr %>%
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif
#' @importFrom stats rnorm runif
#' @importFrom rlang :=
#' @importFrom dplyr filter
#' @importFrom data.table rbindlist
#' @import furrr

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
#'
#' @examples
#' n <- 800
#' p <- 4
#' x <- matrix(rnorm(n * p), n, p)
#' colnames(x) <- c("A1", "A2", "W1", "W2")
#' y_prob <- plogis(3 * sin(x[, 1]) + sin(x[, 2]), sin(x[, 4]))
#' Y <- rbinom(n = n, size = 1, prob = y_prob)
#' data <- as.data.frame(cbind(x,Y))
#'
#' CVtreeMLE_fit <- CVtreeMLE(data = data,
#'                           w = c("W1", "W2"),
#'                           a = c("A1", "A2"),
#'                           y = "Y",
#'                           family = "binomial",
#'                           parallel = FALSE,
#'                           n_folds = 2)
#'
#' @return Object of class \code{CVtreeMLE}, containing a list of table results
#' for: marginal ATEs, mixture ATEs, RMSE of marginal model fits, RMSE of
#' mixture model fits, marginal rules, and mixture rules.
#'
#' \itemize{
#'   \item Model RMSEs: Root mean square error for marginal and interaction
#'   models in the iterative backfitting procedure
#'   \item Pooled TMLE Marginal Results: Data frame of pooled TMLE Marginal
#'   Results: Pooled ATE results using TMLE for thresholds identified for
#'   each mixture component found
#'   \item V-Specific Marg Results: A list of the v-fold marginal results.
#'   These are grouped by variable and direction of the ATE.
#'   \item Pooled TMLE Mixture Results: Data frame of pooled TMLE Mixture
#'   Results
#'   \item V-Specific Mix Results: A list of the v-fold mixture results.
#'   These are grouped by variable and direction of the ATE.
#'   \item Pooled Marginal Refs: A data frame of the reference categories
#'   determined in each of the marginal results.
#'   \item Marginal Rules: A data frame that includes the marginal rules and
#'   details related to fold found and RMSE
#'   \item Mixture Rules: A data frame that includes the mixture rules and
#'   details related to fold found and RMSE
#' }
#'
#' @export

# Start CVtreeMLE ---------------------------

CVtreeMLE <- function(w,
                      a,
                      y,
                      data,
                      w_stack = NULL,
                      aw_stack = NULL,
                      a_stack = NULL,
                      n_folds,
                      seed = 6442,
                      family,
                      parallel = TRUE,
                      parallel_cv = TRUE,
                      parallel_type = "multi_session",
                      num_cores = 2,
                      max_iter = 5,
                      verbose = FALSE) {

  if (any(sapply(data[, a], is.factor))) {
    print("Factor variable detected in exposures, converting to numeric")
    data[, a] <- sapply(data[, a], as.numeric)
  }

  if (any(sapply(data[, w], is.factor))) {
    print("Factor variable detected in the covariate space,
          converting to numeric to avoid issues in different factors
          between training and validation folds")
    data[, w] <- sapply(data[, w], as.numeric)
  }

  # Ensure that Y is numeric; e.g. can't be a factor.
  stopifnot(class(data[, y]) %in% c("numeric", "integer"))

  if (anyNA(data[, y])) {
    stop("NA values cannot be in the outcome variable")
  }

  if (family == "binomial" &&
    (min(data[, y], na.rm = TRUE) < 0 || max(data[, y], na.rm = TRUE) > 1)) {
    stop("With binomial family Y must be bounded by [0, 1].
         Specify family=\"gaussian\" otherwise.")
  }

  if (!family %in% c("binomial", "gaussian")) {
    stop('Family must be either "binomial" or "gaussian".')
  }

  set.seed(seed)

  if (inherits(data, "data.frame")) {
    data <- as.data.frame(data)
  }

  if (is.null(w_stack)) {
    sls <- create_sls()
    w_stack <- sls$W_stack
  }

  if (is.null(a_stack)) {
    sls <- create_sls()
    a_stack <- sls$A_stack
  }

  if (is.null(aw_stack)) {
    sls <- create_sls()
    aw_stack <- sls$AW_stack
  }


  if (family == "binomial") {
    ## create the CV folds
    data[, "y_scaled"] <- data[y]
    data$folds <- create_cv_folds(n_folds, data$y_scaled)
  } else {
    data[, "y_scaled"] <- scale_to_unit(data[y])
    ## create the CV folds
    data$folds <- create_cv_folds(n_folds, data$y_scaled)
  }

  if (parallel == TRUE) {
    if( parallel_type == "multi_session" ) {
    future::plan(future::multisession,
                 workers = num_cores,
                 gc = TRUE)
    }else {
      future::plan(future::multicore,
                   workers = num_cores,
                   gc = TRUE)
    }
  } else {
    future::plan(future::sequential,
                 gc = TRUE)
  }

  data <- as.data.frame(data)

  fold_results_mix_rules <- list()
  fold_results_marg_rules <- list()
  fold_results_marg_directions <- list()

  mix_fold_data <- list()
  mix_fold_directions <- list()

  fold_results_marginal_data <- list()

  fold_results_marginal_combo_data <- list()

  fold_results_mix_combo_data <- list()
  fold_results_mix_Sls <- list()

  # Iterative Back-fitting on Mixture ---------------------------

  if (verbose) {
    print("-----Starting Iterative Backfitting of Exposures on Covariates-----")
  }

  fold_mixture_rules <- furrr::future_map_dfr(unique(data$folds),
                                              function(fold_k) {

    at <- data[data$folds != fold_k, ]
    av <- data[data$folds == fold_k, ]

    rules <-
      fit_mix_rule_backfitting(
        at = at,
        a = a,
        w = w,
        y = y,
        w_stack,
        fold = fold_k,
        max_iter,
        verbose,
        parallel_cv = parallel_cv,
        seed = seed
      )

    rules
  }, .options = furrr::furrr_options(seed = seed, packages = "CVtreeMLE"))


  fold_mixture_rules <- fold_mixture_rules[
    fold_mixture_rules$description != "No Rules Found", ]

  fold_mixture_rules$description <- round_rules(fold_mixture_rules$description)

  mixture_RMSE <- calc_mixture_rule_rmses(fold_mixture_rules)

  if (dim(fold_mixture_rules)[1] == 0) {
    no_mixture_rules <- TRUE
  } else {
    no_mixture_rules <- FALSE
  }

  if (verbose) {
      print("Mixture results found")
      print(fold_mixture_rules)
    }

  # Iterative Back-fitting on Marginals ---------------------------

  if (verbose) {
    print("-----Starting Iterative Backfitting of Each Exposure on Covariates to
          Determine if Partitioning Nodes Exist-----")
  }

  fold_marginal_rules <- furrr::future_map_dfr(unique(data$folds),
                                               function(fold_k) {

    at <- data[data$folds != fold_k, ]
    av <- data[data$folds == fold_k, ]

    marg_decisions <- fit_marg_rule_backfitting(
      mix_comps = a,
      at = at,
      w = w,
      w_stack = w_stack,
      tree_stack = a_stack,
      fold = fold_k,
      max_iter = max_iter,
      verbose = verbose,
      parallel_cv = parallel_cv,
      seed = seed
    )

    marg_decisions
  }, .options = furrr::furrr_options(seed = seed, packages = "CVtreeMLE"))

  filt_fold_marginal_rules <- filter_marginal_rules(fold_marginal_rules, n_folds)
  filt_fold_marginal_rules$rules <- round_rules(filt_fold_marginal_rules$rules)
  marginal_RMSE <- calc_marginal_rule_rmses(data = filt_fold_marginal_rules)

  backfit_model_RMSEs <- rbind(marginal_RMSE, mixture_RMSE)

  if (dim(filt_fold_marginal_rules)[1] == 0) {
    no_marginal_rules <- TRUE
  } else {
    no_marginal_rules <- FALSE
  }

  if (no_marginal_rules == TRUE && no_mixture_rules == TRUE) {
    return("No Mixture or Marginal Rules Found")
  }

  if (verbose) {
    print("Marginal results found")
    print(filt_fold_marginal_rules)
  }

  # Estimate nuisance parameters ---------------------------


  if (verbose) {
    print("-----Starting Fitting Estimators to the Rules Determined to Calculate
          the Average Treatment Effects of Terminal Nodes in Trees-----")
  }

  results <- furrr::future_map(unique(data$folds), function(fold_k) {

    at <- data[data$folds != fold_k, ]
    av <- data[data$folds == fold_k, ]

    rules <- filter_rules(fold_mixture_rules, fold_k = fold_k)
    marg_decisions <- filter_rules(filt_fold_marginal_rules, fold_k = fold_k)

    if (dim(rules)[1] == 0) {
      no_mix_rules <- TRUE
    } else {
      no_mix_rules <- FALSE
    }

    if (dim(marg_decisions)[1] == 0) {
      no_marg_rules <- TRUE
    } else {
      no_marg_rules <- FALSE
    }

    fold_results_mix_rules[[paste("Fold", fold_k)]] <- rules

    mix_nuisance_params <- est_mix_nuisance_params(
      at,
      av,
      w,
      no_mix_rules,
      aw_stack,
      family,
      rules,
      h_aw_trunc_lvl,
      parallel_cv,
      seed
    )

    mix_interaction_data <- mix_nuisance_params$data

    mix_fold_data[[paste("Fold", fold_k)]] <- mix_interaction_data


    marg_nuisance_params <- est_marg_nuisance_params(
      at,
      av,
      w,
      aw_stack,
      family,
      a,
      no_marg_rules,
      marg_decisions,
      h_aw_trunc_lvl,
      parallel_cv,
      seed
    )

    non_ref_rules <- marg_decisions[marg_decisions$quantile > 1, ]

    names(marg_nuisance_params) <- non_ref_rules$var_quant_group

    fold_results_marginal_data[[paste("Fold", fold_k)]] <- marg_nuisance_params
    fold_results_marg_rules[[paste("Fold", fold_k)]] <- marg_decisions

    marginal_data <- marg_nuisance_params

    rules_evaluated_on_training <-
      evaluate_marginal_rules(
        data = at,
        marg_decisions = marg_decisions,
        no_marg_rules,
        mix_comps = a
      )

    rules_evaluated_on_valiation <-
      evaluate_marginal_rules(
        data = av,
        marg_decisions = marg_decisions,
        no_marg_rules,
        mix_comps = a
      )

    marg_rule_train <- rules_evaluated_on_training$marginals
    marg_rule_valid <- rules_evaluated_on_valiation$marginals

    at_w_rule_evals <- rules_evaluated_on_training$data
    av_w_rule_evals <- rules_evaluated_on_valiation$data

    comb_results <- est_comb_exposure(at,
                                      av,
                                      y,
                                      w,
                                      marg_rule_train,
                                      marg_rule_valid,
                                      no_marg_rules,
                                      aw_stack,
                                      family = family,
                                      parallel_cv,
                                      seed)

    fold_results_marginal_combo_data[[
      paste("Fold", fold_k)]] <- comb_results$data


    results_list <- list(
      fold_results_mix_rules,
      mix_fold_data,
      fold_results_marg_rules,
      fold_results_marginal_data,
      fold_results_marginal_combo_data
    )

    names(results_list) <- c(
      "mix rules",
      "mix data",
      "marg rules",
      "marg data",
      "marg combo data")

    results_list
  }, .options = furrr::furrr_options(seed = seed, packages = "CVtreeMLE"))


  # Aggregate results ---------------------------

  marginal_data <- purrr::map(results, c("marg data"))
  marginal_rules <- purrr::map(results, c("marg rules"))

  mix_data <- purrr::map(results, c("mix data"))
  mix_rules <- purrr::map(results, c("mix rules"))

  marg_combo_data <- purrr::map(results, c("marg combo data"))


  if (no_mixture_rules == FALSE) {

    # Calculate mixture ATEs ---------------------------

    mixture_results <- calc_mixtures_ate(
      input_mix_rules = mix_rules,
      input_mix_data = mix_data,
      outcome = y,
      n_folds = n_folds,
      no_mixture_rules = no_mixture_rules
    )

    group_list <- mixture_results$group_list
    mixture_data_list <- mixture_results$mixture_data_list
    mixture_results <- mixture_results$results

    # Calculate v-fold mixture ATEs ---------------------------

    v_fold_mixture_results <- calc_v_fold_mixtures_ate(
      input_mix_rules = mix_rules,
      input_mix_data = mix_data,
      outcome = y,
      n_folds = n_folds
      )

    v_fold_mixture_results <- v_fold_mixture_results
    v_fold_mixture_results_w_pooled <- meta_mix_results(v_fold_mixture_results,
                                                                mix_comps = a,
                                                                n_folds,
                                                                data = data)

    # Create union mixture rules ---------------------------

      mixture_results <- common_mixture_rules(group_list,
                                                 data = data,
                                                 mix_comps = a,
                                                 mixture_results,
                                                 n_folds = n_folds,
                                                 no_mixture_rules
    )

    mix_rules <- unlist(mix_rules, recursive = FALSE, use.names = FALSE)
    mix_rules <- mix_rules[!sapply(mix_rules, is.null)]
    mix_rules <-
      data.table::rbindlist(mix_rules)


  } else {
    mixture_results <- NA
    mixture_rules <- NA
    v_fold_mixture_results_w_pooled <- NA
  }

  if (no_marginal_rules == FALSE) {
    # Calculate marginal ATEs ---------------------------

    marginal_results <- calc_marginal_ate(marginal_data,
                                          mix_comps = a,
                                          marginal_rules,
                                          y = y,
                                          n_folds)

    updated_marginal_data <- marginal_results$data
    marginal_results <- marginal_results$marginal_results

    # Calculate marginal v-fold ATEs ---------------------------

    v_fold_marginal_results <- calc_v_fold_marginal_ate(marginal_data,
                                                        mix_comps = a,
                                                        marginal_rules,
                                                        y,
                                                        n_folds)

    v_fold_marginal_results <- v_fold_marginal_results$marginal_results

    v_fold_marginal_results_w_pooled <- compute_meta_marg_results(
                              v_fold_marginal_results = v_fold_marginal_results,
                              mix_comps = a,
                              data = data,
                              n_folds = n_folds)

    # Create union marginals rules ---------------------------

    total_marginal_rules <- find_common_marginal_rules(
      fold_rules = marginal_rules,
      data = data,
      mix_comps = a,
      marginal_results = marginal_results,
      n_folds = n_folds
    )

    # Parse out ref results ---------------------------
    non_ref_rules <- total_marginal_rules[
      !str_detect(total_marginal_rules$`Variable Quantile`, "_1"), ]

    non_ref_rules <- non_ref_rules[
      order(non_ref_rules$`Variable Quantile`), ]

    ref_rules <- total_marginal_rules[
      str_detect(total_marginal_rules$`Variable Quantile`, "_1"), ]

    marginal_results <- marginal_results[order(rownames(marginal_results)), ]
    marginal_results <- cbind(marginal_results, non_ref_rules)

    marg_vars <- stringr::str_extract(string =
                                        marginal_results$`Variable Quantile`,
                                      pattern = paste(a, collapse = "|"))

    marginal_results$var <- marg_vars

    marginal_rules <- unlist(marginal_rules,
                             recursive = FALSE,
                             use.names = FALSE)

    marginal_rules <- marginal_rules[!sapply(marginal_rules, is.null)]
    marginal_rules <- do.call(rbind, marginal_rules)
  }

  results_list <- list(
    "Model RMSEs" = backfit_model_RMSEs,
    "Pooled TMLE Marginal Results" = marginal_results,
    "V-Specific Marg Results" = v_fold_marginal_results_w_pooled,
    "Pooled TMLE Mixture Results" = mixture_results,
    "V-Specific Mix Results" = v_fold_mixture_results_w_pooled,
    "Pooled Marginal Refs" = ref_rules,
    "Mixture Rules" = mix_rules,
    "Marginal Rules" = marginal_rules
  )

  class(results_list) <- "CVtreeMLE"


  return(results_list)
}
