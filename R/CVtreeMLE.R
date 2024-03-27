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
#' @param family Family ('binomial' or 'continuous').
#' @param region If a predetermined region is of interest, put here: like "A < 0.02"
#' @param parallel Use parallel processing if a backend is registered; enabled
#' by default.
#' @param parallel_cv Use parallel processing on CV procedure vs. parallel
#' processing on Super Learner model fitting
#' @param parallel_type default is `multi_session`, if parallel is true
#' which type of parallelization to do `multi_session` or `multicore`
#' @param num_cores If using parallel, the number of cores to parallelize over
#' @param h_aw_trunc_lvl Level to truncate the clever covariate
#' to control variance, default is 10.
#' @param pooled_rule_type is "average" or "union" how to construct the rule
#' across folds. The average take the average cutpoints and returns an average
#' rule with lower and upper bounds for each cutpoint. The union rule creates
#' a new rule that is the space that contains all the rules found across the
#' fold and is therefore more conservative.
#' @param min_max Which oracle region to go after the one that minimizes or maximizes the outcome.
#' @param min_obs Minimum number of observations to have in a region.
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
#' data <- as.data.frame(cbind(x, Y))
#'
#' CVtreeMLE_fit <- CVtreeMLE(
#'   data = data,
#'   w = c("W1", "W2"),
#'   a = c("A1", "A2"),
#'   y = "Y",
#'   family = "binomial",
#'   parallel = FALSE,
#'   n_folds = 2
#' )
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
                      h_aw_trunc_lvl = 50,
                      pooled_rule_type = "average",
                      min_max = "min",
                      region = NULL,
                      min_obs = 25) {
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

  if (anyNA(data[, a])) {
    stop("NA values cannot be in the exposure variables")
  }

  if (family == "binomial" &&
    (min(data[, y], na.rm = TRUE) < 0 || max(data[, y], na.rm = TRUE) > 1)) {
    stop("With binomial family Y must be bounded by [0, 1].
         Specify family=\"continuous\" otherwise.")
  }

  if (!family %in% c("binomial", "continuous")) {
    stop('Family must be either "binomial" or "continuous".')
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

  data$folds <- create_cv_folds(n_folds, data[, y])

  if (parallel == TRUE) {
    if (parallel_type == "multi_session") {
      future::plan(future::multisession,
        workers = num_cores,
        gc = TRUE
      )
    } else {
      future::plan(future::multicore,
        workers = num_cores,
        gc = TRUE
      )
    }
  } else {
    future::plan(future::sequential,
      gc = TRUE
    )
  }

  data <- as.data.frame(data)

  # Set up lists to store data in the CV

  fold_results_mix_rules <- list()
  mix_fold_data <- list()


  if (is.null(region)) {
    fold_mixture_rules <- furrr::future_map(unique(data$folds),
      function(fold_k) {
        at <- data[data$folds != fold_k, ]

        rules <-
          fit_min_ave_tree_algorithm(
            at = at,
            a = a,
            w = w,
            y = y,
            fold = fold_k,
            parallel_cv = parallel_cv,
            min_max = min_max,
            min_obs = min_obs
          )


        rules
      },
      .options = furrr::furrr_options(seed = seed, packages = c(
        "CVtreeMLE",
        "ranger"
      ))
    )

    mixture_rules <- purrr::map(
      fold_mixture_rules,
      c("rules")
    )
    mixture_rules <- do.call(rbind, mixture_rules)

    fold_mixture_rules <- mixture_rules[
      mixture_rules$description != "No Rules Found",
    ]

    fold_mixture_rules$description <- round_rules(
      fold_mixture_rules$description
    )
  } else {
    fold_mixture_rules <- data.frame(
      rule = rep("rule1", n_folds),
      coefficient = rep(NA, n_folds),
      description = rep(region, n_folds),
      test = sort(paste(
        str_extract_all(region, paste(c(a, w),
          collapse = "|"
        ))[[1]],
        collapse = "-"
      )),
      direction = rep(NA, n_folds),
      fold = 1:n_folds,
      RMSE = rep(NA, n_folds),
      effect_modifiers = as.numeric(str_detect(region, paste(w,
        collapse = "|"
      )))
    )
  }

  # Estimate nuisance parameters ---------------------------


  results <- furrr::future_map(unique(data$folds), function(fold_k) {
    at <- data[data$folds != fold_k, ]
    av <- data[data$folds == fold_k, ]

    rules <- filter_rules(fold_mixture_rules, fold_k = fold_k)

    fold_results_mix_rules[[paste("Fold", fold_k)]] <- rules

    mix_nuisance_params <- est_mix_nuisance_params(
      at = at,
      av = av,
      w = w,
      a = a,
      y = y,
      aw_stack = aw_stack,
      family = family,
      rules = rules,
      parallel_cv = parallel_cv,
      seed = seed,
      h_aw_trunc_lvl = h_aw_trunc_lvl
    )

    mix_interaction_data <- mix_nuisance_params$data

    mix_fold_data[[paste("Fold", fold_k)]] <- mix_interaction_data

    results_list <- list(
      fold_results_mix_rules,
      mix_fold_data
    )

    names(results_list) <- c(
      "mix rules",
      "mix data"
    )

    results_list
  }, .options = furrr::furrr_options(seed = seed, packages = c(
    "CVtreeMLE",
    "sl3"
  )))


  mix_data <- purrr::map(results, c("mix data"))
  mix_rules <- purrr::map(results, c("mix rules"))

  mixture_results <- calc_mixtures_ate(
    input_mix_rules = mix_rules,
    input_mix_data = mix_data,
    y = y,
    n_folds = n_folds
  )

  group_list <- mixture_results$group_list
  inv_var_mixture_results <- mixture_results$inv_var_results
  oracle_region_results <- mixture_results$region_results

  mixture_results <- mixture_results$results

  v_fold_mixture_results <- calc_v_fold_mixtures_ate(
    input_mix_rules = mix_rules,
    input_mix_data = mix_data,
    y = y,
    n_folds = n_folds
  )

  if (pooled_rule_type == "average") {
    mixture_results <- average_mixture_rules(
      group_list = group_list,
      data = data,
      mix_comps = c(a, w),
      n_folds = n_folds,
      mixture_results = mixture_results
    )

    inv_var_mixture_results <- average_mixture_rules(
      group_list = group_list,
      data = data,
      mix_comps = c(a, w),
      n_folds,
      inv_var_mixture_results
    )
  } else {
    mixture_results <- common_mixture_rules(group_list,
      data = data,
      mix_comps = c(a, w),
      mixture_results,
      n_folds = n_folds)

    inv_var_mixture_results <- common_mixture_rules(group_list,
      data = data,
      mix_comps = c(a, w),
      inv_var_mixture_results,
      n_folds = n_folds)
  }

  mix_rules <- unlist(mix_rules, recursive = FALSE, use.names = FALSE)
  mix_rules <- mix_rules[!sapply(mix_rules, is.null)]
  mix_rules <-
    data.table::rbindlist(mix_rules)


  results_list <- list(
    "Pooled TMLE Mixture Results" = mixture_results,
    "Pooled TMLE Inver Variance Results" = inv_var_mixture_results,
    "V-Specific Mix Results" = v_fold_mixture_results,
    "Oracle Region Results" = oracle_region_results,
    "Mixture Rules" = mix_rules
  )

  class(results_list) <- "CVtreeMLE"


  return(results_list)
}
