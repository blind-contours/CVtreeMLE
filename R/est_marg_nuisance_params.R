#' For each marginal mixture component rule found, create a g estimator for the probability of being exposed to the rule thresholds,
#' and a Q estimator for the outcome E(Y| A = a_mix, W). Get estimates of g and Q using the validation data and
#' calculate the clever covariate used in the TMLE fluctuation step.
#'
#' @param At Training data
#' @param Av Validation data
#' @param W Vector of characters denoting covariates
#' @param Q1_stack Super Learner library for fitting Q (outcome mechanism) and g (treatment mechanism)
#' @param family Binomial or gaussian
#' @param A Vector of characters that denote the mixture components
#' @param marg_decisions List of rules found within the fold for each mixture component
#' @param H.AW_trunc_lvl Truncation level of the clever covariate (induces more bias to reduce variance)
#' @param no_marg_rules TRUE/FALSE if no marginal rules were found across the folds
#'
#' @import sl3
#' @importFrom magrittr %>%
#' @importFrom rlang :=

#' @importFrom dplyr group_by filter top_n
#' @return Rules object. TODO: add more detail here.

#'
#' @export

est_marg_nuisance_params <- function(At,
                                     Av,
                                     W,
                                     Q1_stack,
                                     family,
                                     A,
                                     no_marg_rules,
                                     marg_decisions,
                                     H.AW_trunc_lvl) {

  future::plan(future::sequential, gc = TRUE)

  marginal_data <- list()

  marg_decisions_groups <- marg_decisions %>% dplyr::group_by(target_m)
  marg_decisions_groups <- dplyr::group_split(marg_decisions_groups)

  if (no_marg_rules == FALSE) {
    for (i in seq(marg_decisions_groups)) {
      At_c <- At
      Av_c <- Av

      variable_decisions <-marg_decisions_groups[[i]]

      quant_one_row <- variable_decisions[variable_decisions$quantile == 1,]
      quant_one_rule <- quant_one_row$rules

      # base_rule_name <- paste(quant_one_row$target_m, "ref_rule", sep = "_")

      At_c_ref_data <- At_c %>%
        dplyr::mutate("A" := ifelse(eval(parse(text = quant_one_rule)), 1, 0))

      Av_c_ref_data <- Av_c %>%
        dplyr::mutate("A" := ifelse(eval(parse(text = quant_one_rule)), 1, 0))

      At_c_ref_data <- At_c_ref_data[At_c_ref_data[, "A"] == 1,]
      Av_c_ref_data <- Av_c_ref_data[Av_c_ref_data[,"A"] == 1,]

      At_c_ref_data[, "A"] <- 0
      Av_c_ref_data[, "A"] <- 0

      quant_comparisons <- variable_decisions[variable_decisions$quantile > 1,]

      in_group_marg_data <- list()

      for (j in seq(nrow(quant_comparisons))) {

        target_m_row <- quant_comparisons[j,] # marg_decisions[marg_decisions$target_m == A[i],]$rules

        # comp_rule_name <- paste(target_m_row$var_quant_group, "rule", sep = "_")

        At_c_comp_data <- At_c %>%
          dplyr::mutate("A" := ifelse(eval(parse(text = target_m_row$rules)), 1, 0))

        Av_c_comp_data <- Av_c %>%
          dplyr::mutate("A" := ifelse(eval(parse(text = target_m_row$rules)), 1, 0))

        At_c_comp_data <- At_c_comp_data[At_c_comp_data[,"A"] == 1,]
        Av_c_comp_data <- At_c_comp_data[At_c_comp_data[,"A"] == 1,]

        At_data <- rbind(At_c_ref_data, At_c_comp_data)
        Av_data <- rbind(Av_c_ref_data, Av_c_comp_data)

        task_At <- sl3::make_sl3_Task(
          data = At_data,
          covariates = W,
          outcome = "A",
          outcome_type = "binomial"
        )

        task_Av <- sl3::make_sl3_Task(
          data = Av_data,
          covariates = W,
          outcome = "A",
          outcome_type = "binomial"
        )

        discrete_sl_metalrn <- sl3::Lrnr_cv_selector$new()

        discrete_sl <- sl3::Lrnr_sl$new(
          learners = Q1_stack,
          metalearner = discrete_sl_metalrn,
        )

        sl_fit <- discrete_sl$train(task_At)


        gHat1W <- bound_precision(sl_fit$predict(task_Av))

        H.AW <- calc_clever_covariate(gHat1W = gHat1W, data = Av_data, exposure = "A", H.AW_trunc_lvl = 10, type = "reg")

        ## add treatment mechanism results to Av_c dataframe
        Av_data$gHat1W <- gHat1W
        Av_data$H.AW <- H.AW

        task_At <- sl3::make_sl3_Task(
          data = At_data,
          covariates = c(W, "A"),
          outcome = "y_scaled",
          outcome_type = family,
          folds = 2
        )

        X_m1 <- X_m0 <- Av_data
        X_m1$A <- 1 # under exposure
        X_m0$A <- 0 # under control

        task_Av <- sl3::make_sl3_Task(
          data = Av_data,
          covariates = c(W, "A"),
          outcome = "y_scaled",
          outcome_type = family
        )

        task_Av_1 <- sl3::make_sl3_Task(
          data = X_m1,
          covariates = c(W, "A"),
          outcome = "y_scaled",
          outcome_type = family
        )

        task_Av_0 <- sl3::make_sl3_Task(
          data = X_m0,
          covariates = c(W, "A"),
          outcome = "y_scaled",
          outcome_type = family
        )

        sl_fit <- discrete_sl$train(task_At)


        QbarAW <- bound_precision(sl_fit$predict(task_Av))
        Qbar1W <- bound_precision(sl_fit$predict(task_Av_1))
        Qbar0W <- bound_precision(sl_fit$predict(task_Av_0))

        flux_results <- fit_least_fav_submodel(H.AW, data = Av_data, QbarAW, Qbar1W, Qbar0W)

        QbarAW.star <- flux_results$QbarAW.star
        Qbar1W.star <- flux_results$Qbar1W.star
        Qbar0W.star <- flux_results$Qbar0W.star

        Av_data$QbarAW <- QbarAW
        Av_data$Qbar1W <- Qbar1W
        Av_data$Qbar0W <- Qbar0W

        in_group_marg_data[[j]] <- Av_data
      }
      marginal_data[[i]] <- in_group_marg_data
    }

    marginal_data <- unlist(marginal_data,recursive=FALSE, use.names = FALSE)

    }
  else {
      Av_data$gHat1W <- NA
      Av_data$H.AW <- NA
      Av_data$QbarAW <- NA
      Av_data$Qbar1W <- NA
      Av_data$Qbar0W <- NA
      marginal_data[[1]] <- Av_data
      print(
        paste(
          "No ATEs calculated in the validation for",
          marg_decisions$target_m[i],
          "due to no rule found in training set for marginal impact"
        )
      )
    }

    return(marginal_data)
  }
