#' Fit a decision tree on the remaining variance of Y after removing variance attributable to W
#' and other non-target mixture components M
#'
#' Currently uses rpart package but could be generalized
#'
#' @param mix_comps A vector of characters indicating variables for the mixture components
#' @param At Training data
#' @param covariates A vector of characters indicating variables that are covariates
#' on the mixture data to predict remaining variance in Y
#'
#' @param SL.library The Super Learner library to be used to remove the variance in Y|A,W
#' @param family Family ('binomial' or 'gaussian')
#' @param minbucket The minimum number of observations in any terminal <leaf> node.
#' @return Rules object. TODO: add more detail here.
#' @importFrom rpart rpart
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect
#' @importFrom stats na.omit
#' @importFrom purrr discard
#' @importFrom rpart.plot rpart.rules
#'
#' @export

fit_marginal_decision_trees <- function(mix_comps,
                                        At,
                                        covariates,
                                        SL.library,
                                        family,
                                        minbucket) {
  marg_decisions <- list()

  for (i in seq(mix_comps)) {
    target_m <- mix_comps[i]
    covars_m <- mix_comps[-i]

    Y <- At %>% dplyr::select(target_m)
    X_m_part <- At %>% dplyr::select(covariates, covars_m)

    binary_mix_var <- all(na.omit(Y[[1]]) %in% 0:1)

    QbarWSL <- SuperLearner::SuperLearner(
      Y = At$y_scaled,
      X = X_m_part,
      SL.library = SL.library,
      family = family,
      verbose = FALSE
    )


    QbarMi <- predict(QbarWSL, newdata = X_m_part)$pred
    QbarMi_res <- QbarMi - At$y_scaled

    dt_i_At <- cbind(Y, QbarMi_res)

    marg_dt_formula <-
      as.formula(paste("QbarMi_res", target_m, sep = "~"))

    control <- rpart::rpart.control(maxdepth = 1,  minbucket = minbucket, cp = 0.001, xval = 10)

    mixture_rpart_model <- rpart::rpart(
      formula = marg_dt_formula,
      data = dt_i_At,
      method = "anova",
      control = control
    )

    decisions <- rpart.plot::rpart.rules(mixture_rpart_model, roundint = FALSE, digits = 5)
    decisions <- as.data.frame(decisions)
    colnames(decisions)[2] <- "model"

    # TODO: Look into the rpart.rules and see if output could be made cleaner to avoid this coding below

    if (decisions$model[1] != "null model") {
      decisions <- decisions[dim(decisions)[1], 3:dim(decisions)[2]]
      decisions <- decisions %>% purrr::discard(~ all(is.na(.) | . == ""))
      if (any(grepl("to", decisions)) == TRUE) {
        if (decisions[tuple::matchAll("to", decisions) - 1] == decisions[tuple::matchAll("to", decisions) + 1]) {
          decisions[tuple::matchAll("to", decisions) - 1] <- ""
          decisions[tuple::matchAll("to", decisions)] <- "<="
        }
      }

      decisions <- paste(decisions, collapse = " ")
      decisions <- base::strsplit(decisions, "[[:space:]]")[[1]]
      indxs <- match(mix_comps, decisions)

      to_indxs <- tuple::matchAll("to", decisions)

      to_var_inxs <- base::sort(c(indxs, to_indxs))

      for (j in to_indxs) {
        to_in_tovar_str <- base::match(j, to_var_inxs)
        trgt_var_idx <- to_var_inxs[to_in_tovar_str - 1]
        tgt_var <- decisions[trgt_var_idx]

        decisions[j + 1] <-
          paste("&", tgt_var, " <= ", decisions[j + 1])
        decisions[j - 1] <-
          paste("", tgt_var, " >= ", decisions[j - 1])
        decisions[trgt_var_idx] <- ""
      }

      decisions <- base::gsub("^to$", "", decisions)
      decisions <- base::gsub("^is$", "", decisions)


      if (binary_mix_var == FALSE) {
        decisions <- base::gsub("^is$", "", decisions)
      } else {
        decisions <- base::gsub("^is$", " == ", decisions)
      }

      decisions <- base::gsub("^to$", "", decisions)

      decisions <- base::paste(decisions, collapse = "")

      if (stringr::str_detect(decisions, "<=")) {
        decisions <- base::gsub("<=", " <= ", decisions)
      } else {
        decisions <- base::gsub("<", " < ", decisions)
      }

      if (stringr::str_detect(decisions, ">=")) {
        decisions <- base::gsub(">=", " >= ", decisions)
      } else {
        decisions <- base::gsub(">", " < ", decisions)
      }

      marg_decisions[i] <- decisions
    } else {
      marg_decisions[i] <- "1"
    }
  }

  return(marg_decisions)
}
