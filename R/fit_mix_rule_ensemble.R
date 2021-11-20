#' Fit decision tree ensembles on the remaining variance of Y using all mixture components as features
#'
#' Currently uses PRE Predictive Rule Ensembles package but could be generalized
#'
#' @param At Training dataframe
#' @param formula Formula passed to the pre tree fitting algorithm
#' @param rule_boot_n Number of times to resample data with replacement and refit the pre algorithm
#' on the mixture data to predict remaining variance in Y
#'
#' @param n_stable_trees Number of times a set of variables must appear in the same rule across the bootstrap to be considered stable
#' and move to downstream analysis
#' @param A Variable names in the mixture

#' @importFrom pre pre maxdepth_sampler
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by filter top_n
#' @return Rules object. TODO: add more detail here.

#'
#' @export

fit_mix_rule_ensemble <- function(At, formula, rule_boot_n, n_stable_trees, A) {

  pull_out_rule_vars <- function(x) {
    x <- x[3]
    if (x != "1") {
      x_split <- strsplit(x, " ")[[1]]
      hits <- x_split[grep(paste(A, collapse = "|"), x_split)]
      hits <- sort(hits)
      if (length(hits) == 1) {
        hits <- 0
      } else {
        hits <- paste(hits, collapse = "")
      }
      return(hits)
    } else {
      hits <- NA
    }
    return(hits)
  }

  pre_boot_list <- list()

  for (i in seq(rule_boot_n)) {
    At_boot <- At[base::sample(nrow(At), dim(At)[1], replace = FALSE), ]

    pre_model <- pre::pre(formula,
      data = At_boot,
      family = "gaussian",
      use.grad = FALSE,
      tree.unbiased = TRUE,
      removecomplements = TRUE,
      removeduplicates = TRUE,
      maxdepth = pre::maxdepth_sampler(),
      sampfrac = min(1, (11 * sqrt(dim(At)[1]) + 1) / dim(At)[1]),
      nfolds = 10
    )

    pre_model_coefs <- stats::coef(pre_model, penalty.par.val = "lambda.min")

    pre_coefs_no_zero <- base::subset(pre_model_coefs, coefficient != 0)
    pre_coefs_no_zero$test <- apply(pre_coefs_no_zero, 1, function(x) pull_out_rule_vars(x))

    rules <- pre_coefs_no_zero
    rules$boot_num <- i
    pre_boot_list[[i]] <- rules
  }

  pre_boot_df <- do.call(rbind, pre_boot_list)
  pre_boot_df$direction <- ifelse(pre_boot_df$coefficient > 0, 1, 0)
  rules <- pre_boot_df %>%
    dplyr::group_by(test, direction) %>%
    dplyr::filter(n() >= n_stable_trees)

  rules <- rules[!is.na(rules$test), ]
  rules <- rules[!rules$test == 0, ]

  rules <- rules %>%
    dplyr::group_by(test, direction) %>%
    dplyr::top_n(1, abs(coefficient))

  return(rules)
}
