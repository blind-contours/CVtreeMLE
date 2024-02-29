#' @title Calculates the Inverse Variance Pooled Estimate Including Null Folds
#' @description Does a weighted combination estimate for folds with estimates
#' and null folds finding 0 using inverse variance
#' @param results_df Table of results
#' @param n_folds Total number of folds specified
#' @export

calculatePooledEstimate <- function(est, se,  n_folds, n_intxn ) {
  # Identifying the row with the current pooled TMLE
  # pooledRow <- results_df[results_df$Fold == "Pooled TMLE",]

  # Calculate n_0 and n_1
  n_1 <- n_intxn
  n_0 <- n_folds - n_1

  # Variance of the current pooled estimates
  var_pooled <- se^2

  # Estimating the variance of the "null" estimate
  var_null <- var_pooled * n_1 / n_0

  # Weights
  w_1 <- 1 / var_pooled
  w_0 <- 1/ var_null

  # Calculate new pooled estimate and its variance
  pooled_estimate <- est
  psi_null <- 0  # Assuming null value is 0 for additive estimates
  new_pooled_psi <- (w_1 * pooled_estimate + w_0 * psi_null) / (w_1 + w_0)
  var_new_pooled <- 1 / (w_1 + w_0)

  # Calculate standard error and confidence intervals
  se_new_pooled <- sqrt(var_new_pooled)
  lower_ci <- new_pooled_psi - 1.96 * se_new_pooled
  upper_ci <- new_pooled_psi + 1.96 * se_new_pooled

  p_value_pooled <- 2 * stats::pnorm(abs(new_pooled_psi / se_new_pooled), lower.tail = F)


  # Add new row to the table
  new_row <- data.table(
    `Mixture ATE` = new_pooled_psi,
    `Standard Error` = se_new_pooled,
    `Lower CI` = lower_ci,
    `Upper CI` = upper_ci,
    `P-value` = p_value_pooled,
    `P-value Adj` = p_value_pooled,
    `Vars` = results_df$Vars,
    RMSE = results_df$RMSE
  )

  results_df$Type <- "Pooled TMLE"
  new_row$Type <- "Inverse Variance Pooled TMLE"

  bind_rows(results_df, new_row)
}


#' @title  Pull rules out of results table of pre results
#' @details Simply function that looks at variables used in the pre model
#' for matches
#'  with the mixture variables and pulls them out to check for consistent
#'  sets of rules found over the iterative backfitting procedure
#' @param x Row of a dataframe that contains output from the pre fit
#' @param a Vector of characters indicating column names for the mixture
#' variables
#' @return Vector of matches of mixture variables
#' @export

pull_out_rule_vars <- function(x, a) {
  x <- x[3]
  if (x != "1") {
    x_split <- strsplit(x, " ")[[1]]
    hits <- x_split[grep(paste(a, collapse = "|"), x_split)]
    hits <- sort(unique(hits))
    hits <- paste(hits, collapse = "-")
  } else {
    hits <- NA
  }
  return(hits)
}



###################################################################
#' @title Round rules found for easier reading
#' @param rules Vector or rules
#' @importFrom rlang .data
round_rules <- function(rules) {
  rounded_rules <- list()
  split_rules <- strsplit(rules, " ")
  for (i in seq(split_rules)) {
    rule_split <- split_rules[[i]]
    for (j in seq(rule_split)) {
      element <- rule_split[[j]]
      if (grepl("-*\\d+\\.*\\d*", element) == TRUE &&
        grepl("[[:alpha:]]", element) == FALSE) {
        element_round <- round(as.numeric(element), 3)
        rule_split[[j]] <- element_round
      }
    }
    round_temp <- paste(rule_split, collapse = " ")
    rounded_rules[[i]] <- round_temp
  }

  rounded_rules <- unlist(rounded_rules)
  return(rounded_rules)
}


###################################################################
#' @title Evaluate Rules to Binary Indicators
#' @details Takes in a list of rules and outputs a binary matrix
#' @param rules List of rules
#' @param data Data to evaluate rules on
#' @param Y Outcome that is appended to the binary data
#' @importFrom rlang .data
#'
#' @export
evaluate_rules_to_binary <- function(rules, data, Y) {
  # Create a new data frame to store the binary indicators
  binary_data <- data.frame(matrix(ncol = length(rules), nrow = nrow(data)))
  names(binary_data) <- paste("Rule", seq_along(rules), sep = "_")

  # Evaluate each rule and create a binary indicator column
  for (i in seq_along(rules)) {
    # Construct the expression to evaluate
    expression_to_evaluate <- paste0("as.numeric(", rules[[i]], ")")

    # Safely evaluate the expression within the data environment
    binary_data[[i]] <- as.numeric(eval(parse(text = expression_to_evaluate), envir = data))
  }
  colnames(binary_data) <- rules
  binary_data <- cbind(Y, binary_data)
  return(binary_data)
}

###################################################################
#' @title From HAL fit, create string based rules
#' @details Inputs the basis list from HAL then returns a list of rules
#' @param basis_list Basis list from HAL
#' @param col_names Column names
#' @importFrom rlang .data
#'
#' @export

create_rules <- function(basis_list, col_names) {
  rules <- list() # Initialize a list to hold all rules

  for (i in seq_along(basis_list)) {
    basis <- basis_list[[i]]
    cols <- basis$cols
    cutoffs <- basis$cutoffs
    orders <- basis$orders # If you need to use orders in the rules

    # Initialize an empty string for this rule
    rule <- ""

    for (j in seq_along(cols)) {
      # If this is not the first condition in the rule, prepend an " & "
      if (j > 1) {
        rule <- paste0(rule, " & ")
      }

      # Append the condition to the rule
      col_index <- cols[j]
      col_name <- col_names[col_index]
      cutoff <- cutoffs[j]
      rule <- paste0(rule, col_name, " <= ", cutoff)
    }

    # Add this rule to the list of rules
    rules[[i]] <- rule
  }

  return(rules)
}

###################################################################
#' @title Filter data based on fold
#' @details  Filter data to only the training fold of interest during CV
#' @param data Input data
#' @param fold_k Current fold
#' @importFrom rlang .data
#'
#' @export
filter_rules <- function(data, fold_k) {
  data <- data %>%
    dplyr::filter(.data$fold == fold_k)
  return(data)
}

###################################################################
#' @title Group by fold
#' @param data Input data
#' @importFrom rlang .data
#'
#' @export
groupby_fold <- function(data) {
  data <- data %>%
    dplyr::group_by(.data$fold)
  return(data)
}

###################################################################
#' @title Filter mixture rules across the folds for only those that have the
#' same variables and directions across all folds
#' @param data Input data
#' @param n_folds Number of folds in the CV
#' @importFrom rlang .data
#'
#' @export
filter_mixture_rules <- function(data, n_folds) {
  data <- data %>%
    dplyr::group_by(.data$test, .data$direction) %>%
    dplyr::filter(dplyr::n() >= n_folds)
  return(data)
}

###################################################################
#' @title Calculate the mean RMSE in each interaction group
#' @param data Input data
#' @importFrom rlang .data
#'
#' @export
calc_mixture_rule_rmses <- function(data) {
  data <- data %>%
    dplyr::group_by(.data$test) %>%
    dplyr::summarize(RMSE = mean(.data$RMSE, na.rm = TRUE))
  colnames(data) <- c("Var(s)", "RMSE")

  return(data)
}
###################################################################
#' @title Filter marginal rules across the folds for only those that
#' have the same variables
#' @param data Input data
#' @param n_folds Number of folds in the CV
#' @importFrom rlang .data
#'
#' @export
filter_marginal_rules <- function(data, n_folds) {
  data <- as.data.frame(data)
  data$fold <- as.numeric(data$fold)

  data <- data[data$rules != "No Rules Found", ]

  data$var_quant_group <- paste(data$target_m, data$quantile, sep = "_")
  return(data)
}

###################################################################
#' @title Calculate the mean RMSE in each marginal group
#' @param data Input data
#' @importFrom rlang .data
#'
#' @export
calc_marginal_rule_rmses <- function(data) {
  data <- data %>%
    dplyr::group_by(.data$target_m) %>%
    dplyr::summarize(RMSE = mean(.data$RMSE, na.rm = TRUE))

  colnames(data) <- c("Var(s)", "RMSE")

  return(data)
}

###################################################################
#' @title Group by fold
#' @param data Input data
#' @importFrom rlang .data
#'
#' @export
groupby_fold <- function(data) {
  data <- data %>%
    dplyr::group_by(.data$fold)
  return(data)
}

###################################################################
#' @title Group split by marginal variable
#' @param data Input data
#' @importFrom rlang .data
#'
#' @export
marginal_group_split <- function(data) {
  data <- data %>% dplyr::group_by(.data$target_m)
  data <- dplyr::group_split(data)
  return(data)
}

###################################################################
#' @title v-fold marginal group split
#' @param data Input data
#' @importFrom rlang .data
#'
#' @export
v_fold_marginal_qgroup_split <- function(data) {
  data <- data %>% dplyr::group_by(.data$Levels)
  data <- dplyr::group_split(data)
  return(data)
}

###################################################################
#' @title v-fold group split
#' @param data Input data
#' @importFrom rlang .data
#'
#' @export
v_fold_mixture_group_split <- function(data) {
  data <- data %>% dplyr::group_by(.data$variables)
  data <- dplyr::group_split(data)
  return(data)
}

###################################################################
#' @title Get rules from partykit object in rule fitting
#' @param x Partykit glmtree model object
#' @param i null
#' @param ... additional arguments
#' @return List of rules
#'
#' @export
# Copied from internal partykit function
list_rules_party <- function(x, i = NULL, ...) {
  if (is.null(i)) {
    i <- partykit::nodeids(x, terminal = TRUE)
  }
  if (length(i) > 1) {
    ret <- sapply(i, list_rules_party, x = x)
    names(ret) <- if (is.character(i)) {
      i
    } else {
      names(x)[i]
    }
    return(ret)
  }
  if (is.character(i) && !is.null(names(x))) {
    i <- which(names(x) %in% i)
  }
  stopifnot(length(i) == 1 & is.numeric(i))
  stopifnot(i <= length(x) & i >= 1)
  i <- as.integer(i)
  dat <- partykit::data_party(x, i)
  if (!is.null(x$fitted)) {
    findx <- which("(fitted)" == names(dat))[1]
    dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0) {
      dat <- x$data
    }
  } else {
    dat <- x$data
  }
  rule <- c()
  rec_fun <- function(node) {
    if (partykit::id_node(node) == i) {
      return(NULL)
    }
    kid <- sapply(partykit::kids_node(node), partykit::id_node)
    whichkid <- max(which(kid <= i))
    split <- partykit::split_node(node)
    ivar <- partykit::varid_split(split)
    svar <- names(dat)[ivar]
    index <- partykit::index_split(split)
    if (is.factor(dat[, svar])) {
      if (is.null(index)) {
        index <- ((1:nlevels(dat[, svar])) > partykit::breaks_split(split)) +
          1
      }
      slevels <- levels(dat[, svar])[index == whichkid]
      srule <- paste(svar, " %in% c(\"", paste(slevels,
        collapse = "\", \"", sep = ""
      ), "\")", sep = "")
    } else {
      if (is.null(index)) {
        index <- seq_along(kid)
      }
      breaks <- cbind(c(-Inf, partykit::breaks_split(split)), c(
        partykit::breaks_split(split),
        Inf
      ))
      sbreak <- breaks[index == whichkid, ]
      right <- partykit::right_split(split)
      srule <- c()
      if (is.finite(sbreak[1])) {
        srule <- c(srule, paste(svar, ifelse(right, ">",
          ">="
        ), sbreak[1]))
      }
      if (is.finite(sbreak[2])) {
        srule <- c(srule, paste(svar, ifelse(right, "<=",
          "<"
        ), sbreak[2]))
      }
      srule <- paste(srule, collapse = " & ")
    }
    rule <<- c(rule, srule)
    return(rec_fun(node[[whichkid]]))
  }
  node <- rec_fun(partykit::node_party(x))
  paste(rule, collapse = " & ")
}
