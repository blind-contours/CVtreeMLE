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
    if (length(hits) == 1) {
      hits <- 0
    } else {
      hits <- paste(hits, collapse = "-")
    }
    return(hits)
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
#' @title Flip Negative Mixture Rules
#' @param rules Row of absolute max rule results from pre
#' @param at training data to evaluate rule on
#' @param mix_comps mixture components
#' @importFrom rlang .data
flip_rule <- function(rules, at, mix_comps) {
  rule <- rules$description
  rule_split <- strsplit(rule, split = "&")[[1]]

  new_rule <- list()

  for (rule_i in rule_split) {

    var <- mix_comps[str_detect(rule_i, mix_comps)]

    intxn_data <- at %>%
      mutate(rule_i_eva = ifelse(eval(parse(text = rule_i)), 0, 1))

    var_min <-
      intxn_data %>%
      group_by(rule_i_eva) %>%
      summarise(min = min(!!(as.name(var))))

    var_min <- subset(var_min, rule_i_eva == 1, select = min)

    var_max <-
      intxn_data %>%
      group_by(rule_i_eva) %>%
      summarise(max = max(!!(as.name(var))))

    var_max <- subset(var_max, rule_i_eva == 1, select = max)

    augmented_rule <- paste(var, ">=", round(var_min, 3), "&", var,
                            "<=", round(var_max, 3))

    new_rule <- append(new_rule, augmented_rule)
  }

  inv_rule <- paste(new_rule, collapse = " | ")

  return(inv_rule)
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
