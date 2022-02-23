#' Simply function that looks at variables used in the pre model for matches with the mixture variables and pulls them out to check for consistent sets of rules found over the iterative backfitting procedure
#' @param x Row of a dataframe that contains output from the pre fit
#' @param A Vector of characters indicating column names for the mixture variables
#' @return Vector of matches of mixture variables
#' @export

pull_out_rule_vars <- function(x, A) {
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
      if (grepl("-*\\d+\\.*\\d*", element) == TRUE & grepl("[[:alpha:]]", element) == FALSE) {
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
#' @title Filter data based on fold
#' @param data Input data
#' @param fold_k Current fold

#' @importFrom rlang .data
filter_rules <- function(data, fold_k) {
  data <- data %>%
    dplyr::filter(.data$fold == fold_k)
  return(data)
}

###################################################################
#' @title Group by fold
#' @param data Input data

#' @importFrom rlang .data
groupby_fold <- function(data) {
  data <- data %>%
    dplyr::group_by(.data$fold)
  return(data)
}

###################################################################
#' @title Filter mixture rules across the folds for only those that have the same variables and directions across all folds
#' @param data Input data
#' @param n_folds Number of folds in the CV

#' @importFrom rlang .data
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
calc_mixture_rule_RMSEs <- function(data) {
  data <- data %>%
    dplyr::group_by(.data$test) %>%
    dplyr::summarize(RMSE = mean(.data$RMSE, na.rm = TRUE))
  colnames(data) <- c("Var(s)", "RMSE")

  return(data)
}
###################################################################
#' @title Filter marginal rules across the folds for only those that have the same variables
#' @param data Input data
#' @param n_folds Number of folds in the CV

#' @importFrom rlang .data
filter_marginal_rules <- function(data, n_folds) {
  data <- as.data.frame(data)
  data$fold <- as.numeric(data$fold)

  data <- data[data$rules != "No Rules Found", ]

  # data <- data %>%
  #   dplyr::group_by(.data$target_m, .data$quantile) %>%
  #   dplyr::filter(dplyr::n() >= n_folds)

  data$var_quant_group <- paste(data$target_m, data$quantile, sep = "_")
  return(data)
}

###################################################################
#' @title Calculate the mean RMSE in each marginal group
#' @param data Input data
#' @importFrom rlang .data
calc_marginal_rule_RMSEs <- function(data) {
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
groupby_fold <- function(data) {
  data <- data %>%
    dplyr::group_by(.data$fold)
  return(data)
}

###################################################################
#' @title Group split by marginal variable
#' @param data Input data

#' @importFrom rlang .data
marginal_group_split <- function(data) {
  data <- data %>% dplyr::group_by(.data$target_m)
  data <- dplyr::group_split(data)
  return(data)
}

###################################################################
#' @title v-fold marginal group split
#' @param data Input data

#' @importFrom rlang .data
v_fold_marginal_qgroup_split <- function(data) {
  data <- data %>% dplyr::group_by(.data$comparison)
  data <- dplyr::group_split(data)
  return(data)
}

###################################################################
#' @title v-fold group split
#' @param data Input data

#' @importFrom rlang .data
v_fold_mixture_group_split <- function(data) {
  data <- data %>% dplyr::group_by(.data$Variables)
  data <- dplyr::group_split(data)
  return(data)
}

###################################################################
#' @title Get rules from partykit object in rule fitting
#' @param x Partykit glmtree model object
#' @param i null
#' @param ... additional arguments
#' @return List of rules
#' @import partykit

# Copied from internal partykit function
list.rules.party <- function(x, i = NULL, ...) {
  if (is.null(i)) {
    i <- partykit::nodeids(x, terminal = TRUE)
  }
  if (length(i) > 1) {
    ret <- sapply(i, list.rules.party, x = x)
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
    fit <- dat[, findx:ncol(dat), drop = FALSE]
    dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
    if (ncol(dat) == 0) {
      dat <- x$data
    }
  } else {
    fit <- NULL
    dat <- x$data
  }
  rule <- c()
  recFun <- function(node) {
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
        index <- 1:length(kid)
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
    return(recFun(node[[whichkid]]))
  }
  node <- recFun(partykit::node_party(x))
  paste(rule, collapse = " & ")
}
