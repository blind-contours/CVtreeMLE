library(tibble)
library(purrr)
library(Hmisc)
library(readxl)
library(tidyr)
library(MASS)
library(stringr)
library(dplyr)
library(sl3)
library(pre)
library(partykit)


# Generate continuous exposures with cut labels ---------------------------

gen_exposures <- function(mu, sigma, n, n_cuts) {
  exposures <- mvrnorm(n = n, mu = mu, Sigma = sigma)

  exposure_cut <- apply(exposures, 2, cut, breaks = n_cuts,
                        label = FALSE)

  exposure_cut_labels <- apply(exposures, 2, cut, breaks = n_cuts)

  exposures_labeled <- cbind.data.frame(exposures,
                                        apply(exposure_cut, 1, paste,
                                              collapse = " "))

  m1_key <- unique(cbind.data.frame(exposure_cut[,1],
                                    exposure_cut_labels[,1]))

  m2_key <- unique(cbind.data.frame(exposure_cut[,2],
                             exposure_cut_labels[,2]))

  colnames(m1_key) <- c("m1_labels", "m1_levels")
  colnames(m2_key) <- c("m2_labels", "m2_levels")

  colnames(exposures_labeled) <- c("M1", "M2", "Label")

  return(list("exposures_labeled" = exposures_labeled,
              "m1_key" = m1_key,
              "m2_key" = m2_key
              )
         )
}

# Generate covariates ---------------------------
gen_covariates = function(n) {
  tibble(
    age = rnorm(n, 37, 3),
    bmi = rnorm(n, 20, 1),
    sex = as.numeric(rbernoulli(n, 0.5))
  )
}

# Generate exposure cubes ---------------------------
gen_multinom <- function(exposure_grid, covars){

  n_cubes_combs <- nrow(exposure_grid)

  b0i <- round(rnorm(n_cubes_combs, 0.3, 0.01), 2)
  b1i <- round(rnorm(n_cubes_combs, 0.4, 0.01), 2)
  b2i <- round(rnorm(n_cubes_combs, 0.5, 0.01), 2)
  b3i <- round(rnorm(n_cubes_combs, 0.5, 0.01), 2)

  probs_list <- c()

  for (i in seq(nrow(covars))) {
    age <- covars$age[i]
    bmi <- covars$bmi[i]
    sex <- covars$sex[i]

    gen_denominator <- function(index, b0i, b1i, b3i, covars) {
      1 + exp(b0i[index] + (b1i[index] * age) + (b2i[index] * bmi) +
            (b3i[index] * sex))
    }

    gen_probs <- function(index, b0i, b1i, b3i, covars, denominator) {
      exp(b0i[index] + (b1i[index] * age) + (b2i[index] * bmi) +
            (b3i[index] * sex)) / (1 + denominator)
    }

    denominator <- sum(sapply(seq(from = 1, to = n_cubes_combs),
                              FUN = gen_denominator, b0i, b1i, b3i, covars))

    probs <- sapply(seq(from = 1, to = n_cubes_combs),
                        FUN = gen_probs, b0i, b1i, b3i, covars, denominator)

    probs_list[[i]] <- probs
  }

  probs_df <- as.data.frame(do.call(rbind, probs_list))
  colnames(probs_df) <- exposure_grid$labels

  res <- probs_df %>%
    dplyr::mutate(Label = rMultinom(probs_df, 1))

  res <- as.data.frame(res)

  covars_w_exp_probs <- cbind.data.frame(res, covars)

  return(covars_w_exp_probs)
}

# Assign cont. exposures to regions ---------------------------
assign_exposures <- function(covars_w_exp_probs,
                             exposure_grid,
                             exposures){

  ms <- as.data.frame(matrix(data = NA,
                             ncol = 2,
                             nrow = nrow(covars_w_exp_probs)))

  labels <- unique(exposure_grid$labels)

  for (label in labels) {
    cube_region_density <- covars_w_exp_probs[covars_w_exp_probs$Label == label ,]
    region_size <- nrow(cube_region_density)

    exposure_region <- exposures[exposures$Label == label,]
    if (nrow(exposure_region) != 0) {
      exposure_region_sample <- exposure_region[sample(nrow(exposure_region),
                                                       region_size,
                                                       replace =TRUE),
                                                       c("M1", "M2") ]
    }else{
      next
    }

    ms[covars_w_exp_probs$Label == label ,] <- exposure_region_sample
  }

  colnames(ms) <- c("m1", "m2")

  data <- cbind.data.frame(covars_w_exp_probs,
                           ms)
  return(data)

}

# Generate and assign outcomes to each region ---------------------------
assign_outcomes <- function(exposure_grid, c_matrix, data){

  gen_outcome <- function(exposure, c_matrix){
    outcome <- exposure %*% c_matrix %*% exposure
  }

  outcomes <- apply(exposure_grid[,1:2], 1, gen_outcome, c_matrix)

  empty_outcomes <- as.data.frame(matrix(data = NA, ncol = 1, nrow = n))

  cubes_outcomes <- cbind.data.frame(exposure_grid, outcomes)

  colnames(cubes_outcomes) <- c("Region 1", "Region 2", "Label", "Outcome")

  for (label in cubes_outcomes$Label) {
    outcome_val <- cubes_outcomes[cubes_outcomes$Label == label ,]$Outcome
    empty_outcomes[data$Label == label ,] <- outcome_val
  }

  empty_outcomes_confounded <- empty_outcomes + 0.2*data$age + 0.4*data$sex +
    rnorm(nrow(empty_outcomes), mean = 0, sd = 0.01)

  empty_outcomes_no_error <- empty_outcomes + 0.2*data$age + 0.4*data$sex

  data_w_outcomes <- cbind.data.frame(data,
                                      empty_outcomes_confounded,
                                      empty_outcomes_no_error,
                                      empty_outcomes,
                                      t(as.data.frame(
                                        sapply(data$Label,
                                               str_split, pattern = " "))))


  colnames(data_w_outcomes)[(ncol(data_w_outcomes)-4):ncol(data_w_outcomes)] <-
    c("outcome_obs", "outcome_true", "outcome_exposures", "region1", "region2")

  data_w_outcomes$region1 <- as.numeric(data_w_outcomes$region1)
  data_w_outcomes$region2 <- as.numeric(data_w_outcomes$region2)

  rownames(data_w_outcomes) <- NULL

  return(list("data" = data_w_outcomes,
              "cubes_outcomes" = cubes_outcomes
              )
         )

}

# Calc empirical truth treating each region as a threshold -----------------

calc_empir_truth <- function(data, rule){
  data_A <- data %>%
    dplyr::mutate(A = ifelse(eval(parse(text = rule)), 1, 0))

  data_A$A_inv <- 1- data_A$A

  all_regions <- colnames(data_A)[1:25]

  if(nrow(table(data_A$A)) == 2) {

  data_groups <- data_A %>%
    group_by(A)

  data_groups_list <- group_split(data_groups)
  regions <- as.vector(unique(data_groups_list[[2]]$Label))
  regions_inv <- as.vector(unique(data_groups_list[[1]]$Label))

  if (length(regions) > 1) {
    region_probs <- rowSums(data_A[, c(regions)])
  }else{
    region_probs <- data_A[, c(regions)]
  }

  if (length(regions_inv) > 1) {
    region_comp_probs <- rowSums(data_A[, c(regions_inv)])
  }else{
    region_comp_probs <- data_A[, c(regions_inv)]
  }

  ave_region_outcome <-
    mean(data_A$outcome_true * (data_A$A / (region_probs)))

  ave_compl_outcome <-
    mean(data_A$outcome_true * (data_A$A_inv / (region_comp_probs)))

  ate <- ave_region_outcome - ave_compl_outcome

  results_ate <- c(ave_region_outcome, ave_compl_outcome, ate)
  }else{
    results_ate <- c(NA, NA, NA)
  }
  return(results_ate)
}


